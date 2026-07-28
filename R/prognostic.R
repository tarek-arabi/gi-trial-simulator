#' Prognostic covariate adjustment
#'
#' The PROCOVA-style method: train a model on historical control-arm data to
#' predict each patient's expected outcome, freeze it, then enter that single
#' predicted value as one pre-specified covariate in the trial analysis.
#' Randomisation is untouched, so the treatment contrast stays unconfounded and
#' the type I error rate is preserved; what changes is that outcome variation
#' explained by the score no longer inflates the residual.
#'
#' The EMA CHMP Qualification Opinion on PROCOVA (September 2022) concluded that
#' prognostic score adjustment could increase power or precision in randomised
#' trials with continuous outcomes, with reported sample size reductions of
#' roughly 10 to 30 percent. That opinion is about continuous outcomes. This
#' module implements the binary-outcome case, where the arithmetic is different
#' and must not be quoted as if it were the same: see [procova_gain()] for what
#' is and is not claimed.
#'
#' @name prognostic-adjustment
NULL

# Mann-Whitney identity: AUC = P(score of a case exceeds score of a control),
# estimated by the rank-sum statistic with mid-ranks for ties.
auc_rank <- function(score, y) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) {
    return(NA_real_)
  }
  r <- rank(score)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Complete or quasi-complete separation leaves the binomial log likelihood
# monotone along some direction of the coefficient space. IRLS walks out along
# that ridge until the deviance stops moving and then reports success, so
# `glm`'s own `converged` flag is TRUE while the coefficient it returns is set
# by the iteration limit rather than by the data, and its standard error has
# exploded. stats::glm carries no test for this, so detection uses the practical
# signal: a coefficient far outside any plausible effect paired with a standard
# error larger still.
#
# Both quantities are standardised by the standard deviation of the term's own
# column of the model matrix, so the rule is scale free: a covariate measured in
# micrograms is not flagged for its units. The intercept column is constant, so
# it is compared on the linear predictor scale directly. Both conditions must
# hold, so a large but well-estimated effect is not flagged; on 400 healthy
# simulated fits spanning covariate scales from 1e-3 to 1e3 the rule fired zero
# times.
gi_separation_estimate_limit <- 8
gi_separation_se_limit <- 25

separated_terms <- function(model) {
  cf <- stats::coef(summary(model))
  if (is.null(cf) || !nrow(cf)) {
    return(character())
  }
  x <- try(stats::model.matrix(model), silent = TRUE)
  if (inherits(x, "try-error") || is.null(colnames(x))) {
    return(character())
  }
  s <- apply(x, 2L, stats::sd)
  s[is.na(s) | !is.finite(s) | s == 0] <- 1
  s <- s[rownames(cf)]
  s[is.na(s)] <- 1
  est <- abs(cf[, 1L]) * s
  se <- cf[, 2L] * s
  flag <- is.finite(est) & is.finite(se) &
    est > gi_separation_estimate_limit & se > gi_separation_se_limit
  rownames(cf)[flag]
}

# Coefficients glm reports as NA are terms it silently dropped because the
# design matrix was rank deficient. The model that was fitted is then not the
# model that was requested, and any quantity computed from the full coefficient
# vector is undefined.
aliased_terms <- function(model) {
  b <- stats::coef(model)
  if (is.null(b)) {
    return(character())
  }
  nm <- names(b)[is.na(b)]
  if (is.null(nm)) character() else nm
}

# Cross-validation folds assigned within outcome class. With unstratified folds
# a rare-event training set can put every event into one fold, which leaves that
# fold's model fitted to an outcome with no variation; the out-of-fold score is
# then not a prediction at all, and the pooled AUC it produces can sit far below
# 0.5 while looking like an honest discrimination estimate.
stratified_folds <- function(y, folds) {
  id <- integer(length(y))
  for (cls in c(0L, 1L)) {
    idx <- which(y == cls)
    id[idx] <- sample(rep_len(seq_len(folds), length(idx)))
  }
  id
}

as_binary_outcome <- function(x, arg) {
  if (inherits(x, "gi_outcome")) x <- x$y
  if (is.factor(x)) {
    if (nlevels(x) != 2L) {
      stop("`", arg, "` is a factor with ", nlevels(x), " levels; it must have 2.",
        call. = FALSE
      )
    }
    x <- as.integer(x) - 1L
  }
  if (is.logical(x)) x <- as.integer(x)
  if (!is.numeric(x) || anyNA(x) || !all(x %in% c(0, 1))) {
    stop("`", arg, "` must contain only 0 and 1 (no missing values).", call. = FALSE)
  }
  as.integer(x)
}

#' Fit a prognostic model on historical control-arm data
#'
#' @param train_cohort Data frame of covariates for the historical control-arm
#'   patients, as from [simulate_cohort()].
#' @param train_outcome Their binary outcome: a 0/1 vector, a two-level factor,
#'   or a `gi_outcome` from [outcome_model()].
#' @param formula Optional one-sided formula naming the predictors, for example
#'   `~ age + bilirubin`. Defaults to `~ .`, every column of `train_cohort`.
#' @param folds Number of cross-validation folds used for the honest
#'   discrimination estimate. Folds are assigned within outcome class, so the
#'   estimate is only produced when there are at least `folds` events and at
#'   least `folds` non-events; see the `auc_cv` note under Value.
#' @param seed Integer seed for the fold assignment.
#' @return An object of class `gi_prognostic`: a list with `model` (the fitted
#'   [stats::glm()]), `formula`, `predictors`, `n_train`, `event_rate`,
#'   `auc_apparent` (in-sample, optimistic), `auc_cv`, `oof_score`, `fold_id`,
#'   `folds`, `seed`, `cv_performed`, `separation` and `converged`.
#'
#'   `auc_cv` is the out-of-fold AUC. It is an honest discrimination estimate
#'   only when every fold holds both events and non-events, so folds are
#'   assigned within outcome class to guarantee that whenever the class sizes
#'   allow it. That needs at least `folds` events and at least `folds`
#'   non-events. `fold_id` returns the assignment actually used, so the
#'   guarantee can be checked rather than taken on trust. When the training data
#'   is too rare-event for it to hold, no cross-validation is run: `auc_cv` is
#'   `NA_real_`, `cv_performed` is `FALSE`, `oof_score` and `fold_id` are all
#'   `NA`, and a warning reports the two class counts. The alternative would be
#'   an out-of-fold score built partly from folds whose training set contained
#'   no events at all, which yields a number that is frequently below 0.5, is
#'   not a discrimination estimate, and is indistinguishable from one on
#'   inspection.
#'
#'   `converged` is `FALSE` when [stats::glm()] failed to converge **or** when
#'   the fit shows complete or quasi-complete separation, in which case
#'   `separation` names the affected terms and a warning is raised. `glm`'s own
#'   convergence flag is `TRUE` under separation, so it is not sufficient on its
#'   own. Separation is detected from the standardised coefficient and standard
#'   error, which is a heuristic and not an exact test: it is deliberately
#'   conservative, so it will not flag a large but well-estimated effect and can
#'   miss the mildest quasi-separation.
#' @seealso [prognostic_score()], [analyse_with_prognostic()], [procova_gain()]
#' @examples
#' spec <- cohort_spec(list(
#'   covariate_spec("x1", "normal", mean = 0, sd = 1),
#'   covariate_spec("x2", "normal", mean = 0, sd = 1)
#' ))
#' train <- simulate_cohort(spec, n = 1000, seed = 1)
#' y <- outcome_model(train, coefs = c(x1 = 1, x2 = 0.7), intercept = -1, seed = 1)
#' fit <- fit_prognostic(train, y)
#' c(apparent = fit$auc_apparent, cv = fit$auc_cv)
#' @export
fit_prognostic <- function(train_cohort, train_outcome, formula = NULL,
                           folds = 5L, seed = 1) {
  if (!is.data.frame(train_cohort)) {
    stop("`train_cohort` must be a data frame of covariates.", call. = FALSE)
  }
  y <- as_binary_outcome(train_outcome, "train_outcome")
  if (length(y) != nrow(train_cohort)) {
    stop(
      "`train_outcome` has length ", length(y), " but `train_cohort` has ",
      nrow(train_cohort), " rows.",
      call. = FALSE
    )
  }
  if (length(unique(y)) < 2L) {
    stop("`train_outcome` has no variation; a prognostic model cannot be fitted.",
      call. = FALSE
    )
  }
  if (!is.numeric(folds) || length(folds) != 1L || is.na(folds) || folds < 2 ||
    folds != round(folds) || folds > length(y)) {
    stop(
      "`folds` must be a single whole number between 2 and the training n (",
      length(y), ").",
      call. = FALSE
    )
  }
  folds <- as.integer(folds)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single number.", call. = FALSE)
  }
  if (".outcome" %in% names(train_cohort)) {
    stop("`train_cohort` must not contain a column named '.outcome'.", call. = FALSE)
  }

  if (is.null(formula)) formula <- ~.
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("`formula` must be a one-sided formula, for example ~ age + bilirubin.",
      call. = FALSE
    )
  }
  full <- stats::as.formula(
    paste(".outcome", paste(deparse(formula), collapse = " ")),
    env = environment(formula)
  )

  d <- cbind(data.frame(.outcome = y), train_cohort)
  model <- stats::glm(full, data = d, family = stats::binomial())
  predictors <- attr(stats::terms(model), "term.labels")

  n_event <- sum(y == 1L)
  n_nonevent <- length(y) - n_event
  cv_performed <- min(n_event, n_nonevent) >= folds
  oof <- rep(NA_real_, length(y))
  fold_id <- rep(NA_integer_, length(y))
  if (cv_performed) {
    set.seed(as.integer(seed))
    fold_id <- stratified_folds(y, folds)
    for (k in seq_len(folds)) {
      held <- fold_id == k
      mk <- stats::glm(full, data = d[!held, , drop = FALSE], family = stats::binomial())
      oof[held] <- stats::predict(mk, newdata = d[held, , drop = FALSE], type = "link")
    }
  } else {
    hint <- if (min(n_event, n_nonevent) >= 2L) {
      paste0(
        "Set `folds` to at most ", min(n_event, n_nonevent),
        ", or train on more patients."
      )
    } else {
      "Train on more patients; cross-validation needs at least two of each."
    }
    warning(
      "`auc_cv` is NA: cross-validation was skipped. The training data has ",
      n_event, " event(s) and ", n_nonevent, " non-event(s), so ", folds,
      " folds cannot each contain at least one of each. A fold whose training ",
      "set holds no events produces an out-of-fold score that is not a ",
      "prediction, and the AUC it yields is meaningless rather than honest. ",
      hint,
      call. = FALSE
    )
  }
  auc_cv <- if (cv_performed) auc_rank(oof, y) else NA_real_

  separation <- separated_terms(model)
  if (length(separation)) {
    warning(
      "The prognostic model shows complete or quasi-complete separation on ",
      paste(separation, collapse = ", "),
      ": the maximum likelihood estimate does not exist, so the reported ",
      "coefficients and standard errors are artefacts of the iteration limit ",
      "rather than estimates. `converged` is FALSE. Drop or coarsen the ",
      "predictor(s) named, or train on data where the outcome is not perfectly ",
      "predicted.",
      call. = FALSE
    )
  }

  structure(
    list(
      model = model,
      formula = full,
      predictors = predictors,
      n_train = length(y),
      event_rate = mean(y),
      auc_apparent = auc_rank(stats::predict(model, type = "link"), y),
      auc_cv = auc_cv,
      oof_score = oof,
      fold_id = fold_id,
      folds = folds,
      seed = as.integer(seed),
      cv_performed = cv_performed,
      separation = separation,
      converged = isTRUE(model$converged) && length(separation) == 0L
    ),
    class = c("gi_prognostic", "list")
  )
}

#' @export
print.gi_prognostic <- function(x, ...) {
  cat("<gi_prognostic> trained on n = ", x$n_train,
    " historical control patients\n",
    sep = ""
  )
  cat("predictors: ", paste(x$predictors, collapse = ", "), "\n", sep = "")
  cat(sprintf(
    "AUC  apparent %.3f   %d-fold cross-validated %s\n",
    x$auc_apparent, x$folds,
    if (isFALSE(x$cv_performed)) "NA (too few events to fold)" else sprintf("%.3f", x$auc_cv)
  ))
  if (length(x$separation)) {
    cat("WARNING: separation on ", paste(x$separation, collapse = ", "),
      "; those coefficients are not estimable.\n",
      sep = ""
    )
  } else if (!x$converged) {
    cat("WARNING: the fitting algorithm did not converge.\n")
  }
  invisible(x)
}

#' Score new patients with a frozen prognostic model
#'
#' @param fit A `gi_prognostic` from [fit_prognostic()].
#' @param cohort Data frame of covariates for the patients to be scored. Must
#'   contain every predictor the model uses.
#' @return Numeric vector of predicted linear predictors, one per row of
#'   `cohort`.
#' @seealso [fit_prognostic()], [analyse_with_prognostic()]
#' @examples
#' spec <- cohort_spec(list(covariate_spec("x1", "normal", mean = 0, sd = 1)))
#' train <- simulate_cohort(spec, n = 800, seed = 1)
#' y <- outcome_model(train, coefs = c(x1 = 1), intercept = -1, seed = 1)
#' fit <- fit_prognostic(train, y)
#' new <- simulate_cohort(spec, n = 5, seed = 2)
#' round(prognostic_score(fit, new), 3)
#' @export
prognostic_score <- function(fit, cohort) {
  if (!inherits(fit, "gi_prognostic")) {
    stop("`fit` must be a gi_prognostic from fit_prognostic().", call. = FALSE)
  }
  if (!is.data.frame(cohort)) {
    stop("`cohort` must be a data frame of covariates.", call. = FALSE)
  }
  # Read the predictors off the fitted terms, not the formula, so that a `~ .`
  # specification reports the variables it actually expanded to.
  needed <- all.vars(stats::delete.response(stats::terms(fit$model)))
  missing <- setdiff(needed, names(cohort))
  if (length(missing)) {
    stop(
      "`cohort` is missing predictor(s) the prognostic model needs: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  as.numeric(stats::predict(fit$model, newdata = cohort, type = "link"))
}

one_sided_p <- function(z, direction) {
  if (identical(direction, "lower_is_better")) stats::pnorm(z) else stats::pnorm(-z)
}

# Marginal risk difference by g-computation: predict every patient under both
# arms, average, difference. Its standard error is the delta method applied to
# that smooth function of the coefficient vector, holding the covariates fixed.
# Unlike the arm coefficient, this estimand is the same quantity whether or not
# the model contains the prognostic score, so the two standard errors are
# comparable. The delta-method variance is checked against the Monte Carlo
# empirical standard deviation in tests/testthat/test-prognostic.R rather than
# taken on trust.
marginal_rd <- function(model, x1, x0, direction) {
  b <- stats::coef(model)
  if (anyNA(b) || length(b) != ncol(x1)) {
    return(list(rd = NA_real_, rd_se = NA_real_, rd_z = NA_real_, rd_p_one_sided = NA_real_))
  }
  p1 <- stats::plogis(drop(x1 %*% b))
  p0 <- stats::plogis(drop(x0 %*% b))
  grad <- colMeans(p1 * (1 - p1) * x1) - colMeans(p0 * (1 - p0) * x0)
  v <- drop(crossprod(grad, stats::vcov(model) %*% grad))
  rd <- mean(p1) - mean(p0)
  se <- sqrt(max(v, 0))
  z <- rd / se
  list(rd = rd, rd_se = se, rd_z = z, rd_p_one_sided = one_sided_p(z, direction))
}

arm_summary <- function(model, direction, x1, x0) {
  cf <- stats::coef(summary(model))
  wald <- if (!"arm" %in% rownames(cf)) {
    list(estimate = NA_real_, se = NA_real_, z = NA_real_, p_one_sided = NA_real_)
  } else {
    est <- unname(cf["arm", 1L])
    se <- unname(cf["arm", 2L])
    list(
      estimate = est, se = se, z = est / se,
      p_one_sided = one_sided_p(est / se, direction)
    )
  }
  c(wald, marginal_rd(model, x1, x0, direction))
}

#' Analyse a trial with and without prognostic adjustment
#'
#' Fits the pre-specified adjusted model `outcome ~ arm + score` and, on the
#' same data, the unadjusted `outcome ~ arm`, so the two can be compared
#' directly.
#'
#' Each model is summarised on two scales. The **arm log odds ratio** is the
#' conditional effect the model estimates directly. The **marginal risk
#' difference** is obtained by g-computation, predicting every patient under
#' both arms and averaging. The second scale exists because the logistic model
#' is non-collapsible: the adjusted and unadjusted log odds ratios are estimates
#' of two different quantities, so their standard errors cannot be compared,
#' whereas the risk difference is the same estimand either way and its standard
#' errors can.
#'
#' The one-sided p value is computed in the tail implied by `direction`: for
#' `"lower_is_better"` evidence accumulates as the treatment log-odds and the
#' risk difference fall, so `p = pnorm(z)`.
#'
#' @param outcome Binary trial outcome: a 0/1 vector, a two-level factor, or a
#'   `gi_outcome`.
#' @param arm Treatment assignment: 0 for control and 1 for treatment, or a
#'   two-level factor whose second level is the treatment arm.
#' @param score Prognostic score, one value per patient, from
#'   [prognostic_score()].
#' @param direction `"lower_is_better"` (the default, appropriate for mortality
#'   and other harm endpoints) or `"higher_is_better"`.
#' @return An object of class `gi_adjusted_analysis`: a list with `adjusted` and
#'   `unadjusted`, each holding `estimate` (arm log odds ratio), `se`, `z` and
#'   `p_one_sided`, plus `rd` (marginal risk difference), `rd_se`, `rd_z` and
#'   `rd_p_one_sided`; then `se_ratio` (adjusted over unadjusted variance of the
#'   risk difference, the variance-reduction factor), `se_ratio_logor` (the same
#'   ratio for the log odds ratio, which normally exceeds 1), `z_ratio`, `n`,
#'   `n_treatment`, `n_events`, `direction`, `separation`, `rank_deficient` and
#'   `converged`.
#'
#'   `converged` is `FALSE` when either [stats::glm()] fit failed to converge,
#'   when either shows complete or quasi-complete separation, or when either
#'   design matrix was rank deficient. Each of the latter two raises a warning
#'   and fills the corresponding element:
#'
#'   * `separation` names the terms whose maximum likelihood estimate does not
#'     exist because the outcome is perfectly predicted. `glm` reports
#'     convergence in that state, so its own flag is not enough; detection here
#'     is a conservative heuristic on the standardised coefficient and standard
#'     error, not an exact test.
#'   * `rank_deficient` names the terms `glm` dropped because the design matrix
#'     was not of full rank, the usual cause being a `score` that is a linear
#'     function of `arm`. The g-computation risk difference is undefined for
#'     such a fit, so `rd`, `rd_se`, `rd_z`, `rd_p_one_sided` and both
#'     `se_ratio` values come back `NA`.
#' @seealso [procova_gain()]
#' @examples
#' set.seed(1)
#' n <- 400
#' score <- stats::rnorm(n)
#' arm <- rep(0:1, each = n / 2)
#' y <- stats::rbinom(n, 1, stats::plogis(-1 + score - 0.4 * arm))
#' res <- analyse_with_prognostic(y, arm, score)
#' c(unadjusted = res$unadjusted$rd_se, adjusted = res$adjusted$rd_se)
#' @export
analyse_with_prognostic <- function(outcome, arm, score,
                                    direction = c("lower_is_better", "higher_is_better")) {
  direction <- match.arg(direction)
  y <- as_binary_outcome(outcome, "outcome")
  a <- as_binary_outcome(arm, "arm")
  if (!is.numeric(score) || anyNA(score) || !all(is.finite(score))) {
    stop("`score` must be a numeric vector with no missing or infinite values.",
      call. = FALSE
    )
  }
  if (length(y) != length(a) || length(y) != length(score)) {
    stop(
      "`outcome`, `arm` and `score` must have the same length; got ",
      length(y), ", ", length(a), " and ", length(score), ".",
      call. = FALSE
    )
  }
  if (length(unique(a)) < 2L) {
    stop("`arm` has only one distinct value; there is no treatment contrast.",
      call. = FALSE
    )
  }

  d <- data.frame(outcome = y, arm = a, score = score)
  m_adj <- stats::glm(outcome ~ arm + score, data = d, family = stats::binomial())
  m_unadj <- stats::glm(outcome ~ arm, data = d, family = stats::binomial())

  ones <- rep(1, length(y))
  zeros <- rep(0, length(y))
  adj <- arm_summary(
    m_adj, direction,
    x1 = cbind(ones, ones, score), x0 = cbind(ones, zeros, score)
  )
  unadj <- arm_summary(
    m_unadj, direction,
    x1 = cbind(ones, ones), x0 = cbind(ones, zeros)
  )

  rank_deficient <- unique(c(aliased_terms(m_adj), aliased_terms(m_unadj)))
  if (length(rank_deficient)) {
    warning(
      "The design matrix is rank deficient: glm dropped ",
      paste(rank_deficient, collapse = ", "),
      ". The usual cause is a `score` that is a linear function of `arm`, ",
      "which leaves the treatment effect and the score indistinguishable. ",
      "The marginal risk difference is undefined for such a fit, so `rd`, ",
      "`rd_se`, `rd_z`, `rd_p_one_sided` and `se_ratio` are NA and ",
      "`converged` is FALSE.",
      call. = FALSE
    )
  }
  separation <- unique(c(separated_terms(m_adj), separated_terms(m_unadj)))
  if (length(separation)) {
    warning(
      "Complete or quasi-complete separation detected on ",
      paste(separation, collapse = ", "),
      ": the outcome is perfectly predicted, so the maximum likelihood ",
      "estimate does not exist and the reported coefficient and standard ",
      "error are artefacts of the iteration limit rather than estimates. ",
      "`converged` is FALSE. Treat the p values as uninterpretable and use an ",
      "exact or penalised method for this data.",
      call. = FALSE
    )
  }

  structure(
    list(
      adjusted = adj,
      unadjusted = unadj,
      se_ratio = adj$rd_se^2 / unadj$rd_se^2,
      se_ratio_logor = adj$se^2 / unadj$se^2,
      z_ratio = adj$z / unadj$z,
      n = length(y),
      n_treatment = sum(a == 1L),
      n_events = sum(y == 1L),
      direction = direction,
      separation = separation,
      rank_deficient = rank_deficient,
      converged = isTRUE(m_adj$converged) && isTRUE(m_unadj$converged) &&
        length(separation) == 0L && length(rank_deficient) == 0L
    ),
    class = c("gi_adjusted_analysis", "list")
  )
}

#' @export
print.gi_adjusted_analysis <- function(x, ...) {
  cat("<gi_adjusted_analysis> n = ", x$n, " (", x$n_treatment, " treated), ",
    x$n_events, " events\n",
    sep = ""
  )
  cat(sprintf(
    "%-12s %9s %8s %10s %9s %8s %10s\n",
    "", "logOR", "SE", "p (1-sided)", "riskdiff", "SE", "p (1-sided)"
  ))
  for (nm in c("unadjusted", "adjusted")) {
    r <- x[[nm]]
    cat(sprintf(
      "%-12s %9.4f %8.4f %10.4f %9.4f %8.4f %10.4f\n",
      nm, r$estimate, r$se, r$p_one_sided, r$rd, r$rd_se, r$rd_p_one_sided
    ))
  }
  if (length(x$rank_deficient)) {
    cat("WARNING: rank deficient; glm dropped ",
      paste(x$rank_deficient, collapse = ", "), ".\n",
      sep = ""
    )
  }
  if (length(x$separation)) {
    cat("WARNING: separation on ", paste(x$separation, collapse = ", "),
      "; those estimates are not interpretable.\n",
      sep = ""
    )
  }
  if (!x$converged && !length(x$rank_deficient) && !length(x$separation)) {
    cat("WARNING: at least one model did not converge.\n")
  }
  invisible(x)
}

# Monte Carlo standard error of a ratio of two means computed on the same
# replicates, by the delta method. The two numerators are strongly positively
# correlated here (same simulated data), so ignoring the covariance would badly
# overstate the uncertainty.
ratio_mcse <- function(a, b) {
  nsim <- length(a)
  ma <- mean(a)
  mb <- mean(b)
  va <- stats::var(a) / nsim
  vb <- stats::var(b) / nsim
  cab <- stats::cov(a, b) / nsim
  v <- va / mb^2 + ma^2 * vb / mb^4 - 2 * ma * cab / mb^3
  sqrt(max(v, 0))
}

procova_replicate <- function(spec, coefs, link, n_per_arm, a0, delta, fit,
                              direction, rep_seed) {
  n <- 2L * n_per_arm
  cohort <- simulate_cohort(spec, n = n, seed = rep_seed)
  arm <- sample(rep(c(0L, 1L), each = n_per_arm))
  lp <- linear_predictor(cohort, coefs, intercept = a0)
  inv <- link_inverse(link)
  p_alt <- inv(lp + delta * arm)
  p_null <- inv(lp)
  # Common random numbers: the null and alternative replicates differ only
  # through delta, which removes most of the Monte Carlo noise from any
  # comparison between them.
  u <- stats::runif(n)
  y_alt <- as.integer(u < p_alt)
  y_null <- as.integer(u < p_null)
  score <- prognostic_score(fit, cohort)

  alt <- analyse_with_prognostic(y_alt, arm, score, direction = direction)
  nul <- analyse_with_prognostic(y_null, arm, score, direction = direction)

  c(
    est_unadj = alt$unadjusted$estimate, se_unadj = alt$unadjusted$se,
    z_unadj = alt$unadjusted$z, p_unadj = alt$unadjusted$p_one_sided,
    est_adj = alt$adjusted$estimate, se_adj = alt$adjusted$se,
    z_adj = alt$adjusted$z, p_adj = alt$adjusted$p_one_sided,
    rd_unadj = alt$unadjusted$rd, rd_se_unadj = alt$unadjusted$rd_se,
    rd_p_unadj = alt$unadjusted$rd_p_one_sided,
    rd_adj = alt$adjusted$rd, rd_se_adj = alt$adjusted$rd_se,
    rd_p_adj = alt$adjusted$rd_p_one_sided,
    p_unadj_null = nul$unadjusted$p_one_sided,
    p_adj_null = nul$adjusted$p_one_sided,
    rd_p_unadj_null = nul$unadjusted$rd_p_one_sided,
    rd_p_adj_null = nul$adjusted$rd_p_one_sided,
    rate_control = mean(p_alt[arm == 0L]), rate_treatment = mean(p_alt[arm == 1L]),
    converged = as.numeric(alt$converged && nul$converged)
  )
}

#' Efficiency gain from prognostic covariate adjustment
#'
#' Simulates the same trial twice, unadjusted and adjusted for a frozen
#' prognostic score, on identical simulated data, and reports the difference.
#' The prognostic model is trained once on `train_n` simulated historical
#' control-arm patients and then held fixed across every replicate, which
#' reproduces the freezing step PROCOVA requires.
#'
#' **Every efficiency number this function reports is an upper bound, not the
#' gain a real trial should expect.** The prognostic model is fitted with
#' `stats::reformulate(names(coefs))`, that is, on exactly the covariates that
#' generate the outcome, through the same link, with nothing omitted and nothing
#' spurious added. It is correctly specified by construction and differs from
#' the true outcome model only by estimation error on `train_n` patients. That
#' is the best case, not the deployed case.
#'
#' A deployed prognostic model is not in that position. It is trained on
#' historical patients from a different time, place and standard of care; the
#' covariates it can measure are proxies for the ones that actually drive the
#' outcome; and its functional form is chosen rather than known. It is
#' misspecified relative to the trial's outcome model, so its score is less
#' prognostic and the variance reduction it achieves is smaller. Nothing here
#' estimates how much smaller, because that depends on the degree of
#' misspecification, which this simulation does not model: `coefs` fixes the
#' true outcome model and the prognostic model together, so the two cannot be
#' made to disagree through this interface. To study a misspecified score,
#' build the study by hand from [fit_prognostic()] and
#' [analyse_with_prognostic()], giving the prognostic fit a formula that omits
#' or mismeasures a true driver.
#'
#' Read `se_ratio` and `n_reduction` as the ceiling a perfectly specified score
#' could reach, report them with that qualifier attached, and plan against a
#' more conservative figure.
#'
#' Read the efficiency measures carefully, because for a binary outcome two of
#' them point in opposite directions and only one comparison is legitimate.
#'
#' `se_ratio` is the variance-reduction factor: the ratio of squared standard
#' errors of the **marginal risk difference**, adjusted over unadjusted. The
#' risk difference is the same estimand in both analyses, so this is a fair
#' comparison, and a genuinely prognostic score drives it below 1.
#' `n_reduction = 1 - se_ratio` is the implied sample size reduction, since the
#' variance of a trial estimate scales as 1 / n.
#'
#' `se_ratio_logor` is the corresponding ratio for the arm log odds ratio, and
#' for a binary outcome it normally **exceeds 1**. That is not a failure of the
#' method. The logistic model is non-collapsible, so the conditional log odds
#' ratio the adjusted model estimates is a larger quantity than the marginal log
#' odds ratio the unadjusted model estimates, and its standard error rises with
#' it; the estimate rises faster, which is why power still improves. Quoting
#' this number as evidence that adjustment hurts would be comparing the
#' precision of two different estimands.
#'
#' `information_ratio` reads an efficiency gain off the log odds ratio Wald test
#' instead: the squared ratio of the mean test statistics under the alternative,
#' `(mean(z_adjusted) / mean(z_unadjusted))^2`. Required sample size scales
#' inversely with squared noncentrality, so this is a sample size ratio too.
#'
#' It is tempting to read `information_ratio` against `1 / se_ratio` as a
#' cross-check, and the print method shows them side by side, but the agreement
#' is conditional and the condition is easy to violate. `information_ratio` is a
#' log odds ratio quantity while `se_ratio` is a risk difference quantity, and
#' the two estimands carry the same efficiency only in the limit of a vanishing
#' treatment effect, where non-collapsibility is a second-order term. Under a
#' local alternative they match: at the log odds ratio of about -0.53 the
#' package tests use, 4000 replicates put the gap at roughly 2 Monte Carlo
#' standard errors, which is noise. Away from it they do not: at a log odds
#' ratio near -1.4 `information_ratio` sits about 5 percent below `1 / se_ratio`
#' and at -2.2 about 10 percent below, in both cases tens of Monte Carlo
#' standard errors out, systematically and reproducibly. So treat a gap as
#' informative only when the treatment effect is small. At a large effect a gap
#' is the expected behaviour of two different estimands and says nothing about
#' whether either number is right.
#'
#' @param scenario A `gi_scenario` from [scenario()], supplying the control and
#'   treatment event rates, the direction of benefit and the default alpha.
#' @param spec A `gi_cohort_spec` describing the baseline covariates.
#' @param coefs Named numeric vector of true covariate coefficients on the link
#'   scale. Larger values make the score more prognostic.
#' @param n_per_arm Patients per arm in the simulated trial.
#' @param nsim Number of simulated trials.
#' @param train_n Number of historical control-arm patients used to train the
#'   prognostic model.
#' @param seed Master seed. Each replicate is given its own seed drawn from this
#'   one, so results do not depend on `workers`.
#' @param workers Number of parallel processes. Values above 1 use
#'   `parallel::mclapply` where available and fall back to serial evaluation
#'   otherwise.
#' @param alpha One-sided type I error rate. Defaults to the scenario's
#'   `defaults$alpha`, or 0.025.
#' @param link One of `"logit"`, `"probit"` or `"cloglog"`, for the true
#'   data-generating model.
#' @param folds Cross-validation folds used when reporting the prognostic
#'   model's honest discrimination. If the training cohort has fewer than
#'   `folds` events, [fit_prognostic()] warns and the reported CV AUC is `NA`;
#'   the simulation itself is unaffected.
#' @param calibrate_n Size of the reference cohort used to solve for the
#'   intercepts that reproduce the scenario's marginal event rates.
#' @return An object of class `gi_procova`: a list with the inputs, the fitted
#'   `prognostic` model, the calibrated `intercept_control`,
#'   `intercept_treatment` and `arm_effect`, then `se_ratio`, `n_reduction`,
#'   `se_ratio_logor`, `information_ratio`, and, for both the log odds ratio
#'   Wald test and the marginal risk difference test, `power_unadjusted`,
#'   `power_adjusted`, `type1_unadjusted` and `type1_adjusted` (risk difference
#'   versions carry an `rd_` prefix). Every one of those carries a matching
#'   `*_mcse` element. Also returned: mean model-based and empirical Monte Carlo
#'   standard errors of both estimators on both scales, the achieved marginal
#'   event rates, the count of non-converged replicates, and the full
#'   `replicates` matrix for further ADEMP summaries.
#'
#'   `se_ratio`, `n_reduction`, `se_ratio_logor` and `information_ratio` are
#'   upper bounds on the efficiency gain, obtained with a prognostic model that
#'   is correctly specified by construction; see the description above for why a
#'   deployed model does less well. `n_nonconverged` counts replicates whose fit
#'   failed to converge, showed separation, or had a rank deficient design
#'   matrix.
#' @seealso [fit_prognostic()], [analyse_with_prognostic()]
#' @examples
#' spec <- cohort_spec(list(
#'   covariate_spec("x1", "normal", mean = 0, sd = 1),
#'   covariate_spec("x2", "normal", mean = 0, sd = 1)
#' ))
#' gain <- procova_gain(
#'   scenario("ercp_acute_cholangitis"), spec,
#'   coefs = c(x1 = 1, x2 = 0.8),
#'   n_per_arm = 300, nsim = 40, train_n = 800, calibrate_n = 4000
#' )
#' gain$information_ratio
#' \donttest{
#' procova_gain(
#'   scenario("ercp_acute_cholangitis"), spec,
#'   coefs = c(x1 = 1, x2 = 0.8),
#'   n_per_arm = 1500, nsim = 2000
#' )
#' }
#' @export
procova_gain <- function(scenario, spec, coefs, n_per_arm, nsim = 2000,
                         train_n = 5000, seed = 1, workers = 1,
                         alpha = NULL, link = "logit", folds = 5L,
                         calibrate_n = 50000L) {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`scenario` must be a gi_scenario from scenario().", call. = FALSE)
  }
  if (!inherits(spec, "gi_cohort_spec")) {
    stop("`spec` must be a gi_cohort_spec from cohort_spec().", call. = FALSE)
  }
  if (!is.character(link) || length(link) != 1L || !link %in% gi_links) {
    stop("`link` must be one of: ", paste(gi_links, collapse = ", "), ".", call. = FALSE)
  }
  whole <- function(v, arg, min_value) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) || v < min_value) {
      stop("`", arg, "` must be a single whole number of at least ", min_value, ".",
        call. = FALSE
      )
    }
    as.integer(v)
  }
  n_per_arm <- whole(n_per_arm, "n_per_arm", 10)
  nsim <- whole(nsim, "nsim", 2)
  train_n <- whole(train_n, "train_n", 50)
  calibrate_n <- whole(calibrate_n, "calibrate_n", 1000)
  workers <- whole(workers, "workers", 1)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a single number.", call. = FALSE)
  }
  missing <- setdiff(names(coefs %||% stats::setNames(numeric(), character())), spec$names)
  if (!is.numeric(coefs) || is.null(names(coefs)) || length(missing)) {
    stop(
      "`coefs` must be a named numeric vector whose names are covariates in `spec` (",
      paste(spec$names, collapse = ", "), ").",
      if (length(missing)) paste0(" Not found: ", paste(missing, collapse = ", "), ".") else "",
      call. = FALSE
    )
  }
  alpha <- alpha %||% scenario$defaults$alpha %||% 0.025
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number strictly between 0 and 0.5.", call. = FALSE)
  }

  set.seed(as.integer(seed))
  stream <- sample.int(.Machine$integer.max, size = nsim + 2L)
  calib_seed <- stream[1L]
  train_seed <- stream[2L]
  rep_seeds <- stream[-(1:2)]

  reference <- simulate_cohort(spec, n = calibrate_n, seed = calib_seed)
  a0 <- calibrate_intercept(reference, coefs, scenario$control_rate, link = link)
  a1 <- calibrate_intercept(reference, coefs, scenario$treatment_rate, link = link)
  delta <- a1 - a0

  train_cohort <- simulate_cohort(spec, n = train_n, seed = train_seed)
  train_outcome <- outcome_model(train_cohort, coefs, a0, link = link, seed = NULL)
  fit <- fit_prognostic(
    train_cohort, train_outcome,
    formula = stats::reformulate(names(coefs)),
    folds = folds, seed = train_seed
  )

  run_one <- function(i) {
    suppressWarnings(procova_replicate(
      spec, coefs, link, n_per_arm, a0, delta, fit,
      scenario$direction, rep_seeds[i]
    ))
  }
  results <- if (workers > 1L && requireNamespace("parallel", quietly = TRUE) &&
    .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nsim), run_one, mc.cores = workers)
  } else {
    if (workers > 1L) {
      warning("Parallel evaluation is unavailable here; running serially.", call. = FALSE)
    }
    lapply(seq_len(nsim), run_one)
  }
  res <- do.call(rbind, results)

  prop_mcse <- function(p) sqrt(p * (1 - p) / nsim)
  reject <- function(column) {
    p <- mean(res[, column] < alpha)
    c(p, prop_mcse(p))
  }
  power_unadj <- reject("p_unadj")
  power_adj <- reject("p_adj")
  t1_unadj <- reject("p_unadj_null")
  t1_adj <- reject("p_adj_null")
  rd_power_unadj <- reject("rd_p_unadj")
  rd_power_adj <- reject("rd_p_adj")
  rd_t1_unadj <- reject("rd_p_unadj_null")
  rd_t1_adj <- reject("rd_p_adj_null")

  se_ratio <- mean(res[, "rd_se_adj"]^2) / mean(res[, "rd_se_unadj"]^2)
  se_ratio_mcse <- ratio_mcse(res[, "rd_se_adj"]^2, res[, "rd_se_unadj"]^2)
  se_ratio_logor <- mean(res[, "se_adj"]^2) / mean(res[, "se_unadj"]^2)
  se_ratio_logor_mcse <- ratio_mcse(res[, "se_adj"]^2, res[, "se_unadj"]^2)
  z_ratio <- mean(res[, "z_adj"]) / mean(res[, "z_unadj"])
  z_ratio_mcse <- ratio_mcse(res[, "z_adj"], res[, "z_unadj"])
  info_ratio <- z_ratio^2
  info_ratio_mcse <- 2 * abs(z_ratio) * z_ratio_mcse

  structure(
    list(
      scenario = scenario, spec = spec, coefs = coefs, link = link,
      n_per_arm = n_per_arm, n_total = 2L * n_per_arm,
      alpha = alpha, direction = scenario$direction,
      nsim = nsim, train_n = train_n, seed = as.integer(seed), workers = workers,
      prognostic = fit,
      intercept_control = a0, intercept_treatment = a1, arm_effect = delta,
      achieved_control_rate = mean(res[, "rate_control"]),
      achieved_treatment_rate = mean(res[, "rate_treatment"]),
      power_unadjusted = power_unadj[1L], power_unadjusted_mcse = power_unadj[2L],
      power_adjusted = power_adj[1L], power_adjusted_mcse = power_adj[2L],
      type1_unadjusted = t1_unadj[1L], type1_unadjusted_mcse = t1_unadj[2L],
      type1_adjusted = t1_adj[1L], type1_adjusted_mcse = t1_adj[2L],
      rd_power_unadjusted = rd_power_unadj[1L],
      rd_power_unadjusted_mcse = rd_power_unadj[2L],
      rd_power_adjusted = rd_power_adj[1L],
      rd_power_adjusted_mcse = rd_power_adj[2L],
      rd_type1_unadjusted = rd_t1_unadj[1L],
      rd_type1_unadjusted_mcse = rd_t1_unadj[2L],
      rd_type1_adjusted = rd_t1_adj[1L],
      rd_type1_adjusted_mcse = rd_t1_adj[2L],
      se_ratio = se_ratio, se_ratio_mcse = se_ratio_mcse,
      n_reduction = 1 - se_ratio, n_reduction_mcse = se_ratio_mcse,
      se_ratio_logor = se_ratio_logor, se_ratio_logor_mcse = se_ratio_logor_mcse,
      information_ratio = info_ratio, information_ratio_mcse = info_ratio_mcse,
      est_unadjusted = mean(res[, "est_unadj"]),
      est_adjusted = mean(res[, "est_adj"]),
      mean_se_unadjusted = mean(res[, "se_unadj"]),
      mean_se_adjusted = mean(res[, "se_adj"]),
      emp_se_unadjusted = stats::sd(res[, "est_unadj"]),
      emp_se_adjusted = stats::sd(res[, "est_adj"]),
      rd_unadjusted = mean(res[, "rd_unadj"]),
      rd_adjusted = mean(res[, "rd_adj"]),
      mean_rd_se_unadjusted = mean(res[, "rd_se_unadj"]),
      mean_rd_se_adjusted = mean(res[, "rd_se_adj"]),
      emp_rd_se_unadjusted = stats::sd(res[, "rd_unadj"]),
      emp_rd_se_adjusted = stats::sd(res[, "rd_adj"]),
      n_nonconverged = sum(res[, "converged"] == 0),
      replicates = res
    ),
    class = c("gi_procova", "list")
  )
}

#' @export
print.gi_procova <- function(x, ...) {
  cat("<gi_procova> ", x$scenario$pack_id, " / ", x$scenario$endpoint, "\n", sep = "")
  cat(sprintf(
    "%d patients per arm, %d simulations, one-sided alpha %.4g\n",
    x$n_per_arm, x$nsim, x$alpha
  ))
  cat(sprintf(
    "prognostic score: %d-fold CV AUC %s (trained on n = %d)\n",
    x$prognostic$folds,
    if (isFALSE(x$prognostic$cv_performed)) {
      "NA (too few events to fold)"
    } else {
      sprintf("%.3f", x$prognostic$auc_cv)
    },
    x$prognostic$n_train
  ))
  cat("                          unadjusted        adjusted\n")
  cat(sprintf(
    "power        logOR     %6.3f (%.3f)   %6.3f (%.3f)\n",
    x$power_unadjusted, x$power_unadjusted_mcse,
    x$power_adjusted, x$power_adjusted_mcse
  ))
  cat(sprintf(
    "power        riskdiff  %6.3f (%.3f)   %6.3f (%.3f)\n",
    x$rd_power_unadjusted, x$rd_power_unadjusted_mcse,
    x$rd_power_adjusted, x$rd_power_adjusted_mcse
  ))
  cat(sprintf(
    "type I error logOR     %6.4f (%.4f)  %6.4f (%.4f)\n",
    x$type1_unadjusted, x$type1_unadjusted_mcse,
    x$type1_adjusted, x$type1_adjusted_mcse
  ))
  cat(sprintf(
    "type I error riskdiff  %6.4f (%.4f)  %6.4f (%.4f)\n",
    x$rd_type1_unadjusted, x$rd_type1_unadjusted_mcse,
    x$rd_type1_adjusted, x$rd_type1_adjusted_mcse
  ))
  cat(sprintf(
    "variance-reduction factor %.3f (%.3f) on the risk difference,\n",
    x$se_ratio, x$se_ratio_mcse
  ))
  cat(sprintf(
    "  implying a %.1f%% sample size reduction at equal power, an upper bound\n",
    100 * x$n_reduction
  ))
  cat("  obtained with a correctly specified prognostic model.\n")
  cat(sprintf(
    "local-alternative cross-check: information ratio %.3f (%.3f) versus\n",
    x$information_ratio, x$information_ratio_mcse
  ))
  cat(sprintf(
    "  1 / variance ratio %.3f; these agree only for a small treatment effect,\n",
    1 / x$se_ratio
  ))
  cat("  so a gap at a large effect is expected, not a fault.\n")
  cat(sprintf(
    "log odds ratio SE ratio %.3f (%.3f); above 1 is expected for a binary\n",
    x$se_ratio_logor, x$se_ratio_logor_mcse
  ))
  cat("  outcome because the adjusted model targets a different, larger estimand.\n")
  if (x$n_nonconverged > 0L) {
    cat("WARNING: ", x$n_nonconverged, " replicate(s) had a non-converged fit.\n", sep = "")
  }
  invisible(x)
}

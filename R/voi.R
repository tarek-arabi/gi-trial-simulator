#' Value of information for a proposed trial
#'
#' Power is the wrong question to stop at. A design can be well powered and
#' still not be worth running, and a design can be underpowered for the
#' conventional target and still resolve the decision that clinicians actually
#' face. Value of information analysis asks the second question directly: given
#' what is already believed about the treatment effect, how much is it worth to
#' learn more, and how much of that value does a trial of a given size actually
#' deliver.
#'
#' Three quantities appear here. The value under current information is the net
#' benefit of the best decision that could be taken today. The expected value of
#' perfect information is what would be gained if the treatment effect were
#' revealed exactly. The expected value of sample information is what would be
#' gained from a trial of a stated size, which is always less than perfect
#' information and approaches it slowly.
#'
#' All of these are denominated in whatever unit the user's net benefit function
#' returns, per patient exposed to the decision. Deaths averted, quality
#' adjusted life years and currency are all valid, and none of them is supplied
#' by this package: the net benefit function is the user's own judgement about
#' what the decision is worth, and every number below inherits its assumptions.
#' A value of information analysis with a carelessly written net benefit
#' function is a precise answer to a question nobody asked.
#'
#' @name voi
#' @importFrom stats quantile rbinom sd
NULL

# Deliberately duplicated in R/emulator.R so the two modules stay independent of
# each other's internals.
voi_with_seed <- function(seed, expr) {
  has_old <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has_old) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

voi_draws_matrix <- function(prior_draws, arg = "prior_draws") {
  if (is.data.frame(prior_draws)) prior_draws <- data.matrix(prior_draws)
  if (is.numeric(prior_draws) && is.null(dim(prior_draws))) {
    prior_draws <- matrix(prior_draws, ncol = 1L, dimnames = list(NULL, "theta"))
  }
  if (!is.matrix(prior_draws) || !is.numeric(prior_draws)) {
    stop("`", arg, "` must be a numeric vector, matrix or data frame of parameter draws.",
      call. = FALSE
    )
  }
  if (nrow(prior_draws) < 2L) {
    stop("`", arg, "` must hold at least two draws.", call. = FALSE)
  }
  if (!all(is.finite(prior_draws))) {
    stop("`", arg, "` must contain only finite values.", call. = FALSE)
  }
  prior_draws
}

voi_net_benefit_matrix <- function(draws, net_benefit_fn, one_dimensional) {
  if (!is.function(net_benefit_fn)) {
    stop("`net_benefit_fn` must be a function of one parameter draw.", call. = FALSE)
  }
  n <- nrow(draws)
  # A named scalar would push its own name into the returned option names.
  take <- function(i) if (one_dimensional) unname(draws[i, 1L]) else draws[i, ]
  first <- net_benefit_fn(take(1L))
  if (!is.numeric(first) || !is.null(dim(first)) || length(first) < 2L) {
    stop(
      "`net_benefit_fn` must return a plain numeric vector with one net benefit ",
      "per decision option, and there must be at least two options to choose ",
      "between. It returned an object of length ", length(first), ".",
      call. = FALSE
    )
  }
  k <- length(first)
  option_names <- names(first)
  if (is.null(option_names)) {
    option_names <- paste0("option", seq_len(k))
  } else {
    # c(label = value) appends the name of `value` to `label` when the draw was
    # passed in named, which turns a decision option into "label.parameter".
    for (cn in colnames(draws)) {
      option_names <- sub(paste0("\\.", cn, "$"), "", option_names)
    }
  }
  nb <- matrix(NA_real_, nrow = n, ncol = k, dimnames = list(NULL, option_names))
  nb[1L, ] <- first
  for (i in seq_len(n)[-1L]) {
    v <- net_benefit_fn(take(i))
    if (!is.numeric(v) || length(v) != k) {
      stop(
        "`net_benefit_fn` returned ", length(v), " values at draw ", i,
        " but ", k, " at draw 1; it must return the same options every time.",
        call. = FALSE
      )
    }
    nb[i, ] <- v
  }
  if (!all(is.finite(nb))) {
    stop("`net_benefit_fn` returned non-finite net benefits.", call. = FALSE)
  }
  nb
}

#' Expected value of perfect information
#'
#' The gain from resolving all uncertainty about the treatment effect before
#' deciding. Under current information the best that can be done is to take the
#' single decision with the highest expected net benefit. Under perfect
#' information the best decision could be taken separately for every possible
#' truth. EVPI is the difference, and it is the ceiling on what any trial,
#' however large, could be worth.
#'
#' EVPI is estimated draw by draw as
#' `mean(max over options) - max over options of mean`, which is exactly the
#' mean of the per draw opportunity loss of the decision that is optimal under
#' current information. That identity is what makes the estimate non negative by
#' construction and gives the Monte Carlo standard error below.
#'
#' @param prior_draws Draws from the current belief about the parameters that
#'   drive the decision: a numeric vector for a single parameter, or a matrix or
#'   data frame with one row per draw and one column per parameter. These are
#'   assumed to be independent draws from the current distribution, equally
#'   weighted.
#' @param net_benefit_fn Function taking one parameter draw (a scalar when
#'   `prior_draws` is a vector, otherwise the row as a numeric vector) and
#'   returning a numeric vector of net benefits, one per decision option. Name
#'   the returned vector to have the options named in the output.
#'
#' @return An object of class `gi_evpi`, a list with `evpi`, its Monte Carlo
#'   standard error `mcse`, `n_draws`, `options`, `value_current_info`,
#'   `value_perfect_info`, `best_option_current` and `prob_optimal`, the
#'   proportion of draws for which each option is the best one.
#'
#' @note EVPI is per patient exposed to the decision, in the units of
#'   `net_benefit_fn`, and inherits every assumption in it. The Monte Carlo
#'   standard error reflects only the finite number of prior draws, not any
#'   uncertainty about whether the prior itself is right.
#'
#' @seealso [evsi_trial()], [voi_curve()], [population_evsi()]
#' @examples
#' set.seed(1)
#' # A lognormal prior on the risk ratio for 30 day mortality, centred on the
#' # 0.60 the ercp_acute_cholangitis pack implies, turned into risk differences
#' # against a control rate of 0.0658. Negative means lives saved.
#' delta <- 0.0658 * (rlnorm(2000, meanlog = log(0.6), sdlog = 0.35) - 1)
#' nb <- function(d) c(standard_care = 0, early_ercp = -d - 0.005)
#' evpi(delta, nb)
#' @export
evpi <- function(prior_draws, net_benefit_fn) {
  one_d <- is.numeric(prior_draws) && is.null(dim(prior_draws))
  draws <- voi_draws_matrix(prior_draws)
  voi_evpi_from_nb(voi_net_benefit_matrix(draws, net_benefit_fn, one_d))
}

voi_evpi_from_nb <- function(nb) {
  mean_nb <- colMeans(nb)
  best <- which.max(mean_nb)
  best_per_draw <- do.call(pmax, as.data.frame(nb))
  loss <- best_per_draw - nb[, best]
  n <- nrow(nb)

  structure(
    list(
      evpi = mean(loss),
      mcse = stats::sd(loss) / sqrt(n),
      n_draws = n,
      options = colnames(nb),
      value_current_info = mean_nb[[best]],
      value_perfect_info = mean(best_per_draw),
      best_option_current = colnames(nb)[best],
      prob_optimal = prop.table(table(factor(
        colnames(nb)[max.col(nb, ties.method = "first")],
        levels = colnames(nb)
      )))
    ),
    class = c("gi_evpi", "list")
  )
}

#' Expected value of sample information for a proposed trial
#'
#' What a trial of a stated size is worth, in the same units as [evpi()]. For
#' each of `nsim` simulated trials the routine draws one parameter value from
#' the current belief, simulates the two arm binary outcome data that trial would
#' produce under that value, updates the belief by weighting the prior draws by
#' their likelihood of having produced those data, and takes the decision that
#' maximises expected net benefit under the updated belief. The value realised
#' by that trial is the net benefit of the decision it leads to, less the net
#' benefit of the decision that would have been taken without it, both evaluated
#' under the same updated belief. EVSI is the average of those differences.
#'
#' Differencing within each simulated trial rather than differencing two
#' separately estimated averages gives the same quantity, because the expected
#' net benefit of a fixed decision under the posterior averages over datasets to
#' its value under the prior. It is a much better estimator: the large common
#' component of net benefit cancels within each simulated trial, which cuts the
#' Monte Carlo error by roughly a factor of four in the examples here, and the
#' estimate cannot come out negative. The unpaired form is reported under
#' `detail$evsi_unpaired` and the two agree as `nsim` grows.
#'
#' The posterior is formed by likelihood reweighting of the supplied prior draws
#' rather than by a separate conjugate model, so that the belief being updated
#' and the belief being valued are the same object. This keeps EVSI below EVPI
#' as it must be, but it degrades when the trial is large relative to the number
#' of prior draws, because then a handful of draws carry all the posterior
#' weight. The effective sample size of the weights is reported for exactly that
#' reason: if `ess_mean` is small, supply more prior draws.
#'
#' @param prior_draws Either a numeric vector of draws of the absolute risk
#'   difference `p_treatment - p_control` (so a negative value favours the
#'   intervention when the endpoint is a harm such as death), or a two column
#'   matrix or data frame of draws of `p_control` and `p_treatment` in that
#'   order.
#' @param n_per_arm Number of patients per arm in the proposed trial, assuming
#'   one to one allocation.
#' @param control_rate Assumed event rate in the control arm, required when
#'   `prior_draws` is a vector of risk differences and not allowed when
#'   `prior_draws` already carries both arm rates.
#' @param net_benefit_fn Function taking one parameter draw and returning a
#'   numeric vector of net benefits, one per decision option, exactly as in
#'   [evpi()]. It receives the risk difference as a scalar, or the pair of rates
#'   as a numeric vector, matching the form of `prior_draws`.
#' @param nsim Number of simulated trials in the outer loop.
#' @param seed Integer seed, recorded in the return value. The ambient random
#'   number state is restored on exit.
#'
#' @return An object of class `gi_evsi`, a list with `evsi`, its Monte Carlo
#'   standard error `mcse`, the `evpi` for the same prior and net benefit
#'   function, `n_per_arm`, `nsim`, `seed`, `n_draws`, `options`,
#'   `value_current_info`, and `detail`, which holds two alternative estimators
#'   of the same quantity (`evsi_unpaired`, the average value under sample
#'   information minus the value under current information, each estimated
#'   separately; and `evsi_plugin`, using the net benefit at the parameter draw
#'   that generated the data rather than its posterior expectation), their
#'   standard errors, the weight effective sample sizes `ess_mean`, `ess_q10`
#'   and `ess_min`, the count `n_clipped` of prior draws clipped to the probability
#'   scale, and `decisions`, how often each option was chosen.
#'
#' @note Assumptions, stated plainly. One to one allocation and a binary
#'   endpoint observed on every patient, so no dropout and no interim looks.
#'   When `prior_draws` is a vector, the control rate is treated as known, which
#'   makes the control arm uninformative about the effect and so overstates what
#'   is learned per patient a little; supply both arm rates to avoid that. Any
#'   draw implying an arm rate outside (0, 1) is clipped for the purpose of
#'   simulating outcomes, with a warning, and the net benefit function still
#'   sees the unclipped draw. EVSI is per patient exposed to the decision and
#'   inherits every assumption in `net_benefit_fn`.
#'
#' @seealso [evpi()], [voi_curve()], [population_evsi()]
#' @examples
#' set.seed(1)
#' delta <- 0.0658 * (rlnorm(4000, meanlog = log(0.6), sdlog = 0.35) - 1)
#' nb <- function(d) c(standard_care = 0, early_ercp = -d - 0.005)
#' evsi_trial(delta, n_per_arm = 250, control_rate = 0.0658,
#'            net_benefit_fn = nb, nsim = 200)
#' @export
evsi_trial <- function(prior_draws, n_per_arm, control_rate = NULL,
                       net_benefit_fn, nsim = 2000, seed = 1) {
  if (!is.numeric(n_per_arm) || length(n_per_arm) != 1L || !is.finite(n_per_arm) ||
    n_per_arm < 1) {
    stop("`n_per_arm` must be a single positive number of patients per arm.", call. = FALSE)
  }
  n_per_arm <- as.integer(round(n_per_arm))
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) || nsim < 2) {
    stop("`nsim` must be a single number of at least 2.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be a single finite number.", call. = FALSE)
  }

  one_d <- is.numeric(prior_draws) && is.null(dim(prior_draws))
  draws <- voi_draws_matrix(prior_draws)

  if (one_d) {
    if (is.null(control_rate)) {
      stop(
        "`control_rate` is required when `prior_draws` is a vector of risk ",
        "differences, because the sampling model needs the control arm rate.",
        call. = FALSE
      )
    }
    if (!is.numeric(control_rate) || length(control_rate) != 1L ||
      is.na(control_rate) || control_rate <= 0 || control_rate >= 1) {
      stop("`control_rate` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
    p_control <- rep(control_rate, nrow(draws))
    p_treatment <- control_rate + draws[, 1L]
  } else {
    if (ncol(draws) != 2L) {
      stop(
        "`prior_draws` must be a vector of risk differences or a two column ",
        "matrix of control and treatment rates, not ", ncol(draws), " columns.",
        call. = FALSE
      )
    }
    if (!is.null(control_rate)) {
      stop(
        "`control_rate` must be NULL when `prior_draws` already carries the ",
        "control arm rate in its first column.",
        call. = FALSE
      )
    }
    p_control <- draws[, 1L]
    p_treatment <- draws[, 2L]
  }

  eps <- 1e-6
  n_clipped <- sum(p_control < eps | p_control > 1 - eps |
    p_treatment < eps | p_treatment > 1 - eps)
  if (n_clipped > 0L) {
    warning(
      n_clipped, " of ", nrow(draws),
      " prior draws imply an event rate outside (0, 1); they were clipped to ",
      "simulate outcomes. Consider a prior that respects the probability scale.",
      call. = FALSE
    )
  }
  p_control <- pmin(pmax(p_control, eps), 1 - eps)
  p_treatment <- pmin(pmax(p_treatment, eps), 1 - eps)

  nb <- voi_net_benefit_matrix(draws, net_benefit_fn, one_d)
  n_draws <- nrow(nb)
  mean_nb <- colMeans(nb)
  best_now <- which.max(mean_nb)
  value_current <- mean_nb[[best_now]]

  log_pc <- log(p_control)
  log_qc <- log1p(-p_control)
  log_pt <- log(p_treatment)
  log_qt <- log1p(-p_treatment)

  sim <- voi_with_seed(seed, {
    idx <- sample.int(n_draws, nsim, replace = TRUE)
    list(
      idx = idx,
      x_control = stats::rbinom(nsim, n_per_arm, p_control[idx]),
      x_treatment = stats::rbinom(nsim, n_per_arm, p_treatment[idx])
    )
  })

  gain <- numeric(nsim)
  value <- numeric(nsim)
  plugin <- numeric(nsim)
  ess <- numeric(nsim)
  chosen <- integer(nsim)
  for (s in seq_len(nsim)) {
    xc <- sim$x_control[s]
    xt <- sim$x_treatment[s]
    # The binomial coefficient does not depend on the draw, so it cancels in
    # the normalised weights and is left out of the log likelihood.
    loglik <- xc * log_pc + (n_per_arm - xc) * log_qc +
      xt * log_pt + (n_per_arm - xt) * log_qt
    w <- exp(loglik - max(loglik))
    w <- w / sum(w)
    post_nb <- colSums(w * nb)
    d <- which.max(post_nb)
    chosen[s] <- d
    gain[s] <- post_nb[[d]] - post_nb[[best_now]]
    value[s] <- post_nb[[d]]
    plugin[s] <- nb[sim$idx[s], d]
    ess[s] <- 1 / sum(w^2)
  }

  ess_mean <- mean(ess)
  ess_q10 <- unname(stats::quantile(ess, 0.1))
  if (ess_mean < 50 || ess_q10 < 20) {
    warning(
      "The posterior weights have a mean effective sample size of ",
      format(ess_mean, digits = 3), " and a tenth percentile of ",
      format(ess_q10, digits = 3), " out of ", n_draws,
      " prior draws, so this EVSI rests on very few draws. Supply more prior ",
      "draws before quoting it.",
      call. = FALSE
    )
  }

  ref <- voi_evpi_from_nb(nb)

  structure(
    list(
      evsi = mean(gain),
      mcse = stats::sd(gain) / sqrt(nsim),
      evpi = ref$evpi,
      n_per_arm = n_per_arm,
      nsim = nsim,
      seed = seed,
      n_draws = n_draws,
      options = colnames(nb),
      value_current_info = value_current,
      detail = list(
        evsi_unpaired = mean(value) - value_current,
        mcse_unpaired = stats::sd(value) / sqrt(nsim),
        evsi_plugin = mean(plugin) - value_current,
        mcse_plugin = stats::sd(plugin) / sqrt(nsim),
        ess_mean = ess_mean,
        ess_q10 = ess_q10,
        ess_min = min(ess),
        n_clipped = n_clipped,
        decisions = prop.table(table(factor(
          colnames(nb)[chosen],
          levels = colnames(nb)
        )))
      )
    ),
    class = c("gi_evsi", "list")
  )
}

#' Expected value of sample information across trial sizes
#'
#' Runs [evsi_trial()] over a grid of per arm sample sizes and returns the curve
#' of information value against trial size. The shape is the point: EVSI rises
#' steeply while the trial is small, then flattens as the decision becomes
#' settled, and the sample size beyond which the curve is flat is buying
#' precision that changes no decision. The EVPI line is the asymptote the curve
#' can never cross.
#'
#' Every point on the curve uses the same seed, so the simulated trials share
#' their parameter draws. That is a common random numbers design: it makes
#' differences between neighbouring sample sizes far less noisy than the
#' individual points, which is what the curve is read for.
#'
#' @param prior_draws Draws of the treatment effect, as in [evsi_trial()].
#' @param n_grid Numeric vector of per arm sample sizes to evaluate.
#' @param ... Passed to [evsi_trial()]: `control_rate`, `net_benefit_fn`,
#'   `nsim` and `seed`.
#'
#' @return A data frame with one row per sample size and columns `n_per_arm`,
#'   `evsi`, `mcse`, `evsi_unpaired`, `ess_mean`, `evpi` and `fraction_of_evpi`.
#'   The column names are chosen to be readable directly by [plot_evsi()].
#'
#' @note The curve is a Monte Carlo estimate at every point, so it need not be
#'   exactly monotone even though the quantity it estimates is, and once it
#'   approaches saturation `fraction_of_evpi` can print slightly above 1. Both
#'   are simulation noise: read the mcse column before reading a dip, or a
#'   fraction above 1, as real.
#'
#' @seealso [evsi_trial()], [evpi()], [population_evsi()]
#' @examples
#' set.seed(1)
#' delta <- 0.0658 * (rlnorm(3000, meanlog = log(0.6), sdlog = 0.35) - 1)
#' nb <- function(d) c(standard_care = 0, early_ercp = -d - 0.005)
#' voi_curve(delta, n_grid = c(100, 400), control_rate = 0.0658,
#'           net_benefit_fn = nb, nsim = 100)
#' @export
voi_curve <- function(prior_draws, n_grid, ...) {
  if (!is.numeric(n_grid) || length(n_grid) < 1L || any(!is.finite(n_grid)) ||
    any(n_grid < 1)) {
    stop("`n_grid` must be a numeric vector of positive per arm sample sizes.", call. = FALSE)
  }
  n_grid <- sort(unique(as.integer(round(n_grid))))
  rows <- lapply(n_grid, function(n) {
    fit <- evsi_trial(prior_draws, n_per_arm = n, ...)
    data.frame(
      n_per_arm = fit$n_per_arm,
      evsi = fit$evsi,
      mcse = fit$mcse,
      evsi_unpaired = fit$detail$evsi_unpaired,
      ess_mean = fit$detail$ess_mean,
      evpi = fit$evpi,
      fraction_of_evpi = if (fit$evpi > 0) fit$evsi / fit$evpi else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Scale per patient information value to a population
#'
#' A trial is run once and its answer is used on every patient who presents
#' afterwards, so the value of the trial is the per patient value multiplied by
#' the number of patients whose care it will inform, discounted because value
#' arriving in ten years is worth less than value arriving now.
#'
#' The arithmetic is deliberately transparent: incidence is assumed constant
#' over the horizon, the first year is undiscounted, and year `t` is discounted
#' by `(1 + discount_rate)^-t` with `t` running from 0. Nothing here accounts
#' for the years the trial itself takes to report, during which the current
#' decision continues to be made in ignorance, nor for the possibility that the
#' answer is overtaken by a change in practice. Both shorten the effective
#' horizon, so this figure is an upper bound on the population value.
#'
#' @param evsi_per_patient Per patient expected value of sample information, as
#'   returned in the `evsi` element of [evsi_trial()], in whatever unit the net
#'   benefit function used.
#' @param incidence Number of patients per year facing the decision in the
#'   population of interest.
#' @param horizon_years Number of years over which the trial result will inform
#'   care.
#' @param discount_rate Annual discount rate as a proportion, defaulting to 0.03,
#'   the conventional health economic rate. Use 0 for an undiscounted total.
#'
#' @return A list with `total` (the population value), `per_patient`,
#'   `effective_population` (the discounted patient count that `total` divides
#'   by), `undiscounted_population`, the inputs `incidence`, `horizon_years` and
#'   `discount_rate`, and `annual`, a data frame with one row per year giving
#'   `year`, `patients`, `discount_factor`, `discounted_patients` and `value`.
#'
#' @seealso [evsi_trial()], [voi_curve()]
#' @examples
#' population_evsi(evsi_per_patient = 0.0021, incidence = 45000,
#'                 horizon_years = 10)
#' @export
population_evsi <- function(evsi_per_patient, incidence, horizon_years,
                            discount_rate = 0.03) {
  if (!is.numeric(evsi_per_patient) || length(evsi_per_patient) != 1L ||
    !is.finite(evsi_per_patient)) {
    stop("`evsi_per_patient` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(incidence) || length(incidence) != 1L || !is.finite(incidence) ||
    incidence <= 0) {
    stop("`incidence` must be a single positive number of patients per year.", call. = FALSE)
  }
  if (!is.numeric(horizon_years) || length(horizon_years) != 1L ||
    !is.finite(horizon_years) || horizon_years < 1) {
    stop("`horizon_years` must be a single number of years, at least 1.", call. = FALSE)
  }
  horizon_years <- as.integer(round(horizon_years))
  if (!is.numeric(discount_rate) || length(discount_rate) != 1L ||
    !is.finite(discount_rate) || discount_rate < 0 || discount_rate >= 1) {
    stop("`discount_rate` must be a single number in [0, 1).", call. = FALSE)
  }

  year <- seq_len(horizon_years)
  factor_t <- (1 + discount_rate)^-(year - 1)
  discounted <- incidence * factor_t

  list(
    total = evsi_per_patient * sum(discounted),
    per_patient = evsi_per_patient,
    effective_population = sum(discounted),
    undiscounted_population = incidence * horizon_years,
    incidence = incidence,
    horizon_years = horizon_years,
    discount_rate = discount_rate,
    annual = data.frame(
      year = year,
      patients = rep(incidence, horizon_years),
      discount_factor = factor_t,
      discounted_patients = discounted,
      value = evsi_per_patient * discounted
    )
  )
}

#' @export
print.gi_evpi <- function(x, ...) {
  cat("<gi_evpi> expected value of perfect information\n")
  cat("  EVPI per patient: ", format(x$evpi, digits = 4),
    "  (MCSE ", format(x$mcse, digits = 2), ", ", x$n_draws, " draws)\n",
    sep = ""
  )
  cat("  best option under current information: ", x$best_option_current, "\n", sep = "")
  cat("  probability each option is optimal: ",
    paste(sprintf("%s %.3f", names(x$prob_optimal), as.numeric(x$prob_optimal)),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
print.gi_evsi <- function(x, ...) {
  cat("<gi_evsi> expected value of sample information\n")
  cat("  trial size: ", x$n_per_arm, " per arm, ", x$nsim,
    " simulated trials, seed ", x$seed, "\n",
    sep = ""
  )
  cat("  EVSI per patient: ", format(x$evsi, digits = 4),
    "  (MCSE ", format(x$mcse, digits = 2), ")\n",
    sep = ""
  )
  cat("  EVPI per patient: ", format(x$evpi, digits = 4),
    "   EVSI is ", sprintf("%.1f", 100 * x$evsi / x$evpi), " percent of it\n",
    sep = ""
  )
  cat("  posterior weight effective sample size: mean ",
    format(x$detail$ess_mean, digits = 3), ", minimum ",
    format(x$detail$ess_min, digits = 3), " of ", x$n_draws, " draws\n",
    sep = ""
  )
  invisible(x)
}

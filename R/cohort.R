#' Virtual patient cohorts
#'
#' Trial planning frequently needs patient-level data when only aggregate
#' summaries are publishable: marginal distributions of baseline covariates and
#' a correlation matrix between them. These functions build a virtual cohort
#' from exactly that aggregate input, using a Gaussian copula, so no
#' patient-level source data is required or implied.
#'
#' @name virtual-cohorts
NULL

gi_dists <- c("normal", "lognormal", "gamma", "binary", "ordinal")

#' Describe the marginal distribution of one covariate
#'
#' @param name Covariate name. Becomes the column name in the simulated cohort.
#' @param dist One of `"normal"`, `"lognormal"`, `"gamma"`, `"binary"` or
#'   `"ordinal"`.
#' @param ... Distribution parameters. `normal` takes `mean` and `sd`;
#'   `lognormal` takes `meanlog` and `sdlog`; `gamma` takes `shape` and one of
#'   `rate` or `scale`; `binary` takes `prob`; `ordinal` takes `probs` (a vector
#'   of category probabilities summing to 1) and optionally `values` (the
#'   numeric code assigned to each category, defaulting to `1:length(probs)`).
#' @return An object of class `gi_covariate`: a list with elements `name`,
#'   `dist` and `pars` (the validated parameter list).
#' @seealso [cohort_spec()], [simulate_cohort()]
#' @examples
#' covariate_spec("age", "normal", mean = 68, sd = 13)
#' covariate_spec("bilirubin", "lognormal", meanlog = 1.2, sdlog = 0.6)
#' covariate_spec("charlson", "ordinal", probs = c(0.4, 0.35, 0.25))
#' @export
covariate_spec <- function(name, dist, ...) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(dist) || length(dist) != 1L || !dist %in% gi_dists) {
    stop(
      "`dist` must be one of: ", paste(gi_dists, collapse = ", "),
      ". Got: ", paste(as.character(dist), collapse = ", "),
      call. = FALSE
    )
  }
  pars <- list(...)
  pars <- switch(dist,
    normal = validate_normal_pars(name, pars),
    lognormal = validate_lognormal_pars(name, pars),
    gamma = validate_gamma_pars(name, pars),
    binary = validate_binary_pars(name, pars),
    ordinal = validate_ordinal_pars(name, pars)
  )
  structure(
    list(name = name, dist = dist, pars = pars),
    class = c("gi_covariate", "list")
  )
}

positive_scalar <- function(value, arg, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
    stop(
      "Covariate '", name, "': `", arg, "` must be a single positive number.",
      call. = FALSE
    )
  }
  as.numeric(value)
}

finite_scalar <- function(value, arg, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(
      "Covariate '", name, "': `", arg, "` must be a single finite number.",
      call. = FALSE
    )
  }
  as.numeric(value)
}

reject_unknown_pars <- function(pars, allowed, name, dist) {
  extra <- setdiff(names(pars), allowed)
  if (length(extra)) {
    stop(
      "Covariate '", name, "': unknown parameter(s) for dist '", dist, "': ",
      paste(extra, collapse = ", "), ". Allowed: ", paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_normal_pars <- function(name, pars) {
  reject_unknown_pars(pars, c("mean", "sd"), name, "normal")
  list(
    mean = finite_scalar(pars$mean %||% 0, "mean", name),
    sd = positive_scalar(pars$sd %||% 1, "sd", name)
  )
}

validate_lognormal_pars <- function(name, pars) {
  reject_unknown_pars(pars, c("meanlog", "sdlog"), name, "lognormal")
  list(
    meanlog = finite_scalar(pars$meanlog %||% 0, "meanlog", name),
    sdlog = positive_scalar(pars$sdlog %||% 1, "sdlog", name)
  )
}

validate_gamma_pars <- function(name, pars) {
  reject_unknown_pars(pars, c("shape", "rate", "scale"), name, "gamma")
  if (!is.null(pars$rate) && !is.null(pars$scale)) {
    stop(
      "Covariate '", name, "': give only one of `rate` or `scale`, not both.",
      call. = FALSE
    )
  }
  rate <- if (!is.null(pars$scale)) {
    1 / positive_scalar(pars$scale, "scale", name)
  } else {
    positive_scalar(pars$rate %||% 1, "rate", name)
  }
  list(shape = positive_scalar(pars$shape, "shape", name), rate = rate)
}

validate_binary_pars <- function(name, pars) {
  reject_unknown_pars(pars, "prob", name, "binary")
  prob <- pars$prob
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) || prob <= 0 || prob >= 1) {
    stop(
      "Covariate '", name,
      "': `prob` must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }
  list(prob = as.numeric(prob))
}

validate_ordinal_pars <- function(name, pars) {
  reject_unknown_pars(pars, c("probs", "values"), name, "ordinal")
  probs <- pars$probs
  if (!is.numeric(probs) || length(probs) < 2L || anyNA(probs) || any(probs <= 0)) {
    stop(
      "Covariate '", name,
      "': `probs` must be a numeric vector of at least two strictly positive ",
      "category probabilities.",
      call. = FALSE
    )
  }
  if (abs(sum(probs) - 1) > 1e-8) {
    stop(
      "Covariate '", name, "': `probs` must sum to 1; they sum to ",
      format(sum(probs), digits = 8), ".",
      call. = FALSE
    )
  }
  values <- pars$values %||% seq_along(probs)
  if (!is.numeric(values) || length(values) != length(probs) || anyNA(values)) {
    stop(
      "Covariate '", name, "': `values` must be a numeric vector of length ",
      length(probs), " (one code per category).",
      call. = FALSE
    )
  }
  if (is.unsorted(values, strictly = TRUE)) {
    stop(
      "Covariate '", name,
      "': `values` must be strictly increasing so the marginal stays ordinal.",
      call. = FALSE
    )
  }
  list(probs = as.numeric(probs), values = as.numeric(values))
}

#' @export
print.gi_covariate <- function(x, ...) {
  pars <- vapply(
    names(x$pars),
    function(k) paste0(k, "=", paste(format(x$pars[[k]], digits = 4), collapse = ",")),
    character(1)
  )
  cat("<gi_covariate> ", x$name, " ~ ", x$dist,
    "(", paste(pars, collapse = ", "), ")\n",
    sep = ""
  )
  invisible(x)
}

#' Describe a correlated covariate structure
#'
#' @param covariates A list of `gi_covariate` objects (or a single one).
#' @param correlation Target correlation matrix on the latent normal scale.
#'   Either name both dimensions, in which case the matrix may be in any order
#'   and is reordered to match `covariates`, or name neither, in which case it
#'   is read positionally in the order of `covariates`. A matrix carrying only
#'   `rownames` or only `colnames` is rejected: its intended ordering is
#'   ambiguous, and reading it positionally would silently pair covariates with
#'   the wrong correlations. Defaults to the identity, meaning independent
#'   covariates.
#' @return An object of class `gi_cohort_spec`: a list with `covariates`,
#'   `correlation` (named, validated) and `names`.
#' @seealso [simulate_cohort()], [cohort_correlation_check()]
#' @examples
#' spec <- cohort_spec(
#'   list(
#'     covariate_spec("age", "normal", mean = 68, sd = 13),
#'     covariate_spec("bilirubin", "lognormal", meanlog = 1.2, sdlog = 0.6),
#'     covariate_spec("sepsis", "binary", prob = 0.28)
#'   ),
#'   correlation = matrix(
#'     c(
#'       1.0, 0.15, 0.20,
#'       0.15, 1.0, 0.45,
#'       0.20, 0.45, 1.0
#'     ),
#'     nrow = 3, byrow = TRUE
#'   )
#' )
#' spec
#' @export
cohort_spec <- function(covariates, correlation = NULL) {
  if (inherits(covariates, "gi_covariate")) covariates <- list(covariates)
  if (!is.list(covariates) || length(covariates) == 0L) {
    stop("`covariates` must be a non-empty list of gi_covariate objects.", call. = FALSE)
  }
  ok <- vapply(covariates, inherits, logical(1), what = "gi_covariate")
  if (!all(ok)) {
    stop(
      "`covariates` element(s) ", paste(which(!ok), collapse = ", "),
      " are not gi_covariate objects; build them with covariate_spec().",
      call. = FALSE
    )
  }
  nms <- vapply(covariates, function(z) z$name, character(1))
  if (anyDuplicated(nms)) {
    stop(
      "`covariates` has duplicated name(s): ",
      paste(unique(nms[duplicated(nms)]), collapse = ", "),
      call. = FALSE
    )
  }
  names(covariates) <- nms
  p <- length(nms)

  if (is.null(correlation)) {
    correlation <- diag(p)
  }
  correlation <- validate_correlation(correlation, nms)

  structure(
    list(covariates = covariates, correlation = correlation, names = nms),
    class = c("gi_cohort_spec", "list")
  )
}

validate_correlation <- function(correlation, nms) {
  p <- length(nms)
  if (!is.matrix(correlation) || !is.numeric(correlation)) {
    stop("`correlation` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(correlation) != ncol(correlation)) {
    stop(
      "`correlation` must be square; it is ", nrow(correlation), " by ",
      ncol(correlation), ".",
      call. = FALSE
    )
  }
  if (nrow(correlation) != p) {
    stop(
      "`correlation` must be ", p, " by ", p, " to match the ", p,
      " covariates; it is ", nrow(correlation), " by ", ncol(correlation), ".",
      call. = FALSE
    )
  }
  if (anyNA(correlation)) {
    stop("`correlation` contains NA.", call. = FALSE)
  }
  has_rownames <- !is.null(rownames(correlation))
  has_colnames <- !is.null(colnames(correlation))
  # One set of names on its own is refused rather than dropped. Reading such a
  # matrix positionally would pair covariates with the wrong correlations
  # without saying so, and honouring the one set that is present would require
  # assuming the missing dimension follows the same order, which is exactly the
  # assumption that goes wrong.
  if (xor(has_rownames, has_colnames)) {
    stop(
      "`correlation` carries ",
      if (has_rownames) "rownames but no colnames" else "colnames but no rownames",
      ". Supply both sets of dimnames or neither. With only one set the ",
      "intended order of the other dimension is unknown, so the matrix would ",
      "have to be read positionally and could silently pair covariates with ",
      "the wrong correlations. Positional order is: ",
      paste(nms, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (has_rownames && has_colnames) {
    if (!setequal(rownames(correlation), nms) || !setequal(colnames(correlation), nms)) {
      stop(
        "`correlation` dimnames do not match the covariate names (",
        paste(nms, collapse = ", "), ").",
        call. = FALSE
      )
    }
    correlation <- correlation[nms, nms, drop = FALSE]
  }
  dimnames(correlation) <- list(nms, nms)

  asym <- which(abs(correlation - t(correlation)) > 1e-8, arr.ind = TRUE)
  if (nrow(asym)) {
    i <- asym[1L, "row"]
    j <- asym[1L, "col"]
    stop(
      "`correlation` is not symmetric: entry [", nms[i], ", ", nms[j], "] is ",
      format(correlation[i, j], digits = 6), " but [", nms[j], ", ", nms[i],
      "] is ", format(correlation[j, i], digits = 6), ".",
      call. = FALSE
    )
  }
  bad_diag <- which(abs(diag(correlation) - 1) > 1e-8)
  if (length(bad_diag)) {
    stop(
      "`correlation` must have 1 on its diagonal; entry [", nms[bad_diag[1L]],
      ", ", nms[bad_diag[1L]], "] is ",
      format(diag(correlation)[bad_diag[1L]], digits = 6), ".",
      call. = FALSE
    )
  }
  if (any(correlation < -1 | correlation > 1)) {
    stop("`correlation` has entries outside [-1, 1].", call. = FALSE)
  }
  ev <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) {
    stop(
      "`correlation` is not positive definite: its smallest eigenvalue is ",
      format(min(ev), digits = 6),
      ". No joint distribution has this correlation matrix; ",
      "check the pairwise values for mutual inconsistency.",
      call. = FALSE
    )
  }
  correlation
}

#' @export
print.gi_cohort_spec <- function(x, ...) {
  cat("<gi_cohort_spec> ", length(x$names), " covariates\n", sep = "")
  for (cv in x$covariates) print(cv)
  off <- x$correlation[upper.tri(x$correlation)]
  cat("target correlation: ",
    if (length(off) == 0L) "none (single covariate)" else {
      sprintf(
        "%d pair(s), |r| from %.2f to %.2f",
        length(off), min(abs(off)), max(abs(off))
      )
    },
    "\n",
    sep = ""
  )
  invisible(x)
}

marginal_quantile <- function(cv, u) {
  switch(cv$dist,
    normal = stats::qnorm(u, mean = cv$pars$mean, sd = cv$pars$sd),
    lognormal = stats::qlnorm(u, meanlog = cv$pars$meanlog, sdlog = cv$pars$sdlog),
    gamma = stats::qgamma(u, shape = cv$pars$shape, rate = cv$pars$rate),
    binary = as.numeric(u > 1 - cv$pars$prob),
    ordinal = {
      breaks <- cumsum(cv$pars$probs)
      breaks <- breaks[-length(breaks)]
      cv$pars$values[findInterval(u, breaks) + 1L]
    }
  )
}

#' Simulate a virtual cohort by Gaussian copula
#'
#' Draws a latent multivariate normal with the specified correlation using a
#' Cholesky factor, maps it to uniforms with [stats::pnorm()], then through each
#' covariate's quantile function.
#'
#' Known limitation, stated plainly: the correlation supplied in [cohort_spec()]
#' is imposed on the **latent normal** scale. For any non-normal marginal the
#' realised product-moment correlation of the observed variables is attenuated
#' towards zero relative to that target, severely so for binary marginals with
#' rates far from 0.5. This is a property of the Gaussian copula, not a bug, and
#' it is not corrected here. Use [cohort_correlation_check()] to see the size of
#' the attenuation for your own specification before relying on it.
#'
#' @param spec A `gi_cohort_spec` from [cohort_spec()].
#' @param n Number of patients to simulate.
#' @param seed Integer seed passed to [set.seed()]. Use `NULL` to draw from the
#'   current RNG stream instead, which is what you want when calling this inside
#'   a loop that has already set a seed for the replicate.
#' @return A data frame with `n` rows and one column per covariate, carrying
#'   attributes `seed` and `spec`.
#' @seealso [cohort_spec()], [cohort_correlation_check()], [outcome_model()]
#' @examples
#' spec <- cohort_spec(
#'   list(
#'     covariate_spec("age", "normal", mean = 68, sd = 13),
#'     covariate_spec("sepsis", "binary", prob = 0.28)
#'   ),
#'   correlation = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
#' )
#' head(simulate_cohort(spec, n = 5, seed = 42))
#' @export
simulate_cohort <- function(spec, n, seed = 1) {
  if (!inherits(spec, "gi_cohort_spec")) {
    stop("`spec` must be a gi_cohort_spec from cohort_spec().", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != round(n)) {
    stop("`n` must be a single positive whole number.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single number, or NULL to use the current RNG stream.",
        call. = FALSE
      )
    }
    set.seed(as.integer(seed))
  }

  p <- length(spec$names)
  latent <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(spec$correlation)
  u <- stats::pnorm(latent)

  cols <- lapply(seq_len(p), function(j) marginal_quantile(spec$covariates[[j]], u[, j]))
  out <- as.data.frame(cols, stringsAsFactors = FALSE)
  names(out) <- spec$names
  attr(out, "seed") <- seed
  attr(out, "spec") <- spec
  out
}

#' Compare target and realised covariate correlations
#'
#' Reports, pair by pair, the correlation asked for on the latent normal scale
#' against the correlation actually achieved on the observed scale. The gap is
#' the Gaussian copula's attenuation for non-normal marginals.
#'
#' @param spec A `gi_cohort_spec`.
#' @param n Number of patients to simulate for the check. Large values give a
#'   sharper estimate of the attenuation.
#' @param seed Integer seed for the check.
#' @return A data frame with one row per covariate pair: `var1`, `var2`,
#'   `target`, `realised` and `attenuation` (target minus realised).
#'
#'   Errors, naming the covariate, if any covariate comes out constant in the
#'   simulated cohort, since its correlations are then undefined. That happens
#'   when a marginal is extreme relative to `n`, for example a binary covariate
#'   whose `prob` is small enough that `n` draws contain a single category.
#' @seealso [simulate_cohort()]
#' @examples
#' spec <- cohort_spec(
#'   list(
#'     covariate_spec("x", "normal", mean = 0, sd = 1),
#'     covariate_spec("z", "binary", prob = 0.1)
#'   ),
#'   correlation = matrix(c(1, 0.6, 0.6, 1), nrow = 2)
#' )
#' cohort_correlation_check(spec, n = 5000)
#' @export
cohort_correlation_check <- function(spec, n = 20000, seed = 1) {
  if (!inherits(spec, "gi_cohort_spec")) {
    stop("`spec` must be a gi_cohort_spec from cohort_spec().", call. = FALSE)
  }
  if (length(spec$names) < 2L) {
    stop("`spec` must contain at least two covariates to have a correlation.",
      call. = FALSE
    )
  }
  cohort <- simulate_cohort(spec, n = n, seed = seed)
  # A covariate that came out constant has no correlation with anything.
  # stats::cor would return NA for its whole row and column behind base R's bare
  # "the standard deviation is zero" warning, which names nothing and is easy to
  # read past in a loop.
  spread <- vapply(cohort, function(v) stats::sd(as.numeric(v)), numeric(1))
  degenerate <- names(cohort)[is.na(spread) | !is.finite(spread) | spread == 0]
  if (length(degenerate)) {
    stop(
      "Covariate(s) ", paste0("'", degenerate, "'", collapse = ", "),
      " took a single constant value in the simulated cohort of n = ", n,
      ", so every correlation involving them is undefined. This happens when a ",
      "marginal is extreme enough that n draws contain only one category, for ",
      "example a binary covariate with a very small `prob`. Raise `n`, or move ",
      "the marginal away from its boundary.",
      call. = FALSE
    )
  }
  realised <- stats::cor(as.matrix(cohort))
  pairs <- utils::combn(spec$names, 2L)
  target <- spec$correlation[cbind(pairs[1L, ], pairs[2L, ])]
  got <- realised[cbind(pairs[1L, ], pairs[2L, ])]
  data.frame(
    var1 = pairs[1L, ],
    var2 = pairs[2L, ],
    target = as.numeric(target),
    realised = as.numeric(got),
    attenuation = as.numeric(target) - as.numeric(got),
    stringsAsFactors = FALSE
  )
}

# Ceiling of the current implementation: a plain Gaussian copula matches the
# target only on the latent scale. The standard upgrade is NORTA (Cario and
# Nelson 1997), which iteratively perturbs the latent correlation until the
# realised observed correlation hits the target. That is a strictly larger piece
# of machinery and is deliberately not implemented here; cohort_correlation_check
# exists so the shortfall is visible rather than hidden.

gi_links <- c("logit", "probit", "cloglog")

link_inverse <- function(link) {
  switch(link,
    logit = stats::plogis,
    probit = stats::pnorm,
    cloglog = function(eta) -expm1(-exp(eta))
  )
}

linear_predictor <- function(cohort, coefs, intercept) {
  if (!is.data.frame(cohort)) {
    stop("`cohort` must be a data frame, as returned by simulate_cohort().", call. = FALSE)
  }
  if (!is.numeric(coefs) || is.null(names(coefs)) || any(!nzchar(names(coefs)))) {
    stop("`coefs` must be a named numeric vector, one name per covariate used.",
      call. = FALSE
    )
  }
  missing <- setdiff(names(coefs), names(cohort))
  if (length(missing)) {
    stop(
      "`coefs` names not found in `cohort`: ", paste(missing, collapse = ", "),
      ". Cohort columns are: ", paste(names(cohort), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(intercept) || length(intercept) != 1L || !is.finite(intercept)) {
    stop("`intercept` must be a single finite number.", call. = FALSE)
  }
  x <- as.matrix(cohort[, names(coefs), drop = FALSE])
  if (!is.numeric(x)) {
    stop("All covariates named in `coefs` must be numeric.", call. = FALSE)
  }
  as.vector(intercept + x %*% coefs)
}

#' Generate a binary outcome from cohort covariates
#'
#' Gives the virtual cohort an outcome with real, known covariate signal, so a
#' prognostic model has something to learn and the strength of that signal is
#' under the analyst's control.
#'
#' @param cohort A data frame of covariates, from [simulate_cohort()].
#' @param coefs Named numeric vector of covariate coefficients on the link
#'   scale. Names must be columns of `cohort`.
#' @param intercept Intercept on the link scale. Use [calibrate_intercept()] to
#'   choose one that reproduces a target marginal event rate.
#' @param link One of `"logit"`, `"probit"` or `"cloglog"`.
#' @param seed Integer seed for the Bernoulli draws, or `NULL` to use the
#'   current RNG stream.
#' @return An object of class `gi_outcome`: a list with `y` (the simulated 0/1
#'   outcome), `p` (each patient's true event probability), `eta` (the linear
#'   predictor), `event_rate` (the true marginal rate, `mean(p)`),
#'   `observed_rate` (`mean(y)`), plus `intercept`, `coefs`, `link`, `n` and
#'   `seed`.
#' @seealso [calibrate_intercept()], [fit_prognostic()]
#' @examples
#' spec <- cohort_spec(list(covariate_spec("age", "normal", mean = 0, sd = 1)))
#' cohort <- simulate_cohort(spec, n = 500, seed = 3)
#' out <- outcome_model(cohort, coefs = c(age = 0.8), intercept = -2, seed = 3)
#' out$event_rate
#' @export
outcome_model <- function(cohort, coefs, intercept, link = "logit", seed = 1) {
  if (!is.character(link) || length(link) != 1L || !link %in% gi_links) {
    stop(
      "`link` must be one of: ", paste(gi_links, collapse = ", "),
      ". Got: ", paste(as.character(link), collapse = ", "),
      call. = FALSE
    )
  }
  eta <- linear_predictor(cohort, coefs, intercept)
  p <- link_inverse(link)(eta)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single number, or NULL to use the current RNG stream.",
        call. = FALSE
      )
    }
    set.seed(as.integer(seed))
  }
  y <- stats::rbinom(length(p), size = 1L, prob = p)
  structure(
    list(
      y = y, p = p, eta = eta,
      event_rate = mean(p), observed_rate = mean(y),
      intercept = intercept, coefs = coefs, link = link,
      n = length(y), seed = seed
    ),
    class = c("gi_outcome", "list")
  )
}

#' @export
print.gi_outcome <- function(x, ...) {
  cat("<gi_outcome> n = ", x$n, ", link = ", x$link, "\n", sep = "")
  cat(sprintf(
    "  true marginal event rate  %.4f\n  observed event rate       %.4f\n",
    x$event_rate, x$observed_rate
  ))
  invisible(x)
}

#' Find the intercept that yields a target marginal event rate
#'
#' Solves for the intercept, holding the covariate coefficients fixed, so that
#' the average predicted event probability over the cohort equals `target_rate`.
#' This is what lets a virtual cohort be calibrated to the control event rate a
#' parameter pack cites. The mean predicted probability is strictly increasing
#' in the intercept, so the root is unique.
#'
#' @param cohort A data frame of covariates, from [simulate_cohort()].
#' @param coefs Named numeric vector of covariate coefficients on the link
#'   scale.
#' @param target_rate Target marginal event rate, strictly between 0 and 1.
#' @param link One of `"logit"`, `"probit"` or `"cloglog"`.
#' @return A single number: the intercept on the link scale.
#' @seealso [outcome_model()]
#' @examples
#' spec <- cohort_spec(list(covariate_spec("age", "normal", mean = 0, sd = 1)))
#' cohort <- simulate_cohort(spec, n = 5000, seed = 7)
#' a <- calibrate_intercept(cohort, coefs = c(age = 0.8), target_rate = 0.0658)
#' round(a, 3)
#' @export
calibrate_intercept <- function(cohort, coefs, target_rate, link = "logit") {
  if (!is.character(link) || length(link) != 1L || !link %in% gi_links) {
    stop(
      "`link` must be one of: ", paste(gi_links, collapse = ", "),
      ". Got: ", paste(as.character(link), collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(target_rate) || length(target_rate) != 1L || is.na(target_rate) ||
    target_rate <= 0 || target_rate >= 1) {
    stop("`target_rate` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  lp <- linear_predictor(cohort, coefs, intercept = 0)
  inv <- link_inverse(link)
  gap <- function(a) mean(inv(a + lp)) - target_rate
  stats::uniroot(
    gap,
    interval = c(-20, 20), extendInt = "upX",
    tol = .Machine$double.eps^0.5
  )$root
}

#' Gaussian process emulation of expensive simulation output
#'
#' Monte Carlo evaluation of a trial design is expensive: a single design point
#' can cost minutes of compute, and a design space of even three inputs has more
#' points than anyone will ever simulate. An emulator is the standard computer
#' experiment answer to that problem (Sacks and colleagues 1989, Kennedy and
#' O'Hagan 2001): run the expensive simulator at a modest set of design points,
#' fit a Gaussian process to the outputs, and use the fitted process to
#' interpolate the rest of the space together with an honest statement of how
#' uncertain that interpolation is.
#'
#' The emulator here is a surrogate for the simulation layer of this package. It
#' is not a design engine and it does not compete with `rpact` or `gsDesign`. It
#' interpolates numbers those engines, or a Monte Carlo loop over them, already
#' produced.
#'
#' An emulator interpolates simulation output. It does not replace simulation.
#' Any operating characteristic reported in a paper, a protocol or a regulatory
#' document must come from an actual simulation run at that design point. The
#' emulator is for searching, screening and deciding where to spend the next
#' block of compute, and its posterior mean at an unvisited point is a guess
#' with a standard deviation attached, not a result.
#'
#' @name emulator
#' @importFrom stats optim predict runif sd
NULL

emulator_kernels <- c("matern52", "squared_exponential")

# Deliberately duplicated in R/voi.R so the two modules stay independent of each
# other's internals.
emulator_with_seed <- function(seed, expr) {
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

emulator_as_matrix <- function(x, arg, d = NULL) {
  if (is.data.frame(x)) x <- data.matrix(x)
  if (is.numeric(x) && is.null(dim(x))) {
    if (is.null(d) || identical(d, 1L)) {
      x <- matrix(x, ncol = 1L)
    } else if (length(x) == d) {
      x <- matrix(x, nrow = 1L)
    } else {
      stop(
        "`", arg, "` is a plain vector of length ", length(x),
        " but the emulator has ", d, " inputs. Supply a matrix with ", d,
        " columns, one row per point.",
        call. = FALSE
      )
    }
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", arg, "` must be a numeric matrix, data frame or numeric vector.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop("`", arg, "` must contain only finite values.", call. = FALSE)
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  x
}

emulator_sqdist <- function(a, b, lengthscale) {
  za <- sweep(a, 2L, lengthscale, "/")
  zb <- sweep(b, 2L, lengthscale, "/")
  d2 <- outer(rowSums(za^2), rowSums(zb^2), "+") - 2 * tcrossprod(za, zb)
  d2[d2 < 0] <- 0
  d2
}

emulator_correlation <- function(a, b, lengthscale, kernel) {
  d2 <- emulator_sqdist(a, b, lengthscale)
  if (identical(kernel, "squared_exponential")) {
    return(exp(-0.5 * d2))
  }
  r <- sqrt(d2)
  root5 <- sqrt(5)
  (1 + root5 * r + (5 / 3) * d2) * exp(-root5 * r)
}

# Cholesky with a retry ladder: a correlation matrix built from near duplicate
# design points is positive definite in theory and indefinite in floating point,
# so the nugget is raised until the factorisation succeeds.
emulator_chol <- function(R, nugget, max_tries = 10L) {
  eps <- nugget
  for (i in seq_len(max_tries)) {
    A <- R
    diag(A) <- diag(A) + eps
    L <- tryCatch(chol(A), error = function(e) NULL)
    if (!is.null(L)) return(list(chol = L, nugget = eps))
    eps <- if (eps <= 0) 1e-10 else eps * 10
  }
  NULL
}

emulator_profile_fit <- function(lengthscale, nugget, xs, ys, kernel) {
  R <- emulator_correlation(xs, xs, lengthscale, kernel)
  fac <- emulator_chol(R, nugget)
  if (is.null(fac)) return(NULL)
  L <- fac$chol
  alpha <- backsolve(L, backsolve(L, ys, transpose = TRUE))
  n <- length(ys)
  sigma2 <- sum(ys * alpha) / n
  if (!is.finite(sigma2) || sigma2 <= 0) return(NULL)
  logdet <- 2 * sum(log(diag(L)))
  loglik <- -0.5 * (n * log(sigma2) + logdet + n * (1 + log(2 * pi)))
  if (!is.finite(loglik)) return(NULL)
  list(
    chol = L, alpha = alpha, sigma2 = sigma2, loglik = loglik,
    nugget = fac$nugget, lengthscale = lengthscale
  )
}

#' Fit a Gaussian process emulator to simulation output
#'
#' Fits a zero mean Gaussian process to centred and scaled simulation output,
#' with either a Matern 5/2 or a squared exponential correlation function and a
#' separate lengthscale per input dimension. The signal variance is profiled out
#' of the log marginal likelihood analytically; the lengthscales and the nugget,
#' when not supplied, are estimated by maximising the exact profile log marginal
#' likelihood with [stats::optim()] from several starting values.
#'
#' Inputs are standardised internally (each column centred at its mean and
#' scaled by its standard deviation) and the transformation is recorded, so
#' `newdata` in [predict.gi_emulator()] and `bounds` in [suggest_next_point()]
#' are always given on the user's own scale. Lengthscales are reported on the
#' user's scale as well.
#'
#' @param x Design inputs: a numeric matrix with one row per simulated design
#'   point and one column per input (for example control rate, effect size,
#'   sample size per arm). A data frame or, for a single input, a plain numeric
#'   vector is also accepted.
#' @param y Numeric vector of simulator output at those design points, one value
#'   per row of `x` (for example simulated power).
#' @param kernel Correlation function, either `"matern52"` (the default, twice
#'   differentiable and the usual choice for output from a stochastic simulator)
#'   or `"squared_exponential"` (infinitely smooth, which can be optimistic).
#' @param nugget Non negative number added to the diagonal of the correlation
#'   matrix, so it is expressed as a fraction of the signal variance, that is as
#'   a noise to signal ratio. Use a small value such as `1e-10` for an
#'   interpolating emulator of deterministic output, a larger one for output
#'   carrying Monte Carlo noise. `NULL`, the default, estimates it.
#' @param lengthscale Correlation lengthscale on the scale of `x`, either one
#'   number applied to every input or one per input. `NULL`, the default,
#'   estimates them.
#'
#' @return An object of class `gi_emulator`, a list with elements `x`, `y`, `n`,
#'   `d`, `kernel`, `lengthscale` (named, on the scale of `x`), `nugget`,
#'   `sigma2` (signal variance on the standardised output scale), `loglik` (the
#'   profile log marginal likelihood at the fit), `estimated` (which
#'   hyperparameters were estimated rather than supplied), the standardisation
#'   constants `x_center`, `x_scale`, `y_center`, `y_scale`, and the internal
#'   quantities `x_scaled`, `chol` and `alpha` used by the predict method.
#'
#' @note An emulator interpolates simulation output, it does not replace
#'   simulation. Any operating characteristic quoted in a paper must come from
#'   an actual simulation run at that design point, not from the emulator.
#'   Because the process mean is fixed at the sample mean of `y` rather than
#'   estimated jointly, predictive intervals far outside the design region are
#'   mildly optimistic; they are also the region where the emulator should not
#'   be trusted in the first place.
#'
#' @seealso [predict.gi_emulator()], [emulator_loo()], [suggest_next_point()]
#' @examples
#' x <- matrix(seq(0, 1, length.out = 10), ncol = 1)
#' y <- sin(3 * x[, 1]) + 0.2 * x[, 1]^2
#' fit <- fit_emulator(x, y)
#' fit$lengthscale
#' predict(fit, newdata = matrix(c(0.15, 0.45), ncol = 1))
#' @export
fit_emulator <- function(x, y, kernel = "matern52", nugget = NULL, lengthscale = NULL) {
  if (length(kernel) != 1L || !is.character(kernel) || !kernel %in% emulator_kernels) {
    stop(
      "`kernel` must be one of ", paste(emulator_kernels, collapse = ", "),
      ", not '", paste(kernel, collapse = ", "), "'.",
      call. = FALSE
    )
  }
  x <- emulator_as_matrix(x, "x")
  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("`y` must be a plain numeric vector.", call. = FALSE)
  }
  if (length(y) != nrow(x)) {
    stop(
      "`y` has length ", length(y), " but `x` has ", nrow(x),
      " rows; there must be one output value per design point.",
      call. = FALSE
    )
  }
  if (!all(is.finite(y))) stop("`y` must contain only finite values.", call. = FALSE)
  n <- nrow(x)
  d <- ncol(x)
  if (n < d + 2L) {
    stop(
      "`x` has ", n, " design points for ", d,
      " inputs, which is too few to fit an emulator. Simulate at least ",
      d + 2L, " points.",
      call. = FALSE
    )
  }

  x_center <- colMeans(x)
  x_scale <- apply(x, 2L, stats::sd)
  if (any(x_scale <= 0)) {
    bad <- colnames(x)[x_scale <= 0]
    stop(
      "`x` column(s) ", paste(bad, collapse = ", "),
      " are constant, so no lengthscale is identified. Drop them.",
      call. = FALSE
    )
  }
  y_center <- mean(y)
  y_scale <- stats::sd(y)
  if (y_scale <= 0) {
    stop("`y` is constant, so there is nothing to emulate.", call. = FALSE)
  }
  xs <- sweep(sweep(x, 2L, x_center, "-"), 2L, x_scale, "/")
  ys <- (y - y_center) / y_scale

  if (!is.null(nugget)) {
    if (!is.numeric(nugget) || length(nugget) != 1L || !is.finite(nugget) || nugget < 0) {
      stop("`nugget` must be a single non negative finite number, or NULL.", call. = FALSE)
    }
  }
  ls_scaled <- NULL
  if (!is.null(lengthscale)) {
    if (!is.numeric(lengthscale) || !all(is.finite(lengthscale)) || any(lengthscale <= 0)) {
      stop("`lengthscale` must be positive and finite, or NULL.", call. = FALSE)
    }
    if (length(lengthscale) == 1L) lengthscale <- rep(lengthscale, d)
    if (length(lengthscale) != d) {
      stop(
        "`lengthscale` must have length 1 or ", d,
        " (one per input), not ", length(lengthscale), ".",
        call. = FALSE
      )
    }
    ls_scaled <- lengthscale / x_scale
  }

  free_ls <- is.null(ls_scaled)
  free_nug <- is.null(nugget)
  ls_bounds <- log(c(1e-2, 1e3))
  nug_bounds <- log(c(1e-10, 1))

  if (free_ls || free_nug) {
    unpack <- function(par) {
      ls_i <- if (free_ls) exp(par[seq_len(d)]) else ls_scaled
      nug_i <- if (free_nug) exp(par[length(par)]) else nugget
      list(lengthscale = ls_i, nugget = nug_i)
    }
    objective <- function(par) {
      p <- unpack(par)
      f <- emulator_profile_fit(p$lengthscale, p$nugget, xs, ys, kernel)
      if (is.null(f)) 1e10 else -f$loglik
    }
    lower <- c(if (free_ls) rep(ls_bounds[1], d), if (free_nug) nug_bounds[1])
    upper <- c(if (free_ls) rep(ls_bounds[2], d), if (free_nug) nug_bounds[2])
    starts <- lapply(c(0.2, 0.5, 1, 3), function(m) {
      c(if (free_ls) rep(log(m), d), if (free_nug) log(1e-6))
    })
    best <- NULL
    for (s0 in starts) {
      opt <- tryCatch(
        stats::optim(
          par = s0, fn = objective, method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(maxit = 300)
        ),
        error = function(e) NULL
      )
      if (!is.null(opt) && is.finite(opt$value) && (is.null(best) || opt$value < best$value)) {
        best <- opt
      }
    }
    if (is.null(best)) {
      stop(
        "Could not maximise the log marginal likelihood for this design. ",
        "Supply `lengthscale` and `nugget` directly.",
        call. = FALSE
      )
    }
    p <- unpack(best$par)
    ls_scaled <- p$lengthscale
    nugget_used <- p$nugget
  } else {
    nugget_used <- nugget
  }

  fitted <- emulator_profile_fit(ls_scaled, nugget_used, xs, ys, kernel)
  if (is.null(fitted)) {
    stop(
      "The correlation matrix could not be factorised even after raising the nugget. ",
      "Check for duplicated rows in `x`.",
      call. = FALSE
    )
  }
  if (fitted$nugget > nugget_used) {
    warning(
      "Correlation matrix was not positive definite at nugget ",
      format(nugget_used, digits = 3), "; raised it to ",
      format(fitted$nugget, digits = 3),
      " to complete the Cholesky factorisation.",
      call. = FALSE
    )
  }

  ls_original <- ls_scaled * x_scale
  names(ls_original) <- colnames(x)

  structure(
    list(
      x = x,
      y = y,
      n = n,
      d = d,
      kernel = kernel,
      lengthscale = ls_original,
      nugget = fitted$nugget,
      sigma2 = fitted$sigma2,
      loglik = fitted$loglik,
      estimated = c(lengthscale = free_ls, nugget = free_nug),
      x_center = x_center,
      x_scale = x_scale,
      y_center = y_center,
      y_scale = y_scale,
      x_scaled = xs,
      lengthscale_scaled = ls_scaled,
      chol = fitted$chol,
      alpha = fitted$alpha
    ),
    class = c("gi_emulator", "list")
  )
}

#' Predict from a Gaussian process emulator
#'
#' Returns the posterior mean and the posterior standard deviation of the
#' emulated simulator output at new design points. The standard deviation is the
#' point of the exercise: it is how a user sees where the emulator is
#' interpolating between simulations and where it is guessing.
#'
#' @param object A `gi_emulator` from [fit_emulator()].
#' @param newdata Design points to predict at, on the same scale as the `x`
#'   given to [fit_emulator()]: a matrix with one row per point, a data frame,
#'   a plain vector for a one input emulator, or a plain vector of length `d`
#'   for a single point of a `d` input emulator. Columns are matched by name
#'   when `newdata` names every training input, and by position otherwise.
#' @param include_nugget If `TRUE`, the returned standard deviation is that of a
#'   new noisy simulator run at the point (latent uncertainty plus the nugget).
#'   The default `FALSE` returns the uncertainty about the underlying smooth
#'   response surface.
#' @param ... Unused, present for compatibility with [stats::predict()].
#'
#' @return A data frame with one row per row of `newdata` and columns `mean`,
#'   `sd`, `lower` and `upper`, the last two being a 95 percent Gaussian
#'   predictive interval.
#'
#' @note The interval describes the emulator's uncertainty about the simulator,
#'   conditional on the fitted hyperparameters. It does not carry the
#'   uncertainty in those hyperparameters, and it says nothing about whether the
#'   simulator itself is a fair description of the trial.
#'
#' @examples
#' x <- matrix(seq(0, 1, length.out = 12), ncol = 1)
#' fit <- fit_emulator(x, sin(4 * x[, 1]))
#' predict(fit, newdata = c(0.05, 0.5, 0.95))
#' @export
predict.gi_emulator <- function(object, newdata, include_nugget = FALSE, ...) {
  stopifnot(inherits(object, "gi_emulator"))
  if (missing(newdata)) {
    stop("`newdata` is required; give the design points to predict at.", call. = FALSE)
  }
  if (!is.logical(include_nugget) || length(include_nugget) != 1L || is.na(include_nugget)) {
    stop("`include_nugget` must be TRUE or FALSE.", call. = FALSE)
  }
  nd <- emulator_as_matrix(newdata, "newdata", d = object$d)
  input_names <- names(object$lengthscale)
  if (!is.null(colnames(nd)) && all(input_names %in% colnames(nd))) {
    nd <- nd[, input_names, drop = FALSE]
  }
  if (ncol(nd) != object$d) {
    stop(
      "`newdata` has ", ncol(nd), " columns but the emulator was fitted on ",
      object$d, " inputs.",
      call. = FALSE
    )
  }
  zs <- sweep(sweep(nd, 2L, object$x_center, "-"), 2L, object$x_scale, "/")
  r <- emulator_correlation(zs, object$x_scaled, object$lengthscale_scaled, object$kernel)
  mu <- as.vector(r %*% object$alpha)
  v <- backsolve(object$chol, t(r), transpose = TRUE)
  prior_var <- 1 + if (isTRUE(include_nugget)) object$nugget else 0
  var_s <- object$sigma2 * (prior_var - colSums(v^2))
  var_s[var_s < 0] <- 0
  mean_out <- object$y_center + object$y_scale * mu
  sd_out <- object$y_scale * sqrt(var_s)
  data.frame(
    mean = mean_out,
    sd = sd_out,
    lower = mean_out - 1.959964 * sd_out,
    upper = mean_out + 1.959964 * sd_out
  )
}

#' Leave one out diagnostics for an emulator
#'
#' Leave one out cross validation using the exact Gaussian process identity
#' (Rasmussen and Williams 2006, equation 5.12), so no refitting loop is needed.
#' Two things matter: the root mean squared error, which says whether the
#' emulator interpolates well, and the proportion of held out points falling
#' inside their own 95 percent predictive interval, which says whether its
#' uncertainty is honest. A well calibrated emulator has coverage near 0.95. Far
#' below that and the emulator is overconfident and must not be used to make
#' decisions; far above and it is uninformative.
#'
#' @param fit A `gi_emulator` from [fit_emulator()].
#'
#' @return A list with elements `n`, `rmse`, `rmse_relative` (RMSE divided by
#'   the range of the training output), `mae`, `coverage_95`, `y_range` and
#'   `predictions`, a data frame of the held out `y`, the leave one out `mean`
#'   and `sd`, the interval `lower` and `upper`, and the standardised residual
#'   `z`.
#'
#' @note Hyperparameters are estimated once on the full design and are not re
#'   estimated within each fold, which is the usual practice for emulators and
#'   makes these diagnostics mildly optimistic.
#'
#' @examples
#' x <- cbind(rep(seq(0, 1, length.out = 6), each = 6),
#'            rep(seq(0, 1, length.out = 6), times = 6))
#' fit <- fit_emulator(x, sin(3 * x[, 1]) + x[, 2]^2)
#' loo <- emulator_loo(fit)
#' c(rmse = loo$rmse, coverage = loo$coverage_95)
#' @export
emulator_loo <- function(fit) {
  if (!inherits(fit, "gi_emulator")) {
    stop("`fit` must be a gi_emulator, as returned by fit_emulator().", call. = FALSE)
  }
  a_inv <- chol2inv(fit$chol)
  ys <- (fit$y - fit$y_center) / fit$y_scale
  diag_inv <- diag(a_inv)
  mu_s <- ys - fit$alpha / diag_inv
  var_s <- fit$sigma2 / diag_inv
  var_s[var_s < 0] <- 0

  mean_out <- fit$y_center + fit$y_scale * mu_s
  sd_out <- fit$y_scale * sqrt(var_s)
  lower <- mean_out - 1.959964 * sd_out
  upper <- mean_out + 1.959964 * sd_out
  resid <- fit$y - mean_out
  y_range <- diff(range(fit$y))

  list(
    n = fit$n,
    rmse = sqrt(mean(resid^2)),
    rmse_relative = sqrt(mean(resid^2)) / y_range,
    mae = mean(abs(resid)),
    coverage_95 = mean(fit$y >= lower & fit$y <= upper),
    y_range = y_range,
    predictions = data.frame(
      y = fit$y, mean = mean_out, sd = sd_out,
      lower = lower, upper = upper,
      z = ifelse(sd_out > 0, resid / sd_out, 0)
    )
  )
}

#' Suggest the next design point to simulate
#'
#' Active learning by uncertainty sampling: a large random candidate set is
#' drawn inside `bounds`, the emulator's posterior standard deviation is
#' evaluated at every candidate, and the candidate where the emulator is least
#' certain is returned. Spending the next expensive simulation there reduces
#' emulator uncertainty faster than filling in space the emulator already
#' understands.
#'
#' This criterion targets a globally accurate emulator. It is deliberately not
#' an optimiser: it does not chase the best design, it chases the least known
#' part of the design space.
#'
#' @param fit A `gi_emulator` from [fit_emulator()].
#' @param bounds Search region on the scale of the original inputs: a numeric
#'   matrix or data frame with one row per input and two columns holding the
#'   lower and upper limit, or a length two vector for a one input emulator.
#' @param n_candidates Number of random candidate points to screen. Larger is
#'   more thorough and costs nothing but a matrix multiply.
#' @param seed Integer seed for the candidate draw, recorded in the return value
#'   so the suggestion is reproducible. The ambient random number state is
#'   restored on exit.
#'
#' @return A list with `point` (named numeric vector, the suggested design
#'   point), `sd` (posterior standard deviation there), `mean` (posterior mean
#'   there), `n_candidates`, `seed` and `bounds`.
#'
#' @examples
#' x <- matrix(c(seq(0, 0.3, length.out = 6), seq(0.7, 1, length.out = 6)), ncol = 1)
#' fit <- fit_emulator(x, sin(4 * x[, 1]))
#' suggest_next_point(fit, bounds = c(0, 1))$point
#' @export
suggest_next_point <- function(fit, bounds, n_candidates = 2000, seed = 1) {
  if (!inherits(fit, "gi_emulator")) {
    stop("`fit` must be a gi_emulator, as returned by fit_emulator().", call. = FALSE)
  }
  if (is.data.frame(bounds)) bounds <- data.matrix(bounds)
  if (is.numeric(bounds) && is.null(dim(bounds)) && length(bounds) == 2L) {
    bounds <- matrix(bounds, nrow = 1L)
  }
  if (!is.matrix(bounds) || !is.numeric(bounds) || ncol(bounds) != 2L) {
    stop(
      "`bounds` must be a numeric matrix with one row per input and two columns ",
      "(lower, upper).",
      call. = FALSE
    )
  }
  if (nrow(bounds) != fit$d) {
    stop(
      "`bounds` describes ", nrow(bounds), " input(s) but the emulator has ",
      fit$d, ". Give one row of limits per input.",
      call. = FALSE
    )
  }
  if (!all(is.finite(bounds))) stop("`bounds` must be finite.", call. = FALSE)
  if (any(bounds[, 2] <= bounds[, 1])) {
    bad <- which(bounds[, 2] <= bounds[, 1])
    stop(
      "`bounds` row(s) ", paste(bad, collapse = ", "),
      " have an upper limit that is not above the lower limit.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_candidates) || length(n_candidates) != 1L ||
    !is.finite(n_candidates) || n_candidates < 1) {
    stop("`n_candidates` must be a single positive number.", call. = FALSE)
  }
  n_candidates <- as.integer(n_candidates)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be a single finite number.", call. = FALSE)
  }

  cand <- emulator_with_seed(seed, {
    m <- matrix(stats::runif(n_candidates * fit$d), nrow = n_candidates, ncol = fit$d)
    sweep(sweep(m, 2L, bounds[, 2] - bounds[, 1], "*"), 2L, bounds[, 1], "+")
  })
  colnames(cand) <- names(fit$lengthscale)

  pr <- predict(fit, newdata = cand)
  best <- which.max(pr$sd)
  point <- cand[best, ]
  names(point) <- names(fit$lengthscale)

  list(
    point = point,
    sd = pr$sd[best],
    mean = pr$mean[best],
    n_candidates = n_candidates,
    seed = seed,
    bounds = bounds
  )
}

#' @export
print.gi_emulator <- function(x, ...) {
  cat("<gi_emulator> ", x$kernel, " kernel, ", x$d, " input(s), ",
    x$n, " design points\n",
    sep = ""
  )
  cat("  lengthscale: ",
    paste(sprintf("%s %.4g", names(x$lengthscale), x$lengthscale), collapse = ", "),
    if (isTRUE(x$estimated[["lengthscale"]])) "  (estimated)" else "  (supplied)",
    "\n",
    sep = ""
  )
  cat("  nugget: ", format(x$nugget, digits = 3),
    if (isTRUE(x$estimated[["nugget"]])) "  (estimated)" else "  (supplied)",
    "\n",
    sep = ""
  )
  cat("  log marginal likelihood: ", format(x$loglik, digits = 6), "\n", sep = "")
  cat("  emulator output is interpolation, not simulation\n")
  invisible(x)
}

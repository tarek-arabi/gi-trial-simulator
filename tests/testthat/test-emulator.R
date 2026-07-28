smooth_1d <- function(n = 12) {
  x <- matrix(seq(0, 1, length.out = n), ncol = 1)
  list(x = x, y = sin(4 * x[, 1]))
}

smooth_2d <- function(n_side = 7) {
  g <- seq(0, 1, length.out = n_side)
  x <- as.matrix(expand.grid(control_rate = g, effect = g))
  list(x = x, y = sin(3 * x[, 1]) + 0.5 * x[, 2]^2 - 0.4 * x[, 1] * x[, 2])
}

test_that("fit_emulator returns a documented gi_emulator", {
  d <- smooth_1d()
  fit <- fit_emulator(d$x, d$y)

  expect_s3_class(fit, "gi_emulator")
  expect_true(all(
    c("x", "y", "n", "d", "kernel", "lengthscale", "nugget", "sigma2", "loglik") %in%
      names(fit)
  ))
  expect_identical(fit$n, 12L)
  expect_identical(fit$d, 1L)
  expect_identical(fit$kernel, "matern52")
  expect_gt(fit$lengthscale[[1]], 0)
  expect_true(fit$estimated[["lengthscale"]])
  expect_true(fit$estimated[["nugget"]])
  expect_output(print(fit), "gi_emulator")
})

test_that("the emulator interpolates its training data when the nugget is tiny", {
  d <- smooth_1d()
  fit <- fit_emulator(d$x, d$y, nugget = 1e-10)
  p <- predict(fit, newdata = d$x)

  y_range <- diff(range(d$y))
  expect_lt(max(abs(p$mean - d$y)) / y_range, 1e-6)
  expect_lt(max(p$sd) / y_range, 1e-4)
})

test_that("posterior standard deviation grows away from the design points", {
  x <- matrix(c(seq(0, 0.3, length.out = 7), seq(0.7, 1, length.out = 7)), ncol = 1)
  fit <- fit_emulator(x, sin(4 * x[, 1]))

  at_design <- predict(fit, newdata = matrix(0.15, ncol = 1))
  in_gap <- predict(fit, newdata = matrix(0.5, ncol = 1))

  expect_gt(in_gap$sd, at_design$sd)
  expect_true(all(in_gap$lower < in_gap$mean, in_gap$upper > in_gap$mean))
})

test_that("both kernels fit and predict, and include_nugget widens the interval", {
  d <- smooth_2d(6)
  for (k in c("matern52", "squared_exponential")) {
    fit <- fit_emulator(d$x, d$y, kernel = k, nugget = 1e-4)
    p <- predict(fit, newdata = d$x[1:5, , drop = FALSE])
    expect_identical(nrow(p), 5L)
    expect_identical(names(p), c("mean", "sd", "lower", "upper"))
    expect_true(all(is.finite(p$mean)), info = k)

    noisy <- predict(fit, newdata = d$x[1:5, , drop = FALSE], include_nugget = TRUE)
    expect_true(all(noisy$sd > p$sd), info = k)
    # Ordering alone passes for any wrong multiple of the nugget, so the size of
    # the widening is pinned as well: a new noisy run at a point differs from the
    # latent surface there by exactly the nugget variance, sigma2 * nugget on the
    # standardised output scale, and nothing else.
    expect_equal(
      noisy$sd^2 - p$sd^2,
      rep(fit$sigma2 * fit$nugget * fit$y_scale^2, nrow(p)),
      tolerance = 1e-9, info = k
    )
  }
})

test_that("the predictive standard deviation equals the prior one far from the design", {
  # Every correlation with the design is numerically zero out here, so the
  # posterior is the prior and its standard deviation is known in closed form:
  # y_scale * sqrt(sigma2) for the latent surface and y_scale * sqrt(sigma2 *
  # (1 + nugget)) for a new noisy run. This pins the magnitude of the returned
  # standard deviation, not merely its ordering or its sign.
  d <- smooth_2d(6)
  nug <- 0.05
  fit <- fit_emulator(d$x, d$y, nugget = nug, lengthscale = 0.4)
  far <- matrix(c(20, 20), ncol = 2)

  prior_sd <- fit$y_scale * sqrt(fit$sigma2)
  expect_equal(predict(fit, far)$sd, prior_sd, tolerance = 1e-9)
  expect_equal(
    predict(fit, far, include_nugget = TRUE)$sd,
    prior_sd * sqrt(1 + nug),
    tolerance = 1e-9
  )
  expect_equal(predict(fit, far)$mean, fit$y_center, tolerance = 1e-9)
})

test_that("the 95 percent predictive interval covers 95 percent of held out draws", {
  # Data generated from the model the emulator assumes: a Matern 5/2 process with
  # a known lengthscale and a known nugget, drawn jointly at 45 training and 60
  # held out points, forty times. The emulator is given the true hyperparameters,
  # so what is being tested is the predictive standard deviation itself. Coverage
  # of the nominal interval must come out near 0.95; a standard deviation 15
  # percent too small drops it to about 0.90 and one 15 percent too large lifts
  # it to about 0.975, so the band below is a real constraint on the magnitude.
  matern52_cor <- function(a, b, lengthscale) {
    d2 <- outer(a, b, function(u, v) ((u - v) / lengthscale)^2)
    r <- sqrt(d2)
    (1 + sqrt(5) * r + (5 / 3) * d2) * exp(-sqrt(5) * r)
  }
  set.seed(20)
  true_ls <- 0.15
  nug <- 1e-3
  train <- seq(0, 1, length.out = 45)
  held_out <- seq(0.008, 0.992, length.out = 60)
  all_x <- c(train, held_out)
  root <- t(chol(matern52_cor(all_x, all_x, true_ls) + diag(nug, length(all_x))))

  inside <- 0L
  z <- numeric(0)
  for (rep in seq_len(40)) {
    y <- as.vector(root %*% stats::rnorm(length(all_x)))
    fit <- fit_emulator(matrix(train, ncol = 1), y[seq_along(train)],
      lengthscale = true_ls, nugget = nug
    )
    p <- predict(fit, newdata = matrix(held_out, ncol = 1), include_nugget = TRUE)
    y_out <- y[length(train) + seq_along(held_out)]
    inside <- inside + sum(y_out >= p$lower & y_out <= p$upper)
    z <- c(z, (y_out - p$mean) / p$sd)
  }
  coverage <- inside / (40 * length(held_out))

  expect_gt(coverage, 0.92)
  expect_lt(coverage, 0.97)
  expect_lt(abs(stats::sd(z) - 1), 0.1)
})

test_that("newdata columns are matched by name when they are all named", {
  d <- smooth_2d(6)
  fit <- fit_emulator(d$x, d$y, nugget = 1e-8)
  point <- data.frame(control_rate = 0.4, effect = 0.7)
  reversed <- data.frame(effect = 0.7, control_rate = 0.4)

  expect_equal(predict(fit, point)$mean, predict(fit, reversed)$mean)
})

test_that("leave one out error is small and coverage is near nominal on a smooth surface", {
  d <- smooth_2d(7)
  fit <- fit_emulator(d$x, d$y)
  loo <- emulator_loo(fit)

  expect_identical(loo$n, 49L)
  expect_identical(nrow(loo$predictions), 49L)
  expect_lt(loo$rmse_relative, 0.02)
  # An interpolating emulator of a deterministic function is conservative
  # rather than overconfident, so only the lower side of nominal is tested.
  expect_gte(loo$coverage_95, 0.9)
})

test_that("leave one out coverage is close to nominal on noisy simulation output", {
  # With Monte Carlo noise present the estimated nugget carries it, and coverage
  # can be checked on both sides. n = 81 gives a binomial standard error of
  # about 0.024 on the coverage, so the band below is roughly plus or minus 2 SE.
  set.seed(42)
  g <- seq(0, 1, length.out = 9)
  x <- as.matrix(expand.grid(a = g, b = g))
  y <- sin(3 * x[, 1]) + 0.5 * x[, 2]^2 + stats::rnorm(nrow(x), sd = 0.05)
  fit <- fit_emulator(x, y)
  loo <- emulator_loo(fit)

  expect_gt(loo$coverage_95, 0.88)
  expect_lt(loo$coverage_95, 1.0)
  # The seed is fixed, so this is a deterministic quantity and does not need a
  # wide band. The tolerance is about two sampling standard errors of the
  # standard deviation of 81 residuals, 1 / sqrt(2 * 80) = 0.079, which leaves
  # room for platform level differences in the optimiser while still failing if
  # the leave one out standard deviation is off by more than about 12 percent.
  expect_lt(abs(stats::sd(loo$predictions$z) - 1), 0.15)
})

test_that("leave one out predictions never use the held out point", {
  # The exact identity is checked against the thing it is an identity for: refit
  # the emulator on the retained points, with the hyperparameters held at their
  # full design values, and predict the held out point. A fold that centres or
  # scales the output using the point it is predicting does not match this.
  d <- smooth_1d(12)
  # One point well off the surface makes any leakage of it into its own fold
  # large enough to see rather than lost in rounding.
  d$y[4] <- d$y[4] + 3
  fit <- fit_emulator(d$x, d$y, lengthscale = 0.3, nugget = 1e-3)
  loo <- emulator_loo(fit)

  brute <- vapply(seq_len(fit$n), function(i) {
    refit <- fit_emulator(
      d$x[-i, , drop = FALSE], d$y[-i],
      kernel = fit$kernel, lengthscale = fit$lengthscale, nugget = fit$nugget
    )
    p <- predict(refit, newdata = d$x[i, , drop = FALSE], include_nugget = TRUE)
    c(mean = p$mean, sd = p$sd)
  }, numeric(2))

  expect_equal(loo$predictions$mean, unname(brute["mean", ]), tolerance = 1e-8)
  expect_equal(loo$predictions$sd, unname(brute["sd", ]), tolerance = 1e-8)
})

test_that("the estimated nugget recovers a known noise level", {
  set.seed(9)
  x <- matrix(seq(0, 1, length.out = 60), ncol = 1)
  noise_sd <- 0.08
  y <- sin(3 * x[, 1]) + stats::rnorm(60, sd = noise_sd)
  fit <- fit_emulator(x, y)

  implied_noise_sd <- fit$y_scale * sqrt(fit$sigma2 * fit$nugget)
  expect_gt(implied_noise_sd, noise_sd / 2)
  expect_lt(implied_noise_sd, noise_sd * 2)
})

test_that("a known lengthscale is recovered to within a factor of two", {
  # Hyperparameter estimation for a GP is a hard, weakly identified problem even
  # with data generated from the model, so a factor of two is the honest
  # tolerance here, not a placeholder for a tighter one.
  set.seed(7)
  n <- 80
  true_ls <- 0.25
  x <- matrix(seq(0, 1, length.out = n), ncol = 1)
  d2 <- outer(x[, 1], x[, 1], function(a, b) ((a - b) / true_ls)^2)
  r <- sqrt(d2)
  R <- (1 + sqrt(5) * r + (5 / 3) * d2) * exp(-sqrt(5) * r)
  y <- as.vector(t(chol(R + diag(1e-8, n))) %*% stats::rnorm(n))

  fit <- fit_emulator(x, y, kernel = "matern52")

  expect_gt(fit$lengthscale[[1]], true_ls / 2)
  expect_lt(fit$lengthscale[[1]], true_ls * 2)
})

test_that("a supplied lengthscale is used as given and reported on the input scale", {
  x <- matrix(seq(0, 200, length.out = 20), ncol = 1)
  fit <- fit_emulator(x, sin(x[, 1] / 30), lengthscale = 40, nugget = 1e-8)

  expect_equal(unname(fit$lengthscale), 40)
  expect_false(fit$estimated[["lengthscale"]])
  expect_false(fit$estimated[["nugget"]])
  expect_identical(fit$nugget, 1e-8)
})

test_that("hyperparameter estimation says when it stopped on the search boundary", {
  # Binomial counts whose signal across the design, a 0.02 change in the event
  # rate, sits far below the noise, a binomial standard deviation of about 0.07
  # per point. The likelihood then pushes the lengthscale down and the nugget up
  # until both stop on the edge of the search box, and the fit that comes back
  # interpolates nothing. That used to be reported as an ordinary fit.
  set.seed(4)
  x <- matrix(seq(0, 1, length.out = 80), ncol = 1)
  y <- stats::rbinom(80, 50, 0.5 + 0.02 * x[, 1]) / 50

  said <- character(0)
  withCallingHandlers(
    fit <- fit_emulator(x, y),
    warning = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  said <- paste(said, collapse = " || ")

  # The box really did bind, on both hyperparameters at once.
  expect_equal(unname(fit$lengthscale_scaled), 1e-2)
  expect_equal(fit$nugget, 1)

  expect_match(said, "boundary of its search box")
  expect_match(said, "lengthscale for input 'x1' at its lower limit")
  expect_match(said, "nugget at its upper limit")
})

test_that("a design too small for its hyperparameters is refused", {
  # Three points cannot support a lengthscale, a nugget and a signal variance.
  # This used to return a fit whose posterior mean was flat at the mean of `y`
  # between the design points.
  expect_error(fit_emulator(matrix(c(0, 0.5, 1), ncol = 1), c(0, 1, 0.2)), "too few")
  x8 <- matrix(seq(0, 1, length.out = 8), ncol = 1)
  expect_error(fit_emulator(x8, sin(4 * x8[, 1])), "too few")

  x9 <- matrix(seq(0, 1, length.out = 9), ncol = 1)
  expect_silent(fit <- fit_emulator(x9, sin(4 * x9[, 1])))
  expect_identical(fit$n, 9L)
})

test_that("a lengthscale that leaves the design uncorrelated is reported as degenerate", {
  # A supplied lengthscale bypasses the search box entirely, so the degeneracy
  # has to be caught from the fit itself. With a lengthscale far below the
  # spacing of the design every point is independent of every other, the
  # posterior mean is flat at the mean of `y`, and the emulator has learned
  # nothing while looking like a successful fit.
  d <- smooth_1d(12)
  expect_warning(
    fit <- fit_emulator(d$x, d$y, lengthscale = 0.001, nugget = 1e-6),
    "essentially uncorrelated"
  )
  between <- predict(fit, newdata = matrix(0.5 / 11, ncol = 1))
  expect_equal(between$mean, mean(d$y), tolerance = 1e-6)
})

test_that("a fractional seed is refused rather than silently truncated", {
  # set.seed() truncates towards zero, so seed = 1.7 used to be recorded in the
  # result while the candidates actually came from seed = 1, giving a recorded
  # seed that does not reproduce its own result.
  d <- smooth_2d(6)
  fit <- fit_emulator(d$x, d$y, nugget = 1e-8)
  bounds <- matrix(c(0, 0, 1, 1), ncol = 2)

  expect_error(
    suggest_next_point(fit, bounds, n_candidates = 200, seed = 1.7),
    "whole number"
  )
  expect_error(
    suggest_next_point(fit, bounds, n_candidates = 200, seed = -0.5),
    "whole number"
  )
  s <- suggest_next_point(fit, bounds, n_candidates = 200, seed = -3)
  expect_identical(s$seed, -3)
  expect_identical(
    suggest_next_point(fit, bounds, n_candidates = 200, seed = s$seed)$point,
    s$point
  )
})

test_that("suggest_next_point stays inside the bounds and fills an obvious gap", {
  x <- matrix(c(seq(0, 0.3, length.out = 7), seq(0.7, 1, length.out = 7)), ncol = 1)
  fit <- fit_emulator(x, sin(4 * x[, 1]))
  s <- suggest_next_point(fit, bounds = c(0, 1), n_candidates = 4000)

  expect_gte(s$point[[1]], 0)
  expect_lte(s$point[[1]], 1)
  expect_gt(s$point[[1]], 0.35)
  expect_lt(s$point[[1]], 0.65)
  expect_gt(s$sd, 0)
  expect_identical(s$seed, 1)
})

test_that("suggest_next_point is reproducible and leaves the RNG state alone", {
  d <- smooth_2d(6)
  fit <- fit_emulator(d$x, d$y, nugget = 1e-8)
  bounds <- matrix(c(0, 0, 1, 1), ncol = 2)

  set.seed(99)
  interrupted <- c(stats::runif(1), NA)
  a <- suggest_next_point(fit, bounds, n_candidates = 500, seed = 4)
  interrupted[2] <- stats::runif(1)

  set.seed(99)
  undisturbed <- stats::runif(2)
  expect_identical(interrupted, undisturbed)

  b <- suggest_next_point(fit, bounds, n_candidates = 500, seed = 4)
  expect_identical(a$point, b$point)
  expect_identical(names(a$point), c("control_rate", "effect"))
  expect_true(all(a$point >= bounds[, 1] & a$point <= bounds[, 2]))
})

test_that("a singular correlation matrix raises the nugget and says so", {
  # Duplicated design points make the correlation matrix singular, which is the
  # realistic way a simulation grid breaks a Cholesky factorisation.
  x <- matrix(c(seq(0, 1, length.out = 10), 0.5, 0.5), ncol = 1)
  y <- c(sin(4 * x[1:10, 1]), sin(2), sin(2))

  expect_warning(
    fit <- fit_emulator(x, y, nugget = 0, lengthscale = 0.5),
    "not positive definite"
  )
  expect_gt(fit$nugget, 0)
  expect_equal(predict(fit, newdata = matrix(0.5, ncol = 1))$mean, sin(2))
})

test_that("fit_emulator rejects malformed input with a message naming the argument", {
  d <- smooth_1d()

  expect_error(fit_emulator(d$x, d$y[-1]), "`y` has length")
  expect_error(fit_emulator(d$x, d$y, kernel = "linear"), "`kernel` must be one of")
  expect_error(fit_emulator(d$x, rep(1, 12)), "`y` is constant")
  expect_error(
    fit_emulator(cbind(d$x, 1), c(d$y)),
    "constant"
  )
  expect_error(fit_emulator(d$x, d$y, nugget = -1), "`nugget` must be")
  expect_error(fit_emulator(d$x, d$y, lengthscale = c(1, 2)), "`lengthscale` must have length")
  expect_error(fit_emulator(matrix(1:4, ncol = 2), c(1, 2)), "too few")
  # A single training point (n = 1, d = 1) is the plainest case of "too few".
  expect_error(fit_emulator(matrix(0.5, ncol = 1), 5), "too few")
  # The message says how many points the design would need.
  expect_error(fit_emulator(matrix(1:4, ncol = 2), c(1, 2)), "at least 12 points")
  expect_error(fit_emulator(matrix(letters[1:12], ncol = 1), d$y), "numeric matrix")
  expect_error(fit_emulator(d$x, letters[1:12]), "`y` must be a plain numeric vector")
})

test_that("predict and suggest_next_point reject malformed input", {
  d <- smooth_2d(6)
  fit <- fit_emulator(d$x, d$y, nugget = 1e-8)

  expect_error(predict(fit, newdata = matrix(0.5, ncol = 1)), "columns but the emulator")
  expect_error(predict(fit, newdata = c(0.1, 0.2, 0.3)), "plain vector of length 3")
  expect_error(suggest_next_point(fit, bounds = c(0, 1)), "describes 1 input")
  expect_error(
    suggest_next_point(fit, bounds = matrix(c(1, 1, 0, 0), ncol = 2)),
    "upper limit that is not above"
  )
  expect_error(emulator_loo(fit$x), "must be a gi_emulator")
})

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
  }
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
  expect_lt(abs(stats::sd(loo$predictions$z) - 1), 0.35)
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

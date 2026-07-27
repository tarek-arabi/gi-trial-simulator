test_that("covariate_spec validates the distribution and its parameters", {
  expect_error(covariate_spec("x", "weibull", shape = 1), "must be one of")
  expect_error(covariate_spec(c("a", "b"), "normal"), "single non-empty character")
  expect_error(covariate_spec("x", "normal", sd = 0), "`sd` must be a single positive")
  expect_error(covariate_spec("x", "normal", mu = 3), "unknown parameter")
  expect_error(covariate_spec("x", "binary", prob = 1), "strictly between 0 and 1")
  expect_error(covariate_spec("x", "gamma", shape = 2, rate = 1, scale = 1), "only one of")
  expect_error(covariate_spec("x", "ordinal", probs = c(0.5, 0.4)), "must sum to 1")
  expect_error(covariate_spec("x", "ordinal", probs = c(0.5, 0.5, 0)), "strictly positive")
  expect_error(
    covariate_spec("x", "ordinal", probs = c(0.5, 0.5), values = c(3, 1)),
    "strictly increasing"
  )
})

test_that("covariate_spec accepts every supported marginal", {
  cvs <- list(
    covariate_spec("a", "normal", mean = 1, sd = 2),
    covariate_spec("b", "lognormal", meanlog = 0.5, sdlog = 0.3),
    covariate_spec("c", "gamma", shape = 2, scale = 10),
    covariate_spec("d", "binary", prob = 0.3),
    covariate_spec("e", "ordinal", probs = c(0.5, 0.3, 0.2))
  )
  expect_true(all(vapply(cvs, inherits, logical(1), what = "gi_covariate")))
  expect_equal(cvs[[3]]$pars$rate, 0.1)
})

test_that("cohort_spec rejects malformed correlation matrices by name", {
  cvs <- list(
    covariate_spec("a", "normal"),
    covariate_spec("b", "normal")
  )
  expect_error(cohort_spec(cvs, matrix(c(1, 0.3, 0.4, 1), nrow = 2)), "not symmetric")
  expect_error(cohort_spec(cvs, matrix(c(1, 0.3, 0.3, 1), nrow = 2)[, 1, drop = FALSE]), "square")
  expect_error(cohort_spec(cvs, diag(3)), "2 by 2")
  expect_error(cohort_spec(cvs, matrix(c(2, 0, 0, 1), nrow = 2)), "1 on its diagonal")
  three <- list(
    covariate_spec("a", "normal"),
    covariate_spec("b", "normal"),
    covariate_spec("c", "normal")
  )
  singular <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 3)
  expect_error(cohort_spec(three, singular), "not positive definite")
  inconsistent <- matrix(c(1, 0.9, -0.9, 0.9, 1, 0.9, -0.9, 0.9, 1), nrow = 3)
  expect_error(cohort_spec(three, inconsistent), "not positive definite")
})

test_that("cohort_spec rejects duplicate names and non-covariates", {
  expect_error(
    cohort_spec(list(covariate_spec("a", "normal"), covariate_spec("a", "normal"))),
    "duplicated name"
  )
  expect_error(cohort_spec(list(covariate_spec("a", "normal"), 1)), "not gi_covariate")
  expect_error(cohort_spec(list()), "non-empty list")
})

test_that("cohort_spec reorders a named correlation matrix to the covariate order", {
  cvs <- list(covariate_spec("a", "normal"), covariate_spec("b", "normal"))
  r <- matrix(c(1, 0.4, 0.4, 1), nrow = 2, dimnames = list(c("b", "a"), c("b", "a")))
  spec <- cohort_spec(cvs, r)
  expect_equal(rownames(spec$correlation), c("a", "b"))
  expect_equal(spec$correlation["a", "b"], 0.4)
  bad <- matrix(c(1, 0.4, 0.4, 1), nrow = 2, dimnames = list(c("a", "q"), c("a", "q")))
  expect_error(cohort_spec(cvs, bad), "dimnames do not match")
})

test_that("simulate_cohort reproduces each marginal distribution", {
  spec <- cohort_spec(list(
    covariate_spec("norm", "normal", mean = 68, sd = 13),
    covariate_spec("lnorm", "lognormal", meanlog = 1.2, sdlog = 0.6),
    covariate_spec("gam", "gamma", shape = 2, rate = 0.05),
    covariate_spec("bin", "binary", prob = 0.28),
    covariate_spec("ord", "ordinal", probs = c(0.5, 0.3, 0.2))
  ))
  x <- simulate_cohort(spec, n = 40000, seed = 12)

  expect_s3_class(x, "data.frame")
  expect_equal(nrow(x), 40000L)
  expect_named(x, c("norm", "lnorm", "gam", "bin", "ord"))

  expect_equal(mean(x$norm), 68, tolerance = 0.02)
  expect_equal(stats::sd(x$norm), 13, tolerance = 0.02)
  expect_equal(mean(x$lnorm), exp(1.2 + 0.6^2 / 2), tolerance = 0.03)
  expect_equal(mean(x$gam), 2 / 0.05, tolerance = 0.03)
  expect_equal(mean(x$bin), 0.28, tolerance = 0.02)
  expect_true(all(x$bin %in% c(0, 1)))
  expect_equal(as.numeric(prop.table(table(x$ord))), c(0.5, 0.3, 0.2), tolerance = 0.02)
})

test_that("the Gaussian copula reproduces target correlations for normal marginals", {
  target <- matrix(
    c(
      1.0, 0.6, -0.3,
      0.6, 1.0, 0.2,
      -0.3, 0.2, 1.0
    ),
    nrow = 3, byrow = TRUE
  )
  spec <- cohort_spec(
    list(
      covariate_spec("a", "normal", mean = 5, sd = 2),
      covariate_spec("b", "normal", mean = -1, sd = 0.5),
      covariate_spec("c", "normal")
    ),
    correlation = target
  )
  # n = 30000 gives a Monte Carlo SE near 0.004 on each correlation, so a
  # tolerance of 0.02 is roughly five Monte Carlo standard errors.
  chk <- cohort_correlation_check(spec, n = 30000, seed = 4)
  expect_equal(nrow(chk), 3L)
  expect_equal(chk$realised, chk$target, tolerance = 0.02)
  expect_true(max(abs(chk$attenuation)) < 0.02)
})

test_that("cohort_correlation_check reports the attenuation for binary marginals", {
  spec <- cohort_spec(
    list(
      covariate_spec("x", "normal"),
      covariate_spec("z", "binary", prob = 0.1)
    ),
    correlation = matrix(c(1, 0.6, 0.6, 1), nrow = 2)
  )
  chk <- cohort_correlation_check(spec, n = 30000, seed = 4)
  expect_equal(chk$target, 0.6)
  # A rare binary marginal cannot attain a product-moment correlation of 0.6
  # against a normal. The realised value must be materially lower, and the
  # function must say so rather than reporting the target back.
  expect_lt(chk$realised, 0.45)
  expect_gt(chk$attenuation, 0.15)
  expect_equal(chk$attenuation, chk$target - chk$realised)
})

test_that("cohort_correlation_check needs at least two covariates", {
  spec <- cohort_spec(list(covariate_spec("a", "normal")))
  expect_error(cohort_correlation_check(spec, n = 100), "at least two covariates")
})

test_that("simulate_cohort is reproducible and seed-controlled", {
  spec <- cohort_spec(
    list(covariate_spec("a", "normal"), covariate_spec("b", "binary", prob = 0.4)),
    correlation = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  )
  expect_identical(simulate_cohort(spec, 200, seed = 9), simulate_cohort(spec, 200, seed = 9))
  expect_false(identical(
    simulate_cohort(spec, 200, seed = 9)$a,
    simulate_cohort(spec, 200, seed = 10)$a
  ))

  set.seed(3)
  first <- simulate_cohort(spec, 50, seed = NULL)
  set.seed(3)
  second <- simulate_cohort(spec, 50, seed = NULL)
  expect_identical(first, second)

  expect_error(simulate_cohort(spec, 0), "positive whole number")
  expect_error(simulate_cohort(spec, 10, seed = "a"), "single number")
  expect_error(simulate_cohort(list(), 10), "gi_cohort_spec")
})

test_that("outcome_model gives the outcome real covariate signal", {
  spec <- cohort_spec(list(
    covariate_spec("x1", "normal"),
    covariate_spec("x2", "normal")
  ))
  cohort <- simulate_cohort(spec, n = 5000, seed = 21)
  out <- outcome_model(cohort, coefs = c(x1 = 1, x2 = -0.7), intercept = -1, seed = 21)

  expect_s3_class(out, "gi_outcome")
  expect_equal(out$n, 5000L)
  expect_true(all(out$y %in% c(0L, 1L)))
  expect_equal(out$event_rate, mean(out$p))
  expect_equal(out$observed_rate, mean(out$y))
  # Monte Carlo SE of the observed rate at n = 5000 is about 0.007, so allow
  # three of them in absolute terms rather than a relative tolerance.
  expect_lt(abs(out$observed_rate - out$event_rate), 0.021)
  expect_equal(out$p, stats::plogis(out$eta))
  # The signal must be real: cases have systematically higher x1 than controls.
  expect_gt(mean(cohort$x1[out$y == 1]), mean(cohort$x1[out$y == 0]) + 0.5)
})

test_that("outcome_model validates its arguments", {
  spec <- cohort_spec(list(covariate_spec("x1", "normal")))
  cohort <- simulate_cohort(spec, n = 50, seed = 1)
  expect_error(outcome_model(cohort, c(nope = 1), 0), "not found in `cohort`")
  expect_error(outcome_model(cohort, c(1), 0), "named numeric vector")
  expect_error(outcome_model(cohort, c(x1 = 1), NA_real_), "single finite number")
  expect_error(outcome_model(cohort, c(x1 = 1), 0, link = "identity"), "`link` must be one of")
  expect_error(outcome_model(as.matrix(cohort), c(x1 = 1), 0), "must be a data frame")
})

test_that("calibrate_intercept hits the target marginal event rate", {
  spec <- cohort_spec(
    list(
      covariate_spec("x1", "normal"),
      covariate_spec("x2", "binary", prob = 0.3)
    ),
    correlation = matrix(c(1, 0.25, 0.25, 1), nrow = 2)
  )
  cohort <- simulate_cohort(spec, n = 30000, seed = 33)
  coefs <- c(x1 = 1.1, x2 = -0.8)

  for (link in c("logit", "probit", "cloglog")) {
    for (target in c(0.0395, 0.0658, 0.3, 0.75)) {
      a <- calibrate_intercept(cohort, coefs, target, link = link)
      achieved <- outcome_model(cohort, coefs, a, link = link, seed = 1)$event_rate
      expect_equal(achieved, target, tolerance = 1e-6)
    }
  }
})

test_that("calibrate_intercept validates its arguments", {
  spec <- cohort_spec(list(covariate_spec("x1", "normal")))
  cohort <- simulate_cohort(spec, n = 100, seed = 1)
  expect_error(calibrate_intercept(cohort, c(x1 = 1), 0), "strictly between 0 and 1")
  expect_error(calibrate_intercept(cohort, c(x1 = 1), 1), "strictly between 0 and 1")
  expect_error(calibrate_intercept(cohort, c(x1 = 1), 0.2, link = "log"), "`link` must be one of")
})

test_that("print methods run without error", {
  spec <- cohort_spec(
    list(
      covariate_spec("a", "normal"),
      covariate_spec("b", "ordinal", probs = c(0.6, 0.4))
    ),
    correlation = matrix(c(1, 0.2, 0.2, 1), nrow = 2)
  )
  expect_output(print(spec), "gi_cohort_spec")
  expect_output(print(spec$covariates$a), "gi_covariate")
  cohort <- simulate_cohort(spec, 100, seed = 1)
  expect_output(print(outcome_model(cohort, c(a = 1), -1, seed = 1)), "gi_outcome")
})

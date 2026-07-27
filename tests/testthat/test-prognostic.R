# A three-covariate structure with a strongly prognostic signal. Reused across
# the expensive tests below so procova_gain is only run a handful of times.
prognostic_spec <- function() {
  cohort_spec(
    list(
      covariate_spec("x1", "normal"),
      covariate_spec("x2", "normal"),
      covariate_spec("x3", "binary", prob = 0.35)
    ),
    correlation = matrix(
      c(
        1.00, 0.35, 0.20,
        0.35, 1.00, 0.15,
        0.20, 0.15, 1.00
      ),
      nrow = 3, byrow = TRUE
    )
  )
}
strong_coefs <- c(x1 = 1.0, x2 = 0.8, x3 = -0.7)

test_that("the Mann-Whitney AUC matches a brute-force pairwise count", {
  set.seed(1)
  y <- rbinom(60, 1, 0.4)
  s <- rnorm(60) + y
  cases <- s[y == 1]
  controls <- s[y == 0]
  brute <- mean(outer(cases, controls, function(a, b) (a > b) + 0.5 * (a == b)))
  expect_equal(auc_rank(s, y), brute)
  expect_equal(auc_rank(c(1, 2, 3, 4), c(0L, 0L, 1L, 1L)), 1)
  expect_equal(auc_rank(c(4, 3, 2, 1), c(0L, 0L, 1L, 1L)), 0)
  expect_true(is.na(auc_rank(c(1, 2, 3), c(0L, 0L, 0L))))
})

test_that("fit_prognostic recovers real signal with an honest CV AUC", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 3000, seed = 51)
  y <- outcome_model(train, strong_coefs, calibrate_intercept(train, strong_coefs, 0.3),
    seed = 51
  )
  fit <- fit_prognostic(train, y, folds = 5, seed = 51)

  expect_s3_class(fit, "gi_prognostic")
  expect_equal(fit$n_train, 3000L)
  expect_true(fit$converged)
  expect_setequal(fit$predictors, c("x1", "x2", "x3"))
  expect_gt(fit$auc_cv, 0.7)
  # Cross-validation must not flatter the model: the honest estimate cannot
  # exceed the in-sample one by any material margin.
  expect_lte(fit$auc_cv, fit$auc_apparent + 0.005)
  expect_equal(stats::coef(fit$model)[["x1"]], 1.0, tolerance = 0.2)
})

test_that("a pure-noise model has a cross-validated AUC near 0.5", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 3000, seed = 52)
  y <- outcome_model(train, c(x1 = 0, x2 = 0, x3 = 0), stats::qlogis(0.3), seed = 52)
  fit <- fit_prognostic(train, y, folds = 5, seed = 52)
  expect_lt(abs(fit$auc_cv - 0.5), 0.05)
  # In-sample fitting of noise is optimistic; cross-validation removes it.
  expect_lte(fit$auc_cv, fit$auc_apparent)
})

test_that("fit_prognostic honours an explicit one-sided formula", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 800, seed = 53)
  y <- outcome_model(train, strong_coefs, -1, seed = 53)
  fit <- fit_prognostic(train, y, formula = ~ x1 + x3)
  expect_setequal(fit$predictors, c("x1", "x3"))
  expect_error(fit_prognostic(train, y, formula = y ~ x1), "one-sided formula")
})

test_that("fit_prognostic validates its arguments", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 200, seed = 54)
  y <- outcome_model(train, strong_coefs, -1, seed = 54)
  expect_error(fit_prognostic(as.matrix(train), y), "must be a data frame")
  expect_error(fit_prognostic(train, y$y[1:10]), "has length 10")
  expect_error(fit_prognostic(train, rep(0L, 200)), "no variation")
  expect_error(fit_prognostic(train, y, folds = 1), "between 2 and the training n")
  expect_error(fit_prognostic(train, rep(c(0, 2), 100)), "only 0 and 1")
})

test_that("prognostic_score scores new patients and refuses incomplete ones", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 800, seed = 55)
  y <- outcome_model(train, strong_coefs, -1, seed = 55)
  fit <- fit_prognostic(train, y)
  fresh <- simulate_cohort(spec, n = 40, seed = 56)

  s <- prognostic_score(fit, fresh)
  expect_length(s, 40L)
  expect_true(all(is.finite(s)))
  expect_equal(s, unname(stats::predict(fit$model, newdata = fresh, type = "link")))
  # The frozen model has no random component left, so scoring the same data
  # twice must give bit-for-bit identical output, not merely similar output.
  expect_identical(prognostic_score(fit, fresh), s)
  expect_error(prognostic_score(fit, fresh[, c("x1", "x2")]), "missing predictor")
  expect_error(prognostic_score(fit$model, fresh), "must be a gi_prognostic")
})

test_that("the g-computation risk difference reproduces the classic two-proportion result", {
  set.seed(7)
  n <- 400
  arm <- rep(0:1, each = n / 2)
  y <- stats::rbinom(n, 1, 0.3)
  score <- stats::rnorm(n)
  res <- analyse_with_prognostic(y, arm, score)

  p1 <- mean(y[arm == 1])
  p0 <- mean(y[arm == 0])
  expect_equal(res$unadjusted$rd, p1 - p0)
  # With no covariate in the model, the delta-method standard error must reduce
  # to the textbook two-proportion standard error. This is the check that
  # validates the g-computation gradient itself. The agreement is limited by
  # glm's IRLS convergence tolerance rather than by machine precision.
  expect_equal(
    res$unadjusted$rd_se,
    sqrt(p1 * (1 - p1) / (n / 2) + p0 * (1 - p0) / (n / 2)),
    tolerance = 1e-6
  )
})

test_that("analyse_with_prognostic validates its arguments and directions", {
  set.seed(8)
  y <- rbinom(100, 1, 0.4)
  arm <- rep(0:1, each = 50)
  s <- rnorm(100)
  expect_error(analyse_with_prognostic(y, arm, s[1:50]), "same length")
  expect_error(analyse_with_prognostic(y, arm[1:60], s), "same length")
  expect_error(analyse_with_prognostic(y, rep(0L, 100), s), "only one distinct value")
  expect_error(analyse_with_prognostic(y, rep(1L, 100), s), "only one distinct value")
  expect_error(analyse_with_prognostic(y, arm, c(NA, s[-1])), "no missing or infinite")
  expect_error(analyse_with_prognostic(rep(2, 100), arm, s), "only 0 and 1")
  expect_error(analyse_with_prognostic(y, arm, s, direction = "sideways"), "should be one of")

  lo <- analyse_with_prognostic(y, arm, s, direction = "lower_is_better")
  hi <- analyse_with_prognostic(y, arm, s, direction = "higher_is_better")
  expect_equal(lo$adjusted$p_one_sided + hi$adjusted$p_one_sided, 1)
  expect_equal(lo$adjusted$estimate, hi$adjusted$estimate)
  expect_output(print(lo), "gi_adjusted_analysis")
})

test_that("a prognostic score shrinks the risk difference SE and noise does not", {
  spec <- prognostic_spec()
  train <- simulate_cohort(spec, n = 4000, seed = 61)
  a0 <- calibrate_intercept(train, strong_coefs, 0.3)
  fit <- fit_prognostic(train, outcome_model(train, strong_coefs, a0, seed = 61), seed = 61)

  nrep <- 300
  set.seed(62)
  seeds <- sample.int(.Machine$integer.max, nrep)
  ratios <- vapply(seq_len(nrep), function(i) {
    co <- simulate_cohort(spec, n = 600, seed = seeds[i])
    arm <- sample(rep(0:1, each = 300))
    y <- stats::rbinom(600, 1, stats::plogis(
      linear_predictor(co, strong_coefs, a0) - 0.4 * arm
    ))
    real <- analyse_with_prognostic(y, arm, prognostic_score(fit, co))
    noise <- analyse_with_prognostic(y, arm, stats::rnorm(600))
    c(real$se_ratio, noise$se_ratio, real$se_ratio_logor)
  }, numeric(3))

  # A genuinely prognostic score reduces the variance of the marginal risk
  # difference, the estimand both analyses share.
  expect_lt(mean(ratios[1, ]), 0.85)
  # A pure-noise covariate costs one degree of freedom and buys nothing.
  expect_equal(mean(ratios[2, ]), 1, tolerance = 0.02)
  # On the log odds ratio scale the adjusted standard error is LARGER, because
  # the logistic model is non-collapsible and the two analyses are estimating
  # different quantities. This is expected, is documented, and is the reason the
  # risk difference is used for the efficiency comparison above.
  expect_gt(mean(ratios[3, ]), 1)
})

# The expensive run, shared by the tests below. nsim is set at the top of the
# range this module allows, because the type I error test underneath it is the
# one assertion the whole method stands on and it is worth nothing without the
# Monte Carlo resolution to back it. Measured cost of this single call is about
# 13 seconds, which is the entire reason there is no cheaper compromise here.
gain_strong <- procova_gain(
  scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.22),
  prognostic_spec(),
  coefs = strong_coefs,
  n_per_arm = 300, nsim = 2000, train_n = 4000, seed = 101, calibrate_n = 30000
)

test_that("CRITICAL: prognostic adjustment does not inflate type I error", {
  mcse <- gain_strong$type1_adjusted_mcse

  # This test is only meaningful at a Monte Carlo resolution fine enough to see
  # an inflation worth caring about. At nsim = 2000 the three-MCSE band is about
  # alpha +/- 0.0105; at nsim = 400 it widens to alpha +/- 0.023, which would
  # pass even if adjustment nearly doubled the rejection rate. These two
  # assertions exist so that shrinking nsim fails loudly instead of quietly
  # turning the test below into a formality. Do not reduce nsim to speed this up.
  expect_gte(gain_strong$nsim, 2000)
  expect_lt(3 * mcse, 0.6 * gain_strong$alpha)

  expect_gt(gain_strong$prognostic$auc_cv, 0.75)
  expect_lt(abs(gain_strong$type1_adjusted - gain_strong$alpha), 3 * mcse)
  expect_lt(abs(gain_strong$rd_type1_adjusted - gain_strong$alpha), 3 * mcse)
  # The unadjusted analysis is the reference: adjustment must not be worse.
  expect_lt(abs(gain_strong$type1_unadjusted - gain_strong$alpha), 3 * mcse)
  expect_equal(mcse, sqrt(0.025 * 0.975 / 2000), tolerance = 0.15)
})

test_that("procova_gain reports a real efficiency gain on a common estimand", {
  expect_s3_class(gain_strong, "gi_procova")
  expect_gt(gain_strong$power_adjusted, gain_strong$power_unadjusted)
  expect_gt(
    gain_strong$power_adjusted - gain_strong$power_unadjusted,
    3 * gain_strong$power_adjusted_mcse
  )
  expect_lt(gain_strong$se_ratio, 1)
  expect_equal(gain_strong$n_reduction, 1 - gain_strong$se_ratio)
  expect_true(gain_strong$n_reduction > 0.10 && gain_strong$n_reduction < 0.45)
  # Two independent routes to the same efficiency gain must agree: the ratio of
  # squared risk difference standard errors, and the squared ratio of log odds
  # ratio test statistics.
  expect_equal(
    gain_strong$information_ratio, 1 / gain_strong$se_ratio,
    tolerance = 0.03
  )
  expect_gt(gain_strong$se_ratio_logor, 1)
  expect_equal(gain_strong$n_nonconverged, 0)
})

test_that("procova_gain calibrates to the scenario event rates", {
  expect_equal(gain_strong$achieved_control_rate, 0.30, tolerance = 0.01)
  expect_equal(gain_strong$achieved_treatment_rate, 0.22, tolerance = 0.01)
  expect_lt(gain_strong$arm_effect, 0)
  expect_equal(gain_strong$n_total, 600L)
})

test_that("the delta-method standard errors match their Monte Carlo counterparts", {
  # Model-based standard error against the empirical spread of the estimates
  # across replicates. Any error in the g-computation gradient would show up
  # here as a ratio away from 1.
  expect_equal(
    gain_strong$mean_rd_se_unadjusted / gain_strong$emp_rd_se_unadjusted, 1,
    tolerance = 0.06
  )
  expect_equal(
    gain_strong$mean_rd_se_adjusted / gain_strong$emp_rd_se_adjusted, 1,
    tolerance = 0.06
  )
  expect_equal(
    gain_strong$mean_se_unadjusted / gain_strong$emp_se_unadjusted, 1,
    tolerance = 0.06
  )
  # Both estimators must be unbiased for the true marginal risk difference,
  # -0.08 = 0.22 - 0.30. The Monte Carlo error of a mean over nsim replicates
  # is the empirical spread of the per-replicate estimate divided by
  # sqrt(nsim); allow 4 of those, a bit more generous than the 3 MCSE used for
  # the type I error check because this is a two-sided bias check rather than
  # a bound on a single tail.
  mcse_rd_unadj <- gain_strong$emp_rd_se_unadjusted / sqrt(gain_strong$nsim)
  mcse_rd_adj <- gain_strong$emp_rd_se_adjusted / sqrt(gain_strong$nsim)
  expect_lt(abs(gain_strong$rd_unadjusted - (-0.08)), 4 * mcse_rd_unadj)
  expect_lt(abs(gain_strong$rd_adjusted - (-0.08)), 4 * mcse_rd_adj)
})

test_that("a non-prognostic covariate buys nothing but costs nothing", {
  # nsim = 300 again, for the same runtime reason as gain_strong above.
  gain_null <- procova_gain(
    scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.22),
    prognostic_spec(),
    coefs = c(x1 = 0, x2 = 0, x3 = 0),
    n_per_arm = 300, nsim = 300, train_n = 2000, seed = 202, calibrate_n = 10000
  )
  expect_lt(abs(gain_null$prognostic$auc_cv - 0.5), 0.06)
  expect_equal(gain_null$se_ratio, 1, tolerance = 0.02)
  expect_lt(abs(gain_null$n_reduction), 0.02)
  expect_lt(
    abs(gain_null$type1_adjusted - gain_null$alpha),
    3 * gain_null$type1_adjusted_mcse
  )
})

test_that("procova_gain is reproducible and independent of the worker count", {
  args <- list(
    scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.22),
    prognostic_spec(),
    coefs = strong_coefs,
    n_per_arm = 120, nsim = 30, train_n = 600, seed = 7, calibrate_n = 5000
  )
  a <- do.call(procova_gain, args)
  b <- do.call(procova_gain, args)
  expect_identical(a$replicates, b$replicates)

  c2 <- do.call(procova_gain, c(args, list(workers = 2)))
  expect_equal(c2$replicates, a$replicates)
})

test_that("procova_gain validates its arguments", {
  spec <- prognostic_spec()
  sc <- scenario("ercp_acute_cholangitis")
  expect_error(procova_gain(list(), spec, strong_coefs, 100), "must be a gi_scenario")
  expect_error(procova_gain(sc, list(), strong_coefs, 100), "must be a gi_cohort_spec")
  expect_error(procova_gain(sc, spec, c(nope = 1), 100), "Not found: nope")
  expect_error(procova_gain(sc, spec, strong_coefs, 5), "`n_per_arm` must be")
  expect_error(procova_gain(sc, spec, strong_coefs, 100, nsim = 1), "`nsim` must be")
  expect_error(
    procova_gain(sc, spec, strong_coefs, 100, nsim = 5, alpha = 0.9),
    "strictly between 0 and 0.5"
  )
  expect_error(
    procova_gain(sc, spec, strong_coefs, 100, nsim = 5, link = "identity"),
    "`link` must be one of"
  )
})

test_that("print.gi_procova reports the numbers a reviewer would ask for", {
  out <- utils::capture.output(print(gain_strong))
  expect_true(any(grepl("type I error", out)))
  expect_true(any(grepl("variance-reduction factor", out)))
  expect_true(any(grepl("sample size reduction", out)))
  expect_true(any(grepl("cross-check", out)))
  expect_true(any(grepl("CV AUC", out)))
})

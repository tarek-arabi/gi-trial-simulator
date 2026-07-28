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

test_that("fit_prognostic refuses to report a CV AUC it cannot compute", {
  # Reproduction of the defect: 120 patients at a 2 percent event rate leaves 2
  # events. With unstratified folds both could land in one fold, so that fold's
  # model is fitted to an outcome with no variation and its out-of-fold scores
  # are not predictions. The old code pooled them anyway and reported
  # auc_cv = 0.20 against auc_apparent = 0.96, converged = TRUE, silently.
  spec <- cohort_spec(list(
    covariate_spec("x1", "normal"),
    covariate_spec("x2", "normal")
  ))
  train <- simulate_cohort(spec, n = 120, seed = 1)
  coefs <- c(x1 = 1, x2 = 0.8)
  y <- outcome_model(train, coefs, calibrate_intercept(train, coefs, 0.02), seed = 1)
  expect_lt(sum(y$y), 5L)

  expect_warning(
    fit <- fit_prognostic(train, y, folds = 5, seed = 1),
    "cross-validation was skipped"
  )
  expect_true(is.na(fit$auc_cv))
  expect_false(fit$cv_performed)
  expect_true(all(is.na(fit$oof_score)))
  # The apparent AUC is still reported; only the honest one is withheld.
  expect_false(is.na(fit$auc_apparent))
  # The warning has to say what went wrong, not just that something did.
  expect_warning(fit_prognostic(train, y, folds = 5, seed = 1), "event\\(s\\)")
  expect_output(print(fit), "NA \\(too few events to fold\\)")
})

test_that("fit_prognostic stratifies folds so every fold holds both classes", {
  # Rare enough that unstratified folds would routinely leave a fold without an
  # event, but not so rare that cross-validation has to be abandoned.
  spec <- cohort_spec(list(
    covariate_spec("x1", "normal"),
    covariate_spec("x2", "normal")
  ))
  train <- simulate_cohort(spec, n = 400, seed = 71)
  coefs <- c(x1 = 1, x2 = 0.8)
  y <- outcome_model(train, coefs, calibrate_intercept(train, coefs, 0.03), seed = 71)
  yint <- as.integer(y$y)
  expect_equal(sum(yint), 10L)

  empty_fold <- function(id) {
    any(tabulate(id[yint == 1L], 5L) == 0L) || any(tabulate(id[yint == 0L], 5L) == 0L)
  }

  # Fold seed 6 is one of the many at which the unstratified draw the old code
  # used leaves a fold with no events at all. Asserting on `fold_id`, the
  # assignment the fit actually used, is what makes this catch a regression to
  # unstratified sampling; asserting on the helper alone would not.
  set.seed(6L)
  expect_true(empty_fold(sample(rep_len(seq_len(5L), length(yint)))))

  fit <- expect_silent(fit_prognostic(train, y, folds = 5, seed = 6))
  expect_true(fit$cv_performed)
  expect_false(is.na(fit$auc_cv))
  expect_true(all(is.finite(fit$oof_score)))
  expect_false(empty_fold(fit$fold_id))
  expect_equal(sort(unique(fit$fold_id)), 1:5)

  # The guarantee is not "works on this seed", it is "works on every seed", so
  # assert it across 200 of them: stratified assignment must never leave a fold
  # without an event, while the unstratified draw does so on roughly half of
  # them. A run of 200 in which the unstratified draw never failed would itself
  # be evidence that the two schemes had become the same thing.
  strat_bad <- 0L
  naive_bad <- 0L
  for (s in seq_len(200)) {
    set.seed(s)
    strat_bad <- strat_bad + empty_fold(stratified_folds(yint, 5L))
    set.seed(s)
    naive_bad <- naive_bad + empty_fold(sample(rep_len(seq_len(5L), length(yint))))
  }
  expect_equal(strat_bad, 0L)
  expect_gt(naive_bad, 10L)
})

test_that("fit_prognostic reports separation instead of calling it convergence", {
  # A covariate that predicts the outcome perfectly except for one patient.
  # glm walks the coefficient out along a flat ridge, stops when the deviance
  # stops moving, and sets model$converged to TRUE.
  d <- data.frame(z = c(rep(1, 20), rep(0, 60)), w = rep(c(-1, 1), 40))
  yz <- c(rep(1L, 20), rep(c(0L, 1L), c(45, 15)))
  expect_warning(fit <- fit_prognostic(d, yz, folds = 3, seed = 1), "separation")
  expect_true(isTRUE(fit$model$converged))
  expect_true("z" %in% fit$separation)
  expect_false(fit$converged)
  expect_gt(abs(stats::coef(fit$model)[["z"]]), 8)
  expect_output(print(fit), "separation")

  # A strong but estimable signal must not be flagged.
  set.seed(72)
  ok <- data.frame(v = stats::rnorm(400))
  yok <- stats::rbinom(400, 1, stats::plogis(-0.5 + 1.5 * ok$v))
  clean <- expect_silent(fit_prognostic(ok, yok, folds = 5, seed = 72))
  expect_length(clean$separation, 0L)
  expect_true(clean$converged)
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

test_that("analyse_with_prognostic reports separation rather than a fake estimate", {
  # Reproduction of the defect: every treated patient has the event and no
  # control does. The maximum likelihood estimate does not exist, but glm stops
  # on the flat ridge and reports convergence, so the old code returned
  # estimate = 53.13 with se = 92293.68 and converged = TRUE.
  y <- c(rep(1L, 30), rep(0L, 30))
  arm <- c(rep(1L, 30), rep(0L, 30))
  set.seed(2)
  score <- stats::rnorm(60)

  expect_warning(res <- analyse_with_prognostic(y, arm, score), "separation")
  expect_false(res$converged)
  expect_true("arm" %in% res$separation)
  # The underlying glm still claims success; the guard is what catches it.
  expect_gt(abs(res$adjusted$estimate), 8)
  expect_gt(res$adjusted$se, 25)
  expect_output(print(res), "separation")

  # A well-behaved trial of the same size must not be flagged.
  set.seed(73)
  y2 <- stats::rbinom(60, 1, 0.4)
  clean <- expect_silent(analyse_with_prognostic(y2, arm, score))
  expect_true(clean$converged)
  expect_length(clean$separation, 0L)
  expect_length(clean$rank_deficient, 0L)
})

test_that("analyse_with_prognostic surfaces a rank deficient design", {
  # Reproduction of the defect: a score that is a linear function of arm. glm
  # drops it, every g-computation quantity becomes NA, and the old code returned
  # se_ratio = NA with converged = TRUE and no warning at all.
  set.seed(3)
  n <- 200
  arm <- rep(0:1, each = 100)
  y <- stats::rbinom(n, 1, 0.3)
  score <- 2 * arm

  expect_warning(res <- analyse_with_prognostic(y, arm, score), "rank deficient")
  expect_equal(res$rank_deficient, "score")
  expect_false(res$converged)
  expect_true(is.na(res$se_ratio))
  expect_true(is.na(res$adjusted$rd))
  expect_true(is.na(res$adjusted$rd_se))
  expect_output(print(res), "rank deficient")

  # Nudging the score off the arm restores an estimable model.
  set.seed(74)
  ok <- 2 * arm + stats::rnorm(n)
  fine <- expect_silent(analyse_with_prognostic(y, arm, ok))
  expect_length(fine$rank_deficient, 0L)
  expect_true(fine$converged)
  expect_false(is.na(fine$se_ratio))
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

# The expensive run. procova_gain fits two GLMs per replicate under both the
# null and the alternative, so this is the slowest test in the package, and the
# replication count is chosen deliberately rather than for speed.
#
# Preserving the type I error rate is the whole claim of this module, so the
# test has to be able to detect a violation. At alpha = 0.025 the Monte Carlo
# standard error of a rejection rate is sqrt(0.025 * 0.975 / nsim). At nsim =
# 400 that is 0.0078, so a 3 MCSE band admits a true rate near 0.048, which is
# almost double nominal: such a test would pass on a module that was badly
# broken. At nsim = 3000 the MCSE is 0.0028 and the band is 0.017 to 0.033,
# which is tight enough for the assertion to mean something.
#
# Do not lower this to speed the suite up. The runtime is the cost of the test
# being able to fail.
PROCOVA_NSIM <- 3000

gain_strong <- procova_gain(
  scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.22),
  prognostic_spec(),
  coefs = strong_coefs,
  n_per_arm = 300, nsim = PROCOVA_NSIM, train_n = 4000, seed = 101,
  calibrate_n = 30000
)

test_that("CRITICAL: prognostic adjustment does not inflate type I error", {
  mcse <- gain_strong$type1_adjusted_mcse
  expect_gt(gain_strong$prognostic$auc_cv, 0.75)
  expect_lt(abs(gain_strong$type1_adjusted - gain_strong$alpha), 3 * mcse)
  expect_lt(abs(gain_strong$rd_type1_adjusted - gain_strong$alpha), 3 * mcse)
  # The unadjusted analysis is the reference: adjustment must not be worse.
  expect_lt(abs(gain_strong$type1_unadjusted - gain_strong$alpha), 3 * mcse)
  expect_equal(mcse, sqrt(0.025 * 0.975 / PROCOVA_NSIM), tolerance = 0.15)
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

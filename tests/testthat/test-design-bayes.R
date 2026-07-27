bayes_scenario <- function(control_rate = 0.30, treatment_rate = 0.15) {
  scenario("ercp_acute_cholangitis",
    control_rate = control_rate, treatment_rate = treatment_rate
  )
}

# Exact reference for the decision quantity. Pr(p_t < p_c) is the expectation
# over the treatment posterior of the control posterior's upper tail, so it can
# be evaluated by one-dimensional quadrature with no simulation involved.
exact_prob_better <- function(events_t, n_t, events_c, n_c,
                              prior = c(1, 1), direction = "lower_is_better") {
  a_t <- prior[1] + events_t
  b_t <- prior[2] + n_t - events_t
  a_c <- prior[1] + events_c
  b_c <- prior[2] + n_c - events_c
  integrand <- function(x) {
    stats::dbeta(x, a_t, b_t) *
      stats::pbeta(x, a_c, b_c, lower.tail = identical(direction, "higher_is_better"))
  }
  stats::integrate(integrand, 0, 1, rel.tol = 1e-10)$value
}

mc_se <- function(draws) 0.5 / sqrt(draws)


test_that("posterior_prob_better is one half when both arms have identical data", {
  expect_equal(
    posterior_prob_better(10, 100, 10, 100, draws = 2e5, seed = 1),
    0.5,
    tolerance = 4 * mc_se(2e5) / 0.5
  )
  expect_equal(
    posterior_prob_better(0, 40, 0, 40, draws = 2e5, seed = 2),
    0.5,
    tolerance = 4 * mc_se(2e5) / 0.5
  )
  expect_equal(
    posterior_prob_better(7, 15, 7, 15, prior = c(3, 4), draws = 2e5, seed = 3),
    0.5,
    tolerance = 4 * mc_se(2e5) / 0.5
  )
})


test_that("a strongly better treatment arm gives a probability near one, and direction flips it", {
  lower <- posterior_prob_better(2, 200, 60, 200, draws = 2e4, seed = 4)
  higher <- posterior_prob_better(2, 200, 60, 200,
    direction = "higher_is_better", draws = 2e4, seed = 4
  )
  expect_gt(lower, 0.999)
  expect_lt(higher, 0.001)
  expect_equal(lower + higher, 1)

  worse <- posterior_prob_better(60, 200, 2, 200, draws = 2e4, seed = 5)
  expect_lt(worse, 0.001)
})


test_that("the Monte Carlo estimate matches exact numerical integration of the conjugate posteriors", {
  cases <- list(
    list(3, 20, 8, 20, c(1, 1), "lower_is_better"),
    list(0, 10, 4, 10, c(1, 1), "lower_is_better"),
    list(7, 15, 7, 15, c(1, 1), "lower_is_better"),
    list(1, 8, 6, 9, c(1, 1), "lower_is_better"),
    list(12, 30, 5, 30, c(1, 1), "lower_is_better"),
    list(12, 30, 5, 30, c(1, 1), "higher_is_better"),
    list(3, 20, 8, 20, c(2, 5), "lower_is_better"),
    list(9, 25, 4, 18, c(0.5, 0.5), "higher_is_better")
  )
  draws <- 2e5
  for (i in seq_along(cases)) {
    ca <- cases[[i]]
    got <- posterior_prob_better(ca[[1]], ca[[2]], ca[[3]], ca[[4]],
      prior = ca[[5]], direction = ca[[6]], draws = draws, seed = 100 + i
    )
    want <- exact_prob_better(ca[[1]], ca[[2]], ca[[3]], ca[[4]],
      prior = ca[[5]], direction = ca[[6]]
    )
    expect_lt(abs(got - want), 4 * mc_se(draws))
  }
})


test_that("posterior probability accuracy improves with the number of draws", {
  cases <- list(c(3, 20, 8, 20), c(12, 30, 5, 30), c(1, 8, 6, 9), c(9, 25, 4, 18))
  rmse <- function(draws) {
    err <- numeric(0)
    for (i in seq_along(cases)) {
      ca <- cases[[i]]
      want <- exact_prob_better(ca[1], ca[2], ca[3], ca[4])
      for (s in 1:4) {
        got <- posterior_prob_better(ca[1], ca[2], ca[3], ca[4],
          draws = draws, seed = 1000 * s + i
        )
        err <- c(err, got - want)
      }
    }
    sqrt(mean(err^2))
  }
  coarse <- rmse(500)
  fine <- rmse(5e4)
  expect_lt(fine, coarse)
  expect_lt(fine, 4 * mc_se(5e4))
  expect_lt(coarse, 4 * mc_se(500))
})


test_that("posterior_prob_better validates its arguments", {
  expect_error(posterior_prob_better(21, 20, 5, 20), "`events_t` \\(21\\) cannot exceed")
  expect_error(posterior_prob_better(5, 20, 21, 20), "`events_c` \\(21\\) cannot exceed")
  expect_error(posterior_prob_better(-1, 20, 5, 20), "`events_t`")
  expect_error(posterior_prob_better(5, 0, 5, 20), "`n_t`")
  expect_error(posterior_prob_better(5, 20, 5, 20, prior = 1), "`prior`")
  expect_error(posterior_prob_better(5, 20, 5, 20, prior = c(0, 1)), "`prior`")
  expect_error(posterior_prob_better(5, 20, 5, 20, direction = "sideways"), "`direction`")
  expect_error(posterior_prob_better(5, 20, 5, 20, draws = 1), "`draws`")
})


test_that("design_bayesian honours the shared gi_design contract", {
  sc <- bayes_scenario()
  d <- design_bayesian(sc, n_max = 200, looks = 3, nsim = 300, post_draws = 400, seed = 11)

  expect_s3_class(d, "gi_design")
  expect_identical(class(d), c("gi_design", "list"))
  expect_identical(d$type, "bayesian_adaptive")
  expect_identical(d$engine, "gitrialsim conjugate beta-binomial (simulation)")
  expect_identical(d$scenario, sc)
  expect_identical(d$n_total, 200L)
  expect_identical(d$n_per_arm, 100L)
  expect_true(is.numeric(d$power) && d$power >= 0 && d$power <= 1)
  expect_true(is.numeric(d$alpha) && d$alpha >= 0 && d$alpha <= 1)
  expect_true(all(c(
    "type", "scenario", "alpha", "power", "n_total", "n_per_arm", "engine", "detail"
  ) %in% names(d)))

  expect_identical(d$alpha, d$detail$alpha_null)
  expect_identical(d$power, d$detail$power)
  expect_equal(d$detail$power_mcse, sqrt(d$power * (1 - d$power) / 300))
  expect_equal(d$detail$alpha_mcse, sqrt(d$alpha * (1 - d$alpha) / 300))
  expect_identical(d$detail$look_n_total, c(66L, 134L, 200L))
  expect_equal(d$detail$information_rates, c(33, 67, 100) / 100)
  expect_identical(d$detail$rng_kind, "L'Ecuyer-CMRG")
  expect_null(d$detail$calibration)

  st <- d$detail$stopping
  expect_s3_class(st, "data.frame")
  expect_identical(nrow(st), 3L)
  expect_equal(sum(st$efficacy_alt), d$power)
  expect_equal(sum(st$efficacy_null), d$alpha)
  expect_equal(sum(st$efficacy_alt) + sum(st$futility_alt) + d$detail$inconclusive_alt, 1)
  expect_equal(sum(st$efficacy_null) + sum(st$futility_null) + d$detail$inconclusive_null, 1)

  expect_gte(d$detail$expected_n_alt, min(d$detail$look_n_total))
  expect_lte(d$detail$expected_n_alt, d$n_total)
  expect_gte(d$detail$expected_n_null, min(d$detail$look_n_total))
  expect_lte(d$detail$expected_n_null, d$n_total)
})


test_that("an odd n_max is rounded down to an even total and reported honestly", {
  d <- design_bayesian(bayes_scenario(), n_max = 201, looks = 1, nsim = 50, post_draws = 200, seed = 1)
  expect_identical(d$n_per_arm, 100L)
  expect_identical(d$n_total, 200L)
})


test_that("with a single look the trial always runs to the maximum sample size", {
  d <- design_bayesian(bayes_scenario(), n_max = 200, looks = 1, nsim = 200, post_draws = 300, seed = 6)
  expect_identical(nrow(d$detail$stopping), 1L)
  expect_equal(d$detail$expected_n_alt, 200)
  expect_equal(d$detail$expected_n_null, 200)
  expect_equal(d$detail$expected_n_alt_mcse, 0)
})


test_that("interim looks reduce the expected sample size below the maximum", {
  d <- design_bayesian(bayes_scenario(), n_max = 200, looks = 4, nsim = 400, post_draws = 300, seed = 8)
  expect_lt(d$detail$expected_n_alt, d$n_total)
  expect_lt(d$detail$expected_n_null, d$n_total)
})


test_that("design_bayesian validates its arguments", {
  sc <- bayes_scenario()
  expect_error(design_bayesian("ercp_acute_cholangitis", 100), "`scenario` must be a gi_scenario")
  expect_error(design_bayesian(sc, 3), "`n_max`")
  expect_error(design_bayesian(sc, 100.5), "`n_max`")
  expect_error(design_bayesian(sc, 100, looks = 200), "too many for n_max")
  expect_error(design_bayesian(sc, 100, prior = c(0, 1)), "`prior`")
  expect_error(
    design_bayesian(sc, 100, efficacy_threshold = 0.05, futility_threshold = 0.5),
    "must be below"
  )
  expect_error(design_bayesian(sc, 100, nsim = 1), "`nsim`")
  expect_error(design_bayesian(sc, 100, workers = 0), "`workers`")
})


test_that("a calibrated design controls simulated type I error at the target", {
  sc <- bayes_scenario()
  target <- 0.025
  nsim <- 2500
  cal <- calibrate_bayesian(sc,
    n_max = 200, target_alpha = target, looks = 2,
    nsim = nsim, post_draws = 1200, seed = 101
  )
  info <- cal$detail$calibration

  expect_true(info$bracketed)
  expect_gt(info$threshold, 0.5)
  expect_lt(info$threshold, 0.9999)
  expect_lte(info$achieved_alpha, target)
  expect_gt(info$achieved_alpha, target - 3 * info$achieved_alpha_mcse)
  expect_equal(info$achieved_alpha_mcse, sqrt(info$achieved_alpha * (1 - info$achieved_alpha) / nsim))

  # The design returned must be the design that was calibrated, not a rerun
  # with different numbers.
  expect_equal(cal$alpha, info$achieved_alpha)
  expect_equal(cal$detail$efficacy_threshold, info$threshold)

  # Out-of-sample check. Applying the calibrated threshold to null data from an
  # independent seed must land on the target within Monte Carlo error, which is
  # the claim that matters: in-sample agreement is guaranteed by construction.
  fresh <- design_bayesian(sc,
    n_max = 200, looks = 2, efficacy_threshold = info$threshold,
    nsim = nsim, post_draws = 1200, seed = 202
  )
  se <- sqrt(target * (1 - target) / nsim)
  expect_lt(abs(fresh$alpha - target), 3 * se)
})


test_that("more interim looks require a stricter calibrated threshold", {
  sc <- bayes_scenario()
  one <- calibrate_bayesian(sc,
    n_max = 200, target_alpha = 0.025, looks = 1,
    nsim = 1500, post_draws = 800, seed = 55
  )
  three <- calibrate_bayesian(sc,
    n_max = 200, target_alpha = 0.025, looks = 3,
    nsim = 1500, post_draws = 800, seed = 55
  )
  expect_gt(
    three$detail$calibration$threshold,
    one$detail$calibration$threshold
  )
})


test_that("a calibrated single-look design matches the rpact fixed design power", {
  skip_if_not_installed("rpact")
  sc <- bayes_scenario(control_rate = 0.30, treatment_rate = 0.15)
  nsim <- 4000
  d <- calibrate_bayesian(sc,
    n_max = 200, target_alpha = 0.025, looks = 1,
    nsim = nsim, post_draws = 2500, seed = 21
  )

  fixed <- rpact::getDesignGroupSequential(kMax = 1, alpha = 0.025, sided = 1)
  ref <- rpact::getPowerRates(fixed,
    maxNumberOfSubjects = 200, pi1 = sc$treatment_rate, pi2 = sc$control_rate,
    allocationRatioPlanned = 1, directionUpper = FALSE
  )

  # Tolerance, and why it is not tighter. Two things separate the procedures.
  # First, Monte Carlo error: the standard error of the simulated power here is
  # about 0.007, so three standard errors is already 0.021. Second, and more
  # fundamentally, the two tests are not the same procedure. The Bayesian rule
  # rejects when the posterior probability of benefit under independent uniform
  # priors exceeds a simulation-calibrated threshold; rpact rejects when a
  # normal-approximation z statistic exceeds a fixed critical value. They agree
  # asymptotically but not exactly at finite n, so agreement can only ever be
  # approximate. A tolerance of 0.05 on the power scale is a meaningful check
  # that the simulation machinery is right without pretending the procedures
  # are identical.
  expect_lt(abs(d$power - ref$overallReject), 0.05)

  # The calibrated threshold should also land near the nominal 0.975, which is
  # what a single-look Bayesian rule with uniform priors approximates.
  expect_lt(abs(d$detail$calibration$threshold - 0.975), 0.02)
})


test_that("results are reproducible for a given seed and independent of the number of workers", {
  sc <- bayes_scenario()
  args <- list(
    scenario = sc, n_max = 200, looks = 3, nsim = 300,
    post_draws = 400, seed = 11
  )
  a <- do.call(design_bayesian, args)
  b <- do.call(design_bayesian, args)
  expect_identical(a, b)

  different <- do.call(design_bayesian, utils::modifyList(args, list(seed = 12)))
  expect_false(isTRUE(all.equal(a$power, different$power)))

  skip_on_os("windows")
  parallel_run <- do.call(design_bayesian, utils::modifyList(args, list(workers = 2)))
  expect_identical(parallel_run$detail$workers, 2L)

  strip <- function(d) {
    d$detail$workers <- NULL
    d
  }
  expect_identical(strip(a), strip(parallel_run))
})


test_that("calibration is reproducible and independent of the number of workers", {
  skip_on_os("windows")
  sc <- bayes_scenario()
  args <- list(
    scenario = sc, n_max = 200, target_alpha = 0.025, looks = 2,
    nsim = 400, post_draws = 400, seed = 77
  )
  serial <- do.call(calibrate_bayesian, args)
  forked <- do.call(calibrate_bayesian, utils::modifyList(args, list(workers = 2)))
  expect_identical(
    serial$detail$calibration$threshold,
    forked$detail$calibration$threshold
  )
  expect_identical(serial$power, forked$power)
})


test_that("simulation leaves the caller's random number state untouched", {
  sc <- bayes_scenario()

  set.seed(4242)
  before <- runif(3)
  set.seed(4242)
  invisible(design_bayesian(sc, n_max = 100, looks = 2, nsim = 20, post_draws = 100, seed = 5))
  after <- runif(3)
  expect_identical(before, after)
  expect_identical(RNGkind()[1], "Mersenne-Twister")

  set.seed(9)
  before <- runif(2)
  set.seed(9)
  invisible(posterior_prob_better(3, 20, 8, 20, draws = 1000, seed = 123))
  after <- runif(2)
  expect_identical(before, after)
})


test_that("a seeded posterior probability does not depend on the ambient stream", {
  set.seed(1)
  first <- posterior_prob_better(3, 20, 8, 20, draws = 5000, seed = 7)
  set.seed(999)
  runif(50)
  second <- posterior_prob_better(3, 20, 8, 20, draws = 5000, seed = 7)
  expect_identical(first, second)
})


test_that("calibration warns and reports when the target cannot be bracketed", {
  sc <- bayes_scenario()
  expect_warning(
    conservative <- calibrate_bayesian(sc,
      n_max = 100, target_alpha = 0.9, looks = 1,
      nsim = 200, post_draws = 200, seed = 1
    ),
    "already"
  )
  expect_false(conservative$detail$calibration$bracketed)
  expect_equal(conservative$detail$calibration$threshold, 0.5)

  expect_warning(
    strict <- calibrate_bayesian(sc,
      n_max = 600, target_alpha = 1e-5, looks = 5,
      nsim = 800, post_draws = 300, seed = 2
    ),
    "still"
  )
  expect_false(strict$detail$calibration$bracketed)
  expect_equal(strict$detail$calibration$threshold, 0.9999)
})


test_that("calibrate_bayesian validates its arguments", {
  sc <- bayes_scenario()
  expect_error(calibrate_bayesian("nope", 200), "`scenario` must be a gi_scenario")
  expect_error(calibrate_bayesian(sc, 200, target_alpha = 0), "`target_alpha`")
  expect_error(calibrate_bayesian(sc, 200, target_alpha = 1), "`target_alpha`")
  expect_error(calibrate_bayesian(sc, 200, futility_threshold = 0.6), "lower end of the efficacy search bracket")
  expect_error(calibrate_bayesian(sc, 200, tol = 0), "`tol`")
  expect_error(calibrate_bayesian(sc, 200, max_iter = 0), "`max_iter`")
})


test_that("direction is taken from the scenario and reverses the decision", {
  lower <- bayes_scenario(control_rate = 0.30, treatment_rate = 0.15)
  higher <- lower
  higher$direction <- "higher_is_better"
  higher$control_rate <- 0.15
  higher$treatment_rate <- 0.30

  a <- design_bayesian(lower,
    n_max = 200, looks = 1, efficacy_threshold = 0.975,
    nsim = 1000, post_draws = 1000, seed = 21
  )
  b <- design_bayesian(higher,
    n_max = 200, looks = 1, efficacy_threshold = 0.975,
    nsim = 1000, post_draws = 1000, seed = 21
  )
  expect_identical(a$detail$direction, "lower_is_better")
  expect_identical(b$detail$direction, "higher_is_better")
  expect_lt(abs(a$power - b$power), 0.05)

  # A treatment that is worse on a lower-is-better endpoint must have power
  # far below the type I error level, not above it.
  harmful <- bayes_scenario(control_rate = 0.15, treatment_rate = 0.30)
  c_design <- design_bayesian(harmful,
    n_max = 200, looks = 1, efficacy_threshold = 0.975,
    nsim = 1000, post_draws = 1000, seed = 21
  )
  expect_lt(c_design$power, 0.005)
})


test_that("the flagship ERCP scenario runs end to end at its published rates", {
  sc <- scenario("ercp_acute_cholangitis")
  expect_equal(sc$control_rate, 0.0658)
  expect_equal(sc$treatment_rate, 0.0395)
  d <- design_bayesian(sc, n_max = 3000, looks = 3, nsim = 300, post_draws = 500, seed = 3)
  expect_identical(d$n_total, 3000L)
  expect_identical(d$detail$look_n_total, c(1000L, 2000L, 3000L))
  expect_true(d$power > 0 && d$power < 1)
})

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
  d <- design_bayesian(sc, n_max = 200, looks = 3, nsim = 300, post_draws = 500, seed = 11)

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
  expect_identical(d$detail$look_n_total, c(66L, 134L, 200L))
  expect_equal(d$detail$information_rates, c(33, 67, 100) / 100)
  expect_identical(d$detail$rng_kind, "L'Ecuyer-CMRG")
  expect_null(d$detail$calibration)

  # The simulation settings must be reported as they were used, because every
  # error statement below is only interpretable against them.
  expect_identical(d$detail$nsim, 300L)
  expect_identical(d$detail$post_draws, 500L)
  expect_identical(d$detail$seed, 11L)

  # Deliberately not a restatement of the formula in the implementation, which
  # would pass whatever that formula computed. Three independent facts are
  # checked instead: the rates really are averages over the `nsim` replications
  # that were reported, so they are whole multiples of 1 / nsim; the standard
  # errors respect the bound the documentation claims for them; and they are
  # strictly positive for a rate that is neither 0 nor 1. The statistical
  # calibration of these standard errors against the spread of independent
  # replicates is checked in its own test below.
  expect_equal(d$power * 300, round(d$power * 300))
  expect_equal(d$alpha * 300, round(d$alpha * 300))
  expect_gt(d$detail$power_mcse, 0)
  expect_gt(d$detail$alpha_mcse, 0)
  expect_lte(d$detail$power_mcse, 0.5 / sqrt(300))
  expect_lte(d$detail$alpha_mcse, 0.5 / sqrt(300))
  expect_gt(d$detail$expected_n_alt_mcse, 0)

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


test_that("the reported Monte Carlo standard errors match the spread of independent replicates", {
  # What a Monte Carlo standard error claims is that repeating the simulation
  # under a different seed moves the estimate by about that much. That claim is
  # checked here directly, against sixteen independent replicates, rather than
  # by rewriting the implementation's own expression. An error of a factor of
  # two or more in either direction fails.
  sc <- bayes_scenario()
  runs <- lapply(601:616, function(s) {
    design_bayesian(sc, n_max = 200, looks = 2, nsim = 200, post_draws = 500, seed = s)
  })
  pull <- function(f) vapply(runs, f, numeric(1))

  observed_power_sd <- stats::sd(pull(function(d) d$power))
  observed_alpha_sd <- stats::sd(pull(function(d) d$alpha))
  observed_n_sd <- stats::sd(pull(function(d) d$detail$expected_n_alt))
  reported_power_mcse <- mean(pull(function(d) d$detail$power_mcse))
  reported_alpha_mcse <- mean(pull(function(d) d$detail$alpha_mcse))
  reported_n_mcse <- mean(pull(function(d) d$detail$expected_n_alt_mcse))

  expect_gt(observed_power_sd / reported_power_mcse, 0.5)
  expect_lt(observed_power_sd / reported_power_mcse, 2)
  expect_gt(observed_alpha_sd / reported_alpha_mcse, 0.5)
  expect_lt(observed_alpha_sd / reported_alpha_mcse, 2)
  expect_gt(observed_n_sd / reported_n_mcse, 0.5)
  expect_lt(observed_n_sd / reported_n_mcse, 2)

  # Quadrupling the number of replications must halve the standard error. A
  # standard error built on the wrong power of nsim survives the checks above
  # only by accident and fails here.
  small <- design_bayesian(sc, n_max = 200, looks = 2, nsim = 250, post_draws = 500, seed = 91)
  large <- design_bayesian(sc, n_max = 200, looks = 2, nsim = 1000, post_draws = 500, seed = 91)
  expect_equal(small$detail$power_mcse / large$detail$power_mcse, 2, tolerance = 0.15)
  expect_equal(
    small$detail$expected_n_alt_mcse / large$detail$expected_n_alt_mcse, 2,
    tolerance = 0.15
  )
})


test_that("an odd n_max is rounded down to an even total and reported honestly", {
  d <- design_bayesian(bayes_scenario(), n_max = 201, looks = 1, nsim = 50, post_draws = 500, seed = 1)
  expect_identical(d$n_per_arm, 100L)
  expect_identical(d$n_total, 200L)
})


test_that("with a single look the trial always runs to the maximum sample size", {
  d <- design_bayesian(bayes_scenario(), n_max = 200, looks = 1, nsim = 200, post_draws = 500, seed = 6)
  expect_identical(nrow(d$detail$stopping), 1L)
  expect_equal(d$detail$expected_n_alt, 200)
  expect_equal(d$detail$expected_n_null, 200)
  expect_equal(d$detail$expected_n_alt_mcse, 0)
})


test_that("interim looks reduce the expected sample size below the maximum", {
  d <- design_bayesian(bayes_scenario(), n_max = 200, looks = 4, nsim = 400, post_draws = 500, seed = 8)
  expect_lt(d$detail$expected_n_alt, d$n_total)
  expect_lt(d$detail$expected_n_null, d$n_total)
})


test_that("a more aggressive futility rule cuts the expected sample size under the null", {
  # Common random numbers: seed, nsim and post_draws are held fixed, so the
  # simulated posterior trajectories are the same in every arm of this
  # comparison and the only thing that moves is the stopping rule applied to
  # them. Raising the futility threshold can then only stop trials earlier, so
  # expected sample size must fall and the rejection rates cannot rise.
  sc <- bayes_scenario()
  run <- function(fut) {
    design_bayesian(sc,
      n_max = 240, looks = 4, futility_threshold = fut,
      nsim = 500, post_draws = 500, seed = 31
    )
  }
  ladder <- lapply(c(0.02, 0.10, 0.20, 0.35), run)
  en_null <- vapply(ladder, function(d) d$detail$expected_n_null, numeric(1))
  en_alt <- vapply(ladder, function(d) d$detail$expected_n_alt, numeric(1))
  power <- vapply(ladder, function(d) d$power, numeric(1))
  alpha <- vapply(ladder, function(d) d$alpha, numeric(1))

  expect_true(all(diff(en_null) < 0))
  expect_lt(en_null[4], en_null[1])
  expect_true(all(diff(en_alt) <= 0))
  expect_true(all(diff(power) <= 0))
  expect_true(all(diff(alpha) <= 0))

  # The saving has to be real, not a rounding artefact: the harshest rule here
  # stops a large share of null trials at the first look.
  expect_lt(en_null[4], en_null[1] - 20)
  expect_gt(ladder[[4]]$detail$stopping$futility_null[1], 0.10)

  # And it must show up where it comes from, in the futility column of the
  # stopping table rather than anywhere else.
  fut_null <- vapply(ladder, function(d) sum(d$detail$stopping$futility_null), numeric(1))
  expect_true(all(diff(fut_null) > 0))
})


test_that("the stopping data frame is complete and reconciles with the reported summaries", {
  d <- design_bayesian(bayes_scenario(),
    n_max = 240, looks = 3, nsim = 500, post_draws = 500, seed = 17
  )
  st <- d$detail$stopping
  nsim <- 500

  expect_identical(names(st), c(
    "look", "n_total", "information_rate",
    "efficacy_alt", "futility_alt", "efficacy_null", "futility_null"
  ))
  expect_identical(st$look, 1:3)
  expect_identical(st$n_total, d$detail$look_n_total)
  expect_identical(st$n_total, c(80L, 160L, 240L))
  expect_equal(st$information_rate, st$n_total / d$n_total)

  probs <- as.matrix(st[, c("efficacy_alt", "futility_alt", "efficacy_null", "futility_null")])
  expect_true(all(probs >= 0 & probs <= 1))
  # Every entry is a count of simulated trials divided by nsim, so it must be a
  # whole multiple of 1 / nsim.
  expect_equal(probs * nsim, round(probs * nsim))

  # Each simulated trial ends exactly once, under either hypothesis.
  expect_equal(sum(st$efficacy_alt) + sum(st$futility_alt) + d$detail$inconclusive_alt, 1)
  expect_equal(sum(st$efficacy_null) + sum(st$futility_null) + d$detail$inconclusive_null, 1)
  expect_equal(sum(st$efficacy_alt), d$power)
  expect_equal(sum(st$efficacy_null), d$alpha)

  # The alternative and the null are genuinely different simulations, not the
  # same column reported twice, and the real effect produces more efficacy
  # stops than the null does.
  expect_false(isTRUE(all.equal(st$efficacy_alt, st$efficacy_null)))
  expect_gt(sum(st$efficacy_alt), sum(st$efficacy_null))

  # Expected sample size rebuilt from the stopping probabilities. A trial stops
  # at look k for efficacy or futility, and anything still open at the last
  # look ends there inconclusively, so this reproduces the reported mean
  # exactly. It fails if the stopping table and the expected sample size
  # disagree about what happened.
  rebuild <- function(eff, fut, inconclusive) {
    p_stop <- eff + fut
    p_stop[length(p_stop)] <- p_stop[length(p_stop)] + inconclusive
    sum(st$n_total * p_stop)
  }
  expect_equal(
    rebuild(st$efficacy_alt, st$futility_alt, d$detail$inconclusive_alt),
    d$detail$expected_n_alt
  )
  expect_equal(
    rebuild(st$efficacy_null, st$futility_null, d$detail$inconclusive_null),
    d$detail$expected_n_null
  )

  # A single-look design puts everything in one row and nothing stops early.
  one <- design_bayesian(bayes_scenario(),
    n_max = 240, looks = 1, nsim = 200, post_draws = 500, seed = 17
  )
  st1 <- one$detail$stopping
  expect_identical(nrow(st1), 1L)
  expect_identical(st1$n_total, 240L)
  expect_equal(st1$information_rate, 1)
  expect_equal(st1$efficacy_alt + st1$futility_alt + one$detail$inconclusive_alt, 1)
})


test_that("too few posterior draws bias the design, and no amount of nsim repairs it", {
  # The stopping rule compares an *estimated* posterior probability against a
  # fixed threshold, so Monte Carlo error in that estimate does not average out
  # across replications: it pushes trials across the boundary and biases the
  # decision. This test drives the simulation machinery directly, below the
  # guardrail that design_bayesian now applies, so the failure mode stays
  # visible in CI instead of only being described in a comment.
  sc <- bayes_scenario()
  look_n <- 100L
  snap <- rng_snapshot()
  on.exit(rng_restore(snap), add = TRUE)

  alpha_at <- function(post_draws, nsim, seed = 5) {
    streams <- bayes_stream_sets(seed, nsim)$null
    traj <- bayes_trajectories(
      p_control = sc$control_rate, p_treatment = sc$control_rate,
      look_n = look_n, prior = c(1, 1), direction = "lower_is_better",
      post_draws = post_draws, streams = streams, workers = 1L
    )
    mean(bayes_decide(traj, look_n, 0.975, 0.10)$reason == "efficacy")
  }

  coarse <- alpha_at(2, 4000)
  middling <- alpha_at(25, 4000)
  resolved <- alpha_at(1000, 4000)

  # Two posterior draws turn a 0.025 test into a 0.33 test.
  expect_gt(coarse, 0.25)
  expect_gt(coarse, 10 * resolved)
  # Twenty five draws are still visibly inflated.
  expect_gt(middling, 1.2 * resolved)
  # A thousand draws land on the nominal level.
  expect_lt(abs(resolved - 0.025), 0.01)

  # The bias is a property of post_draws alone. Multiplying nsim by eight
  # shrinks the Monte Carlo standard error by a factor of nearly three and
  # leaves the answer just as wrong, which is the whole point: nsim buys
  # precision, post_draws buys accuracy.
  small_nsim <- alpha_at(2, 1000)
  large_nsim <- alpha_at(2, 8000)
  expect_gt(small_nsim, 0.25)
  expect_gt(large_nsim, 0.25)
  expect_lt(abs(large_nsim - small_nsim), 0.05)
})


test_that("post_draws must resolve the efficacy threshold in use", {
  sc <- bayes_scenario()
  # The distance from 0.5 to the threshold has to cover ten standard errors of
  # the estimate, whose largest value is 0.5 / sqrt(post_draws). At 0.975 that
  # is 111 draws, and twenty standard errors, the recommended value, is 444.
  expect_identical(bayes_post_draws_needed(0.975, 10), 111L)
  expect_identical(bayes_post_draws_needed(0.975, 20), 444L)
  expect_identical(bayes_post_draws_floor(), 100L)
  expect_identical(bayes_post_draws_needed(0.5, 10), NA_integer_)

  expect_error(
    design_bayesian(sc, n_max = 200, looks = 1, nsim = 20, post_draws = 2, seed = 1),
    "`post_draws` = 2 is too few"
  )
  expect_error(
    design_bayesian(sc, n_max = 200, looks = 1, nsim = 20, post_draws = 110, seed = 1),
    "at least 111 are needed"
  )
  # A threshold nearer 0.5 is harder to resolve and demands more draws, so the
  # requirement is not a constant bolted on to the function.
  expect_error(
    design_bayesian(sc,
      n_max = 200, looks = 1, efficacy_threshold = 0.6,
      nsim = 20, post_draws = 500, seed = 1
    ),
    "at least 2500 are needed"
  )
  expect_warning(
    design_bayesian(sc, n_max = 200, looks = 1, nsim = 20, post_draws = 200, seed = 1),
    "near the resolution limit"
  )
  expect_silent(
    design_bayesian(sc, n_max = 200, looks = 1, nsim = 20, post_draws = 444, seed = 1)
  )

  # calibrate_bayesian applies the same rule against the least favourable
  # threshold its search can return, and fails before simulating anything.
  expect_error(
    calibrate_bayesian(sc, n_max = 200, target_alpha = 0.025, looks = 1, nsim = 20, post_draws = 50),
    "`post_draws` = 50 is too few"
  )
  seen <- character()
  withCallingHandlers(
    calibrate_bayesian(sc,
      n_max = 200, target_alpha = 0.025, looks = 1,
      nsim = 100, post_draws = 200, seed = 1
    ),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("near the resolution limit", seen)))
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

  # A count argument beyond 32-bit integer range must produce the message the
  # function intends, not base R's "missing value where TRUE/FALSE needed" from
  # comparing against an NA produced by silent coercion.
  for (bad in list(3e9, 2^31, -3e9, .Machine$integer.max + 1)) {
    expect_error(design_bayesian(sc, 200, nsim = bad), "`nsim` must be a single whole number")
    expect_error(design_bayesian(sc, bad), "`n_max` must be a single whole number")
    expect_error(design_bayesian(sc, 200, seed = bad), "`seed` must be a single whole number")
    expect_error(design_bayesian(sc, 200, looks = bad), "`looks` must be a single whole number")
    expect_error(design_bayesian(sc, 200, workers = bad), "`workers` must be a single whole number")
  }
  expect_error(posterior_prob_better(5, 20, 5, 20, draws = 3e9), "`draws` must be a single whole number")
  expect_error(calibrate_bayesian(sc, 200, nsim = 3e9), "`nsim` must be a single whole number")
  expect_error(calibrate_bayesian(sc, 200, max_iter = 3e9), "`max_iter` must be a single whole number")

  # The largest representable value is still accepted, so the range check has
  # not been tightened past what an integer can hold.
  expect_identical(check_count(.Machine$integer.max, "x"), .Machine$integer.max)
  expect_identical(check_count(-.Machine$integer.max, "x", min = -.Machine$integer.max),
                   -.Machine$integer.max)
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
    post_draws = 500, seed = 11
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
    nsim = 400, post_draws = 500, seed = 77
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
  invisible(design_bayesian(sc, n_max = 100, looks = 2, nsim = 20, post_draws = 500, seed = 5))
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
      nsim = 800, post_draws = 500, seed = 2
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

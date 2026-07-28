# Draws of the absolute risk difference in 30 day mortality for early versus
# delayed ERCP, negative meaning the intervention saves lives. The prior is
# lognormal on the risk ratio, centred on the 0.60 implied by the
# ercp_acute_cholangitis pack, so no draw implies an impossible event rate.
control_rate <- 0.0658

effect_draws <- function(n = 4000, log_sd = 0.35, seed = 11) {
  set.seed(seed)
  control_rate * (stats::rlnorm(n, meanlog = log(0.6), sdlog = log_sd) - 1)
}

# Net benefit in deaths averted per patient, charging the intervention a fixed
# harm of 0.005 deaths equivalent per patient for the extra procedure.
deaths_averted <- function(d) c(standard_care = 0, early_ercp = -d - 0.005)

test_that("evpi returns the documented structure and a non-negative estimate", {
  e <- evpi(effect_draws(), deaths_averted)

  expect_s3_class(e, "gi_evpi")
  expect_true(all(
    c("evpi", "mcse", "n_draws", "options", "value_current_info", "prob_optimal") %in%
      names(e)
  ))
  expect_gte(e$evpi, 0)
  expect_gt(e$mcse, 0)
  expect_identical(e$n_draws, 4000L)
  expect_identical(e$options, c("standard_care", "early_ercp"))
  expect_identical(e$best_option_current, "early_ercp")
  expect_output(print(e), "expected value of perfect information")
})

test_that("evpi is exactly zero when one option dominates for every draw", {
  dominant <- function(d) c(always_better = 1 + 0 * d, always_worse = 0)
  e <- evpi(effect_draws(), dominant)

  expect_identical(e$evpi, 0)
  expect_identical(e$mcse, 0)
  expect_identical(e$best_option_current, "always_better")
  expect_identical(unname(e$prob_optimal[["always_better"]]), 1)
})

test_that("evpi is positive when the optimal decision differs across draws", {
  e <- evpi(effect_draws(), deaths_averted)

  expect_gt(e$evpi, 0)
  expect_gt(e$mcse, 0)
  expect_true(all(e$prob_optimal > 0))
  expect_equal(e$evpi, e$value_perfect_info - e$value_current_info)
})

test_that("evpi widens as the prior widens and vanishes as it narrows", {
  wide <- evpi(effect_draws(log_sd = 0.8), deaths_averted)
  narrow <- evpi(effect_draws(log_sd = 0.02), deaths_averted)

  expect_gt(wide$evpi, narrow$evpi)
  expect_lt(narrow$evpi, 1e-6)
})

test_that("evsi is non-negative, below evpi, and rises with trial size", {
  # EVSI cannot exceed EVPI because perfect information is the limit of sample
  # information as the trial grows without bound: no dataset can be worth more
  # than being told the truth. The comparison allows Monte Carlo slack because
  # both sides are simulation estimates.
  draws <- effect_draws()
  sizes <- c(50, 200, 800, 3000)
  fits <- lapply(sizes, function(n) {
    suppressWarnings(evsi_trial(
      draws,
      n_per_arm = n, control_rate = 0.0658,
      net_benefit_fn = deaths_averted, nsim = 400
    ))
  })
  est <- vapply(fits, function(f) f$evsi, numeric(1))
  mcse <- vapply(fits, function(f) f$mcse, numeric(1))
  ceiling_evpi <- fits[[1]]$evpi

  expect_true(all(est >= 0))
  expect_true(all(est <= ceiling_evpi + 3 * mcse))
  expect_true(all(mcse > 0))
  # EVSI is monotone non-decreasing in trial size as an estimand, but each
  # point on the grid is a separate Monte Carlo estimate, so a middle point
  # could dip below its neighbour by chance without the code being wrong.
  # Only the smallest-versus-largest comparison is tested, against a margin
  # of several combined Monte Carlo standard errors rather than zero.
  se_gap <- sqrt(mcse[1]^2 + mcse[length(mcse)]^2)
  expect_gt(est[length(est)], est[1] + 3 * se_gap)
})

test_that("evsi returns the documented structure", {
  f <- suppressWarnings(evsi_trial(
    effect_draws(), 400, 0.0658,
    net_benefit_fn = deaths_averted, nsim = 200, seed = 3
  ))

  expect_s3_class(f, "gi_evsi")
  expect_true(all(
    c("evsi", "mcse", "evpi", "n_per_arm", "nsim", "seed", "n_draws", "detail") %in%
      names(f)
  ))
  expect_identical(f$n_per_arm, 400L)
  expect_identical(f$nsim, 200L)
  expect_identical(f$seed, 3)
  expect_true(all(
    c("evsi_unpaired", "evsi_plugin", "ess_mean", "ess_min", "decisions") %in%
      names(f$detail)
  ))
  expect_gt(f$detail$ess_mean, 0)
  expect_output(print(f), "expected value of sample information")
})

test_that("the paired and unpaired estimators agree to within their Monte Carlo error", {
  f <- suppressWarnings(evsi_trial(
    effect_draws(n = 8000), 600, 0.0658,
    net_benefit_fn = deaths_averted, nsim = 3000, seed = 5
  ))

  gap <- abs(f$evsi - f$detail$evsi_unpaired)
  expect_lt(gap, 3 * f$detail$mcse_unpaired)
  # The paired form is the reason it is the headline estimate.
  expect_lt(f$mcse, f$detail$mcse_unpaired)
})

test_that("evsi is reproducible from its seed and leaves the RNG state alone", {
  draws <- effect_draws(n = 1000)
  args <- list(
    prior_draws = draws, n_per_arm = 300, control_rate = 0.0658,
    net_benefit_fn = deaths_averted, nsim = 100, seed = 8
  )

  set.seed(77)
  interrupted <- c(stats::runif(1), NA)
  a <- suppressWarnings(do.call(evsi_trial, args))
  interrupted[2] <- stats::runif(1)

  set.seed(77)
  expect_identical(interrupted, stats::runif(2))

  b <- suppressWarnings(do.call(evsi_trial, args))
  expect_identical(a$evsi, b$evsi)

  args$seed <- 9
  c9 <- suppressWarnings(do.call(evsi_trial, args))
  expect_false(identical(a$evsi, c9$evsi))
})

test_that("evsi accepts draws of both arm rates directly", {
  set.seed(4)
  p_control <- stats::rbeta(3000, 66, 934)
  p_treatment <- stats::rbeta(3000, 45, 955)
  nb <- function(p) c(standard_care = 0, early_ercp = (p[["p_control"]] - p[["p_treatment"]]) - 0.02)

  f <- evsi_trial(
    cbind(p_control = p_control, p_treatment = p_treatment),
    n_per_arm = 400, net_benefit_fn = nb, nsim = 200
  )

  expect_identical(f$options, c("standard_care", "early_ercp"))
  expect_gte(f$evsi, 0)
  expect_lte(f$evsi, f$evpi + 3 * f$mcse)
})

test_that("evsi warns when prior draws imply impossible event rates", {
  expect_warning(
    evsi_trial(
      stats::rnorm(2000, -0.03, 0.06), 200, 0.0658,
      net_benefit_fn = deaths_averted, nsim = 50
    ),
    "outside \\(0, 1\\)"
  )
})

test_that("evsi warns when the posterior rests on too few prior draws", {
  expect_warning(
    evsi_trial(
      stats::rnorm(200, mean = -0.026, sd = 0.008), 20000, 0.0658,
      net_benefit_fn = deaths_averted, nsim = 50
    ),
    "effective sample size"
  )
})

test_that("voi_curve returns one row per size with columns plot_evsi can read", {
  curve <- suppressWarnings(voi_curve(
    effect_draws(), n_grid = c(1600, 100, 400),
    control_rate = 0.0658, net_benefit_fn = deaths_averted, nsim = 300
  ))

  expect_s3_class(curve, "data.frame")
  expect_identical(nrow(curve), 3L)
  expect_identical(curve$n_per_arm, c(100L, 400L, 1600L))
  expect_identical(
    names(curve),
    c("n_per_arm", "evsi", "mcse", "evsi_unpaired", "ess_mean", "evpi", "fraction_of_evpi")
  )
  expect_true(all(diff(curve$evsi) > 0))
  expect_true(all(curve$fraction_of_evpi > 0))
  # Near saturation the estimated fraction can sit slightly above 1: the
  # inequality binds on the estimands, and the two sides are separate Monte
  # Carlo estimates, so only the mcse-adjusted comparison is guaranteed.
  expect_true(all(curve$evsi <= curve$evpi + 3 * curve$mcse))
  expect_identical(length(unique(curve$evpi)), 1L)
})

test_that("printing an evsi says there is nothing to buy rather than printing NaN", {
  # One option dominates every prior draw, so EVPI is exactly 0 and EVSI is a
  # percentage of nothing. The print method used to divide anyway and report
  # "EVSI is NaN percent of it".
  dominant <- function(d) c(always_better = 1 + 0 * d, always_worse = 0)
  f <- suppressWarnings(evsi_trial(
    effect_draws(n = 500), 100, 0.0658,
    net_benefit_fn = dominant, nsim = 50
  ))

  expect_identical(f$evpi, 0)
  out <- utils::capture.output(print(f))
  expect_false(any(grepl("NaN", out, fixed = TRUE)))
  expect_true(any(grepl("no information for a trial to buy", out, fixed = TRUE)))
})

test_that("voi_curve reports NA for the fraction of EVPI when EVPI is zero", {
  # The documented behaviour of fraction_of_evpi in the one case where the ratio
  # does not exist.
  dominant <- function(d) c(always_better = 1 + 0 * d, always_worse = 0)
  curve <- suppressWarnings(voi_curve(
    effect_draws(n = 500), n_grid = c(100, 200),
    control_rate = 0.0658, net_benefit_fn = dominant, nsim = 50
  ))

  expect_true(all(curve$evpi == 0))
  expect_identical(curve$fraction_of_evpi, rep(NA_real_, 2))
})

test_that("a fractional seed is refused rather than silently truncated", {
  # set.seed() truncates towards zero, so seed = 1.7 used to be recorded in the
  # result while the simulated trials actually came from seed = 1, giving a
  # recorded seed that does not reproduce its own result.
  draws <- effect_draws(n = 500)

  expect_error(
    evsi_trial(draws, 200, 0.0658, deaths_averted, nsim = 20, seed = 1.7),
    "whole number"
  )
  expect_error(
    voi_curve(draws, n_grid = 200, control_rate = 0.0658,
              net_benefit_fn = deaths_averted, nsim = 20, seed = 0.5),
    "whole number"
  )

  f <- suppressWarnings(evsi_trial(draws, 200, 0.0658, deaths_averted,
    nsim = 20, seed = 2
  ))
  again <- suppressWarnings(evsi_trial(draws, 200, 0.0658, deaths_averted,
    nsim = 20, seed = f$seed
  ))
  expect_identical(f$seed, 2)
  expect_identical(again$evsi, f$evsi)
})

test_that("population_evsi discounts a constant incidence stream correctly", {
  p <- population_evsi(0.002, incidence = 45000, horizon_years = 10, discount_rate = 0.03)

  expected_population <- sum(45000 / 1.03^(0:9))
  expect_equal(p$effective_population, expected_population)
  expect_equal(p$total, 0.002 * expected_population)
  expect_equal(p$undiscounted_population, 450000)
  expect_identical(nrow(p$annual), 10L)
  expect_identical(p$annual$discount_factor[1], 1)
  expect_equal(sum(p$annual$value), p$total)
})

test_that("population_evsi with no discounting is a plain multiplication", {
  p <- population_evsi(0.002, incidence = 1000, horizon_years = 5, discount_rate = 0)

  expect_equal(p$total, 0.002 * 1000 * 5)
  expect_equal(p$effective_population, p$undiscounted_population)
})

test_that("population_evsi is monotone in the horizon and decreasing in the rate", {
  short <- population_evsi(0.002, 1000, 5)
  long <- population_evsi(0.002, 1000, 20)
  steep <- population_evsi(0.002, 1000, 20, discount_rate = 0.1)

  expect_gt(long$total, short$total)
  expect_lt(steep$total, long$total)
})

test_that("voi functions reject malformed input with a message naming the argument", {
  draws <- effect_draws(n = 100)

  expect_error(evpi(draws, "not a function"), "`net_benefit_fn` must be a function")
  expect_error(evpi(draws, function(d) 1), "at least two options")
  expect_error(evpi(letters, deaths_averted), "`prior_draws` must be")
  expect_error(evpi(1, deaths_averted), "at least two draws")
  expect_error(
    evpi(draws, function(d) if (d < 0) c(1, 2) else c(1, 2, 3)),
    "must return the same options every time"
  )

  expect_error(
    evsi_trial(draws, 0, 0.0658, deaths_averted, nsim = 10),
    "`n_per_arm` must be"
  )
  expect_error(
    evsi_trial(draws, 100, net_benefit_fn = deaths_averted, nsim = 10),
    "`control_rate` is required"
  )
  expect_error(
    evsi_trial(draws, 100, 1.5, deaths_averted, nsim = 10),
    "`control_rate` must be a single number"
  )
  expect_error(
    evsi_trial(cbind(0.06, 0.04)[rep(1, 50), ], 100, 0.0658, deaths_averted, nsim = 10),
    "`control_rate` must be NULL"
  )
  expect_error(
    evsi_trial(draws, 100, 0.0658, deaths_averted, nsim = 1),
    "`nsim` must be"
  )
  expect_error(
    voi_curve(draws, n_grid = c(-1, 100), control_rate = 0.0658,
              net_benefit_fn = deaths_averted),
    "`n_grid` must be"
  )

  expect_error(population_evsi(0.002, incidence = 0, horizon_years = 5), "`incidence` must be")
  expect_error(population_evsi(0.002, 1000, horizon_years = 0), "`horizon_years` must be")
  expect_error(population_evsi(0.002, 1000, 5, discount_rate = 1), "`discount_rate` must be")
  expect_error(population_evsi(NA_real_, 1000, 5), "`evsi_per_patient` must be")
})

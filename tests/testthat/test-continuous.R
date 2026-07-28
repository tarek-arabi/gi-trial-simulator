# Synthetic continuous endpoints, added to a clone of the shipped
# ercp_acute_cholangitis pack rather than to the shipped YAML, since these
# values are illustrative (sedative dose, sedation quality) and not backed by
# a published source the way the pack's real endpoints are. "jagtap2026" is
# used as the `source` only because it is already declared in the pack's
# `sources`, which is what validate_pack() requires; it is not a claim that
# the anchoring trial reported these numbers.
continuous_pack <- function() {
  pack <- load_pack("ercp_acute_cholangitis")
  pack$endpoints$sedative_dose_mg <- list(
    label = "Sedative dose (mg midazolam-equivalent)",
    type = "continuous",
    direction = "lower_is_better",
    role = "secondary",
    control_mean = 50,
    treatment_mean = 40,
    sd = 20,
    source = "jagtap2026"
  )
  pack$endpoints$sedation_quality_score <- list(
    label = "Sedation quality score (1-10, higher is better)",
    type = "continuous",
    direction = "higher_is_better",
    role = "secondary",
    control_mean = 6,
    treatment_mean = 8,
    sd = 4,
    source = "jagtap2026"
  )
  validate_pack(pack)
}

# Default: lower_is_better, control_mean 50, treatment_mean 40, sd 20, so
# standardised effect |40 - 50| / 20 = 0.5, the same effect size as the
# textbook check below, which the quadrupling checks lean on.
cts <- function(...) {
  scenario(continuous_pack(), endpoint = "sedative_dose_mg", ...)
}

cts_higher <- function(...) {
  scenario(continuous_pack(), endpoint = "sedation_quality_score", ...)
}

# ---- scenario() and validate_pack() ---------------------------------------

test_that("scenario() carries a continuous endpoint through with endpoint_type set", {
  sc <- cts()

  expect_s3_class(sc, "gi_scenario")
  expect_identical(sc$endpoint_type, "continuous")
  expect_equal(sc$control_mean, 50)
  expect_equal(sc$treatment_mean, 40)
  expect_equal(sc$sd, 20)
  expect_identical(sc$direction, "lower_is_better")
  expect_null(sc$control_rate)
  expect_null(sc$treatment_rate)
  expect_false(sc$overridden)
})

test_that("a binary scenario from the same pack leaves the continuous fields NULL", {
  sc <- scenario(continuous_pack())

  expect_identical(sc$endpoint_type, "binary")
  expect_null(sc$control_mean)
  expect_null(sc$treatment_mean)
  expect_null(sc$sd)
})

test_that("control_mean/treatment_mean/sd override the pack for sensitivity analysis", {
  sc <- cts(treatment_mean = 35, sd = 15)

  expect_equal(sc$treatment_mean, 35)
  expect_equal(sc$control_mean, 50)
  expect_equal(sc$sd, 15)
  expect_true(sc$overridden)
})

test_that("scenario() refuses rate overrides on a continuous endpoint and vice versa", {
  expect_error(cts(control_rate = 0.1), "continuous")
  expect_error(
    scenario(continuous_pack(), control_mean = 10),
    "binary"
  )
})

test_that("validate_pack rejects a continuous endpoint with a non-positive sd", {
  for (bad_sd in c(0, -5)) {
    pack <- continuous_pack()
    pack$endpoints$sedative_dose_mg$sd <- bad_sd
    expect_error(validate_pack(pack), "sd must be a single positive finite number")
  }
})

test_that("validate_pack rejects a continuous endpoint with a missing mean", {
  pack <- continuous_pack()
  pack$endpoints$sedative_dose_mg$treatment_mean <- NULL
  expect_error(validate_pack(pack), "treatment_mean must be a single finite number")

  pack2 <- continuous_pack()
  pack2$endpoints$sedative_dose_mg$control_mean <- NULL
  expect_error(validate_pack(pack2), "control_mean must be a single finite number")
})

test_that("validate_pack rejects a binary endpoint that carries a mean", {
  pack <- continuous_pack()
  pack$endpoints$mortality_30d$control_mean <- 10
  expect_error(validate_pack(pack), "binary endpoint must not carry 'control_mean'")
})

test_that("validate_pack rejects a continuous endpoint that carries a rate", {
  pack <- continuous_pack()
  pack$endpoints$sedative_dose_mg$control_rate <- 0.1
  expect_error(validate_pack(pack), "continuous endpoint must not carry 'control_rate'")
})

# ---- design_fixed(): dispatch to rpact::getSampleSizeMeans -----------------

test_that("a continuous sample size matches a direct rpact::getSampleSizeMeans() call exactly", {
  d <- design_fixed(cts())

  rp <- rpact::getDesignGroupSequential(kMax = 1L, alpha = 0.025, beta = 0.1, sided = 1L)
  ss <- rpact::getSampleSizeMeans(
    rp,
    alternative = abs(40 - 50), stDev = 20, allocationRatioPlanned = 1
  )

  expect_identical(d$engine, "rpact::getSampleSizeMeans")
  expect_equal(d$detail$n_fixed, as.numeric(ss$nFixed))
  expect_identical(d$n_per_arm, ceiling(as.numeric(ss$nFixed1)))
  expect_equal(d$detail$standardised_effect, abs(40 - 50) / 20)
})

test_that("allocation ratio breaks symmetry for continuous designs the same way it does for rates", {
  d <- design_fixed(cts(), allocation_ratio = 2)

  rp <- rpact::getDesignGroupSequential(kMax = 1L, alpha = 0.025, beta = 0.1, sided = 1L)
  ss <- rpact::getSampleSizeMeans(
    rp,
    alternative = 10, stDev = 20, allocationRatioPlanned = 2
  )

  expect_identical(d$detail$n_treatment, ceiling(as.numeric(ss$nFixed1)))
  expect_identical(d$detail$n_control, ceiling(as.numeric(ss$nFixed2)))
  expect_gt(d$detail$n_treatment, d$detail$n_control)
  expect_true(is.na(d$n_per_arm))
})

test_that("halving the mean difference roughly quadruples n", {
  full <- design_fixed(cts())
  halved <- design_fixed(cts(treatment_mean = 45))

  ratio <- halved$n_total / full$n_total
  expect_gt(ratio, 3.5)
  expect_lt(ratio, 4.5)
})

test_that("doubling the sd roughly quadruples n", {
  full <- design_fixed(cts())
  doubled_sd <- design_fixed(cts(sd = 40))

  ratio <- doubled_sd$n_total / full$n_total
  expect_gt(ratio, 3.5)
  expect_lt(ratio, 4.5)
})

test_that("standardised effect 0.5 at alpha 0.025 one-sided and 90 percent power gives n per arm about 86", {
  # rpact's exact (non-central t) calculation gives nFixed1 = 85.03 for this
  # design, which ceiling() rounds up to 86; this is the textbook value quoted
  # for a two-sample t-test at this effect size, alpha and power. A tolerance
  # of 1 patient is used rather than exact equality because the underlying
  # rpact computation is a numerical root-find, not a closed-form formula, and
  # a one-patient difference either way would not indicate a broken dispatch,
  # only a different rpact minor version's numerics.
  # sedative_dose_mg is lower_is_better, so the treatment mean must be the
  # smaller of the two for this to be a genuine (not backwards) effect.
  d <- design_fixed(cts(control_mean = 0.5, treatment_mean = 0, sd = 1))

  expect_equal(d$detail$standardised_effect, 0.5)
  expect_lte(abs(d$n_per_arm - 86), 1)
})

test_that("direction is handled correctly for lower_is_better and higher_is_better continuous endpoints", {
  lower <- design_fixed(cts())
  higher <- design_fixed(cts_higher(control_mean = 40, treatment_mean = 50, sd = 20))

  expect_identical(lower$n_total, higher$n_total)
  expect_false(lower$detail$direction_upper)
  expect_true(higher$detail$direction_upper)
})

test_that("a continuous scenario whose means contradict its direction warns", {
  backwards <- cts(treatment_mean = 60)
  expect_warning(design_fixed(backwards), "assumes no treatment benefit")
})

# ---- design_group_sequential(): dispatch -----------------------------------

test_that("design_group_sequential dispatches to getSampleSizeMeans for a continuous scenario", {
  g <- design_group_sequential(cts(), k = 3)

  expect_identical(g$type, "group_sequential")
  expect_identical(g$engine, "rpact::getSampleSizeMeans")
  expect_equal(g$detail$standardised_effect, 0.5)
  expect_false(is.null(g$detail$expected_n_h1))
})

test_that("a group-sequential continuous design has a larger max n and smaller expected n than fixed", {
  f <- design_fixed(cts())
  g <- design_group_sequential(cts(), k = 3)

  expect_gt(g$n_total, f$n_total)
  expect_lt(g$detail$expected_n_h1, f$n_total)
  expect_lt(g$detail$expected_n_h1, g$n_total)
})

# ---- power_at(): dispatch to rpact::getPowerMeans --------------------------

test_that("power_at recovers the target power for a continuous design", {
  d <- design_fixed(cts())

  # Sample size is rounded up to a whole even number, so the realised power must
  # be at least the target and only slightly above it. Asserting equality to the
  # target would be asserting that the rounding does not happen. The upper bound
  # is what makes this a real check: it fails if the design is oversized.
  achieved <- power_at(d, control_mean = 50, treatment_mean = 40, sd = 20)
  expect_gte(achieved, 0.9)
  expect_lt(achieved, 0.91)
  expect_equal(
    power_at(d, control_mean = 50, treatment_mean = 50, sd = 20),
    0.025,
    tolerance = 1e-3
  )
})

test_that("power_at refuses mismatched arguments for the design's endpoint type", {
  binary_d <- design_fixed(scenario("ercp_acute_cholangitis"))
  continuous_d <- design_fixed(cts())

  expect_error(
    power_at(binary_d, control_mean = 50, treatment_mean = 40, sd = 20),
    "binary endpoint"
  )
  expect_error(
    power_at(continuous_d, control_rate = 0.1, treatment_rate = 0.05),
    "continuous endpoint"
  )
})

# ---- print.gi_design and print.gi_scenario for continuous ------------------

test_that("print() on a continuous design shows the means and the standardised effect", {
  out <- capture.output(print(design_fixed(cts())))
  flat <- paste(out, collapse = "\n")

  expect_true(grepl("50.0000", flat, fixed = TRUE))
  expect_true(grepl("40.0000", flat, fixed = TRUE))
  expect_true(grepl("standardised effect", flat, fixed = TRUE))
  expect_true(grepl("0.5000", flat, fixed = TRUE))
})

test_that("print() on a continuous scenario shows the means, not event rates", {
  out <- capture.output(print(cts()))
  flat <- paste(out, collapse = "\n")

  expect_true(grepl("50.0000", flat, fixed = TRUE))
  expect_true(grepl("40.0000", flat, fixed = TRUE))
  expect_true(grepl("common SD", flat, fixed = TRUE))
})

# ---- guards: functions that cannot handle a continuous scenario -----------

test_that("simulate_fixed refuses a continuous scenario with a clear message", {
  expect_error(
    simulate_fixed(cts(), n_per_arm = 40, nsim = 10, seed = 1),
    "simulate_fixed.*continuous|continuous.*simulate_fixed"
  )
})

test_that("simulate_group_sequential refuses a design built on a continuous scenario", {
  g <- design_group_sequential(cts(), k = 2)
  expect_error(
    simulate_group_sequential(g, nsim = 10, seed = 1),
    "simulate_group_sequential.*continuous|continuous.*simulate_group_sequential"
  )
})

test_that("simulate_grid refuses a continuous scenario in the grid", {
  expect_error(
    simulate_grid(list(cts()), design_fn = function(s) 40, nsim = 10, seed = 1),
    "simulate_grid.*continuous|continuous.*simulate_grid"
  )
})

test_that("simulate_grid refuses a design_fn that returns a continuous-scenario design", {
  binary_sc <- scenario("ercp_acute_cholangitis")
  bad_design_fn <- function(s) design_group_sequential(cts(), k = 2)

  expect_error(
    simulate_grid(list(binary_sc), design_fn = bad_design_fn, nsim = 10, seed = 1),
    "simulate_grid.*continuous|continuous.*simulate_grid"
  )
})

test_that("design_bayesian refuses a continuous scenario with a clear message", {
  expect_error(
    design_bayesian(cts(), n_max = 100, nsim = 20, post_draws = 500, seed = 1),
    "design_bayesian.*continuous|continuous.*design_bayesian"
  )
})

test_that("calibrate_bayesian refuses a continuous scenario with a clear message", {
  expect_error(
    calibrate_bayesian(cts(), n_max = 100, nsim = 20, post_draws = 500, seed = 1),
    "calibrate_bayesian.*continuous|continuous.*calibrate_bayesian"
  )
})

test_that("procova_gain refuses a continuous scenario with a clear message", {
  spec <- cohort_spec(list(covariate_spec("x1", "normal", mean = 0, sd = 1)))
  expect_error(
    procova_gain(cts(), spec, coefs = c(x1 = 1), n_per_arm = 50, nsim = 5),
    "procova_gain.*continuous|continuous.*procova_gain"
  )
})

test_that("ademp_skeleton refuses a continuous scenario with a clear message", {
  expect_error(
    ademp_skeleton(cts()),
    "ademp_skeleton.*continuous|continuous.*ademp_skeleton"
  )
})

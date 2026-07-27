ercp <- function(...) scenario("ercp_acute_cholangitis", ...)

test_that("design_fixed returns the shared gi_design contract", {
  d <- design_fixed(ercp())

  expect_s3_class(d, "gi_design")
  expect_true(all(
    c("type", "scenario", "alpha", "power", "n_total", "n_per_arm", "engine", "detail") %in%
      names(d)
  ))
  expect_identical(d$type, "fixed")
  expect_identical(d$engine, "rpact::getSampleSizeRates")
  expect_s3_class(d$scenario, "gi_scenario")
  expect_identical(d$alpha, 0.025)
  expect_identical(d$power, 0.9)
})

test_that("fixed sample size is a positive even number at 1:1 allocation", {
  d <- design_fixed(ercp())

  expect_gt(d$n_total, 0)
  expect_identical(d$n_total %% 2, 0)
  expect_identical(d$n_total, d$n_per_arm * 2L)
  expect_identical(d$detail$n_control, d$detail$n_treatment)
})

test_that("fixed sample size matches a direct rpact call exactly", {
  d <- design_fixed(ercp())

  rp <- rpact::getDesignGroupSequential(kMax = 1L, alpha = 0.025, beta = 0.1, sided = 1L)
  ss <- rpact::getSampleSizeRates(
    rp,
    pi1 = 0.0395, pi2 = 0.0658, allocationRatioPlanned = 1
  )

  expect_equal(d$detail$n_fixed, as.numeric(ss$nFixed))
  expect_identical(d$n_per_arm, ceiling(as.numeric(ss$nFixed1)))
  expect_equal(d$detail$critical_value, as.numeric(rp$criticalValues))
})

test_that("pi1 is the treatment arm, so allocation is applied to the right arm", {
  d <- design_fixed(ercp(), allocation_ratio = 2)

  ss <- rpact::getSampleSizeRates(
    rpact::getDesignGroupSequential(kMax = 1L, alpha = 0.025, beta = 0.1, sided = 1L),
    pi1 = 0.0395, pi2 = 0.0658, allocationRatioPlanned = 2
  )

  expect_identical(d$detail$n_treatment, ceiling(as.numeric(ss$nFixed1)))
  expect_identical(d$detail$n_control, ceiling(as.numeric(ss$nFixed2)))
  expect_gt(d$detail$n_treatment, d$detail$n_control)
})

test_that("halving the treatment effect increases the required sample size", {
  full <- design_fixed(ercp())
  half_rate <- 0.0658 - (0.0658 - 0.0395) / 2
  halved <- design_fixed(ercp(treatment_rate = half_rate))

  expect_gt(halved$n_total, full$n_total)
})

test_that("required sample size falls as power falls and as alpha rises", {
  base <- design_fixed(ercp())

  expect_lt(design_fixed(ercp(), power = 0.8)$n_total, base$n_total)
  expect_lt(design_fixed(ercp(), alpha = 0.05)$n_total, base$n_total)
})

test_that("a higher_is_better endpoint is sized like its mirror image", {
  lower <- design_fixed(ercp())

  mirrored <- lower$scenario
  mirrored$direction <- "higher_is_better"
  mirrored$control_rate <- 0.0395
  mirrored$treatment_rate <- 0.0658
  higher <- design_fixed(mirrored)

  expect_identical(higher$n_total, lower$n_total)
  expect_true(higher$detail$direction_upper)
  expect_false(lower$detail$direction_upper)
})

test_that("a scenario whose rates contradict its direction warns", {
  harm <- ercp(endpoint = "post_ercp_adverse_events")
  expect_warning(design_fixed(harm), "assumes no treatment benefit")
})

test_that("design_fixed validates its arguments", {
  expect_error(design_fixed("not a scenario"), "gi_scenario")
  expect_error(design_fixed(ercp(), alpha = 0), "alpha")
  expect_error(design_fixed(ercp(), alpha = 0.9), "alpha")
  expect_error(design_fixed(ercp(), power = 1), "power")
  expect_error(design_fixed(ercp(), power = 0.01), "must exceed")
  expect_error(design_fixed(ercp(), allocation_ratio = -1), "allocation_ratio")
})

test_that("power_at recovers the target power at the design rates", {
  d <- design_fixed(ercp())

  expect_equal(power_at(d, 0.0658, 0.0395), 0.9, tolerance = 1e-3)
  expect_equal(power_at(d, 0.0658, 0.0658), 0.025, tolerance = 1e-3)
  expect_lt(power_at(d, 0.0658, 0.055), power_at(d, 0.0658, 0.0395))
})

test_that("power_at validates its arguments", {
  d <- design_fixed(ercp())

  expect_error(power_at("nope", 0.1, 0.05), "gi_design")
  expect_error(power_at(d, 0, 0.05), "control_rate")
  expect_error(power_at(d, 0.1, 1), "treatment_rate")
})

test_that("design_group_sequential returns the shared gi_design contract", {
  g <- design_group_sequential(ercp(), k = 3)

  expect_s3_class(g, "gi_design")
  expect_identical(g$type, "group_sequential")
  expect_identical(g$detail$k, 3L)
  expect_identical(g$detail$type_of_design, "asOF")
  expect_length(g$detail$efficacy_z, 3L)
  expect_length(g$detail$information_rates, 3L)
  expect_false(is.null(g$detail$expected_n_h0))
  expect_false(is.null(g$detail$expected_n_h1))
})

test_that("group-sequential max n exceeds fixed n while expected n under H1 is smaller", {
  f <- design_fixed(ercp())
  g <- design_group_sequential(ercp(), k = 3)

  expect_gt(g$n_total, f$n_total)
  expect_lt(g$detail$expected_n_h1, f$n_total)
  expect_lt(g$detail$expected_n_h1, g$n_total)
})

test_that("O'Brien-Fleming boundaries decrease monotonically in z", {
  for (type in c("asOF", "OF")) {
    g <- design_group_sequential(ercp(), k = 4, type_of_design = type)
    expect_true(all(diff(g$detail$efficacy_z) < 0), info = type)
  }
})

test_that("Pocock boundaries are near constant, unlike O'Brien-Fleming", {
  pocock <- design_group_sequential(ercp(), k = 3, type_of_design = "P")
  obf <- design_group_sequential(ercp(), k = 3, type_of_design = "OF")

  expect_equal(diff(range(pocock$detail$efficacy_z)), 0, tolerance = 1e-6)
  expect_gt(diff(range(obf$detail$efficacy_z)), 1)
})

test_that("cumulative alpha spending is increasing and lands on alpha", {
  g <- design_group_sequential(ercp(), k = 3)

  expect_true(all(diff(g$detail$cumulative_alpha_spent) > 0))
  expect_equal(g$detail$cumulative_alpha_spent[3], 0.025, tolerance = 1e-6)
})

test_that("futility is wired through rpact rather than invented", {
  none <- design_group_sequential(ercp(), k = 3, futility = "none")
  nonbinding <- design_group_sequential(ercp(), k = 3, futility = "nonbinding_obf")
  binding <- design_group_sequential(ercp(), k = 3, futility = "binding_obf")

  expect_true(all(is.na(none$detail$futility_z)))
  expect_false(any(is.na(nonbinding$detail$futility_z[1:2])))
  expect_true(is.na(nonbinding$detail$futility_z[3]))

  expect_false(nonbinding$detail$binding_futility)
  expect_true(binding$detail$binding_futility)

  # A binding bound may be relied on, so it buys back alpha at the last look.
  expect_lt(
    binding$detail$efficacy_z[3],
    nonbinding$detail$efficacy_z[3]
  )
  # Non-binding futility does not change the efficacy bounds.
  expect_equal(nonbinding$detail$efficacy_z, none$detail$efficacy_z)
  # Adding futility costs sample size.
  expect_gt(nonbinding$n_total, none$n_total)
})

test_that("group-sequential boundaries match a direct rpact call exactly", {
  g <- design_group_sequential(ercp(), k = 3, type_of_design = "asOF")

  rp <- rpact::getDesignGroupSequential(
    kMax = 3L, alpha = 0.025, beta = 0.1, sided = 1L, typeOfDesign = "asOF"
  )
  expect_equal(g$detail$efficacy_z, as.numeric(rp$criticalValues))
  expect_equal(g$detail$cumulative_alpha_spent, as.numeric(rp$alphaSpent))
})

test_that("custom information rates are honoured", {
  ir <- c(0.4, 0.7, 1)
  g <- design_group_sequential(ercp(), k = 3, information_rates = ir)

  expect_equal(g$detail$information_rates, ir)
  equal_spaced <- design_group_sequential(ercp(), k = 3)
  expect_false(isTRUE(all.equal(g$detail$efficacy_z, equal_spaced$detail$efficacy_z)))
})

test_that("design_group_sequential validates its arguments", {
  expect_error(design_group_sequential(ercp(), k = 1), "`k`")
  expect_error(design_group_sequential(ercp(), k = 2.5), "`k`")
  expect_error(design_group_sequential(ercp(), type_of_design = "nope"), "type_of_design")
  expect_error(design_group_sequential(ercp(), futility = "maybe"), "arg")
  expect_error(
    design_group_sequential(ercp(), k = 3, information_rates = c(0.5, 1)),
    "length"
  )
  expect_error(
    design_group_sequential(ercp(), k = 3, information_rates = c(0.5, 0.4, 1)),
    "increasing"
  )
  expect_error(
    design_group_sequential(ercp(), k = 3, information_rates = c(0.3, 0.6, 0.9)),
    "end at 1"
  )
})

test_that("gs_boundaries returns one tidy row per analysis", {
  g <- design_group_sequential(ercp(), k = 3, futility = "nonbinding_obf")
  b <- gs_boundaries(g)

  expect_s3_class(b, "data.frame")
  expect_identical(nrow(b), 3L)
  expect_identical(
    names(b),
    c(
      "analysis", "information_rate", "n_cumulative", "efficacy_z",
      "futility_z", "cumulative_alpha_spent"
    )
  )
  expect_identical(b$analysis, 1:3)
  expect_true(all(diff(b$n_cumulative) > 0))
  expect_identical(b$n_cumulative[3], g$n_total)
  expect_equal(b$efficacy_z, g$detail$efficacy_z)
})

test_that("gs_boundaries also describes a fixed design as a single analysis", {
  b <- gs_boundaries(design_fixed(ercp()))

  expect_identical(nrow(b), 1L)
  expect_equal(b$information_rate, 1)
  expect_equal(b$efficacy_z, stats::qnorm(1 - 0.025), tolerance = 1e-6)
  expect_true(is.na(b$futility_z))
})

test_that("group-sequential power is close to target when simulated analytically", {
  g <- design_group_sequential(ercp(), k = 3)
  expect_equal(power_at(g, 0.0658, 0.0395), 0.9, tolerance = 5e-3)
})

test_that("print.gi_design is defined in exactly one source file", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package sources not available")

  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  defines <- vapply(files, function(f) {
    any(grepl("^\\s*print\\.gi_design\\s*<-", readLines(f, warn = FALSE)))
  }, logical(1))

  expect_identical(sum(defines), 1L)
  expect_identical(basename(names(defines)[defines]), "design_fixed.R")
})

test_that("print.gi_design works for every design type", {
  f <- design_fixed(ercp())
  g <- design_group_sequential(ercp(), k = 3, futility = "nonbinding_obf")

  expect_output(print(f), "<gi_design> fixed")
  expect_output(print(f), "rpact::getSampleSizeRates")
  expect_output(print(g), "<gi_design> group_sequential")
  expect_output(print(g), "efficacy_z")
  expect_output(print(g), "expected n")

  invisible(capture.output(returned <- print(f)))
  expect_identical(returned, f)
})

test_that("print.gi_design survives a design with no rpact boundaries", {
  bare <- structure(
    list(
      type = "bayesian_adaptive", scenario = ercp(), alpha = 0.025,
      power = NA_real_, n_total = 3000L, n_per_arm = 1500L,
      engine = "test", detail = list(seed = 42, nsim = 1000)
    ),
    class = c("gi_design", "list")
  )

  expect_output(print(bare), "bayesian_adaptive")
  expect_output(print(bare), "by simulation")
  expect_output(print(bare), "seed: 42")
})

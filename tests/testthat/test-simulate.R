sc_effect <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
sc_null <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.30)

test_that("simulate_fixed returns the documented gi_simulation structure", {
  sim <- simulate_fixed(sc_effect, n_per_arm = 40, nsim = 200, seed = 1)

  expect_s3_class(sim, "gi_simulation")
  expect_identical(sim$design_type, "fixed")
  expect_identical(sim$nsim, 200L)
  expect_identical(sim$n_per_arm, 40L)
  expect_identical(sim$n_total, 80L)
  expect_identical(sim$scenario, sc_effect)
  expect_equal(unname(sim$rates), c(0.30, 0.15))

  expect_true(all(
    c("z", "p", "reject", "risk_difference", "risk_difference_se",
      "ci_lower", "ci_upper", "n_total") %in% names(sim$results)
  ))
  expect_identical(nrow(sim$results), 200L)
  expect_true(is.logical(sim$results$reject))
  expect_true(all(sim$results$p >= 0 & sim$results$p <= 1))
  expect_identical(sim$results$reject, sim$results$p < sim$alpha)
  expect_equal(
    sim$results$risk_difference,
    sim$results$events_treatment / 40 - sim$results$events_control / 40
  )
})

test_that("simulate_fixed validates its arguments by name", {
  expect_error(simulate_fixed(list(), 40), "`scenario` must be a gi_scenario")
  expect_error(simulate_fixed(sc_effect, 0), "`n_per_arm`")
  expect_error(simulate_fixed(sc_effect, 40.5), "`n_per_arm`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 0), "`nsim`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, alpha = 1), "`alpha`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = NA), "`seed`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, workers = 0), "`workers`")
})

test_that("the vectorised score statistic matches stats::prop.test(correct = FALSE)", {
  n <- 60
  sim <- simulate_fixed(sc_effect, n_per_arm = n, nsim = 300, seed = 11)
  res <- sim$results
  total <- res$events_control + res$events_treatment
  usable <- utils::head(which(total > 0 & total < 2 * n), 60)
  expect_gt(length(usable), 40)

  for (i in usable) {
    row <- res[i, ]
    pt <- stats::prop.test(
      x = c(row$events_treatment, row$events_control),
      n = c(n, n),
      correct = FALSE,
      alternative = "less"
    )
    expect_equal(abs(row$z), unname(sqrt(pt$statistic)), tolerance = 1e-10)
    expect_equal(row$p, unname(pt$p.value), tolerance = 1e-10)
  }
})

test_that("the score statistic is oriented by the endpoint direction", {
  sc_up <- scenario("hrs_terlipressin")
  expect_identical(sc_up$direction, "higher_is_better")

  n <- 50
  sim <- simulate_fixed(sc_up, n_per_arm = n, nsim = 200, seed = 5)
  res <- sim$results
  total <- res$events_control + res$events_treatment
  usable <- utils::head(which(total > 0 & total < 2 * n), 40)

  for (i in usable) {
    row <- res[i, ]
    pt <- stats::prop.test(
      x = c(row$events_treatment, row$events_control),
      n = c(n, n),
      correct = FALSE,
      alternative = "greater"
    )
    expect_equal(row$p, unname(pt$p.value), tolerance = 1e-10)
  }
  # A treatment that raises the event rate must look good, not bad, here.
  expect_gt(mean(res$reject), 0.5)
})

test_that("simulated power agrees with rpact analytic power at the design sample size", {
  design <- design_fixed(sc_effect, alpha = 0.025, power = 0.9)
  analytic <- power_at(design, control_rate = 0.30, treatment_rate = 0.15)

  nsim <- 40000
  sim <- simulate_fixed(
    sc_effect,
    n_per_arm = design$n_per_arm, nsim = nsim, alpha = design$alpha, seed = 20260727
  )
  simulated <- mean(sim$results$reject)
  mcse <- sqrt(simulated * (1 - simulated) / nsim)

  expect_gt(analytic, 0.89)
  expect_lt(abs(simulated - analytic), 3 * mcse)
})

test_that("simulating under the null recovers the nominal alpha", {
  design <- design_fixed(sc_effect, alpha = 0.025, power = 0.9)

  nsim <- 40000
  sim <- simulate_fixed(
    sc_null,
    n_per_arm = design$n_per_arm, nsim = nsim, alpha = 0.025, seed = 8675309
  )
  simulated <- mean(sim$results$reject)
  mcse <- sqrt(simulated * (1 - simulated) / nsim)

  expect_lt(abs(simulated - 0.025), 3 * mcse)
})

test_that("simulation is reproducible and leaves the caller's RNG untouched", {
  a <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 42)
  b <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 42)
  c_diff <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 43)

  expect_identical(a$results, b$results)
  expect_false(identical(a$results$z, c_diff$results$z))

  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  invisible(simulate_fixed(sc_effect, 40, nsim = 100, seed = 7))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("workers does not change simulate_fixed or simulate_group_sequential results", {
  one <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 3, workers = 1)
  two <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 3, workers = 2)
  expect_identical(one$results, two$results)

  design <- design_group_sequential(sc_effect, k = 2, futility = "nonbinding_obf")
  gs_one <- simulate_group_sequential(design, nsim = 300, seed = 3, workers = 1)
  gs_two <- simulate_group_sequential(design, nsim = 300, seed = 3, workers = 2)
  expect_identical(gs_one$results, gs_two$results)
})

test_that("sim_seeds derives distinct, reproducible L'Ecuyer streams", {
  s1 <- sim_seeds(2026, 4)
  s2 <- sim_seeds(2026, 4)

  expect_length(s1, 4L)
  expect_identical(s1, s2)
  expect_length(unique(vapply(s1, function(x) paste(x, collapse = ","), character(1))), 4L)
  expect_true(all(vapply(s1, function(x) x[1L] %% 100L == 7L, logical(1))))
  expect_false(identical(sim_seeds(2027, 4), s1))

  set.seed(1234)
  before <- get(".Random.seed", envir = globalenv())
  invisible(sim_seeds(11, 3))
  expect_identical(get(".Random.seed", envir = globalenv()), before)

  expect_error(sim_seeds(1, 0), "`n`")
})

test_that("a seed stream is an accepted seed and reproduces itself", {
  stream <- sim_seeds(5, 2)[[2L]]
  a <- simulate_fixed(sc_effect, 40, nsim = 200, seed = stream)
  b <- simulate_fixed(sc_effect, 40, nsim = 200, seed = stream)
  expect_identical(a$results, b$results)
  expect_false(identical(a$results$z, simulate_fixed(sc_effect, 40, nsim = 200,
    seed = sim_seeds(5, 2)[[1L]]
  )$results$z))
})

test_that("simulate_grid returns one tidy row per scenario", {
  grid <- list(
    optimistic = sc_effect,
    null = sc_null
  )
  out <- simulate_grid(grid, design_fn = function(s) 161, nsim = 500, seed = 7)

  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 2L)
  expect_identical(out$scenario, c("optimistic", "null"))
  expect_true(all(
    c("rejection_rate", "rejection_rate_mcse", "bias", "coverage", "n_total_max") %in% names(out)
  ))
  expect_gt(out$rejection_rate[1L], 0.8)
  expect_lt(out$rejection_rate[2L], 0.1)
  expect_length(attr(out, "simulations"), 2L)

  expect_error(simulate_grid(list(), function(s) 100), "`scenario_grid`")
  expect_error(simulate_grid(sc_effect, "not a function"), "`design_fn`")
  expect_error(
    simulate_grid(sc_effect, function(s) "no", nsim = 10),
    "must return a gi_design"
  )
})

test_that("simulate_grid gives identical results with 1 and 2 workers", {
  grid <- list(a = sc_effect, b = sc_null, c = sc_effect)

  one <- simulate_grid(grid, design_fn = function(s) 80, nsim = 400, seed = 99, workers = 1)
  two <- simulate_grid(grid, design_fn = function(s) 80, nsim = 400, seed = 99, workers = 2)

  expect_equal(one, two)
  expect_identical(
    attr(one, "simulations")$b$results,
    attr(two, "simulations")$b$results
  )
})

test_that("simulate_grid accepts a design function returning a gi_design", {
  out <- simulate_grid(
    list(effect = sc_effect),
    design_fn = function(s) design_fixed(s, alpha = 0.025, power = 0.9),
    nsim = 500, seed = 4
  )
  expect_equal(out$n_per_arm_max, design_fixed(sc_effect, alpha = 0.025, power = 0.9)$n_per_arm)
  expect_identical(out$design_type, "fixed")
})

test_that("simulate_group_sequential reads its boundaries from the rpact design", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")
  sim <- simulate_group_sequential(design, nsim = 500, seed = 2)

  expect_s3_class(sim, "gi_simulation")
  expect_identical(sim$design_type, "group_sequential")
  expect_identical(
    sim$boundaries$critical_values,
    as.numeric(design$detail$rpact_design$criticalValues)
  )
  expect_identical(
    sim$boundaries$information_rates,
    as.numeric(design$detail$rpact_design$informationRates)
  )
  expect_identical(sim$boundaries$read_from, "rpact design object")

  expect_true(all(sim$results$stop_stage %in% seq_len(3L)))
  expect_true(all(sim$results$n_total <= design$n_total))
  expect_true(all(sim$results$stopped_for %in% c("efficacy", "futility", "final_analysis")))
  expect_identical(sim$results$reject, sim$results$stopped_for == "efficacy" |
    (sim$results$stop_stage == 3L & sim$results$reject))
  # Every efficacy stop must sit at or above the boundary it crossed.
  eff <- sim$results[sim$results$stopped_for == "efficacy", ]
  expect_true(all(eff$z >= sim$boundaries$critical_values[eff$stop_stage]))
})

test_that("a design without futility bounds never stops for futility", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "none")
  sim <- simulate_group_sequential(design, nsim = 1000, seed = 6)

  expect_true(all(is.infinite(sim$boundaries$futility_bounds)))
  expect_false(any(sim$results$stopped_for == "futility"))
})

test_that("group-sequential simulation reproduces rpact power and expected sample size", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")

  nsim <- 20000
  sim <- simulate_group_sequential(design, nsim = nsim, seed = 3)
  simulated <- mean(sim$results$reject)
  mcse <- sqrt(simulated * (1 - simulated) / nsim)
  analytic <- power_at(design, control_rate = 0.30, treatment_rate = 0.15)

  expect_lt(abs(simulated - analytic), 4 * mcse)
  # Expected sample size under the alternative, as rpact reports it. The
  # tolerance absorbs the rounding of interim looks to whole patients.
  expect_lt(
    abs(mean(sim$results$n_total) - design$detail$expected_n_h1),
    0.02 * design$n_total
  )
})

test_that("group-sequential type I error is controlled under the null", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")

  nsim <- 20000
  sim <- simulate_group_sequential(design, nsim = nsim, seed = 4, rates = 0.30)
  simulated <- mean(sim$results$reject)

  expect_equal(unname(sim$rates), c(0.30, 0.30))
  # The interim analyses are discrete, so agreement with the nominal 0.025 is
  # to within simulation noise plus a small discreteness allowance.
  expect_lt(abs(simulated - design$alpha), 0.005)
})

test_that("simulate_group_sequential validates design and rates", {
  fixed <- design_fixed(sc_effect)
  expect_error(simulate_group_sequential(fixed), "must have type 'group_sequential'")
  expect_error(simulate_group_sequential(list()), "`design` must be a gi_design")

  design <- design_group_sequential(sc_effect, k = 2, futility = "none")
  expect_error(simulate_group_sequential(design, nsim = 10, rates = c(0.3, 1.2)), "`rates`")
  expect_error(simulate_group_sequential(design, nsim = 10, rates = c(1, 2, 3)), "`rates`")

  bare <- design
  bare$detail <- list()
  expect_error(simulate_group_sequential(bare, nsim = 10), "will not recompute them")
})

test_that("boundaries can also be read from detail when no rpact object is carried", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")
  stripped <- design
  stripped$detail <- list(
    efficacy_z = design$detail$efficacy_z,
    futility_z = design$detail$futility_z,
    information_rates = design$detail$information_rates
  )

  from_rpact <- simulate_group_sequential(design, nsim = 500, seed = 21)
  from_detail <- simulate_group_sequential(stripped, nsim = 500, seed = 21)

  expect_identical(from_detail$boundaries$read_from, "design$detail")
  expect_identical(from_rpact$results, from_detail$results)
})

test_that("print.gi_simulation reports the headline operating characteristic", {
  sim <- simulate_fixed(sc_effect, 40, nsim = 200, seed = 1)
  expect_output(print(sim), "gi_simulation")
  expect_output(print(sim), "rejection rate")
})

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

test_that("simulated power sits just above rpact analytic power, by the characterised amount", {
  design <- design_fixed(sc_effect, alpha = 0.025, power = 0.9)
  analytic <- power_at(design, control_rate = 0.30, treatment_rate = 0.15)

  nsim <- 200000
  sim <- simulate_fixed(
    sc_effect,
    n_per_arm = design$n_per_arm, nsim = nsim, alpha = design$alpha, seed = 20260727
  )
  simulated <- mean(sim$results$reject)

  expect_gt(analytic, 0.89)

  # The simulator and rpact do not agree exactly, and are not supposed to.
  # rpact evaluates a normal-approximation formula for two rates; the simulator
  # runs the discrete binomial test a real trial would run. Characterised at
  # 200,000 replications across five scenarios (control rates 0.066 to 0.50,
  # per-arm sizes 161 to 1514), the simulated rejection rate ran between 0.003
  # and 0.007 above the rpact value, with the sign positive in every scenario.
  # The assertion pins that: the same side, and inside one percentage point.
  # A negative difference, or one above 0.01, means something has changed and
  # is worth investigating rather than absorbing into a wider tolerance.
  expect_gt(simulated, analytic)
  expect_lt(simulated - analytic, 0.01)
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

test_that("a recorded seed reproduces its results whatever generator the session was using", {
  ambient <- RNGkind()
  on.exit(
    RNGkind(
      kind = ambient[1L], normal.kind = ambient[2L], sample.kind = ambient[3L]
    ),
    add = TRUE
  )

  RNGkind("Mersenne-Twister")
  default <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 42)

  RNGkind("L'Ecuyer-CMRG")
  lecuyer <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 42)
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")

  RNGkind("Wichmann-Hill")
  wichmann <- simulate_fixed(sc_effect, 40, nsim = 300, seed = 42)
  expect_identical(RNGkind()[1L], "Wichmann-Hill")

  # sim$seed is only a reproduction recipe if it also fixes the generator.
  expect_identical(default$results, lecuyer$results)
  expect_identical(default$results, wichmann$results)
  expect_identical(default$rng_kind, c("Mersenne-Twister", "Inversion", "Rejection"))
  expect_identical(lecuyer$rng_kind, default$rng_kind)

  design <- design_group_sequential(sc_effect, k = 2, futility = "none")
  RNGkind("Mersenne-Twister")
  gs_default <- simulate_group_sequential(design, nsim = 300, seed = 42)
  RNGkind("L'Ecuyer-CMRG")
  gs_lecuyer <- simulate_group_sequential(design, nsim = 300, seed = 42)
  expect_identical(gs_default$results, gs_lecuyer$results)
  expect_identical(gs_default$rng_kind, c("Mersenne-Twister", "Inversion", "Rejection"))

  # A seed stream carries its own generator, so it records L'Ecuyer-CMRG.
  stream <- simulate_fixed(sc_effect, 40, nsim = 50, seed = sim_seeds(5, 1)[[1L]])
  expect_identical(stream$rng_kind[1L], "L'Ecuyer-CMRG")
})

test_that("simulate_fixed rejects seeds it could not reproduce results from", {
  expect_error(
    simulate_fixed(sc_effect, 40, nsim = 10, seed = c(1, 2)),
    "not a valid .Random.seed"
  )
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = 1.5), "whole number")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = 3e9), "whole number")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = Inf), "`seed`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = numeric(0)), "`seed`")
  expect_error(simulate_fixed(sc_effect, 40, nsim = 10, seed = "1"), "`seed`")

  design <- design_group_sequential(sc_effect, k = 2, futility = "none")
  expect_error(
    simulate_group_sequential(design, nsim = 10, seed = c(1, 2)),
    "not a valid .Random.seed"
  )

  # Rejecting a seed must not leave the session's generator in a broken state.
  expect_length(stats::runif(1), 1L)
  expect_length(sim_seeds(1, 2), 2L)
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

test_that("simulate_grid simulates a group-sequential design at each scenario's own rates", {
  design <- design_group_sequential(sc_effect, k = 2, futility = "none")
  out <- simulate_grid(
    list(effect = sc_effect, null = sc_null),
    design_fn = function(s) design,
    nsim = 2000, seed = 31
  )

  expect_identical(out$design_type, c("group_sequential", "group_sequential"))
  expect_equal(out$control_rate, c(0.30, 0.30))
  expect_equal(out$treatment_rate, c(0.15, 0.30))

  # The rates a row is simulated at are the row's own, not the ones the design
  # happened to be powered on. A design built under a 0.15 treatment rate and
  # handed a null scenario must produce a null rejection rate.
  expect_gt(out$rejection_rate[1L], 0.8)
  expect_lt(out$rejection_rate[2L], 0.06)
  expect_gt(out$rejection_rate[1L] - out$rejection_rate[2L], 0.5)

  sims <- attr(out, "simulations")
  expect_equal(unname(sims$effect$rates), c(0.30, 0.15))
  expect_equal(unname(sims$null$rates), c(0.30, 0.30))
  # The boundaries and the maximum sample size still come from the design.
  expect_identical(sims$null$boundaries$critical_values, sims$effect$boundaries$critical_values)
  expect_equal(out$n_total_max, rep(design$n_total, 2L))
})

test_that("simulate_grid refuses a design whose endpoint direction fights the scenario", {
  expect_error(
    simulate_grid(
      list(down = sc_effect),
      design_fn = function(s) design_group_sequential(scenario("hrs_terlipressin"), k = 2),
      nsim = 20, seed = 5
    ),
    "oriented against the wrong tail"
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

  res <- sim$results
  # Rejection is exactly an efficacy crossing. A replicate that reaches the
  # final analysis without crossing the boundary does not reject, and a
  # replicate stopped for futility never does.
  expect_identical(res$reject, res$stopped_for == "efficacy")
  # The recorded sample size is the one belonging to the stage it stopped at.
  expect_equal(res$n_total, 2 * sim$boundaries$n_per_arm_by_stage[res$stop_stage])

  # Every efficacy stop must sit at or above the boundary it crossed.
  eff <- res[res$stopped_for == "efficacy", ]
  expect_gt(nrow(eff), 0L)
  expect_true(all(eff$z >= sim$boundaries$critical_values[eff$stop_stage]))

  # Every futility stop must sit at or below its boundary, and only interim
  # analyses have one.
  fut <- res[res$stopped_for == "futility", ]
  expect_gt(nrow(fut), 0L)
  expect_true(all(fut$stop_stage < 3L))
  expect_true(all(fut$z <= sim$boundaries$futility_bounds[fut$stop_stage]))
  expect_false(any(fut$reject))

  # Everything else runs to the end and finishes below the final boundary.
  fin <- res[res$stopped_for == "final_analysis", ]
  expect_gt(nrow(fin), 0L)
  expect_true(all(fin$stop_stage == 3L))
  expect_true(all(fin$z < sim$boundaries$critical_values[3L]))
  expect_true(all(fin$n_total == design$n_total))
})

test_that("a design without futility bounds never stops for futility", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "none")
  sim <- simulate_group_sequential(design, nsim = 1000, seed = 6)

  expect_true(all(is.infinite(sim$boundaries$futility_bounds)))
  expect_false(any(sim$results$stopped_for == "futility"))
})

test_that("group-sequential simulation tracks rpact power and expected sample size", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")

  nsim <- 100000
  sim <- simulate_group_sequential(design, nsim = nsim, seed = 3)
  simulated <- mean(sim$results$reject)
  analytic <- power_at(design, control_rate = 0.30, treatment_rate = 0.15)

  # The same discrete-versus-normal-approximation gap that simulate_fixed()
  # shows: the simulated rejection rate sits a little above rpact's, never
  # below, and well inside one percentage point.
  expect_gt(simulated, analytic)
  expect_lt(simulated - analytic, 0.01)

  # No rounding of interim looks happens for this design, because its
  # information rates divide the per-arm sample size exactly. Anything that
  # follows is therefore not a rounding artefact.
  expect_equal(
    sim$boundaries$n_per_arm_by_stage,
    design$detail$information_rates * design$n_per_arm
  )

  # rpact's expected sample size under H1 comes from the same normal
  # approximation, so the realised mean sits slightly above it: about 2
  # patients on a maximum of 342, again on the upper side.
  mean_n <- mean(sim$results$n_total)
  expect_gt(mean_n, design$detail$expected_n_h1)
  expect_lt(mean_n - design$detail$expected_n_h1, 0.01 * design$n_total)
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

test_that("degenerate replicates are counted, not silently absorbed", {
  rare <- scenario("ercp_acute_cholangitis", control_rate = 0.01, treatment_rate = 0.005)
  sim <- simulate_fixed(rare, n_per_arm = 12, nsim = 500, seed = 2)
  events <- sim$results$events_control + sim$results$events_treatment
  undefined <- events == 0L | events == 24L

  expect_gt(sim$degenerate$any, 0L)
  expect_identical(sim$degenerate$score_undefined, sum(undefined))
  expect_identical(sim$degenerate$any, sum(sim$results$degenerate))
  expect_identical(sim$results$degenerate, undefined)
  # The substituted statistic is what makes them look like ordinary failures
  # to reject, which is exactly why the count has to be reported.
  expect_true(all(sim$results$z[undefined] == 0))
  expect_false(any(sim$results$reject[undefined]))
  expect_output(print(sim), "replications degenerate")
})

test_that("a zero-width interval counts as degenerate even when the score test is defined", {
  extreme <- scenario("ercp_acute_cholangitis", control_rate = 0.99, treatment_rate = 0.01)
  sim <- simulate_fixed(extreme, n_per_arm = 5, nsim = 300, seed = 4)
  split <- sim$results$events_control == 5L & sim$results$events_treatment == 0L

  expect_gt(sum(split), 0L)
  expect_true(all(sim$results$risk_difference_se[split] == 0))
  expect_true(all(sim$results$ci_lower[split] == sim$results$ci_upper[split]))
  expect_identical(
    sim$degenerate$interval_undefined,
    sum(sim$results$risk_difference_se == 0)
  )
  # Here the pooled score test is perfectly well defined, so the two counts
  # have to be kept apart rather than reported as one number.
  expect_gt(sim$degenerate$interval_undefined, sim$degenerate$score_undefined)
  expect_identical(sim$degenerate$any, sum(sim$results$degenerate))
})

test_that("group-sequential degeneracy is counted at the stopping analysis", {
  rare <- scenario("ercp_acute_cholangitis", control_rate = 0.01, treatment_rate = 0.005)
  design <- design_group_sequential(rare, k = 2, futility = "none")
  small <- design
  small$n_per_arm <- 20L
  small$n_total <- 40L

  sim <- simulate_group_sequential(small, nsim = 500, seed = 8)
  events <- sim$results$events_control + sim$results$events_treatment
  undefined <- events == 0L | events == (2L * sim$results$n_per_arm)

  expect_gt(sim$degenerate$any, 0L)
  expect_identical(sim$degenerate$score_undefined, sum(undefined))
  expect_identical(sim$degenerate$any, sum(sim$results$degenerate))
  expect_true(all(sim$results$z[undefined] == 0))
})

test_that("print.gi_simulation reports the headline operating characteristic", {
  sim <- simulate_fixed(sc_effect, 40, nsim = 200, seed = 1)
  expect_output(print(sim), "gi_simulation")
  expect_output(print(sim), "rejection rate")
})

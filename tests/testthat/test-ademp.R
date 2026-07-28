sc_effect <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
sc_null <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.30)

test_that("ademp_summary reports every measure with a Monte Carlo standard error", {
  sim <- simulate_fixed(sc_effect, n_per_arm = 161, nsim = 2000, seed = 1)
  out <- ademp_summary(sim)

  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("measure", "estimate", "mcse", "definition"))
  expect_identical(
    out$measure,
    c(
      "rejection_rate", "bias", "empirical_se", "mse", "coverage",
      "mean_sample_size", "degenerate_rate"
    )
  )
  expect_false(anyNA(out$estimate))
  expect_false(anyNA(out$mcse))
  expect_true(all(out$mcse >= 0))
  expect_true(all(nzchar(out$definition)))
})

test_that("ademp_summary Monte Carlo standard errors follow Morris, White and Crowther", {
  sim <- simulate_fixed(sc_effect, n_per_arm = 100, nsim = 3000, seed = 12)
  out <- ademp_summary(sim)
  res <- sim$results
  nsim <- sim$nsim
  truth <- 0.15 - 0.30
  theta <- res$risk_difference
  error <- theta - truth

  get_row <- function(m) out[out$measure == m, ]

  rejection <- mean(res$reject)
  expect_equal(get_row("rejection_rate")$estimate, rejection)
  expect_equal(get_row("rejection_rate")$mcse, sqrt(rejection * (1 - rejection) / nsim))

  expect_equal(get_row("bias")$estimate, mean(theta) - truth)
  expect_equal(get_row("bias")$mcse, stats::sd(theta) / sqrt(nsim))

  emp_se <- stats::sd(theta)
  expect_equal(get_row("empirical_se")$estimate, emp_se)
  expect_equal(get_row("empirical_se")$mcse, emp_se / sqrt(2 * (nsim - 1)))

  expect_equal(get_row("mse")$estimate, mean(error^2))
  expect_equal(get_row("mse")$mcse, stats::sd(error^2) / sqrt(nsim))

  coverage <- mean(res$ci_lower <= truth & truth <= res$ci_upper)
  expect_equal(get_row("coverage")$estimate, coverage)
  expect_equal(get_row("coverage")$mcse, sqrt(coverage * (1 - coverage) / nsim))

  expect_equal(get_row("mean_sample_size")$estimate, 200)
  expect_equal(get_row("mean_sample_size")$mcse, 0)

  degenerate <- sum(res$degenerate)
  expect_equal(get_row("degenerate_rate")$estimate, degenerate / nsim)
  expect_equal(
    get_row("degenerate_rate")$mcse,
    sqrt((degenerate / nsim) * (1 - degenerate / nsim) / nsim)
  )
})

test_that("ademp_summary reports how many replicates were degenerate", {
  rare <- scenario("ercp_acute_cholangitis", control_rate = 0.01, treatment_rate = 0.005)
  sim <- simulate_fixed(rare, n_per_arm = 12, nsim = 500, seed = 2)
  row <- ademp_summary(sim)[ademp_summary(sim)$measure == "degenerate_rate", ]

  expect_gt(sim$degenerate$any, 0L)
  expect_equal(row$estimate, sim$degenerate$any / sim$nsim)
  # The count itself has to be readable, not just the proportion, and the
  # definition has to say what happened to those replicates.
  expect_match(row$definition, paste0(sim$degenerate$any, " of ", sim$nsim), fixed = TRUE)
  expect_match(row$definition, "non-rejections")

  # The rejection rate is computed over every replicate, degenerate ones
  # included, which is what makes the count worth reporting alongside it.
  full <- ademp_summary(sim)
  expect_equal(full$estimate[1L], mean(sim$results$reject))

  # A simulation carrying no record of degeneracy is refused rather than
  # reported as having none.
  stripped <- sim
  stripped$degenerate <- NULL
  stripped$results$degenerate <- NULL
  expect_error(ademp_summary(stripped), "no record of degenerate replications")
})

test_that("ademp_summary recovers the truth it was generated from", {
  nsim <- 20000
  sim <- simulate_fixed(sc_effect, n_per_arm = 161, nsim = nsim, seed = 314159)
  out <- ademp_summary(sim)

  bias <- out[out$measure == "bias", ]
  expect_lt(abs(bias$estimate), 3 * bias$mcse)

  coverage <- out[out$measure == "coverage", ]
  expect_lt(abs(coverage$estimate - 0.95), 0.01)

  # MSE is bias squared plus the empirical variance, to Monte Carlo accuracy.
  emp_se <- out[out$measure == "empirical_se", ]$estimate
  mse <- out[out$measure == "mse", ]$estimate
  expect_equal(mse, bias$estimate^2 + emp_se^2 * (nsim - 1) / nsim, tolerance = 1e-8)
})

test_that("ademp_summary labels rejection as power or type I error", {
  power_sim <- simulate_fixed(sc_effect, n_per_arm = 161, nsim = 500, seed = 2)
  null_sim <- simulate_fixed(sc_null, n_per_arm = 161, nsim = 500, seed = 2)

  power_def <- ademp_summary(power_sim)$definition[1L]
  null_def <- ademp_summary(null_sim)$definition[1L]

  expect_match(power_def, "power")
  expect_match(null_def, "type I error")
  expect_equal(ademp_summary(null_sim)$estimate[1L], mean(null_sim$results$reject))
})

test_that("ademp_summary honours an explicit truth", {
  sim <- simulate_fixed(sc_effect, n_per_arm = 100, nsim = 1000, seed = 3)
  default <- ademp_summary(sim)
  shifted <- ademp_summary(sim, truth = 0)

  expect_equal(shifted$estimate[2L], default$estimate[2L] - 0.15)
  expect_lt(shifted$estimate[5L], default$estimate[5L])
  expect_error(ademp_summary(sim, truth = c(0, 1)), "`truth`")
  expect_error(ademp_summary(list()), "`sim` must be a gi_simulation")
})

test_that("ademp_summary adds an early stopping measure for group-sequential designs", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")
  sim <- simulate_group_sequential(design, nsim = 2000, seed = 5)
  out <- ademp_summary(sim)

  expect_true("early_stopping_rate" %in% out$measure)
  early <- out[out$measure == "early_stopping_rate", ]
  expect_equal(early$estimate, mean(sim$results$stop_stage < 3L))
  expect_equal(early$mcse, sqrt(early$estimate * (1 - early$estimate) / sim$nsim))

  mean_n <- out[out$measure == "mean_sample_size", ]
  expect_equal(mean_n$estimate, mean(sim$results$n_total))
  expect_gt(mean_n$mcse, 0)
  expect_lt(mean_n$estimate, design$n_total)
})

test_that("group-sequential estimates at the stopping analysis are biased away from the null", {
  design <- design_group_sequential(sc_effect, k = 3, futility = "nonbinding_obf")
  sim <- simulate_group_sequential(design, nsim = 20000, seed = 9)
  gs_bias <- ademp_summary(sim)$estimate[2L]

  fixed_sim <- simulate_fixed(sc_effect, n_per_arm = design$n_per_arm, nsim = 20000, seed = 9)
  fixed_bias <- ademp_summary(fixed_sim)$estimate[2L]

  # The risk difference is negative under benefit, so bias away from the null
  # is a more negative value than the unbiased fixed-design estimate.
  expect_lt(gs_bias, fixed_bias)
})

test_that("nsim_required inverts the Monte Carlo standard error of a proportion", {
  expect_identical(nsim_required(0.005), 10000L)
  expect_identical(nsim_required(0.005, expected_proportion = 0.9), 3600L)
  expect_identical(nsim_required(0.01), 2500L)

  n <- nsim_required(0.004, expected_proportion = 0.8)
  expect_lte(sqrt(0.8 * 0.2 / n), 0.004)

  expect_error(nsim_required(0), "`target_mcse`")
  expect_error(nsim_required(-1), "`target_mcse`")
  expect_error(nsim_required(0.01, expected_proportion = 1.5), "`expected_proportion`")
})

test_that("nsim_required never returns an unusable replication count", {
  # A degenerate proportion inverts to zero replications, which is never an
  # answer anyone can run, so it is refused and the caller is pointed at 0.5.
  expect_error(
    nsim_required(0.005, expected_proportion = 0),
    "strictly between 0 and 1"
  )
  expect_error(
    nsim_required(0.005, expected_proportion = 1),
    "strictly between 0 and 1"
  )

  # A target this fine needs more replications than R can hold in an integer,
  # so it is refused instead of overflowing to NA.
  expect_error(nsim_required(1e-12), "more than R can hold in an integer")
  expect_error(nsim_required(1e-12), "2.5e\\+23")

  # A Monte Carlo standard error cannot be formed from one replication.
  expect_identical(nsim_required(0.5), 2L)
  expect_identical(nsim_required(1, expected_proportion = 0.1), 2L)

  for (p in c(0.001, 0.1, 0.5, 0.9, 0.999)) {
    n <- nsim_required(0.005, expected_proportion = p)
    expect_gte(n, 2L)
    expect_lte(sqrt(p * (1 - p) / n), 0.005)
  }
})

test_that("the achieved Monte Carlo standard error matches what nsim_required promised", {
  n <- nsim_required(0.01, expected_proportion = 0.9)
  sim <- simulate_fixed(sc_effect, n_per_arm = 161, nsim = n, seed = 77)
  achieved <- ademp_summary(sim)$mcse[1L]
  expect_lt(achieved, 0.011)
})

test_that("ademp_skeleton pre-fills the five ADEMP headings from the scenario", {
  sk <- ademp_skeleton(sc_effect)

  expect_type(sk, "character")
  for (heading in c(
    "## Aims", "## Data-generating mechanisms", "## Estimands",
    "## Methods", "## Performance measures"
  )) {
    expect_true(heading %in% sk)
  }
  body <- paste(sk, collapse = "\n")
  expect_match(body, "ercp_acute_cholangitis")
  expect_match(body, "0.3")
  expect_match(body, "0.15")
  expect_match(body, "Monte Carlo standard error")
  expect_match(body, "TODO")
  # House style bans em dashes in generated prose. Built from its UTF-8 bytes
  # so this file stays pure ASCII.
  em_dash <- rawToChar(as.raw(c(0xe2, 0x80, 0x94)))
  Encoding(em_dash) <- "UTF-8"
  expect_false(grepl(em_dash, body, fixed = TRUE))
})

test_that("ademp_skeleton reflects the endpoint direction and writes to file", {
  up <- paste(ademp_skeleton(scenario("hrs_terlipressin")), collapse = "\n")
  expect_match(up, "A higher event rate is the better outcome")

  down <- paste(ademp_skeleton(scenario("ercp_acute_cholangitis")), collapse = "\n")
  expect_match(down, "A lower event rate is the better outcome")

  path <- tempfile(fileext = ".md")
  on.exit(unlink(path), add = TRUE)
  written <- ademp_skeleton(sc_effect, file = path)
  expect_true(file.exists(path))
  expect_identical(readLines(path), written)

  expect_error(ademp_skeleton(list()), "`scenario` must be a gi_scenario")
  expect_error(ademp_skeleton(sc_effect, file = ""), "`file`")
})

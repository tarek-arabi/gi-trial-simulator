# Flagship simulation study: designing the definitive trial of ERCP timing in
# acute cholangitis. The protocol for this study is in ADEMP.md and was written
# before the study was run.
#
# Run from the repository root:
#   Rscript analysis/ercp_flagship/run_flagship.R
#
# Runtime is roughly 10 to 25 minutes on a laptop. Set GITRIALSIM_WORKERS to
# parallelise; the L'Ecuyer-CMRG streams make results independent of that choice.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

OUT <- "analysis/ercp_flagship"
RESULTS <- file.path(OUT, "results")
FIGURES <- file.path(OUT, "figures")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

WORKERS <- as.integer(Sys.getenv("GITRIALSIM_WORKERS", unset = "4"))
SEED <- 20260727
NSIM <- 10000
NSIM_BAYES <- 5000
NSIM_PROCOVA <- 1000

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

s <- scenario("ercp_acute_cholangitis")
control <- s$control_rate
say("Scenario: ", s$label, " | control ", control, " | treatment ", s$treatment_rate)

# ---------------------------------------------------------------------------
# A5 first: validation. If the engines disagree, nothing below is worth running.
# ---------------------------------------------------------------------------

say("A5: cross-engine benchmark")
bench <- benchmark_report()
write.csv(bench, file.path(RESULTS, "benchmark_rpact_vs_gsdesign.csv"), row.names = FALSE)
stopifnot(all(bench$all_agree))
say("     rpact and gsDesign agree; max abs difference ",
    signif(max(bench$max_abs_diff), 3))

# ---------------------------------------------------------------------------
# A1: required sample size across the plausible effect range
# ---------------------------------------------------------------------------

say("A1: sample size against treatment effect")

rrr_grid <- seq(0.10, 0.60, by = 0.025)
a1 <- do.call(rbind, lapply(rrr_grid, function(rrr) {
  treat <- control * (1 - rrr)
  sc <- scenario("ercp_acute_cholangitis", control_rate = control, treatment_rate = treat)
  fx <- design_fixed(sc, alpha = 0.025, power = 0.90)
  gs <- design_group_sequential(sc, alpha = 0.025, power = 0.90, k = 3)
  data.frame(
    relative_risk_reduction = rrr,
    treatment_rate = treat,
    risk_difference = control - treat,
    n_fixed = fx$n_total,
    n_gs_max = gs$n_total,
    n_gs_expected_h1 = gs$detail$expected_n_h1,
    n_gs_expected_h0 = gs$detail$expected_n_h0
  )
}))
write.csv(a1, file.path(RESULTS, "a1_sample_size_by_effect.csv"), row.names = FALSE)

observed_rrr <- 1 - s$treatment_rate / control
say("     observed effect in the anchoring trial is a ",
    round(observed_rrr * 100, 1), "% relative reduction, needing n = ",
    format(design_fixed(s)$n_total, big.mark = ","))

# Sensitivity to the control rate, which one trial of 304 patients estimates poorly.
control_grid <- seq(0.04, 0.10, by = 0.01)
a1b <- do.call(rbind, lapply(control_grid, function(cr) {
  do.call(rbind, lapply(c(0.20, 0.40, 0.60), function(rrr) {
    sc <- scenario("ercp_acute_cholangitis",
                   control_rate = cr, treatment_rate = cr * (1 - rrr))
    data.frame(
      control_rate = cr, relative_risk_reduction = rrr,
      n_fixed = design_fixed(sc, alpha = 0.025, power = 0.90)$n_total
    )
  }))
}))
write.csv(a1b, file.path(RESULTS, "a1b_sensitivity_control_rate.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# A2: design comparison at the design alternative and under the null
# ---------------------------------------------------------------------------

say("A2: design comparison")

fixed <- design_fixed(s, alpha = 0.025, power = 0.90)
gs <- design_group_sequential(s, alpha = 0.025, power = 0.90, k = 3,
                              futility = "nonbinding_obf")

say("     simulating fixed design under H1 (nsim = ", NSIM, ")")
sim_fixed_h1 <- simulate_fixed(s, n_per_arm = fixed$n_per_arm, nsim = NSIM,
                               alpha = 0.025, seed = SEED, workers = WORKERS)
say("     simulating fixed design under H0")
s_null <- scenario("ercp_acute_cholangitis",
                   control_rate = control, treatment_rate = control - 1e-9)
sim_fixed_h0 <- simulate_fixed(s_null, n_per_arm = fixed$n_per_arm, nsim = NSIM,
                               alpha = 0.025, seed = SEED + 1, workers = WORKERS)

say("     simulating group-sequential design under H1 and H0")
sim_gs_h1 <- simulate_group_sequential(gs, nsim = NSIM, seed = SEED + 2, workers = WORKERS)
sim_gs_h0 <- simulate_group_sequential(gs, nsim = NSIM, seed = SEED + 3, workers = WORKERS,
                                       rates = c(control, control))

say("     calibrating and simulating the Bayesian adaptive design (nsim = ", NSIM_BAYES, ")")
bayes <- calibrate_bayesian(s, n_max = gs$n_total, target_alpha = 0.025, looks = 3,
                            nsim = NSIM_BAYES, seed = SEED + 4, workers = WORKERS)

# The estimator is the risk difference p_treatment minus p_control, so the truth
# it is judged against is negative when urgent ERCP averts deaths. ademp_summary
# defaults to the rates the data were generated from, which is exactly right
# here, so it is left to do that rather than passing a hand-built value.
truth_rd <- s$treatment_rate - s$control_rate
ademp_fixed_h1 <- ademp_summary(sim_fixed_h1)
ademp_fixed_h0 <- ademp_summary(sim_fixed_h0)
ademp_gs_h1 <- ademp_summary(sim_gs_h1)
ademp_gs_h0 <- ademp_summary(sim_gs_h0)

for (nm in c("ademp_fixed_h1", "ademp_fixed_h0", "ademp_gs_h1", "ademp_gs_h0")) {
  df <- get(nm)
  df$scenario <- nm
  assign(nm, df)
}
ademp_all <- rbind(ademp_fixed_h1, ademp_fixed_h0, ademp_gs_h1, ademp_gs_h0)
write.csv(ademp_all, file.path(RESULTS, "a2_ademp_performance.csv"), row.names = FALSE)

grab <- function(tbl, measure) {
  row <- tbl[tbl$measure == measure, ]
  if (nrow(row) == 0) return(c(NA_real_, NA_real_))
  c(row$estimate[1], row$mcse[1])
}

rej_f_h1 <- grab(ademp_fixed_h1, "rejection_rate")
rej_f_h0 <- grab(ademp_fixed_h0, "rejection_rate")
rej_g_h1 <- grab(ademp_gs_h1, "rejection_rate")
rej_g_h0 <- grab(ademp_gs_h0, "rejection_rate")
en_g_h1 <- grab(ademp_gs_h1, "mean_sample_size")
en_g_h0 <- grab(ademp_gs_h0, "mean_sample_size")

num1 <- function(x) {
  if (is.null(x) || length(x) != 1L || !is.numeric(x)) NA_real_ else as.numeric(x)
}

a2 <- data.frame(
  design = c("Fixed", "Group sequential", "Bayesian adaptive"),
  max_n = c(num1(fixed$n_total), num1(gs$n_total), num1(bayes$n_total)),
  power = c(num1(rej_f_h1[1]), num1(rej_g_h1[1]), num1(bayes$detail$power)),
  power_mcse = c(num1(rej_f_h1[2]), num1(rej_g_h1[2]), num1(bayes$detail$power_mcse)),
  type_i_error = c(num1(rej_f_h0[1]), num1(rej_g_h0[1]), num1(bayes$detail$alpha_null)),
  type_i_mcse = c(num1(rej_f_h0[2]), num1(rej_g_h0[2]), num1(bayes$detail$alpha_mcse)),
  expected_n_h1 = c(num1(fixed$n_total), num1(en_g_h1[1]), num1(bayes$detail$expected_n_alt)),
  expected_n_h0 = c(num1(fixed$n_total), num1(en_g_h0[1]), num1(bayes$detail$expected_n_null)),
  analytic_power = c(num1(fixed$power), num1(gs$power), NA_real_),
  efficacy_threshold = c(NA_real_, NA_real_, num1(bayes$detail$efficacy_threshold)),
  stringsAsFactors = FALSE
)

write.csv(
  data.frame(
    look = seq_along(bayes$detail$look_n_total),
    n_cumulative = bayes$detail$look_n_total,
    information_rate = bayes$detail$information_rates
  ),
  file.path(RESULTS, "a2c_bayesian_looks.csv"), row.names = FALSE
)
write.csv(a2, file.path(RESULTS, "a2_design_comparison.csv"), row.names = FALSE)
print(a2, row.names = FALSE)

# The simulator must agree with rpact's analytic power, otherwise the simulated
# quantities that have no analytic counterpart cannot be trusted either.
analytic <- power_at(fixed, control, s$treatment_rate)
discrepancy <- abs(rej_f_h1[1] - analytic)
say("     simulated power ", round(rej_f_h1[1], 4), " vs rpact analytic ",
    round(analytic, 4), " (", round(discrepancy / rej_f_h1[2], 2), " MCSE)")
write.csv(
  data.frame(simulated = rej_f_h1[1], mcse = rej_f_h1[2], analytic = analytic,
             abs_difference = discrepancy, mcse_units = discrepancy / rej_f_h1[2]),
  file.path(RESULTS, "a2b_simulator_vs_rpact.csv"), row.names = FALSE
)

# ---------------------------------------------------------------------------
# A3: prognostic covariate adjustment
# ---------------------------------------------------------------------------

say("A3: prognostic covariate adjustment (nsim = ", NSIM_PROCOVA, ")")

# A plausible baseline covariate structure for acute cholangitis: age, bilirubin,
# a severity score, and the presence of organ dysfunction at presentation. The
# correlations are illustrative rather than estimated, which is stated in the
# report; replacing them with real aggregate estimates is exactly what an
# external parameter pack is for.
spec <- cohort_spec(
  covariates = list(
    covariate_spec("age", "normal", mean = 66, sd = 14),
    covariate_spec("bilirubin", "lognormal", meanlog = 1.2, sdlog = 0.6),
    covariate_spec("severity", "normal", mean = 0, sd = 1),
    covariate_spec("organ_dysfunction", "binary", prob = 0.22)
  ),
  correlation = matrix(
    c(
      1.00, 0.10, 0.25, 0.20,
      0.10, 1.00, 0.35, 0.30,
      0.25, 0.35, 1.00, 0.45,
      0.20, 0.30, 0.45, 1.00
    ),
    nrow = 4, byrow = TRUE
  )
)
coefs <- c(age = 0.030, bilirubin = 0.060, severity = 0.550, organ_dysfunction = 0.900)

# Demonstrated at a smaller per-arm size than the definitive trial because each
# replicate fits two generalised linear models. The variance reduction factor is
# a property of the covariate's prognostic strength, not of n, so it transfers.
procova_n <- 500
procova <- procova_gain(s, spec = spec, coefs = coefs, n_per_arm = procova_n,
                        nsim = NSIM_PROCOVA, train_n = 5000, seed = SEED + 5,
                        workers = WORKERS)
procova_df <- as.data.frame(procova[vapply(procova, function(z) is.numeric(z) && length(z) == 1, logical(1))])
procova_df$n_per_arm <- procova_n
procova_df$nsim <- NSIM_PROCOVA
write.csv(procova_df, file.path(RESULTS, "a3_prognostic_adjustment.csv"), row.names = FALSE)
print(procova_df, row.names = FALSE)

# ---------------------------------------------------------------------------
# A4: value of information
# ---------------------------------------------------------------------------

say("A4: value of information")

# Prior on the log relative risk, taken directly from the anchoring trial's
# reported 30-day mortality hazard ratio of 0.70 with 95% CI 0.25 to 1.93. The
# prior is placed on the log scale rather than on the risk difference so that
# every draw implies an event rate that is actually a probability.
set.seed(SEED + 6)
n_draws <- 20000
log_rr_mean <- log(0.70)
log_rr_sd <- (log(1.93) - log(0.25)) / (2 * stats::qnorm(0.975))
treatment_draws <- control * exp(stats::rnorm(n_draws, log_rr_mean, log_rr_sd))
treatment_draws <- pmin(treatment_draws, 0.99)

# evsi_trial() expects draws of p_treatment minus p_control, so a negative draw
# favours urgent ERCP. Deriving the draws this way rather than sampling the risk
# difference directly guarantees every implied arm rate is a valid probability,
# so no draw is ever clipped.
prior_draws <- treatment_draws - control
say("     prior on log RR: mean ", round(log_rr_mean, 3), ", sd ", round(log_rr_sd, 3),
    "; implied P(urgent ERCP reduces mortality) = ",
    round(mean(prior_draws < 0), 3))

# Net benefit of adopting urgent ERCP, per patient, on a scale where one death
# averted is worth 1 and one additional procedure-related adverse event costs
# 0.05. The excess adverse-event rate is the anchoring trial's observed
# difference. Standard care scores zero by construction.
# delta is the risk difference p_treatment minus p_control, so a negative draw
# means urgent ERCP averts deaths. Early ERCP is the reference and scores zero.
excess_ae <- 0.171 - 0.092
net_benefit <- function(delta) {
  c(early_ercp = 0, urgent_ercp = -delta - 0.05 * excess_ae)
}

evpi_result <- evpi(prior_draws, net_benefit)
say("     EVPI per patient: ", signif(evpi_result$evpi, 4),
    " (MCSE ", signif(evpi_result$mcse, 3), "); probability urgent ERCP is optimal ",
    signif(evpi_result$prob_optimal[["urgent_ercp"]], 3))
write.csv(
  data.frame(
    evpi = evpi_result$evpi, mcse = evpi_result$mcse,
    value_current_info = evpi_result$value_current_info,
    value_perfect_info = evpi_result$value_perfect_info,
    best_option_current = evpi_result$best_option_current,
    prob_urgent_optimal = evpi_result$prob_optimal[["urgent_ercp"]]
  ),
  file.path(RESULTS, "a4a_evpi.csv"), row.names = FALSE
)

n_grid <- c(100, 250, 500, 750, 1000, 1500, 2000)
voi <- voi_curve(prior_draws, n_grid = n_grid, control_rate = control,
                 net_benefit_fn = net_benefit, nsim = 2000, seed = SEED + 7)
voi$evpi <- evpi_result$evpi
write.csv(voi, file.path(RESULTS, "a4_voi_curve.csv"), row.names = FALSE)
print(voi, row.names = FALSE)

best_evsi <- max(voi$evsi, na.rm = TRUE)
pop <- population_evsi(
  evsi_per_patient = best_evsi,
  incidence = 100000, horizon_years = 10, discount_rate = 0.03
)
write.csv(
  data.frame(
    evsi_per_patient = pop$per_patient,
    total = pop$total,
    effective_population = pop$effective_population,
    undiscounted_population = pop$undiscounted_population,
    incidence = pop$incidence,
    horizon_years = pop$horizon_years,
    discount_rate = pop$discount_rate
  ),
  file.path(RESULTS, "a4b_population_evsi.csv"), row.names = FALSE
)
say("     population EVSI over 10 years at 100,000 cases a year: ",
    format(round(pop$total), big.mark = ","), " deaths-equivalent")

# ---------------------------------------------------------------------------
# Figures. No titles are baked in; captions live in REPORT.md.
# ---------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE)) {
  say("Figures")
  library(ggplot2)

  long <- rbind(
    data.frame(rrr = a1$relative_risk_reduction, n = a1$n_fixed, design = "Fixed"),
    data.frame(rrr = a1$relative_risk_reduction, n = a1$n_gs_max,
               design = "Group sequential, maximum"),
    data.frame(rrr = a1$relative_risk_reduction, n = a1$n_gs_expected_h1,
               design = "Group sequential, expected")
  )
  p1 <- ggplot(long, aes(rrr * 100, n, colour = design, linetype = design)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = observed_rrr * 100, linetype = "dotted", colour = "grey40") +
    scale_y_log10(labels = function(z) format(z, big.mark = ",", scientific = FALSE)) +
    scale_colour_gi(3) +
    labs(x = "Relative risk reduction in 30-day mortality (%)",
         y = "Total sample size", colour = NULL, linetype = NULL) +
    gi_theme()
  save_figure(p1, file.path(FIGURES, "fig1_sample_size_by_effect.png"), width = 7, height = 4.5)

  p2 <- ggplot(a1b, aes(control_rate * 100, n_fixed,
                        colour = factor(relative_risk_reduction * 100))) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
    scale_y_log10(labels = function(z) format(z, big.mark = ",", scientific = FALSE)) +
    scale_colour_gi(3) +
    labs(x = "Control arm 30-day mortality (%)", y = "Total sample size",
         colour = "Relative risk\nreduction (%)") +
    gi_theme()
  save_figure(p2, file.path(FIGURES, "fig2_sensitivity_control_rate.png"), width = 7, height = 4.5)

  oc <- data.frame(
    design = rep(a2$design, 2),
    quantity = rep(c("Maximum", "Expected under H1"), each = nrow(a2)),
    n = c(a2$max_n, a2$expected_n_h1)
  )
  oc$design <- factor(oc$design, levels = a2$design)
  p3 <- ggplot(oc, aes(design, n, fill = quantity)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_y_continuous(labels = function(z) format(z, big.mark = ",", scientific = FALSE)) +
    scale_fill_gi(2) +
    labs(x = NULL, y = "Total sample size", fill = NULL) +
    gi_theme()
  save_figure(p3, file.path(FIGURES, "fig3_design_comparison.png"), width = 7, height = 4.5)

  p4 <- plot_boundaries(gs)
  save_figure(p4, file.path(FIGURES, "fig4_gs_boundaries.png"), width = 7, height = 4.5)

  p5 <- plot_evsi(voi)
  save_figure(p5, file.path(FIGURES, "fig5_evsi.png"), width = 7, height = 4.5)
} else {
  say("ggplot2 not available; figures skipped")
}

# ---------------------------------------------------------------------------

sessioninfo <- utils::capture.output(utils::sessionInfo())
writeLines(c(paste("seed:", SEED), paste("nsim:", NSIM), paste("workers:", WORKERS), "",
             sessioninfo),
           file.path(RESULTS, "session_info.txt"))

say("Done. Results in ", RESULTS, ", figures in ", FIGURES)

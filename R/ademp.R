#' ADEMP reporting of simulation studies
#'
#' Simulation results are reported following the ADEMP framework of Morris,
#' White and Crowther (2019) <doi:10.1002/sim.8086>: Aims, Data-generating
#' mechanisms, Estimands, Methods, Performance measures. The central discipline
#' the framework imposes is that no performance measure is reported without its
#' Monte Carlo standard error, so a reader can tell a real difference between
#' methods from simulation noise.
#'
#' @name ademp
NULL

# Monte Carlo standard errors follow Morris, White and Crowther (2019),
# Table 6. Each returns the standard error of the corresponding estimate over
# repeated runs of the whole simulation study.
.mcse_proportion <- function(p, nsim) sqrt(p * (1 - p) / nsim)
.mcse_mean <- function(x) stats::sd(x) / sqrt(length(x))
.mcse_emp_se <- function(emp_se, nsim) emp_se / sqrt(2 * (nsim - 1))

# Degenerate replicates, those with no defined score statistic or a zero-width
# interval, are counted by the simulator. Older objects carry only the
# per-replicate flag; objects carrying neither are refused rather than reported
# as having none, because "none" would be a claim this function cannot support.
.degenerate_counts <- function(sim) {
  counts <- sim$degenerate
  if (is.list(counts) && all(c("score_undefined", "interval_undefined", "any") %in% names(counts))) {
    return(counts)
  }
  flag <- sim$results$degenerate
  if (is.logical(flag)) {
    return(list(
      score_undefined = NA_integer_,
      interval_undefined = NA_integer_,
      any = sum(flag)
    ))
  }
  stop(
    "`sim` carries no record of degenerate replications; it predates ",
    "degeneracy accounting and must be re-run before it can be summarised.",
    call. = FALSE
  )
}

#' Performance measures for a simulation, with Monte Carlo standard errors
#'
#' Summarises a `gi_simulation` into the performance measures an ADEMP-compliant
#' report requires. Every measure carries a Monte Carlo standard error computed
#' with the formulas of Morris, White and Crowther (2019)
#' <doi:10.1002/sim.8086>, Table 6:
#' proportions use `sqrt(p (1 - p) / nsim)`; the mean and hence bias use the
#' empirical standard deviation divided by `sqrt(nsim)`; the empirical standard
#' error uses `empSE / sqrt(2 (nsim - 1))`; the mean squared error uses the
#' standard deviation of the squared errors divided by `sqrt(nsim)`.
#'
#' The estimand is the risk difference, treatment event rate minus control event
#' rate, on its natural scale. For a group-sequential simulation the estimate is
#' the one available at the stopping analysis, so the bias row reports the
#' well-known bias of an estimate taken at a boundary crossing rather than an
#' error in the simulation.
#'
#' A `degenerate_rate` row reports the replications whose score statistic or
#' confidence interval was undefined because an arm had no events or every
#' patient had one. Those replications are carried through the other measures
#' as non-rejections with a zero-width interval, which is what a real analysis
#' would face, so the row exists to make their weight visible: a rejection rate
#' resting on a large degenerate fraction is describing an infeasible trial, not
#' a design.
#'
#' @param sim A `gi_simulation`, as returned by [simulate_fixed()] or
#'   [simulate_group_sequential()].
#' @param truth True risk difference the estimator is judged against. Defaults
#'   to the difference in the event rates the data were generated from.
#' @return A data frame with one row per performance measure and columns
#'   `measure`, `estimate`, `mcse` and `definition`. Rows are
#'   `rejection_rate`, `bias`, `empirical_se`, `mse`, `coverage`,
#'   `mean_sample_size` and `degenerate_rate`, plus `early_stopping_rate` for a
#'   group-sequential simulation.
#' @references Morris TP, White IR, Crowther MJ. Using simulation studies to
#'   evaluate statistical methods. Statistics in Medicine. 2019;38(11):2074-2102.
#'   \doi{10.1002/sim.8086}
#' @seealso [nsim_required()], [ademp_skeleton()]
#' @examples
#' sc <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
#' sim <- simulate_fixed(sc, n_per_arm = 161, nsim = 500, seed = 1)
#' ademp_summary(sim)
#' @export
ademp_summary <- function(sim, truth = NULL) {
  if (!inherits(sim, "gi_simulation")) {
    stop("`sim` must be a gi_simulation, as returned by simulate_fixed().", call. = FALSE)
  }
  if (is.null(truth)) {
    truth <- unname(sim$rates[["treatment"]] - sim$rates[["control"]])
  }
  if (!is.numeric(truth) || length(truth) != 1L || is.na(truth)) {
    stop("`truth` must be a single number, the true risk difference.", call. = FALSE)
  }

  res <- sim$results
  nsim <- sim$nsim
  if (nsim < 2L) {
    stop("`sim` must hold at least 2 replications for Monte Carlo standard errors.", call. = FALSE)
  }

  rejection <- mean(res$reject)
  theta <- res$risk_difference
  error <- theta - truth
  emp_se <- stats::sd(theta)
  mse <- mean(error^2)
  covered <- res$ci_lower <= truth & truth <= res$ci_upper
  coverage <- mean(covered)

  degenerate <- .degenerate_counts(sim)
  degenerate_rate <- degenerate$any / nsim
  degenerate_split <- if (is.na(degenerate$score_undefined)) {
    "the split between the two is not recorded on this simulation"
  } else {
    paste0(
      degenerate$score_undefined,
      " with an undefined score statistic, scored as z = 0 and so as non-rejections; ",
      degenerate$interval_undefined, " with a zero-width confidence interval"
    )
  }
  degenerate_definition <- paste0(
    "Proportion of replications with no defined test statistic or a zero-width ",
    "interval, because an arm had no events or every patient had one: ",
    degenerate$any, " of ", nsim, " replications (", degenerate_split,
    "). These replications are included in every measure above."
  )

  rejection_meaning <- if (isTRUE(all.equal(truth, 0))) {
    "type I error, since the data are generated under the null"
  } else {
    "power, since the data are generated under a non-null effect"
  }

  out <- data.frame(
    measure = c(
      "rejection_rate", "bias", "empirical_se", "mse", "coverage",
      "mean_sample_size", "degenerate_rate"
    ),
    estimate = c(
      rejection,
      mean(theta) - truth,
      emp_se,
      mse,
      coverage,
      mean(res$n_total),
      degenerate_rate
    ),
    mcse = c(
      .mcse_proportion(rejection, nsim),
      .mcse_mean(theta),
      .mcse_emp_se(emp_se, nsim),
      .mcse_mean(error^2),
      .mcse_proportion(coverage, nsim),
      .mcse_mean(res$n_total),
      .mcse_proportion(degenerate_rate, nsim)
    ),
    definition = c(
      paste0(
        "Proportion of replications rejecting the null at one-sided alpha ",
        format(sim$alpha), ". Read as ", rejection_meaning, "."
      ),
      paste0(
        "Mean estimated risk difference minus the true value of ",
        format(round(truth, 5)), "."
      ),
      "Standard deviation of the estimated risk difference across replications.",
      "Mean squared difference between the estimated and true risk difference.",
      "Proportion of replications whose nominal 95 percent Wald interval for the risk difference contains the true value.",
      "Mean total sample size actually used, both arms combined. Constant by construction for a fixed design.",
      degenerate_definition
    ),
    stringsAsFactors = FALSE
  )

  if (identical(sim$design_type, "group_sequential")) {
    stop_early <- mean(res$stop_stage < sim$boundaries$kmax)
    out <- rbind(out, data.frame(
      measure = "early_stopping_rate",
      estimate = stop_early,
      mcse = .mcse_proportion(stop_early, nsim),
      definition = "Proportion of replications stopping before the final analysis, for either efficacy or futility.",
      stringsAsFactors = FALSE
    ))
  }

  rownames(out) <- NULL
  out
}

#' Replications needed for a target Monte Carlo standard error
#'
#' Inverts the Monte Carlo standard error of a proportion,
#' `sqrt(p (1 - p) / nsim)`, to give the number of replications needed to
#' estimate a rejection rate or coverage to a stated precision. This is the
#' justification for the repetition count that Morris, White and Crowther (2019)
#' <doi:10.1002/sim.8086> require a simulation study to report.
#'
#' Two boundaries are handled explicitly rather than left to the arithmetic. An
#' `expected_proportion` of exactly 0 or 1 is rejected: the formula returns zero
#' replications there, which is never a usable answer, and a proportion known to
#' be degenerate needs no simulation to resolve it. A `target_mcse` so small
#' that the required count exceeds the largest representable integer is also
#' rejected, rather than returned as an overflowed or unrunnable number. The
#' result is never below 2, since a Monte Carlo standard error cannot be formed
#' from a single replication.
#'
#' @param target_mcse Target Monte Carlo standard error, on the proportion
#'   scale. For example 0.005 to resolve power to within half a percentage
#'   point.
#' @param expected_proportion The proportion being estimated, strictly between 0
#'   and 1. Defaults to 0.5, which maximises the variance and is therefore the
#'   conservative choice; supply the anticipated power to get a smaller and
#'   still sufficient number.
#' @return A single integer of at least 2, the number of replications required.
#' @references Morris TP, White IR, Crowther MJ. Using simulation studies to
#'   evaluate statistical methods. Statistics in Medicine. 2019;38(11):2074-2102.
#'   \doi{10.1002/sim.8086}
#' @examples
#' nsim_required(0.005)
#' nsim_required(0.005, expected_proportion = 0.9)
#' @export
nsim_required <- function(target_mcse, expected_proportion = 0.5) {
  if (!is.numeric(target_mcse) || length(target_mcse) != 1L || is.na(target_mcse) ||
    !is.finite(target_mcse) || target_mcse <= 0) {
    stop("`target_mcse` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(expected_proportion) || length(expected_proportion) != 1L ||
    is.na(expected_proportion) || expected_proportion < 0 || expected_proportion > 1) {
    stop("`expected_proportion` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (expected_proportion == 0 || expected_proportion == 1) {
    stop(
      "`expected_proportion` must be strictly between 0 and 1; a proportion of ",
      format(expected_proportion), " has no Monte Carlo standard error to target ",
      "and would ask for 0 replications. Use 0.5 for the conservative answer.",
      call. = FALSE
    )
  }

  required <- ceiling(expected_proportion * (1 - expected_proportion) / target_mcse^2)
  if (required > .Machine$integer.max) {
    stop(
      "A target Monte Carlo standard error of ", format(target_mcse), " at a ",
      "proportion of ", format(expected_proportion), " needs about ",
      format(required, scientific = TRUE, digits = 3), " replications, more ",
      "than R can hold in an integer. Ask for a larger `target_mcse`.",
      call. = FALSE
    )
  }
  as.integer(max(required, 2))
}

#' Markdown ADEMP skeleton for a scenario
#'
#' Writes the five ADEMP headings pre-filled from a scenario, so that the aims,
#' data-generating mechanism, estimand, methods and performance measures are
#' specified in writing before any replication is run. Placeholders marked TODO
#' are the decisions a scenario cannot make on the analyst's behalf.
#'
#' @param scenario A `gi_scenario` with a binary endpoint, as returned by
#'   [scenario()]. A continuous scenario is refused with an error naming its
#'   endpoint type, since the Bernoulli data-generating mechanism written below
#'   is specific to a binary outcome.
#' @param file Optional path to write the markdown to. When supplied the file is
#'   written and the text returned invisibly.
#' @return A character vector of markdown lines.
#' @references Morris TP, White IR, Crowther MJ. Using simulation studies to
#'   evaluate statistical methods. Statistics in Medicine. 2019;38(11):2074-2102.
#'   \doi{10.1002/sim.8086}
#' @seealso [ademp_summary()], [nsim_required()]
#' @examples
#' sk <- ademp_skeleton(scenario("ercp_acute_cholangitis"))
#' cat(head(sk, 5), sep = "\n")
#' @export
ademp_skeleton <- function(scenario, file = NULL) {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`scenario` must be a gi_scenario, as returned by scenario().", call. = FALSE)
  }
  gi_require_endpoint_type(scenario, "binary", "ademp_skeleton")
  if (!is.null(file) && (!is.character(file) || length(file) != 1L || !nzchar(file))) {
    stop("`file` must be NULL or a single non-empty path.", call. = FALSE)
  }

  defaults <- scenario$defaults %||% list()
  alpha <- defaults$alpha %||% 0.025
  power <- defaults$power %||% 0.9
  ratio <- defaults$allocation_ratio %||% 1
  rd <- scenario$treatment_rate - scenario$control_rate
  better <- if (identical(scenario$direction, "lower_is_better")) "lower" else "higher"
  reps <- nsim_required(0.005, expected_proportion = power)

  lines <- c(
    paste0("# ADEMP protocol: ", scenario$label),
    "",
    paste0("Scenario source: parameter pack `", scenario$pack_id, "` v",
      scenario$pack_version, ", endpoint `", scenario$endpoint, "`."
    ),
    paste0("Rates cited from: ", trimws(scenario$source %||% "TODO: cite the source")),
    if (isTRUE(scenario$overridden)) {
      "Rates have been overridden from the pack defaults for this analysis."
    } else {
      "Rates are the pack defaults, unmodified."
    },
    "",
    "## Aims",
    "",
    paste0(
      "TODO: state the question. A typical aim is to establish the sample size and ",
      "operating characteristics of a randomised comparison of ", scenario$treatment_arm,
      " against ", scenario$control_arm, " on ", scenario$label, ", and to establish ",
      "how those characteristics degrade when the assumed event rates are wrong."
    ),
    "",
    "## Data-generating mechanisms",
    "",
    paste0(
      "Two parallel arms with independent Bernoulli outcomes. Control (",
      scenario$control_arm, ") event rate ", format(scenario$control_rate),
      "; treatment (", scenario$treatment_arm, ") event rate ",
      format(scenario$treatment_rate), "."
    ),
    paste0("Allocation ratio ", format(ratio), ", treatment to control."),
    paste0(
      "A ", better, " event rate is the better outcome, so superiority is tested ",
      "one-sided in that direction."
    ),
    "TODO: list the sensitivity rates the grid will sweep, and any intercurrent events modelled.",
    "TODO: record the root seed and the number of independent streams used.",
    "",
    "## Estimands",
    "",
    paste0(
      "Risk difference, treatment minus control, on the natural scale. True value ",
      "under the assumed rates: ", format(round(rd, 5)), "."
    ),
    "TODO: state the intercurrent-event strategy per ICH E9(R1) that this estimand assumes.",
    "",
    "## Methods",
    "",
    paste0(
      "One-sided pooled-variance score test of two proportions at alpha ",
      format(alpha), ", targeting ", format(power), " power."
    ),
    "Sample size from rpact; group-sequential boundaries, where used, from rpact alpha spending.",
    "TODO: list every competing design being compared, for example fixed against group sequential against Bayesian adaptive.",
    "",
    "## Performance measures",
    "",
    paste0(
      "Rejection rate (power under the alternative, type I error under the null), ",
      "bias of the estimated risk difference, empirical standard error, mean squared ",
      "error, coverage of the nominal 95 percent interval, and mean realised sample size."
    ),
    "Every measure is reported with its Monte Carlo standard error per Morris, White and Crowther (2019).",
    paste0(
      "Replications: ", format(reps), " gives a Monte Carlo standard error of about 0.005 ",
      "on a rejection rate near ", format(power), " (see nsim_required())."
    ),
    ""
  )
  lines <- lines[!vapply(lines, is.null, logical(1))]

  if (!is.null(file)) {
    writeLines(lines, file)
    return(invisible(lines))
  }
  lines
}

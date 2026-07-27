#' Fixed-sample designs
#'
#' Sizing a two-arm superiority trial with a binary endpoint. Every number is
#' produced by `rpact`; this file contributes argument marshalling, the
#' direction convention, and integer rounding, nothing statistical.
#'
#' @name fixed-designs
NULL

# rpact's getSampleSizeRates places pi1 in the group whose size is the
# numerator of allocationRatioPlanned. Verified empirically: swapping pi1 and
# pi2 at allocationRatioPlanned = 2 changes nFixed (3343.185 vs 3455.004),
# so pi1 is the treatment arm and pi2 the control arm.
gi_rates <- function(scenario) {
  list(pi1 = scenario$treatment_rate, pi2 = scenario$control_rate)
}

# Superiority is defined by the endpoint's direction, not by the assumed
# rates. rpact infers directionUpper from sign(pi1 - pi2), which is the same
# thing only when the scenario actually assumes a benefit.
gi_direction_upper <- function(scenario) {
  identical(scenario$direction, "higher_is_better")
}

gi_check_scenario <- function(scenario) {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`scenario` must be a gi_scenario, as returned by scenario().", call. = FALSE)
  }
  if (!scenario$direction %in% c("lower_is_better", "higher_is_better")) {
    stop(
      "`scenario` has direction '", scenario$direction,
      "'; expected lower_is_better or higher_is_better.",
      call. = FALSE
    )
  }
  favours_treatment <- if (gi_direction_upper(scenario)) {
    scenario$treatment_rate > scenario$control_rate
  } else {
    scenario$treatment_rate < scenario$control_rate
  }
  if (!favours_treatment) {
    warning(
      "`scenario` assumes no treatment benefit on '", scenario$endpoint,
      "' (", scenario$direction, ", control ", signif(scenario$control_rate, 4),
      " vs treatment ", signif(scenario$treatment_rate, 4),
      "). The design is sized for an effect of this magnitude, but the ",
      "one-sided superiority test points the other way.",
      call. = FALSE
    )
  }
  invisible(scenario)
}

gi_resolve_defaults <- function(scenario, alpha, power, allocation_ratio) {
  d <- scenario$defaults %||% list()

  alpha <- alpha %||% d$alpha %||% 0.025
  power <- power %||% d$power %||% 0.9
  allocation_ratio <- allocation_ratio %||% d$allocation_ratio %||% 1

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
    alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number strictly between 0 and 0.5.", call. = FALSE)
  }
  if (!is.numeric(power) || length(power) != 1L || is.na(power) ||
    power <= 0 || power >= 1) {
    stop("`power` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (power <= 1 - alpha) {
    stop(
      "`power` (", power, ") must exceed 1 - alpha (", 1 - alpha,
      ") for the design to be meaningful.",
      call. = FALSE
    )
  }
  if (!is.numeric(allocation_ratio) || length(allocation_ratio) != 1L ||
    is.na(allocation_ratio) || allocation_ratio <= 0) {
    stop("`allocation_ratio` must be a single positive number.", call. = FALSE)
  }

  sided <- d$sided %||% 1
  if (!identical(as.numeric(sided), 1)) {
    stop(
      "Pack default `sided` is ", sided,
      "; gitrialsim builds one-sided superiority designs only.",
      call. = FALSE
    )
  }

  list(alpha = alpha, power = power, allocation_ratio = allocation_ratio)
}

# Per-arm sizes are rounded up so that no arm is under-recruited. With 1:1
# allocation this makes the total an even integer by construction.
gi_round_arms <- function(n_treatment_raw, n_control_raw) {
  n_treatment <- ceiling(n_treatment_raw)
  n_control <- ceiling(n_control_raw)
  list(
    n_treatment = n_treatment,
    n_control = n_control,
    n_total = n_treatment + n_control,
    n_per_arm = max(n_treatment, n_control)
  )
}

#' Size a fixed-sample two-arm superiority trial
#'
#' Delegates entirely to `rpact`: a one-stage group-sequential design supplies
#' the critical value, and [rpact::getSampleSizeRates()] supplies the sample
#' size for the two binomial rates carried by the scenario.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param alpha One-sided type I error rate. Defaults to
#'   `scenario$defaults$alpha`.
#' @param power Target power. Defaults to `scenario$defaults$power`.
#' @param allocation_ratio Planned ratio of treatment-arm to control-arm
#'   sample size. Defaults to `scenario$defaults$allocation_ratio`.
#' @return An object of class `gi_design` with `type = "fixed"`. Beyond the
#'   shared contract (`type`, `scenario`, `alpha`, `power`, `n_total`,
#'   `n_per_arm`, `engine`), `detail` holds `rpact_design`,
#'   `rpact_sample_size`, the unrounded `n_fixed`, the rounded `n_control` and
#'   `n_treatment`, `allocation_ratio`, `direction_upper` and the z-scale
#'   `critical_value`.
#' @seealso [design_group_sequential()], [power_at()]
#' @examples
#' d <- design_fixed(scenario("ercp_acute_cholangitis"))
#' d$n_total
#' print(d)
#' @export
design_fixed <- function(scenario, alpha = NULL, power = NULL,
                         allocation_ratio = NULL) {
  gi_check_scenario(scenario)
  opts <- gi_resolve_defaults(scenario, alpha, power, allocation_ratio)
  rates <- gi_rates(scenario)

  rpact_design <- rpact::getDesignGroupSequential(
    kMax = 1L,
    alpha = opts$alpha,
    beta = 1 - opts$power,
    sided = 1L
  )
  ss <- rpact::getSampleSizeRates(
    design = rpact_design,
    pi1 = rates$pi1,
    pi2 = rates$pi2,
    allocationRatioPlanned = opts$allocation_ratio
  )

  n <- gi_round_arms(ss$nFixed1, ss$nFixed2)

  structure(
    list(
      type = "fixed",
      scenario = scenario,
      alpha = opts$alpha,
      power = opts$power,
      n_total = n$n_total,
      n_per_arm = n$n_per_arm,
      engine = "rpact::getSampleSizeRates",
      detail = list(
        rpact_design = rpact_design,
        rpact_sample_size = ss,
        n_fixed = as.numeric(ss$nFixed),
        n_control = n$n_control,
        n_treatment = n$n_treatment,
        allocation_ratio = opts$allocation_ratio,
        direction_upper = gi_direction_upper(scenario),
        critical_value = as.numeric(rpact_design$criticalValues)
      )
    ),
    class = c("gi_design", "list")
  )
}

#' Analytic power of a sized design at given event rates
#'
#' Re-evaluates an already-sized design at rates that may differ from the ones
#' it was built on, which is how the Monte Carlo simulator in this package is
#' checked against an analytic reference. Computed by
#' [rpact::getPowerRates()], holding the design's maximum sample size fixed.
#'
#' @param design A `gi_design` of type `"fixed"` or `"group_sequential"`.
#' @param control_rate Control-arm event probability, in (0, 1).
#' @param treatment_rate Treatment-arm event probability, in (0, 1).
#' @return A single number: the probability of rejecting the null in the
#'   direction of treatment superiority, given the design's total sample size.
#' @seealso [design_fixed()]
#' @examples
#' d <- design_fixed(scenario("ercp_acute_cholangitis"))
#' power_at(d, control_rate = 0.0658, treatment_rate = 0.0395)
#' power_at(d, control_rate = 0.0658, treatment_rate = 0.0658)
#' @export
power_at <- function(design, control_rate, treatment_rate) {
  if (!inherits(design, "gi_design")) {
    stop("`design` must be a gi_design.", call. = FALSE)
  }
  if (is.null(design$detail$rpact_design)) {
    stop(
      "`design` of type '", design$type,
      "' carries no rpact design; power_at() supports fixed and ",
      "group_sequential designs.",
      call. = FALSE
    )
  }
  for (nm in c("control_rate", "treatment_rate")) {
    r <- get(nm)
    if (!is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0 || r >= 1) {
      stop("`", nm, "` must be a single number strictly between 0 and 1.", call. = FALSE)
    }
  }

  pw <- rpact::getPowerRates(
    design = design$detail$rpact_design,
    pi1 = treatment_rate,
    pi2 = control_rate,
    maxNumberOfSubjects = design$n_total,
    directionUpper = design$detail$direction_upper,
    allocationRatioPlanned = design$detail$allocation_ratio
  )
  as.numeric(pw$overallReject)
}

#' @param x A `gi_design`.
#' @param ... Ignored, present for S3 consistency.
#' @return `x`, invisibly.
#' @rdname design_fixed
#' @export
print.gi_design <- function(x, ...) {
  cat("<gi_design> ", x$type, "\n", sep = "")

  s <- x$scenario
  if (!is.null(s)) {
    cat("  ", s$pack_id, " / ", s$endpoint, "\n", sep = "")
    cat("  ", s$label, "\n", sep = "")
    cat(sprintf(
      "    %-30s %.4f\n    %-30s %.4f\n",
      s$control_arm, s$control_rate, s$treatment_arm, s$treatment_rate
    ))
  }

  cat(sprintf(
    "  alpha %.4g (1-sided)   power %s\n",
    x$alpha,
    if (is.null(x$power) || is.na(x$power)) "by simulation" else sprintf("%.3f", x$power)
  ))
  cat(sprintf(
    "  maximum n %s total, %s per arm\n",
    format(x$n_total, big.mark = ","), format(x$n_per_arm, big.mark = ",")
  ))

  d <- x$detail %||% list()
  if (!is.null(d$efficacy_z)) {
    tab <- data.frame(
      analysis = seq_along(d$efficacy_z),
      info = round(d$information_rates, 4),
      n_cum = d$n_cumulative,
      efficacy_z = round(d$efficacy_z, 4)
    )
    if (!is.null(d$futility_z) && any(!is.na(d$futility_z))) {
      tab$futility_z <- round(d$futility_z, 4)
    }
    cat("\n")
    print(tab, row.names = FALSE)
  }
  if (!is.null(d$expected_n_h1)) {
    cat(sprintf(
      "\n  expected n: %.0f under H1, %.0f under H0\n",
      d$expected_n_h1, d$expected_n_h0
    ))
  }
  if (!is.null(d$n_sim)) {
    cat(sprintf("  replicates: %s\n", format(d$n_sim, big.mark = ",")))
  }
  if (!is.null(d$seed)) {
    cat(sprintf("  seed: %s\n", as.character(d$seed)))
  }

  cat("  engine: ", x$engine, "\n", sep = "")
  invisible(x)
}

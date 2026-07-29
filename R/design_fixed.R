#' Fixed-sample designs
#'
#' Sizing a two-arm superiority trial with a binary or a continuous endpoint.
#' Every number is produced by `rpact`; this file contributes argument
#' marshalling, the direction convention, and integer rounding, nothing
#' statistical.
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

# getSampleSizeMeans has no per-arm mean arguments: it takes a single
# `alternative`, the assumed mean difference, and requires it to be strictly
# positive (thetaH0 defaults to 0 and rpact errors if alternative <= thetaH0).
# Sample size is symmetric in which arm is "larger", so the magnitude is all a
# sizing call needs; direction is applied later, in power_at(), the same way
# gi_direction_upper() is applied for getPowerRates rather than for
# getSampleSizeRates. Verified empirically: getSampleSizeMeans(alternative = 5,
# stDev = 10, allocationRatioPlanned = 2) returns nFixed1 = 127.38,
# nFixed2 = 63.69, ratio exactly 2, so nFixed1 is the numerator arm of
# allocationRatioPlanned exactly as pi1 is for getSampleSizeRates. Treatment is
# therefore assigned to "arm 1" here for the same reason pi1 is the treatment
# rate above: allocation_ratio is documented as treatment-to-control.
gi_means <- function(scenario) {
  list(
    effect = abs(scenario$treatment_mean - scenario$control_mean),
    sd = scenario$sd
  )
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
  if (identical(scenario$endpoint_type, "continuous")) {
    favours_treatment <- if (gi_direction_upper(scenario)) {
      scenario$treatment_mean > scenario$control_mean
    } else {
      scenario$treatment_mean < scenario$control_mean
    }
    if (!favours_treatment) {
      warning(
        "`scenario` assumes no treatment benefit on '", scenario$endpoint,
        "' (", scenario$direction, ", control ", signif(scenario$control_mean, 4),
        " vs treatment ", signif(scenario$treatment_mean, 4),
        "). The design is sized for an effect of this magnitude, but the ",
        "one-sided superiority test points the other way.",
        call. = FALSE
      )
    }
  } else {
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
  }
  invisible(scenario)
}

gi_resolve_defaults <- function(scenario, alpha, power, allocation_ratio) {
  d <- scenario$defaults %||% list()

  # A pack records the convention its source paper used. Published trials
  # overwhelmingly state a two-sided alpha, and a pack that cannot say so has to
  # silently rewrite its source, which is exactly what a validation dataset must
  # not do. Designs here are still built one-sided: for a two-arm superiority
  # test the two conventions give the same critical value and the same sample
  # size, since z(1 - 0.05/2) = z(1 - 0.025). Verified against rpact 4.2.0 for
  # both getSampleSizeRates and getSampleSizeMeans.
  #
  # An `alpha` passed by the caller is one-sided by definition (see the roxygen
  # for design_fixed), so only a pack-supplied alpha is halved.
  alpha_supplied <- !is.null(alpha)
  sided <- as.numeric(d$sided %||% 1)
  if (length(sided) != 1L || is.na(sided) || !sided %in% c(1, 2)) {
    stop(
      "Pack default `sided` is ", format(d$sided),
      "; it must be 1 or 2.",
      call. = FALSE
    )
  }

  alpha <- alpha %||% d$alpha %||% 0.025
  power <- power %||% d$power %||% 0.9
  allocation_ratio <- allocation_ratio %||% d$allocation_ratio %||% 1

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }
  if (!alpha_supplied && sided == 2) {
    alpha <- alpha / 2
  }
  if (alpha >= 0.5) {
    stop("`alpha` must be a single number strictly between 0 and 0.5.", call. = FALSE)
  }
  if (!is.numeric(power) || length(power) != 1L || is.na(power) ||
    power <= 0 || power >= 1) {
    stop("`power` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (power <= alpha) {
    stop(
      "`power` (", power, ") must exceed `alpha` (", alpha,
      ") for the design to be meaningful.",
      call. = FALSE
    )
  }
  if (!is.numeric(allocation_ratio) || length(allocation_ratio) != 1L ||
    is.na(allocation_ratio) || allocation_ratio <= 0) {
    stop("`allocation_ratio` must be a single positive number.", call. = FALSE)
  }

  list(
    alpha = alpha, power = power, allocation_ratio = allocation_ratio,
    sided_source = sided
  )
}

# Per-arm sizes are rounded up so that no arm is under-recruited. With 1:1
# allocation this makes the total an even integer by construction.
#
# `n_per_arm` is only meaningful when the two arms are the same size, so it is
# NA under unequal allocation. It used to be max(n_treatment, n_control),
# which reads as a plausible number and is silently wrong: a consumer that
# recruits n_per_arm into both arms runs a larger, balanced trial than the one
# that was sized, and overstates power by several percentage points. NA makes
# every such consumer fail loudly instead. The arm sizes themselves are always
# carried in `n_treatment` and `n_control`.
gi_round_arms <- function(n_treatment_raw, n_control_raw) {
  n_treatment <- ceiling(n_treatment_raw)
  n_control <- ceiling(n_control_raw)
  list(
    n_treatment = n_treatment,
    n_control = n_control,
    n_total = n_treatment + n_control,
    n_per_arm = if (isTRUE(n_treatment == n_control)) n_treatment else NA_real_
  )
}

#' Size a fixed-sample two-arm superiority trial
#'
#' Delegates entirely to `rpact`: a one-stage group-sequential design supplies
#' the critical value, and [rpact::getSampleSizeRates()] supplies the sample
#' size for a binary endpoint's two rates, or [rpact::getSampleSizeMeans()]
#' the sample size for a continuous endpoint's two means and common SD,
#' whichever the scenario carries.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param alpha One-sided type I error rate. Defaults to
#'   `scenario$defaults$alpha`, which is interpreted on the pack's own `sided`
#'   convention: a pack recording `sided: 2` and `alpha: 0.05`, as a published
#'   trial usually states it, yields a one-sided 0.025 and an identical design.
#'   An `alpha` supplied here is always one-sided and is never halved.
#' @param power Target power. Defaults to `scenario$defaults$power`.
#' @param allocation_ratio Planned ratio of treatment-arm to control-arm
#'   sample size. Defaults to `scenario$defaults$allocation_ratio`.
#' @return An object of class `gi_design` with `type = "fixed"`. Beyond the
#'   shared contract (`type`, `scenario`, `alpha`, `power`, `n_total`,
#'   `n_per_arm`, `engine`), `detail` holds `rpact_design`,
#'   `rpact_sample_size`, the unrounded `n_fixed`, the rounded `n_control` and
#'   `n_treatment`, `allocation_ratio`, `direction_upper` and the z-scale
#'   `critical_value`. For a continuous scenario, `detail` additionally holds
#'   `standardised_effect`, the Cohen's d equivalent
#'   `abs(treatment_mean - control_mean) / sd`, and `engine` is
#'   `"rpact::getSampleSizeMeans"` rather than `"rpact::getSampleSizeRates"`.
#'
#'   `n_total` is always the sum of the two arm sizes. `n_per_arm` is the
#'   common arm size and is defined only under 1:1 allocation; when
#'   `allocation_ratio` is anything else the arms differ and `n_per_arm` is
#'   `NA_real_`, so that code which assumes balanced arms fails rather than
#'   quietly using the wrong number. Read `detail$n_treatment` and
#'   `detail$n_control` when allocation may be unequal.
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
  is_continuous <- identical(scenario$endpoint_type, "continuous")

  rpact_design <- rpact::getDesignGroupSequential(
    kMax = 1L,
    alpha = opts$alpha,
    beta = 1 - opts$power,
    sided = 1L
  )

  if (is_continuous) {
    means <- gi_means(scenario)
    ss <- rpact::getSampleSizeMeans(
      design = rpact_design,
      alternative = means$effect,
      stDev = means$sd,
      allocationRatioPlanned = opts$allocation_ratio
    )
    engine <- "rpact::getSampleSizeMeans"
  } else {
    rates <- gi_rates(scenario)
    ss <- rpact::getSampleSizeRates(
      design = rpact_design,
      pi1 = rates$pi1,
      pi2 = rates$pi2,
      allocationRatioPlanned = opts$allocation_ratio
    )
    engine <- "rpact::getSampleSizeRates"
  }

  n <- gi_round_arms(ss$nFixed1, ss$nFixed2)

  detail <- list(
    rpact_design = rpact_design,
    rpact_sample_size = ss,
    n_fixed = as.numeric(ss$nFixed),
    n_control = n$n_control,
    n_treatment = n$n_treatment,
    allocation_ratio = opts$allocation_ratio,
    direction_upper = gi_direction_upper(scenario),
    critical_value = as.numeric(rpact_design$criticalValues),
    sided_source = opts$sided_source
  )
  if (is_continuous) {
    detail$standardised_effect <- means$effect / means$sd
  }

  structure(
    list(
      type = "fixed",
      scenario = scenario,
      alpha = opts$alpha,
      power = opts$power,
      n_total = n$n_total,
      n_per_arm = n$n_per_arm,
      engine = engine,
      detail = detail
    ),
    class = c("gi_design", "list")
  )
}

#' Analytic power of a sized design at given event rates or means
#'
#' Re-evaluates an already-sized design at rates or means that may differ from
#' the ones it was built on, which is how the Monte Carlo simulator in this
#' package is checked against an analytic reference. Computed by
#' [rpact::getPowerRates()] for a binary design, or [rpact::getPowerMeans()]
#' for a continuous one, holding the design's maximum sample size fixed.
#'
#' @param design A `gi_design` of type `"fixed"` or `"group_sequential"`.
#' @param control_rate,treatment_rate Control- and treatment-arm event
#'   probability, in (0, 1). Used when `design` was built on a binary
#'   endpoint; leave `NULL` (the default) for a continuous design and supply
#'   `control_mean`/`treatment_mean`/`sd` instead.
#' @param control_mean,treatment_mean,sd Control- and treatment-arm means, and
#'   the common within-arm standard deviation to evaluate power at. Used when
#'   `design` was built on a continuous endpoint; leave `NULL` (the default)
#'   for a binary design and supply `control_rate`/`treatment_rate` instead.
#' @return A single number, rpact's `overallReject`: the probability that the
#'   trial crosses an efficacy boundary in the direction of treatment
#'   superiority, given the design's total sample size and allocation ratio.
#'
#'   For a design built with `futility = "nonbinding_obf"` or
#'   `"binding_obf"` this is a probability of reaching efficacy *without having
#'   stopped for futility first*, because rpact applies the futility bounds
#'   when it propagates the trial forward. It is therefore not the type I error
#'   rate of the design when the futility bounds are non-binding: evaluated at
#'   equal rates, a non-binding O'Brien-Fleming futility design returns roughly
#'   0.023 rather than its 0.025 alpha, since the replicates that would have
#'   gone on to reject after crossing a futility bound have been removed. A
#'   sponsor who declines to act on a non-binding bound retains the full alpha,
#'   which this function does not report. Use `futility = "none"` to see the
#'   efficacy-only quantity.
#' @seealso [design_fixed()]
#' @examples
#' d <- design_fixed(scenario("ercp_acute_cholangitis"))
#' power_at(d, control_rate = 0.0658, treatment_rate = 0.0395)
#' power_at(d, control_rate = 0.0658, treatment_rate = 0.0658)
#' @export
power_at <- function(design, control_rate = NULL, treatment_rate = NULL,
                     control_mean = NULL, treatment_mean = NULL, sd = NULL) {
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

  endpoint_type <- design$scenario$endpoint_type %||% "binary"

  if (identical(endpoint_type, "continuous")) {
    if (!is.null(control_rate) || !is.null(treatment_rate)) {
      stop(
        "`design` has a continuous endpoint; power_at() takes `control_mean`, ",
        "`treatment_mean` and `sd` for this design, not `control_rate`/",
        "`treatment_rate`.",
        call. = FALSE
      )
    }
    if (!is.numeric(control_mean) || length(control_mean) != 1L ||
      is.na(control_mean) || !is.finite(control_mean)) {
      stop("`control_mean` must be a single finite number.", call. = FALSE)
    }
    if (!is.numeric(treatment_mean) || length(treatment_mean) != 1L ||
      is.na(treatment_mean) || !is.finite(treatment_mean)) {
      stop("`treatment_mean` must be a single finite number.", call. = FALSE)
    }
    if (!is.numeric(sd) || length(sd) != 1L || is.na(sd) || !is.finite(sd) || sd <= 0) {
      stop("`sd` must be a single positive finite number.", call. = FALSE)
    }

    pw <- rpact::getPowerMeans(
      design = design$detail$rpact_design,
      alternative = treatment_mean - control_mean,
      stDev = sd,
      directionUpper = design$detail$direction_upper,
      maxNumberOfSubjects = design$n_total,
      allocationRatioPlanned = design$detail$allocation_ratio
    )
    return(as.numeric(pw$overallReject))
  }

  if (!is.null(control_mean) || !is.null(treatment_mean) || !is.null(sd)) {
    stop(
      "`design` has a binary endpoint; power_at() takes `control_rate`/",
      "`treatment_rate` for this design, not `control_mean`/`treatment_mean`/`sd`.",
      call. = FALSE
    )
  }
  if (!is.numeric(control_rate) || length(control_rate) != 1L ||
    is.na(control_rate) || control_rate <= 0 || control_rate >= 1) {
    stop("`control_rate` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(treatment_rate) || length(treatment_rate) != 1L ||
    is.na(treatment_rate) || treatment_rate <= 0 || treatment_rate >= 1) {
    stop("`treatment_rate` must be a single number strictly between 0 and 1.", call. = FALSE)
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
  d <- x$detail %||% list()
  if (!is.null(s)) {
    cat("  ", s$pack_id, " / ", s$endpoint, "\n", sep = "")
    cat("  ", s$label, "\n", sep = "")
    if (identical(s$endpoint_type, "continuous")) {
      cat(sprintf(
        "    %-30s %.4f\n    %-30s %.4f\n",
        s$control_arm, s$control_mean, s$treatment_arm, s$treatment_mean
      ))
      cat(sprintf("    %-30s %.4f\n", "common SD", s$sd))
      if (!is.null(d$standardised_effect)) {
        cat(sprintf(
          "    %-30s %.4f\n", "standardised effect (d)", d$standardised_effect
        ))
      }
    } else {
      cat(sprintf(
        "    %-30s %.4f\n    %-30s %.4f\n",
        s$control_arm, s$control_rate, s$treatment_arm, s$treatment_rate
      ))
    }
  }

  cat(sprintf(
    "  alpha %.4g (1-sided)   power %s\n",
    x$alpha,
    if (is.null(x$power) || is.na(x$power)) "by simulation" else sprintf("%.3f", x$power)
  ))
  n_per_arm <- x$n_per_arm
  balanced <- length(n_per_arm) == 1L && !is.na(n_per_arm)
  if (balanced) {
    cat(sprintf(
      "  maximum n %s total, %s per arm\n",
      format(x$n_total, big.mark = ","), format(n_per_arm, big.mark = ",")
    ))
  } else if (!is.null(d$n_treatment) && !is.null(d$n_control)) {
    cat(sprintf(
      "  maximum n %s total, %s treatment / %s control\n",
      format(x$n_total, big.mark = ","),
      format(d$n_treatment, big.mark = ","),
      format(d$n_control, big.mark = ",")
    ))
  } else {
    cat(sprintf(
      "  maximum n %s total, arms unequal\n",
      format(x$n_total, big.mark = ",")
    ))
  }

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
  n_sim <- d$n_sim %||% d$nsim
  if (!is.null(n_sim)) {
    cat(sprintf("  replicates: %s\n", format(n_sim, big.mark = ",")))
  }
  if (!is.null(d$seed)) {
    cat(sprintf("  seed: %s\n", as.character(d$seed)))
  }

  cat("  engine: ", x$engine, "\n", sep = "")
  invisible(x)
}

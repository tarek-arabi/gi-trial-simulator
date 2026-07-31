#' Reconciling published sample sizes against the engine
#'
#' A published trial states the inputs to its sample size calculation and the
#' number it arrived at. Feeding those inputs to [design_fixed()] and comparing
#' is not, on its own, informative: trials are sized with different
#' approximations, and rpact implements one of them. A bare disagreement
#' therefore says nothing about whether the engine is right or the trial's
#' arithmetic was.
#'
#' These functions close that gap by computing the sample size under each
#' standard formulation and reporting which one the published number
#' corresponds to. A discrepancy is then attributable rather than merely
#' observed.
#'
#' The reference formulas here are **diagnostic only**. Every design this
#' package produces still comes from rpact; nothing in this file is used to size
#' a trial. Each formula is a published historical one, cited at its definition
#' and pinned by a test to a worked example. One of them, the pooled-variance
#' normal approximation, is available independently from
#' [stats::power.prop.test()], and the three-way agreement between it, rpact and
#' the closed form is itself checked in the test suite.
#'
#' @name reconciliation
NULL

# The methods a published sample size might have been computed under. The
# pooled-variance normal approximation is what rpact::getSampleSizeRates
# implements, verified to all printed digits against both the closed form and
# stats::power.prop.test.
#
# Continuity-corrected: Casagrande JT, Pike MC, Smith PG. An improved
# approximate formula for calculating sample sizes for comparing two binomial
# distributions. Biometrics 1978;34(3):483-486. Also given as Fleiss's formula.
# This is the standard analytic approximation to Fisher's exact test and is the
# default in several commercial sample size packages, so a trial reporting
# "Fisher's exact" typically lands here rather than on the uncorrected formula.
#
# Arcsine: Fleiss JL, Levin B, Paik MC. Statistical Methods for Rates and
# Proportions, 3rd ed., Wiley 2003, section 4.
GI_RECONCILE_METHODS <- c(
  "rpact", "pooled", "power.prop.test", "unpooled", "continuity", "arcsine"
)

#' Sample size per arm under one named formulation
#'
#' @param p_control,p_treat Assumed control- and treatment-arm event
#'   probabilities, each strictly between 0 and 1.
#' @param alpha_1sided One-sided type I error rate. A trial reporting a
#'   two-sided alpha of 0.05 is reproduced with `alpha_1sided = 0.025`.
#' @param power Target power.
#' @param method One of `"rpact"`, `"pooled"`, `"power.prop.test"`,
#'   `"unpooled"`, `"continuity"` or `"arcsine"`.
#' @return A single number, the unrounded sample size per arm under 1:1
#'   allocation.
#' @seealso [reconcile_sample_size()]
#' @examples
#' # Elmunzer 2012 sized 10% versus 5% at 80% power on Fisher's exact test
#' gi_n_reference(0.10, 0.05, 0.025, 0.80, "rpact")
#' gi_n_reference(0.10, 0.05, 0.025, 0.80, "continuity")
#' @export
gi_n_reference <- function(p_control, p_treat, alpha_1sided, power,
                           method = "rpact") {
  gi_check_rate(p_control, "p_control")
  gi_check_rate(p_treat, "p_treat")
  if (!is.numeric(alpha_1sided) || length(alpha_1sided) != 1L ||
    is.na(alpha_1sided) || alpha_1sided <= 0 || alpha_1sided >= 0.5) {
    stop("`alpha_1sided` must be a single number strictly between 0 and 0.5.",
      call. = FALSE
    )
  }
  if (!is.numeric(power) || length(power) != 1L || is.na(power) ||
    power <= 0 || power >= 1) {
    stop("`power` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  method <- match.arg(method, GI_RECONCILE_METHODS)
  if (isTRUE(all.equal(p_control, p_treat))) {
    stop("`p_control` and `p_treat` must differ; no finite sample size otherwise.",
      call. = FALSE
    )
  }

  za <- stats::qnorm(1 - alpha_1sided)
  zb <- stats::qnorm(power)
  delta <- abs(p_control - p_treat)
  pbar <- (p_control + p_treat) / 2
  var_sep <- p_control * (1 - p_control) + p_treat * (1 - p_treat)

  pooled <- (za * sqrt(2 * pbar * (1 - pbar)) + zb * sqrt(var_sep))^2 / delta^2

  switch(method,
    rpact = {
      design <- rpact::getDesignGroupSequential(
        kMax = 1L, alpha = alpha_1sided, beta = 1 - power, sided = 1L
      )
      ss <- rpact::getSampleSizeRates(
        design = design, pi1 = p_treat, pi2 = p_control,
        allocationRatioPlanned = 1
      )
      as.numeric(ss$nFixed1)
    },
    pooled = pooled,
    power.prop.test = stats::power.prop.test(
      p1 = p_control, p2 = p_treat,
      sig.level = 2 * alpha_1sided, power = power
    )$n,
    unpooled = (za + zb)^2 * var_sep / delta^2,
    continuity = pooled * (1 + sqrt(1 + 4 / (pooled * delta)))^2 / 4,
    arcsine = (za + zb)^2 /
      (2 * (asin(sqrt(p_treat)) - asin(sqrt(p_control)))^2)
  )
}

gi_check_rate <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", nm, "` must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }
  invisible(x)
}

# Time-to-event reference methods. For a survival design the quantity a formula
# returns is the required number of EVENTS; patients follow only after accrual
# and follow-up assumptions. So this is a separate reconciliation, against a
# separate published quantity, and the two strata are not pooled.
#
# Schoenfeld D. The asymptotic properties of nonparametric tests for comparing
# survival distributions. Biometrika 1981;68(1):316-319. Verified to equal
# rpact::getSampleSizeSurvival()$eventsFixed to all printed digits at hazard
# ratios 0.5 through 0.8 and at 80% and 90% power. rpact implements this one.
#
# Freedman LS. Tables of the number of patients required in clinical trials
# using the logrank test. Stat Med 1982;1(2):121-129. Verified to run above
# Schoenfeld by 8.1% at HR 0.5, falling to 0.8% at HR 0.8.
GI_EVENT_METHODS <- c("rpact", "schoenfeld", "freedman")

#' Required events under one named survival formulation
#'
#' @param hazard_ratio Assumed hazard ratio, strictly positive and not 1.
#' @param alpha_1sided One-sided type I error rate.
#' @param power Target power.
#' @param allocation_ratio Treatment-to-control allocation. Only 1 has been
#'   verified against rpact; other values use the same published formulas with
#'   their allocation term but are not pinned.
#' @param method One of `"rpact"`, `"schoenfeld"` or `"freedman"`.
#' @return A single number, the required number of events.
#' @seealso [reconcile_events()]
#' @examples
#' gi_events_reference(0.7, 0.025, 0.8, method = "schoenfeld")
#' gi_events_reference(0.7, 0.025, 0.8, method = "freedman")
#' @export
gi_events_reference <- function(hazard_ratio, alpha_1sided, power,
                                allocation_ratio = 1, method = "rpact") {
  if (!is.numeric(hazard_ratio) || length(hazard_ratio) != 1L ||
    is.na(hazard_ratio) || hazard_ratio <= 0) {
    stop("`hazard_ratio` must be a single positive number.", call. = FALSE)
  }
  if (isTRUE(all.equal(hazard_ratio, 1))) {
    stop("`hazard_ratio` of 1 implies no effect and no finite event count.",
      call. = FALSE
    )
  }
  if (!is.numeric(alpha_1sided) || length(alpha_1sided) != 1L ||
    is.na(alpha_1sided) || alpha_1sided <= 0 || alpha_1sided >= 0.5) {
    stop("`alpha_1sided` must be a single number strictly between 0 and 0.5.",
      call. = FALSE
    )
  }
  if (!is.numeric(power) || length(power) != 1L || is.na(power) ||
    power <= 0 || power >= 1) {
    stop("`power` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(allocation_ratio) || length(allocation_ratio) != 1L ||
    is.na(allocation_ratio) || allocation_ratio <= 0) {
    stop("`allocation_ratio` must be a single positive number.", call. = FALSE)
  }
  method <- match.arg(method, GI_EVENT_METHODS)

  za <- stats::qnorm(1 - alpha_1sided)
  zb <- stats::qnorm(power)
  r <- allocation_ratio

  switch(method,
    rpact = {
      design <- rpact::getDesignGroupSequential(
        kMax = 1L, alpha = alpha_1sided, beta = 1 - power, sided = 1L
      )
      ss <- rpact::getSampleSizeSurvival(
        design = design, hazardRatio = hazard_ratio,
        allocationRatioPlanned = r
      )
      as.numeric(ss$eventsFixed)
    },
    schoenfeld = (za + zb)^2 * (1 + r)^2 / (r * log(hazard_ratio)^2),
    freedman = (za + zb)^2 * (1 + r * hazard_ratio)^2 /
      (r * (1 - hazard_ratio)^2)
  )
}

#' Reconcile a published required-events figure against the survival formulas
#'
#' The time-to-event counterpart of [reconcile_sample_size()]. The published
#' quantity compared against is the required number of **events**, because that
#' is what a survival sample size formula returns.
#'
#' @param published_events Required events stated in the paper.
#' @param hazard_ratio,alpha_1sided,power,allocation_ratio Design inputs as the
#'   paper states them; halve a reported two-sided alpha.
#' @param tolerance Relative agreement within which a method counts as
#'   reproducing the published figure. Defaults to 0.02.
#' @return A data frame with one row per method, ordered by absolute relative
#'   difference, carrying `attribution` and `reproduced` as attributes.
#' @seealso [gi_events_reference()], [reconcile_sample_size()]
#' @examples
#' r <- reconcile_events(250, hazard_ratio = 0.7, alpha_1sided = 0.025, power = 0.8)
#' attr(r, "attribution")
#' @export
reconcile_events <- function(published_events, hazard_ratio, alpha_1sided,
                             power, allocation_ratio = 1, tolerance = 0.02) {
  if (!is.numeric(published_events) || length(published_events) != 1L ||
    is.na(published_events) || published_events <= 0) {
    stop("`published_events` must be a single positive number.", call. = FALSE)
  }
  n <- vapply(
    GI_EVENT_METHODS,
    function(m) {
      gi_events_reference(hazard_ratio, alpha_1sided, power, allocation_ratio, m)
    },
    numeric(1)
  )
  pct <- n / published_events - 1
  out <- data.frame(
    method = GI_EVENT_METHODS,
    events = as.numeric(n),
    pct_vs_published = as.numeric(pct),
    agrees = abs(as.numeric(pct)) <= tolerance,
    stringsAsFactors = FALSE
  )
  out <- out[order(abs(out$pct_vs_published)), ]
  row.names(out) <- NULL
  attr(out, "attribution") <- if (out$agrees[1]) out$method[1] else NA_character_
  attr(out, "reproduced") <- out$agrees[1]
  out
}

#' Reconcile a published sample size against every standard formulation
#'
#' Computes the per-arm sample size under each method in
#' [gi_n_reference()] and ranks them by closeness to the number the trial
#' published, so that a discrepancy between the engine and the publication can
#' be attributed to a named computational choice rather than left unexplained.
#'
#' The attribution is descriptive. A published number matching the
#' continuity-corrected formula is evidence about how that trial was sized, not
#' evidence that either the trial or the engine is wrong.
#'
#' @param published_n_per_arm The sample size per arm stated in the paper,
#'   before any inflation for dropout.
#' @param p_control,p_treat The assumed event probabilities stated in the paper.
#' @param alpha_1sided One-sided type I error rate; halve a reported two-sided
#'   alpha.
#' @param power Target power stated in the paper.
#' @param tolerance Relative agreement, as a proportion, within which a method
#'   counts as reproducing the published number. Defaults to 0.02.
#' @return A data frame with one row per method, ordered by absolute relative
#'   difference: `method`, `n_per_arm`, `pct_vs_published` and `agrees`. The
#'   attributed method, the achieved power at the published sample size, and
#'   whether any method agreed are attached as the attributes `attribution`,
#'   `power_at_published` and `reproduced`.
#' @seealso [gi_n_reference()], [benchmark_report()]
#' @examples
#' # Elmunzer 2012: published 474 per arm, reported as Fisher's exact
#' r <- reconcile_sample_size(474, 0.10, 0.05, 0.025, 0.80)
#' attr(r, "attribution")
#' @export
reconcile_sample_size <- function(published_n_per_arm, p_control, p_treat,
                                  alpha_1sided, power, tolerance = 0.02) {
  if (!is.numeric(published_n_per_arm) || length(published_n_per_arm) != 1L ||
    is.na(published_n_per_arm) || published_n_per_arm <= 0) {
    stop("`published_n_per_arm` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
    tolerance <= 0) {
    stop("`tolerance` must be a single positive number.", call. = FALSE)
  }

  n <- vapply(
    GI_RECONCILE_METHODS,
    function(m) gi_n_reference(p_control, p_treat, alpha_1sided, power, m),
    numeric(1)
  )
  pct <- n / published_n_per_arm - 1

  out <- data.frame(
    method = GI_RECONCILE_METHODS,
    n_per_arm = as.numeric(n),
    pct_vs_published = as.numeric(pct),
    agrees = abs(as.numeric(pct)) <= tolerance,
    stringsAsFactors = FALSE
  )
  out <- out[order(abs(out$pct_vs_published)), ]
  row.names(out) <- NULL

  # Power actually delivered by the number the trial published, computed by
  # rpact at the trial's own assumed rates. This is a design-stage quantity: it
  # uses the assumed effect, never the observed one, so it is not observed power.
  design <- rpact::getDesignGroupSequential(
    kMax = 1L, alpha = alpha_1sided, beta = 1 - power, sided = 1L
  )
  achieved <- rpact::getPowerRates(
    design = design, pi1 = p_treat, pi2 = p_control,
    maxNumberOfSubjects = 2 * published_n_per_arm,
    directionUpper = p_treat > p_control, allocationRatioPlanned = 1
  )

  attr(out, "attribution") <- if (out$agrees[1]) out$method[1] else NA_character_
  attr(out, "reproduced") <- out$agrees[1]
  attr(out, "power_at_published") <- as.numeric(achieved$overallReject)
  attr(out, "power_stated") <- power
  out
}

#' Reconcile a set of trials and reduce to one row each
#'
#' The set-level companion to [reconcile_sample_size()], in the same shape as
#' [benchmark_report()]: one row per trial, suitable for a results table.
#'
#' @param trials A data frame, or a list of lists, with columns `label`,
#'   `published_n_per_arm`, `p_control`, `p_treat`, `alpha_1sided` and `power`.
#'   An optional `stated_method` column records what the paper said it used, for
#'   comparison against what it evidently used.
#' @param tolerance Passed to [reconcile_sample_size()].
#' @return A data frame with `label`, `published_n_per_arm`, `rpact_n_per_arm`,
#'   `pct_rpact_vs_published`, `attribution`, `stated_method`, `reproduced`,
#'   `power_stated` and `power_at_published`.
#' @seealso [reconcile_sample_size()]
#' @export
reconcile_report <- function(trials, tolerance = 0.02) {
  trials <- as.data.frame(trials, stringsAsFactors = FALSE)
  required <- c(
    "label", "published_n_per_arm", "p_control", "p_treat",
    "alpha_1sided", "power"
  )
  missing <- setdiff(required, names(trials))
  if (length(missing)) {
    stop("`trials` is missing column(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(trials)) {
    stop("`trials` has no rows.", call. = FALSE)
  }

  rows <- lapply(seq_len(nrow(trials)), function(i) {
    tr <- trials[i, ]
    r <- reconcile_sample_size(
      published_n_per_arm = tr$published_n_per_arm,
      p_control = tr$p_control, p_treat = tr$p_treat,
      alpha_1sided = tr$alpha_1sided, power = tr$power,
      tolerance = tolerance
    )
    rp <- r[r$method == "rpact", ]
    data.frame(
      label = as.character(tr$label),
      published_n_per_arm = as.numeric(tr$published_n_per_arm),
      rpact_n_per_arm = rp$n_per_arm,
      pct_rpact_vs_published = rp$pct_vs_published,
      attribution = attr(r, "attribution"),
      stated_method = if (is.null(tr$stated_method)) {
        NA_character_
      } else {
        as.character(tr$stated_method)
      },
      reproduced = attr(r, "reproduced"),
      power_stated = attr(r, "power_stated"),
      power_at_published = attr(r, "power_at_published"),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

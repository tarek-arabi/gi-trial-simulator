#' Group-sequential designs
#'
#' Boundaries, information rates and expected sample sizes all come from
#' `rpact`. Futility is requested through rpact's beta-spending machinery
#' rather than by supplying bounds of our own.
#'
#' @name group-sequential-designs
NULL

# rpact reports -6 in futilityBounds when no futility boundary was requested.
GI_RPACT_NO_FUTILITY <- -6

gi_check_type_of_design <- function(type_of_design) {
  allowed <- c("asOF", "asP", "OF", "P")
  if (!is.character(type_of_design) || length(type_of_design) != 1L ||
    !type_of_design %in% allowed) {
    stop(
      "`type_of_design` must be one of ", paste(allowed, collapse = ", "),
      "; got '", paste(type_of_design, collapse = ", "), "'.",
      call. = FALSE
    )
  }
  type_of_design
}

gi_check_information_rates <- function(information_rates, k) {
  if (is.null(information_rates)) {
    return(NA_real_)
  }
  if (!is.numeric(information_rates) || length(information_rates) != k) {
    stop(
      "`information_rates` must be a numeric vector of length k (", k, "); got length ",
      length(information_rates), ".",
      call. = FALSE
    )
  }
  if (anyNA(information_rates) || any(information_rates <= 0) ||
    any(information_rates > 1)) {
    stop("`information_rates` must lie in (0, 1].", call. = FALSE)
  }
  if (any(diff(information_rates) <= 0)) {
    stop("`information_rates` must be strictly increasing.", call. = FALSE)
  }
  if (!isTRUE(all.equal(information_rates[k], 1))) {
    stop("`information_rates` must end at 1.", call. = FALSE)
  }
  information_rates
}

#' Size a group-sequential two-arm superiority trial
#'
#' Builds the boundaries with [rpact::getDesignGroupSequential()] and the
#' sample size with [rpact::getSampleSizeRates()]. Nothing about the
#' boundaries, the alpha spending or the expected sample sizes is computed
#' here.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param alpha One-sided type I error rate. Defaults to
#'   `scenario$defaults$alpha`.
#' @param power Target power. Defaults to `scenario$defaults$power`.
#' @param k Number of analyses, including the final one. A single integer of
#'   at least 2.
#' @param type_of_design One of `"asOF"` (O'Brien-Fleming-type alpha
#'   spending, the default), `"asP"` (Pocock-type alpha spending), `"OF"`
#'   (classical O'Brien-Fleming) or `"P"` (classical Pocock).
#' @param futility `"none"` for efficacy-only monitoring, `"nonbinding_obf"`
#'   or `"binding_obf"` for O'Brien-Fleming-type beta spending
#'   (rpact's `typeBetaSpending = "bsOF"`) with the corresponding
#'   `bindingFutility` setting.
#' @param information_rates Optional numeric vector of length `k`, strictly
#'   increasing and ending at 1. Defaults to equally spaced analyses.
#' @param allocation_ratio Planned ratio of treatment-arm to control-arm
#'   sample size. Defaults to `scenario$defaults$allocation_ratio`.
#' @return An object of class `gi_design` with `type = "group_sequential"`.
#'   Beyond the shared contract, `detail` holds `efficacy_z`, `futility_z`,
#'   `information_rates`, `n_cumulative`, `cumulative_alpha_spent`,
#'   `expected_n_h0`, `expected_n_h1`, `expected_n_h01`, `k`,
#'   `type_of_design`, `futility`, `binding_futility`, the unrounded
#'   `n_fixed`, and the underlying `rpact_design` and `rpact_sample_size`.
#' @seealso [gs_boundaries()], [design_fixed()]
#' @examples
#' d <- design_group_sequential(scenario("ercp_acute_cholangitis"), k = 3)
#' d$n_total
#' gs_boundaries(d)
#'
#' design_group_sequential(
#'   scenario("ercp_acute_cholangitis"),
#'   k = 2, futility = "nonbinding_obf"
#' )
#' @export
design_group_sequential <- function(scenario, alpha = NULL, power = NULL, k = 3,
                                    type_of_design = "asOF",
                                    futility = c("none", "nonbinding_obf", "binding_obf"),
                                    information_rates = NULL,
                                    allocation_ratio = NULL) {
  gi_check_scenario(scenario)
  opts <- gi_resolve_defaults(scenario, alpha, power, allocation_ratio)
  type_of_design <- gi_check_type_of_design(type_of_design)
  futility <- match.arg(futility)

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != round(k) || k < 2) {
    stop("`k` must be a single integer of at least 2.", call. = FALSE)
  }
  k <- as.integer(k)
  rates_arg <- gi_check_information_rates(information_rates, k)

  binding <- switch(futility,
    none = NA,
    nonbinding_obf = FALSE,
    binding_obf = TRUE
  )
  beta_spending <- if (identical(futility, "none")) "none" else "bsOF"

  rpact_design <- rpact::getDesignGroupSequential(
    kMax = k,
    alpha = opts$alpha,
    beta = 1 - opts$power,
    sided = 1L,
    typeOfDesign = type_of_design,
    informationRates = rates_arg,
    typeBetaSpending = beta_spending,
    bindingFutility = binding
  )
  ss <- rpact::getSampleSizeRates(
    design = rpact_design,
    pi1 = scenario$treatment_rate,
    pi2 = scenario$control_rate,
    allocationRatioPlanned = opts$allocation_ratio
  )

  n <- gi_round_arms(
    as.numeric(ss$maxNumberOfSubjects1),
    as.numeric(ss$maxNumberOfSubjects2)
  )

  info <- as.numeric(rpact_design$informationRates)
  n_cumulative <- ceiling(info * n$n_total)
  n_cumulative[k] <- n$n_total

  futility_z <- rep(NA_real_, k)
  if (!identical(futility, "none")) {
    bounds <- as.numeric(rpact_design$futilityBounds)
    bounds[bounds <= GI_RPACT_NO_FUTILITY] <- NA_real_
    futility_z[seq_len(k - 1L)] <- bounds
  }

  structure(
    list(
      type = "group_sequential",
      scenario = scenario,
      alpha = opts$alpha,
      power = opts$power,
      n_total = n$n_total,
      n_per_arm = n$n_per_arm,
      engine = "rpact::getSampleSizeRates",
      detail = list(
        rpact_design = rpact_design,
        rpact_sample_size = ss,
        k = k,
        type_of_design = type_of_design,
        futility = futility,
        binding_futility = isTRUE(binding),
        information_rates = info,
        efficacy_z = as.numeric(rpact_design$criticalValues),
        futility_z = futility_z,
        cumulative_alpha_spent = as.numeric(rpact_design$alphaSpent),
        n_cumulative = n_cumulative,
        n_fixed = as.numeric(ss$nFixed),
        n_control = n$n_control,
        n_treatment = n$n_treatment,
        allocation_ratio = opts$allocation_ratio,
        direction_upper = gi_direction_upper(scenario),
        expected_n_h0 = as.numeric(ss$expectedNumberOfSubjectsH0),
        expected_n_h01 = as.numeric(ss$expectedNumberOfSubjectsH01),
        expected_n_h1 = as.numeric(ss$expectedNumberOfSubjectsH1)
      )
    ),
    class = c("gi_design", "list")
  )
}

#' Tabulate the monitoring boundaries of a design
#'
#' @param design A `gi_design` carrying rpact boundaries, that is one built by
#'   [design_group_sequential()] or [design_fixed()].
#' @return A data frame with one row per analysis and columns `analysis`,
#'   `information_rate`, `n_cumulative`, `efficacy_z`, `futility_z` and
#'   `cumulative_alpha_spent`. `futility_z` is `NA` where no futility bound
#'   applies, including at the final analysis, where rpact reports none.
#' @seealso [design_group_sequential()]
#' @examples
#' gs_boundaries(design_group_sequential(scenario("ercp_acute_cholangitis")))
#' @export
gs_boundaries <- function(design) {
  if (!inherits(design, "gi_design")) {
    stop("`design` must be a gi_design.", call. = FALSE)
  }
  d <- design$detail %||% list()
  if (is.null(d$rpact_design)) {
    stop(
      "`design` of type '", design$type,
      "' carries no rpact boundaries; gs_boundaries() supports fixed and ",
      "group_sequential designs.",
      call. = FALSE
    )
  }

  efficacy <- d$efficacy_z %||% as.numeric(d$critical_value)
  info <- d$information_rates %||% as.numeric(d$rpact_design$informationRates)
  k <- length(efficacy)
  n_cum <- d$n_cumulative %||% design$n_total
  futility <- d$futility_z %||% rep(NA_real_, k)

  data.frame(
    analysis = seq_len(k),
    information_rate = info,
    n_cumulative = n_cum,
    efficacy_z = efficacy,
    futility_z = futility,
    cumulative_alpha_spent = as.numeric(d$rpact_design$alphaSpent),
    stringsAsFactors = FALSE
  )
}

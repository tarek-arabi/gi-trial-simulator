#' Cross-engine benchmark
#'
#' The claim this package makes is that it delegates its statistics rather
#' than reimplementing them. This benchmark is one narrow check on that claim,
#' and it is worth being exact about its scope. It computes the one-sided
#' efficacy boundaries of a group-sequential design twice, once through
#' `rpact` and once through `gsDesign`, and reports the disagreement. What it
#' establishes is that the two engines agree on group-sequential efficacy
#' boundaries for the four design families it covers (`OF`, `P`, `asOF`,
#' `asP`), at the analysis counts and information rates it is run at.
#'
#' What it does not establish: it says nothing about sample sizes, nothing
#' about futility or beta-spending bounds, nothing about the Bayesian designs,
#' and nothing about any quantity this package obtains by simulation. Boundary
#' agreement is evidence that the boundaries are not homegrown. It is not
#' evidence that anything else in the package is right.
#'
#' Classical O'Brien-Fleming and Pocock boundaries are defined by a single
#' shape constraint and a root search, so the two engines should agree to
#' within their root-finding tolerances. The alpha-spending variants (`asOF`,
#' `asP`) are approximations to those shapes and each engine integrates the
#' multivariate normal distribution its own way, so a slightly larger
#' disagreement there is expected rather than alarming. That is why the
#' tolerance is an argument.
#'
#' @name cross-engine-benchmark
NULL

# rpact typeOfDesign to the gsDesign sfu that defines the same boundary
# family. Determined empirically, not from the documentation: classical
# families are named characters in gsDesign, spending families are functions.
gi_gsdesign_sfu <- function(type) {
  switch(type,
    OF = "OF",
    P = "Pocock",
    asOF = gsDesign::sfLDOF,
    asP = gsDesign::sfLDPocock,
    stop("Unsupported design type '", type, "'.", call. = FALSE)
  )
}

#' Compare rpact and gsDesign efficacy boundaries
#'
#' Computes the one-sided efficacy boundaries of a group-sequential design
#' with both engines and reports them side by side. `gsDesign::gsDesign()` is
#' called with `test.type = 1`, which is its pure one-sided efficacy-only
#' design and the counterpart of an rpact one-sided design without beta
#' spending.
#'
#' There is deliberately no type II error argument. In both engines these
#' boundaries are a function of `alpha`, the number of analyses, the
#' information rates and the spending family alone: `test.type = 1` bounds in
#' `gsDesign` and `criticalValues` in an rpact design without beta spending do
#' not move when beta moves. The function used to take a validated `beta` that
#' changed none of its output, which implied a check that was not happening.
#' Each engine is now left on its own default and `tests/testthat` pins the
#' invariance, so a future engine release that made the bounds depend on beta
#' would fail the suite rather than pass silently.
#'
#' @param k Number of analyses, an integer of at least 2.
#' @param alpha One-sided type I error rate.
#' @param type Boundary family: `"OF"` or `"P"` for the classical
#'   O'Brien-Fleming and Pocock boundaries, `"asOF"` or `"asP"` for the
#'   corresponding Lan-DeMets alpha-spending approximations.
#' @param tolerance Absolute z-scale difference at or below which the two
#'   engines are treated as agreeing.
#' @param information_rates Optional numeric vector of length `k`, strictly
#'   increasing and ending at 1. Defaults to equally spaced analyses.
#' @return A data frame with one row per analysis and columns `type`, `k`,
#'   `analysis`, `information_rate`, `rpact_z`, `gsdesign_z`, `abs_diff`,
#'   `tolerance` and `agrees`.
#' @seealso [benchmark_report()]
#' @examples
#' benchmark_gs_boundaries(k = 3, type = "OF")
#' benchmark_gs_boundaries(k = 4, type = "asOF", tolerance = 1e-4)
#' @export
benchmark_gs_boundaries <- function(k = 3, alpha = 0.025,
                                    type = "OF", tolerance = 1e-4,
                                    information_rates = NULL) {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != round(k) || k < 2) {
    stop("`k` must be a single integer of at least 2.", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
    alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number strictly between 0 and 0.5.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) ||
    tolerance <= 0) {
    stop("`tolerance` must be a single positive number.", call. = FALSE)
  }
  type <- gi_check_type_of_design(type)
  rates_arg <- gi_check_information_rates(information_rates, k)

  # Neither call is given a type II error rate: the quantities compared below
  # do not depend on one. See the note in this function's documentation.
  rp <- rpact::getDesignGroupSequential(
    kMax = k,
    alpha = alpha,
    sided = 1L,
    typeOfDesign = type,
    informationRates = rates_arg
  )
  info <- as.numeric(rp$informationRates)

  gs <- gsDesign::gsDesign(
    k = k,
    test.type = 1,
    alpha = alpha,
    sfu = gi_gsdesign_sfu(type),
    timing = if (is.null(information_rates)) 1 else info
  )

  rpact_z <- as.numeric(rp$criticalValues)
  gsdesign_z <- as.numeric(gs$upper$bound)
  abs_diff <- abs(rpact_z - gsdesign_z)

  data.frame(
    type = type,
    k = k,
    analysis = seq_len(k),
    information_rate = info,
    rpact_z = rpact_z,
    gsdesign_z = gsdesign_z,
    abs_diff = abs_diff,
    tolerance = tolerance,
    agrees = abs_diff <= tolerance,
    stringsAsFactors = FALSE
  )
}

#' Summarise cross-engine agreement across a grid of designs
#'
#' Runs [benchmark_gs_boundaries()] over a grid of analysis counts and
#' boundary families and reduces each design to one row. This is the table the
#' README and the validation section of the manuscript cite.
#'
#' @param k Integer vector of analysis counts to benchmark.
#' @param types Character vector of boundary families, drawn from `"OF"`,
#'   `"P"`, `"asOF"` and `"asP"`.
#' @param alpha One-sided type I error rate.
#' @param tolerance Absolute z-scale difference at or below which the two
#'   engines are treated as agreeing.
#' @return A data frame with one row per (`type`, `k`) combination and columns
#'   `type`, `k`, `n_analyses`, `max_abs_diff`, `tolerance` and `all_agree`.
#' @seealso [benchmark_gs_boundaries()]
#' @examples
#' benchmark_report()
#' benchmark_report(k = 2:3, types = c("OF", "P"))
#' @export
benchmark_report <- function(k = 2:5,
                             types = c("OF", "P", "asOF", "asP"),
                             alpha = 0.025, tolerance = 1e-4) {
  if (!is.numeric(k) || length(k) == 0L || anyNA(k) ||
    any(k != round(k)) || any(k < 2)) {
    stop("`k` must be a vector of integers, each at least 2.", call. = FALSE)
  }
  if (!is.character(types) || length(types) == 0L) {
    stop("`types` must be a non-empty character vector.", call. = FALSE)
  }
  for (type in types) gi_check_type_of_design(type)

  grid <- expand.grid(k = as.integer(k), type = types, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    b <- benchmark_gs_boundaries(
      k = grid$k[i], alpha = alpha,
      type = grid$type[i], tolerance = tolerance
    )
    data.frame(
      type = grid$type[i],
      k = grid$k[i],
      n_analyses = nrow(b),
      max_abs_diff = max(b$abs_diff),
      tolerance = tolerance,
      all_agree = all(b$agrees),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(out$type, out$k), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Monte Carlo simulation of trial designs
#'
#' The design engines in this package delegate closed-form calculations to
#' rpact and gsDesign. Simulation exists for the questions those engines cannot
#' answer analytically: realised sample size distributions, estimator bias at a
#' stopping boundary, coverage, and operating characteristics under rates other
#' than the ones a design was powered on. Nothing here recomputes a boundary or
#' a power curve that a validated engine already provides.
#'
#' Every replication is generated from an explicitly recorded random number
#' stream, and every simulation object carries the seed it was produced from.
#'
#' @name simulation
NULL

.rng_snapshot <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

# The caller's complete random number state: the seed vector when one exists,
# and the generator kinds. The kinds have to be captured separately because
# set.seed(kind = ) changes them globally and deleting .Random.seed does not
# put them back.
.rng_capture <- function() {
  list(seed = .rng_snapshot(), kind = RNGkind())
}

.rng_forget <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(NULL)
}

.rng_reset <- function(captured) {
  if (!is.null(captured$seed)) {
    assign(".Random.seed", captured$seed, envir = globalenv())
    return(invisible(NULL))
  }
  # The caller had no .Random.seed, so the kinds are all there is to put back.
  # Whatever is there now is dropped first: RNGkind() reads .Random.seed and
  # refuses to run on a malformed one, which .is_random_seed() can leave behind.
  .rng_forget()
  RNGkind(
    kind = captured$kind[1L],
    normal.kind = captured$kind[2L],
    sample.kind = captured$kind[3L]
  )
  .rng_forget()
  invisible(NULL)
}

# The generator a scalar seed is pinned to. Without pinning, set.seed() runs
# under whatever kind the calling session happened to be using, so a recorded
# seed would not be enough to reproduce a recorded result.
GI_SCALAR_SEED_RNG <- c(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

# A seed is either a single whole number, or a full .Random.seed vector as
# produced by sim_seeds(). Assigning the vector also restores its generator
# kind, which is how a task run under mclapply reproduces a task run under
# lapply. Returns the three generator kinds actually in force, which every
# simulation records alongside its seed.
.rng_set <- function(seed) {
  if (length(seed) > 1L) {
    assign(".Random.seed", as.integer(seed), envir = globalenv())
  } else {
    set.seed(
      as.integer(seed),
      kind = GI_SCALAR_SEED_RNG[["kind"]],
      normal.kind = GI_SCALAR_SEED_RNG[["normal.kind"]],
      sample.kind = GI_SCALAR_SEED_RNG[["sample.kind"]]
    )
  }
  RNGkind()
}

.check_count <- function(x, arg, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != as.integer(x)) {
    stop("`", arg, "` must be a single whole number of at least ", min, ".", call. = FALSE)
  }
  as.integer(x)
}

.check_root_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != as.integer(seed)) {
    stop("`seed` must be a single whole number.", call. = FALSE)
  }
  as.integer(seed)
}

# The only dependable test that a vector is a usable .Random.seed is to install
# it and draw from it, which is what the caller is about to do anyway. An
# unusable vector either errors or triggers R's "ignored, not an integer
# vector" warning, and both must fail here rather than several frames deeper.
.is_random_seed <- function(seed) {
  captured <- .rng_capture()
  on.exit(
    {
      .rng_forget()
      .rng_reset(captured)
    },
    add = TRUE
  )
  isTRUE(tryCatch(
    {
      assign(".Random.seed", seed, envir = globalenv())
      stats::runif(1L)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  ))
}

.check_seed <- function(seed) {
  msg <- paste0(
    "`seed` must be a single whole number, or a .Random.seed vector as ",
    "returned by sim_seeds()."
  )
  if (!is.numeric(seed) || length(seed) < 1L || anyNA(seed) || !all(is.finite(seed))) {
    stop(msg, call. = FALSE)
  }
  if (any(seed != trunc(seed)) || any(abs(seed) > .Machine$integer.max)) {
    stop(
      msg, " Every element must be a whole number no larger than ",
      .Machine$integer.max, " in absolute value.",
      call. = FALSE
    )
  }
  seed <- as.integer(seed)
  if (length(seed) > 1L && !.is_random_seed(seed)) {
    stop(
      msg, " A vector of length ", length(seed),
      " starting at ", seed[1L], " is not a valid .Random.seed for any ",
      "generator this R session knows.",
      call. = FALSE
    )
  }
  seed
}

.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  alpha
}

.check_scenario <- function(scenario, arg = "scenario") {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`", arg, "` must be a gi_scenario, as returned by scenario().", call. = FALSE)
  }
  scenario
}

# The binomial simulators in this file draw event counts, not continuous
# outcomes, so a continuous scenario must be refused here rather than read as
# a NULL control_rate and silently simulated as nonsense (or, with a stray
# non-NULL rate left over from a mistaken override, as a wrong number that
# looks plausible). This is the one guard behind simulate_fixed(),
# simulate_group_sequential() (via design$scenario) and simulate_grid() (via
# every scenario in the grid and every design_fn(scenario)$scenario).
.check_binary_scenario <- function(scenario, arg = "scenario", fn_name) {
  .check_scenario(scenario, arg)
  gi_require_endpoint_type(scenario, "binary", fn_name)
  scenario
}

# +1 when a higher event rate is the benefit, -1 when a lower one is.
.benefit_sign <- function(direction) {
  switch(direction,
    lower_is_better = -1,
    higher_is_better = 1,
    stop(
      "Unknown endpoint direction '", direction,
      "'; expected lower_is_better or higher_is_better.",
      call. = FALSE
    )
  )
}

.resolve_rates <- function(rates, scenario) {
  if (is.null(rates)) {
    return(c(control = scenario$control_rate, treatment = scenario$treatment_rate))
  }
  if (!is.numeric(rates) || anyNA(rates) || !length(rates) %in% c(1L, 2L)) {
    stop(
      "`rates` must be NULL, a single event rate applied to both arms, or ",
      "a numeric vector of length 2 giving the control and treatment rates.",
      call. = FALSE
    )
  }
  if (length(rates) == 1L) rates <- c(rates, rates)
  if (!is.null(names(rates)) && all(c("control", "treatment") %in% names(rates))) {
    rates <- c(rates[["control"]], rates[["treatment"]])
  }
  if (any(rates <= 0) || any(rates >= 1)) {
    stop("`rates` must lie strictly between 0 and 1.", call. = FALSE)
  }
  c(control = unname(rates[1L]), treatment = unname(rates[2L]))
}

# Pooled-variance two-proportion score statistic, signed so that larger values
# always mean more evidence of treatment benefit. This is the statistic
# stats::prop.test(correct = FALSE) reports, and test-simulate.R pins it to
# prop.test on a sample of replicates.
.score_z <- function(x_t, n_t, x_c, n_c, sign) {
  p_t <- x_t / n_t
  p_c <- x_c / n_c
  p_bar <- (x_t + x_c) / (n_t + n_c)
  se <- sqrt(p_bar * (1 - p_bar) * (1 / n_t + 1 / n_c))
  z <- sign * (p_t - p_c) / se
  # With no events, or events in every patient, the score test is undefined and
  # carries no evidence either way. Substituting 0 keeps the replicate in the
  # denominator as a non-rejection; .score_undefined() flags the same replicates
  # so the substitution is counted and reported rather than absorbed silently.
  z[!is.finite(z)] <- 0
  z
}

# TRUE for replicates whose pooled score statistic is undefined, that is where
# the pooled event proportion is exactly 0 or exactly 1.
.score_undefined <- function(x_t, n_t, x_c, n_c) {
  events <- x_t + x_c
  events == 0 | events == (n_t + n_c)
}

.risk_difference <- function(x_t, n_t, x_c, n_c, conf_level = 0.95) {
  p_t <- x_t / n_t
  p_c <- x_c / n_c
  rd <- p_t - p_c
  se <- sqrt(p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c)
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  # se is exactly 0 when both arm proportions sit on 0 or 1, which collapses the
  # Wald interval to a point and makes it cover nothing. Those replicates are
  # returned as they are and counted by the caller, not repaired here.
  list(rd = rd, se = se, lower = rd - crit * se, upper = rd + crit * se)
}

# One degeneracy record per simulation: how many replicates had an undefined
# score statistic, how many had a zero-width interval, and how many had either.
.degeneracy <- function(score_undefined, interval_undefined) {
  list(
    score_undefined = sum(score_undefined),
    interval_undefined = sum(interval_undefined),
    any = sum(score_undefined | interval_undefined)
  )
}

#' Independent random number streams for reproducible simulation
#'
#' Derives `n` non-overlapping L'Ecuyer-CMRG streams from a single integer
#' seed. Because each stream fully determines its task, a set of tasks gives
#' the same answers whether it is evaluated serially or across any number of
#' workers.
#'
#' @param seed Single integer, the root seed.
#' @param n Number of streams to derive.
#' @return A list of `n` integer vectors, each a valid `.Random.seed` for the
#'   L'Ecuyer-CMRG generator.
#' @examples
#' streams <- sim_seeds(2026, 3)
#' length(streams)
#' @export
sim_seeds <- function(seed, n) {
  seed <- .check_root_seed(seed)
  n <- .check_count(n, "n")

  state <- .rng_capture()
  on.exit(.rng_reset(state), add = TRUE)

  set.seed(seed, kind = "L'Ecuyer-CMRG")
  stream <- .rng_snapshot()
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  out
}

.run_tasks <- function(tasks, fun, workers) {
  workers <- .check_count(workers, "workers")
  parallel_ok <- workers > 1L && .Platform$OS.type != "windows"
  out <- if (parallel_ok) {
    parallel::mclapply(tasks, fun, mc.cores = workers, mc.set.seed = FALSE)
  } else {
    lapply(tasks, fun)
  }
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(
      "Simulation task ", which(failed)[1L], " failed: ",
      as.character(out[[which(failed)[1L]]]),
      call. = FALSE
    )
  }
  out
}

#' Simulate a fixed-sample two-arm binary-endpoint trial
#'
#' Draws event counts for both arms directly from the binomial distribution,
#' vectorised over replications, and analyses each replicate with the one-sided
#' pooled-variance score test of treatment superiority. The direction of
#' superiority is taken from the scenario, so a mortality endpoint is tested for
#' a reduction and a response endpoint for an increase.
#'
#' The statistic is the one [stats::prop.test()] reports with
#' `correct = FALSE`; the package test suite pins it to `prop.test` on a sample
#' of replicates rather than trusting the arithmetic.
#'
#' @section Agreement with rpact:
#' Simulated power reproduces the analytic power [power_at()] reads from rpact
#' to within about half a percentage point, and is consistently slightly
#' higher. Characterised at 200,000 replications across five scenarios spanning
#' control rates from 0.066 to 0.50, the simulated rejection rate ran between
#' 0.003 and 0.007 above the rpact value, with the sign positive every time.
#' Neither number is wrong: rpact evaluates a normal-approximation formula for
#' two rates, while this function simulates the discrete binomial test that a
#' real trial would perform, and the two differ by about that much. Treat a gap
#' of this size and sign as expected, and a gap above 0.01, or one of the
#' opposite sign, as a reason to investigate.
#'
#' @section Degenerate replicates:
#' A replicate in which nobody has an event, or everybody does, has no defined
#' score statistic, and one in which both arm proportions sit on 0 or 1 has a
#' zero-width confidence interval. Such replicates are recorded with `z = 0`,
#' which counts them as non-rejections, and are flagged in the `degenerate`
#' column and counted in the `degenerate` element so that they are visible
#' rather than silently absorbed into the rejection rate. They are common only
#' at small sample sizes or very low event rates.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param n_per_arm Number of patients randomised to each arm.
#' @param nsim Number of simulated trials.
#' @param alpha One-sided type I error rate the test is judged against.
#' @param seed Root seed. Either a single whole number or a `.Random.seed`
#'   vector from [sim_seeds()]. A single number also pins the generator kind to
#'   Mersenne-Twister with inversion and rejection sampling, so the recorded
#'   seed reproduces the recorded results whatever generator the calling session
#'   was using. The caller's random number state is left unchanged.
#' @param workers Retained for interface symmetry with [simulate_grid()].
#'   Replication is vectorised rather than looped, so values above 1 have no
#'   effect here and results never depend on it. Parallelism is applied across
#'   scenarios in [simulate_grid()].
#' @return An object of class `gi_simulation`: a list with `results` (a data
#'   frame with one row per replicate holding `z`, `p`, `reject`,
#'   `risk_difference`, its standard error and 95 percent confidence limits,
#'   the two arm event counts, `n_total` and the logical `degenerate`), plus
#'   `design_type`, `n_per_arm`, `n_total`, `nsim`, `alpha`, `seed`,
#'   `rng_kind` (the three generator kinds the replicates were drawn under),
#'   `rates`, `scenario` and `degenerate` (counts of the replicates with an
#'   undefined score statistic, a zero-width interval, and either).
#' @seealso [ademp_summary()], [simulate_group_sequential()], [simulate_grid()]
#' @examples
#' sc <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
#' sim <- simulate_fixed(sc, n_per_arm = 161, nsim = 500, seed = 1)
#' mean(sim$results$reject)
#' @export
simulate_fixed <- function(scenario, n_per_arm, nsim = 10000, alpha = 0.025,
                           seed = 1, workers = 1) {
  .check_binary_scenario(scenario, fn_name = "simulate_fixed")
  n_per_arm <- .check_count(n_per_arm, "n_per_arm")
  nsim <- .check_count(nsim, "nsim")
  alpha <- .check_alpha(alpha)
  seed <- .check_seed(seed)
  .check_count(workers, "workers")

  rates <- .resolve_rates(NULL, scenario)
  sign <- .benefit_sign(scenario$direction)

  state <- .rng_capture()
  on.exit(.rng_reset(state), add = TRUE)
  rng_kind <- .rng_set(seed)

  x_c <- stats::rbinom(nsim, n_per_arm, rates[["control"]])
  x_t <- stats::rbinom(nsim, n_per_arm, rates[["treatment"]])

  z <- .score_z(x_t, n_per_arm, x_c, n_per_arm, sign)
  p <- stats::pnorm(z, lower.tail = FALSE)
  rd <- .risk_difference(x_t, n_per_arm, x_c, n_per_arm)
  score_undefined <- .score_undefined(x_t, n_per_arm, x_c, n_per_arm)
  interval_undefined <- rd$se == 0

  results <- data.frame(
    replicate = seq_len(nsim),
    events_control = x_c,
    events_treatment = x_t,
    z = z,
    p = p,
    reject = p < alpha,
    risk_difference = rd$rd,
    risk_difference_se = rd$se,
    ci_lower = rd$lower,
    ci_upper = rd$upper,
    n_total = 2L * n_per_arm,
    degenerate = score_undefined | interval_undefined,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      design_type = "fixed",
      results = results,
      n_per_arm = n_per_arm,
      n_total = 2L * n_per_arm,
      nsim = nsim,
      alpha = alpha,
      seed = seed,
      rng_kind = rng_kind,
      rates = rates,
      scenario = scenario,
      degenerate = .degeneracy(score_undefined, interval_undefined),
      engine = "stats::rbinom + pooled score test"
    ),
    class = c("gi_simulation", "list")
  )
}

# Locate the rpact design a gi_design was built from, wherever the design
# module chose to store it. Boundaries are always read from that object.
.find_rpact_design <- function(x) {
  if (inherits(x, "TrialDesign")) {
    return(x)
  }
  if (!is.list(x) || length(x) == 0L) {
    return(NULL)
  }
  for (element in x) {
    hit <- .find_rpact_design(element)
    if (!is.null(hit)) {
      return(hit)
    }
  }
  NULL
}

.gs_boundaries <- function(design) {
  rp <- .find_rpact_design(design$detail)
  if (is.null(rp)) rp <- .find_rpact_design(design)

  if (!is.null(rp)) {
    if (!is.null(rp$sided) && identical(as.integer(rp$sided), 2L)) {
      stop(
        "`design` holds a two-sided rpact design; simulation supports one-sided ",
        "superiority designs only.",
        call. = FALSE
      )
    }
    bounds <- list(
      critical = as.numeric(rp$criticalValues),
      futility = as.numeric(rp$futilityBounds %||% numeric()),
      information = as.numeric(rp$informationRates),
      binding = isTRUE(rp$bindingFutility),
      source = "rpact design object"
    )
  } else {
    detail <- design$detail %||% list()
    critical <- detail$efficacy_z %||% detail$critical_values
    information <- detail$information_rates
    if (is.null(critical) || is.null(information)) {
      stop(
        "`design` carries no rpact design object and no `efficacy_z`/",
        "`information_rates` in `detail`; boundaries cannot be read and this ",
        "function will not recompute them.",
        call. = FALSE
      )
    }
    bounds <- list(
      critical = as.numeric(critical),
      futility = as.numeric(detail$futility_z %||% detail$futility_bounds %||% numeric()),
      information = as.numeric(information),
      binding = isTRUE(detail$binding_futility),
      source = "design$detail"
    )
  }

  kmax <- length(bounds$critical)
  if (kmax < 1L || length(bounds$information) != kmax) {
    stop(
      "`design` has ", kmax, " critical values but ",
      length(bounds$information), " information rates.",
      call. = FALSE
    )
  }
  futility <- rep(-Inf, max(kmax - 1L, 0L))
  if (length(bounds$futility)) {
    take <- min(length(bounds$futility), length(futility))
    if (take > 0L) futility[seq_len(take)] <- bounds$futility[seq_len(take)]
  }
  # rpact reports -6 (not -Inf) for a stage with no futility bound, and
  # design$detail uses NA for the same thing. Both mean "never stop here".
  no_futility <- -6
  futility[is.na(futility) | !is.finite(futility) | futility <= no_futility] <- -Inf
  bounds$futility <- futility
  bounds$kmax <- kmax
  bounds
}

#' Simulate a group-sequential trial
#'
#' Accumulates events across the design's information rates, compares the
#' interim score statistics to the efficacy and futility boundaries held on the
#' design, and records where each replicate stopped, whether it rejected, and
#' how many patients it used. Boundaries are read from the rpact design object
#' the design was built with. They are never recomputed here.
#'
#' Futility bounds are applied as given. For a non-binding futility design the
#' simulated type I error is therefore the rate under a trial that does stop for
#' futility, which is below the rate rpact reports for the binding-free
#' boundary.
#'
#' Interim looks are placed at `round(information_rate * n_per_arm)` patients
#' per arm, so a design whose information rates do not divide its sample size
#' exactly is simulated at the nearest whole patient. As with [simulate_fixed()]
#' the simulated operating characteristics are those of the discrete binomial
#' test rather than of rpact's normal approximation, which leaves simulated
#' power slightly above the analytic value and the realised mean sample size
#' slightly above rpact's expected sample size.
#'
#' @param design A `gi_design` of type `group_sequential`.
#' @param nsim Number of simulated trials.
#' @param seed Root seed. Either a single whole number or a `.Random.seed`
#'   vector from [sim_seeds()]. A single number also pins the generator kind, as
#'   documented in [simulate_fixed()]. The caller's random number state is left
#'   unchanged.
#' @param workers Retained for interface symmetry with [simulate_grid()].
#'   Replication is vectorised within each stage, so values above 1 have no
#'   effect here and results never depend on it.
#' @param rates Optional numeric override of the scenario event rates, either a
#'   single rate applied to both arms (which simulates under the null) or
#'   `c(control, treatment)`. Boundaries and sample size still come from the
#'   design; only the data-generating rates change.
#' @return An object of class `gi_simulation`: a list with `results` (one row
#'   per replicate holding `stop_stage`, `stopped_for`, `reject`, `z`, `p`, the
#'   risk difference at the stopping analysis with its 95 percent confidence
#'   limits, the realised `n_total`, and the logical `degenerate` flagging a
#'   stopping analysis whose score statistic or interval was undefined), plus
#'   `design_type`, `n_per_arm`, `n_total`, `nsim`, `alpha`, `seed`,
#'   `rng_kind`, `rates`, `scenario`, `design`, `degenerate` (the counts behind
#'   that flag) and a `boundaries` element recording exactly what was compared
#'   against.
#' @seealso [ademp_summary()], [simulate_fixed()]
#' @examples
#' sc <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
#' rp <- rpact::getDesignGroupSequential(
#'   kMax = 2, alpha = 0.025, sided = 1, typeOfDesign = "asOF"
#' )
#' design <- structure(
#'   list(
#'     type = "group_sequential", scenario = sc, alpha = 0.025, power = 0.9,
#'     n_total = 400, n_per_arm = 200, engine = "rpact::getDesignGroupSequential",
#'     detail = list(rpact_design = rp)
#'   ),
#'   class = c("gi_design", "list")
#' )
#' sim <- simulate_group_sequential(design, nsim = 200, seed = 1)
#' mean(sim$results$reject)
#' @export
simulate_group_sequential <- function(design, nsim = 10000, seed = 1,
                                      workers = 1, rates = NULL) {
  if (!inherits(design, "gi_design")) {
    stop("`design` must be a gi_design.", call. = FALSE)
  }
  if (!identical(design$type, "group_sequential")) {
    stop(
      "`design` must have type 'group_sequential', not '",
      design$type %||% "<missing>", "'.",
      call. = FALSE
    )
  }
  scenario <- .check_binary_scenario(
    design$scenario, "design$scenario",
    fn_name = "simulate_group_sequential"
  )
  nsim <- .check_count(nsim, "nsim")
  seed <- .check_seed(seed)
  .check_count(workers, "workers")

  n_per_arm <- .check_count(design$n_per_arm, "design$n_per_arm")
  alpha <- .check_alpha(design$alpha %||% 0.025)
  rates <- .resolve_rates(rates, scenario)
  sign <- .benefit_sign(scenario$direction)

  bounds <- .gs_boundaries(design)
  kmax <- bounds$kmax

  n_stage <- round(bounds$information * n_per_arm)
  n_stage[kmax] <- n_per_arm
  n_stage <- pmax(n_stage, 1L)
  n_stage <- cummax(n_stage)
  increment <- diff(c(0, n_stage))

  state <- .rng_capture()
  on.exit(.rng_reset(state), add = TRUE)
  rng_kind <- .rng_set(seed)

  cum_c <- integer(nsim)
  cum_t <- integer(nsim)
  active <- rep(TRUE, nsim)
  stop_stage <- integer(nsim)
  stop_z <- numeric(nsim)
  stop_x_c <- integer(nsim)
  stop_x_t <- integer(nsim)
  stop_n <- integer(nsim)
  reject <- logical(nsim)
  stopped_for <- character(nsim)
  stop_undefined <- logical(nsim)

  for (k in seq_len(kmax)) {
    if (increment[k] > 0L) {
      cum_c <- cum_c + stats::rbinom(nsim, increment[k], rates[["control"]])
      cum_t <- cum_t + stats::rbinom(nsim, increment[k], rates[["treatment"]])
    }
    z_k <- .score_z(cum_t, n_stage[k], cum_c, n_stage[k], sign)
    undefined_k <- .score_undefined(cum_t, n_stage[k], cum_c, n_stage[k])

    cross_efficacy <- active & z_k >= bounds$critical[k]
    cross_futility <- if (k < kmax) {
      active & !cross_efficacy & z_k <= bounds$futility[k]
    } else {
      rep(FALSE, nsim)
    }
    final <- if (k == kmax) active & !cross_efficacy else rep(FALSE, nsim)
    ending <- cross_efficacy | cross_futility | final

    if (any(ending)) {
      stop_stage[ending] <- k
      stop_z[ending] <- z_k[ending]
      stop_x_c[ending] <- cum_c[ending]
      stop_x_t[ending] <- cum_t[ending]
      stop_n[ending] <- n_stage[k]
      stop_undefined[ending] <- undefined_k[ending]
      reject[ending] <- cross_efficacy[ending]
      stopped_for[cross_efficacy] <- "efficacy"
      stopped_for[cross_futility] <- "futility"
      stopped_for[final] <- "final_analysis"
      active[ending] <- FALSE
    }
    if (!any(active)) break
  }

  rd <- .risk_difference(stop_x_t, stop_n, stop_x_c, stop_n)
  interval_undefined <- rd$se == 0

  results <- data.frame(
    replicate = seq_len(nsim),
    stop_stage = stop_stage,
    stopped_for = stopped_for,
    events_control = stop_x_c,
    events_treatment = stop_x_t,
    z = stop_z,
    p = stats::pnorm(stop_z, lower.tail = FALSE),
    reject = reject,
    risk_difference = rd$rd,
    risk_difference_se = rd$se,
    ci_lower = rd$lower,
    ci_upper = rd$upper,
    n_per_arm = stop_n,
    n_total = 2L * stop_n,
    degenerate = stop_undefined | interval_undefined,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      design_type = "group_sequential",
      results = results,
      n_per_arm = n_per_arm,
      n_total = 2L * n_per_arm,
      nsim = nsim,
      alpha = alpha,
      seed = seed,
      rng_kind = rng_kind,
      rates = rates,
      scenario = scenario,
      design = design,
      degenerate = .degeneracy(stop_undefined, interval_undefined),
      boundaries = list(
        kmax = kmax,
        information_rates = bounds$information,
        critical_values = bounds$critical,
        futility_bounds = bounds$futility,
        binding_futility = bounds$binding,
        n_per_arm_by_stage = n_stage,
        read_from = bounds$source
      ),
      engine = "stats::rbinom + boundaries from rpact design"
    ),
    class = c("gi_simulation", "list")
  )
}

#' Simulate a grid of scenarios
#'
#' Builds a design for each scenario with `design_fn`, simulates it, and
#' returns one tidy row per scenario. Each scenario is given its own
#' L'Ecuyer-CMRG stream from [sim_seeds()], so the grid returns identical
#' numbers however many workers evaluate it.
#'
#' @param scenario_grid A `gi_scenario`, or a list of them. Names, if present,
#'   are used to label rows.
#' @param design_fn Function of one argument (a `gi_scenario`) returning either
#'   a `gi_design` or a single number giving the per-arm sample size of a fixed
#'   design. Whatever it returns, every scenario is simulated at that
#'   scenario's own event rates: a design built on other rates contributes its
#'   sample size and, for a group-sequential design, its boundaries, but never
#'   the rates the data are drawn from. Its endpoint direction must match the
#'   scenario's, since that is what orients the one-sided test.
#' @param nsim Number of simulated trials per scenario.
#' @param seed Root seed from which the per-scenario streams are derived.
#' @param workers Number of parallel workers. Values above 1 use
#'   [parallel::mclapply()] where the platform supports forking, and fall back
#'   to serial evaluation on Windows.
#' @return A data frame with one row per scenario: identifying columns, the
#'   design used, and the ADEMP performance measures with their Monte Carlo
#'   standard errors. The individual `gi_simulation` objects are attached as the
#'   `"simulations"` attribute.
#' @seealso [ademp_summary()], [sim_seeds()]
#' @examples
#' grid <- list(
#'   optimistic = scenario("ercp_acute_cholangitis",
#'     control_rate = 0.30, treatment_rate = 0.15
#'   ),
#'   null = scenario("ercp_acute_cholangitis",
#'     control_rate = 0.30, treatment_rate = 0.30
#'   )
#' )
#' simulate_grid(grid, design_fn = function(s) 161, nsim = 200, seed = 7)
#' @export
simulate_grid <- function(scenario_grid, design_fn, nsim = 10000, seed = 1,
                          workers = 1) {
  if (inherits(scenario_grid, "gi_scenario")) scenario_grid <- list(scenario_grid)
  if (!is.list(scenario_grid) || length(scenario_grid) == 0L) {
    stop("`scenario_grid` must be a gi_scenario or a non-empty list of them.", call. = FALSE)
  }
  for (i in seq_along(scenario_grid)) {
    .check_binary_scenario(
      scenario_grid[[i]], paste0("scenario_grid[[", i, "]]"),
      fn_name = "simulate_grid"
    )
  }
  if (!is.function(design_fn)) {
    stop("`design_fn` must be a function taking a gi_scenario.", call. = FALSE)
  }
  nsim <- .check_count(nsim, "nsim")
  seed <- .check_root_seed(seed)
  workers <- .check_count(workers, "workers")

  labels <- names(scenario_grid)
  if (is.null(labels)) labels <- rep("", length(scenario_grid))
  blank <- !nzchar(labels)
  labels[blank] <- paste0("scenario_", seq_along(labels))[blank]

  streams <- sim_seeds(seed, length(scenario_grid))
  tasks <- seq_along(scenario_grid)

  run_one <- function(i) {
    sc <- scenario_grid[[i]]
    design <- design_fn(sc)
    if (is.numeric(design) && length(design) == 1L) {
      simulate_fixed(sc, n_per_arm = design, nsim = nsim, seed = streams[[i]])
    } else if (inherits(design, "gi_design") && identical(design$type, "group_sequential")) {
      # The grid scenario, not the scenario the design was powered on, supplies
      # the data-generating rates. Without this a grid of scenarios would be
      # simulated at whatever rates design_fn happened to build its design from.
      .check_binary_scenario(
        design$scenario, "design_fn(scenario)$scenario",
        fn_name = "simulate_grid"
      )
      if (!identical(design$scenario$direction, sc$direction)) {
        stop(
          "`design_fn` returned a group-sequential design whose endpoint ",
          "direction ('", design$scenario$direction, "') differs from that of ",
          "grid scenario ", i, " ('", sc$direction, "'), so the simulated test ",
          "would be oriented against the wrong tail.",
          call. = FALSE
        )
      }
      simulate_group_sequential(design,
        nsim = nsim, seed = streams[[i]],
        rates = c(control = sc$control_rate, treatment = sc$treatment_rate)
      )
    } else if (inherits(design, "gi_design")) {
      simulate_fixed(sc, n_per_arm = design$n_per_arm, nsim = nsim,
        alpha = design$alpha %||% 0.025, seed = streams[[i]]
      )
    } else {
      stop(
        "`design_fn` must return a gi_design or a single per-arm sample size; ",
        "it returned an object of class ", paste(class(design), collapse = "/"), ".",
        call. = FALSE
      )
    }
  }

  sims <- .run_tasks(tasks, run_one, workers)

  rows <- lapply(seq_along(sims), function(i) {
    sim <- sims[[i]]
    ademp <- ademp_summary(sim)
    wide <- as.list(stats::setNames(ademp$estimate, ademp$measure))
    wide_mcse <- as.list(stats::setNames(ademp$mcse, paste0(ademp$measure, "_mcse")))
    # The row is labelled by the grid scenario, so it is identified by the grid
    # scenario too, whatever scenario design_fn built its design from.
    sc <- scenario_grid[[i]]
    base <- data.frame(
      scenario = labels[i],
      pack_id = sc$pack_id,
      endpoint = sc$endpoint,
      control_rate = sim$rates[["control"]],
      treatment_rate = sim$rates[["treatment"]],
      true_risk_difference = sim$rates[["treatment"]] - sim$rates[["control"]],
      design_type = sim$design_type,
      n_per_arm_max = sim$n_per_arm,
      n_total_max = sim$n_total,
      nsim = sim$nsim,
      alpha = sim$alpha,
      stringsAsFactors = FALSE
    )
    cbind(base, as.data.frame(c(wide, wide_mcse), stringsAsFactors = FALSE))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "simulations") <- stats::setNames(sims, labels)
  out
}

#' @export
print.gi_simulation <- function(x, ...) {
  cat("<gi_simulation> ", x$design_type, "\n", sep = "")
  cat("  scenario: ", x$scenario$pack_id, " / ", x$scenario$endpoint, "\n", sep = "")
  cat(sprintf(
    "  rates:    control %.4f, treatment %.4f\n",
    x$rates[["control"]], x$rates[["treatment"]]
  ))
  cat(sprintf(
    "  %d replications, up to %d per arm, one-sided alpha %.4g, seed %s\n",
    x$nsim, x$n_per_arm, x$alpha,
    if (length(x$seed) > 1L) "<stream>" else as.character(x$seed)
  ))
  rejection <- mean(x$results$reject)
  cat(sprintf(
    "  rejection rate %.4f (MCSE %.4f)\n",
    rejection, sqrt(rejection * (1 - rejection) / x$nsim)
  ))
  if (identical(x$design_type, "group_sequential")) {
    cat(sprintf("  mean total sample size %.1f\n", mean(x$results$n_total)))
  }
  degenerate <- x$degenerate$any %||% 0L
  if (degenerate > 0L) {
    cat(sprintf(
      "  %d of %d replications degenerate (no events, or every patient an event)\n",
      degenerate, x$nsim
    ))
  }
  invisible(x)
}

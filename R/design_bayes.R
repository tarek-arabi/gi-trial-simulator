#' Bayesian adaptive designs evaluated by simulation
#'
#' A two-arm binary-endpoint trial with an independent conjugate Beta prior on
#' each arm's event rate. At every interim look the posteriors are
#' `Beta(a0 + events, b0 + n - events)` by conjugacy, and the decision quantity
#' is `Pr(treatment better than control | data)`. Unlike a group-sequential
#' design there is no closed form for the operating characteristics of this
#' procedure, so type I error, power and expected sample size are established
#' by Monte Carlo simulation. This is the pattern the FDA Complex Innovative
#' Trial Design programme expects: a design whose behaviour is demonstrated by
#' simulation rather than asserted from a formula.
#'
#' Nothing here duplicates `rpact` or `gsDesign`. Beta-binomial conjugate
#' arithmetic and Monte Carlo evaluation of a design with no closed form are
#' outside what those packages provide.
#'
#' @name bayesian-designs
NULL

# The whole module runs under L'Ecuyer-CMRG so that every replication can be
# given its own independent stream. That is what makes results identical
# regardless of how the work is split across workers.
gi_rng_kind <- "L'Ecuyer-CMRG"

rng_snapshot <- function() {
  list(
    kind = RNGkind(),
    seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
  )
}

rng_restore <- function(snap) {
  RNGkind(snap$kind[1], snap$kind[2], snap$kind[3])
  if (is.null(snap$seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", snap$seed, envir = globalenv())
  }
  invisible(NULL)
}

check_count <- function(x, name, min = 1L) {
  # The comparison is against trunc(), not as.integer(). as.integer() returns NA
  # with a coercion warning for anything outside 32-bit integer range, so a
  # guard written with it fails on out-of-range input with base R's "missing
  # value where TRUE/FALSE needed" instead of the message intended here. The
  # representable range is therefore checked explicitly and reported.
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == trunc(x) && x >= min && abs(x) <= .Machine$integer.max
  if (!ok) {
    stop(
      "`", name, "` must be a single whole number of at least ", min,
      " and at most ", .Machine$integer.max, ".",
      call. = FALSE
    )
  }
  as.integer(x)
}

# How many standard errors of an estimated posterior probability must fit
# between 0.5 and a decision threshold. Below `gi_post_draws_min_sd` the
# simulated design is not the design that was asked for, so the call is an
# error; below `gi_post_draws_rec_sd` there is residual Monte Carlo bias worth
# warning about. See the `post_draws` entry of [design_bayesian()] for the
# derivation and for the magnitude of the damage.
gi_post_draws_min_sd <- 10
gi_post_draws_rec_sd <- 20

# Draws needed for the largest standard error of an estimated posterior
# probability, 0.5 / sqrt(draws), to be `sd_factor` times smaller than the
# distance from `threshold` to 0.5. That standard error is attained at a true
# probability of 0.5, which is where the null distribution of the decision
# quantity is centred, so this is the ratio that decides whether the bulk of
# the null can cross the threshold on estimation noise alone. Returns NA for a
# threshold at or below 0.5, where no number of draws makes the comparison
# meaningful and only the absolute floor applies.
bayes_post_draws_needed <- function(threshold, sd_factor) {
  gap <- threshold - 0.5
  if (!is.finite(gap) || gap <= 0) {
    return(NA_integer_)
  }
  # The rounding before ceiling() is not cosmetic: (0.5 * 10 / (0.6 - 0.5))^2 is
  # 2500.0000000000005 in binary floating point, and ceiling() alone would turn
  # a requirement of 2500 draws into one of 2501.
  needed <- ceiling(round((0.5 * sd_factor / gap)^2, 8))
  as.integer(min(needed, .Machine$integer.max))
}

# The gap between 0.5 and a threshold can never exceed 0.5, so this is the
# smallest number of draws that can resolve any threshold at all.
bayes_post_draws_floor <- function() {
  bayes_post_draws_needed(1, gi_post_draws_min_sd)
}

check_post_draws <- function(post_draws, threshold) {
  post_draws <- check_count(post_draws, "post_draws", min = 2L)
  needed <- bayes_post_draws_needed(threshold, gi_post_draws_min_sd)
  recommended <- bayes_post_draws_needed(threshold, gi_post_draws_rec_sd)
  required <- if (is.na(needed)) bayes_post_draws_floor() else max(needed, bayes_post_draws_floor())
  if (post_draws < required) {
    stop(
      "`post_draws` = ", post_draws, " is too few to resolve a decision threshold",
      if (is.na(needed)) "" else paste0(" of ", signif(threshold, 4)),
      "; at least ", required, " are needed. Each posterior probability is a mean of ",
      "`post_draws` draws, so its standard error is up to ", signif(0.5 / sqrt(post_draws), 3),
      ", while the null distribution of that probability is centred on 0.5. Trials then cross ",
      "the threshold on estimation noise alone, which biases the simulated type I error ",
      "upwards by more than tenfold at `post_draws` = 2, and no value of `nsim` removes that ",
      "bias. See the `post_draws` entry in ?design_bayesian.",
      call. = FALSE
    )
  }
  if (!is.na(recommended) && post_draws < recommended) {
    warning(
      "`post_draws` = ", post_draws, " is near the resolution limit for a decision threshold of ",
      signif(threshold, 4), ": the standard error of each posterior probability is up to ",
      signif(0.5 / sqrt(post_draws), 3), ", more than a ", gi_post_draws_rec_sd, "th of the ",
      "distance from 0.5 to that threshold. The reported operating characteristics carry ",
      "residual Monte Carlo bias that `nsim` cannot reduce. Use at least ", recommended,
      " draws; the default of 4000 leaves more room again.",
      call. = FALSE
    )
  }
  post_draws
}

check_unit <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", name, "` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

check_prior <- function(prior) {
  if (!is.numeric(prior) || length(prior) != 2L || any(is.na(prior)) ||
    any(!is.finite(prior)) || any(prior <= 0)) {
    stop(
      "`prior` must be two finite positive numbers, the Beta shape parameters c(a0, b0).",
      call. = FALSE
    )
  }
  as.numeric(prior)
}

check_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
    !direction %in% c("lower_is_better", "higher_is_better")) {
    stop(
      "`direction` must be \"lower_is_better\" or \"higher_is_better\", not ",
      deparse(direction), ".",
      call. = FALSE
    )
  }
  direction
}

#' Posterior probability that treatment beats control
#'
#' With independent `Beta(a0, b0)` priors, the posterior for an arm observing
#' `events` responses out of `n` is `Beta(a0 + events, b0 + n - events)`. This
#' function reports the posterior probability that the treatment arm is the
#' better of the two, estimated by Monte Carlo draws from the two posteriors.
#' Which tail counts as better is set by `direction`, so an endpoint like
#' mortality is handled without the caller flipping any signs.
#'
#' The estimate carries Monte Carlo error of at most `0.5 / sqrt(draws)`, so
#' 4000 draws gives a standard error under 0.008 and 100000 draws under 0.0016.
#' Choose `draws` against the decision thresholds the value will be compared to.
#' When the value drives a stopping rule that error becomes a bias in the
#' decision rather than noise around it, which is why [design_bayesian()]
#' enforces a minimum: see the `post_draws` entry there for the rule and for
#' what happens below it.
#'
#' @param events_t Number of events observed in the treatment arm.
#' @param n_t Number of participants observed in the treatment arm.
#' @param events_c Number of events observed in the control arm.
#' @param n_c Number of participants observed in the control arm.
#' @param prior Length-2 numeric vector `c(a0, b0)` giving the Beta prior shape
#'   parameters, applied independently to both arms. Defaults to the uniform
#'   prior `c(1, 1)`.
#' @param direction Either `"lower_is_better"` (treatment wins when its event
#'   rate is lower, the usual case for mortality) or `"higher_is_better"`.
#' @param draws Number of Monte Carlo draws from each posterior.
#' @param seed Optional integer seed. When supplied, the generator kind is also
#'   fixed, so the returned value depends on nothing but the arguments, and the
#'   caller's random number state and generator kind are both restored on exit.
#'   When `NULL` the ambient stream is used, which is what the internal
#'   simulation loop needs in order to drive each replication from its own
#'   independent stream.
#' @return A single number in `[0, 1]`: the posterior probability that the
#'   treatment arm is better, in the sense given by `direction`.
#' @examples
#' posterior_prob_better(events_t = 5, n_t = 50, events_c = 12, n_c = 50, seed = 42)
#' posterior_prob_better(
#'   events_t = 5, n_t = 50, events_c = 12, n_c = 50,
#'   direction = "higher_is_better", seed = 42
#' )
#' @export
posterior_prob_better <- function(events_t, n_t, events_c, n_c,
                                  prior = c(1, 1),
                                  direction = "lower_is_better",
                                  draws = 4000,
                                  seed = NULL) {
  n_t <- check_count(n_t, "n_t", min = 1L)
  n_c <- check_count(n_c, "n_c", min = 1L)
  events_t <- check_count(events_t, "events_t", min = 0L)
  events_c <- check_count(events_c, "events_c", min = 0L)
  if (events_t > n_t) {
    stop("`events_t` (", events_t, ") cannot exceed `n_t` (", n_t, ").", call. = FALSE)
  }
  if (events_c > n_c) {
    stop("`events_c` (", events_c, ") cannot exceed `n_c` (", n_c, ").", call. = FALSE)
  }
  prior <- check_prior(prior)
  direction <- check_direction(direction)
  draws <- check_count(draws, "draws", min = 2L)
  if (!is.null(seed)) {
    seed <- check_count(seed, "seed", min = -.Machine$integer.max)
    snap <- rng_snapshot()
    on.exit(rng_restore(snap), add = TRUE)
    set.seed(seed, kind = gi_rng_kind)
  }

  p_t <- stats::rbeta(draws, prior[1] + events_t, prior[2] + n_t - events_t)
  p_c <- stats::rbeta(draws, prior[1] + events_c, prior[2] + n_c - events_c)
  if (identical(direction, "lower_is_better")) mean(p_t < p_c) else mean(p_t > p_c)
}

bayes_look_sizes <- function(n_per_arm, looks) {
  look_n <- round(n_per_arm * seq_len(looks) / looks)
  look_n[looks] <- n_per_arm
  if (any(diff(c(0L, look_n)) < 1L)) {
    stop(
      "`looks` = ", looks, " is too many for n_max = ", 2L * n_per_arm,
      "; each look must add at least one participant per arm.",
      call. = FALSE
    )
  }
  as.integer(look_n)
}

# One independent stream per replication, split into a set for the alternative
# and a set for the null. Because the stream is a property of the replication
# index and not of the worker, any chunking gives bit-identical results.
bayes_stream_sets <- function(seed, nsim) {
  snap <- rng_snapshot()
  on.exit(rng_restore(snap), add = TRUE)
  set.seed(seed, kind = gi_rng_kind)
  s <- get(".Random.seed", envir = globalenv())
  out <- vector("list", 2L * nsim)
  for (i in seq_len(2L * nsim)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  list(alt = out[seq_len(nsim)], null = out[nsim + seq_len(nsim)])
}

# Returns an nsim x looks matrix of posterior probabilities. Every look is
# evaluated for every replication, whether or not the trial would already have
# stopped, so that the same matrix can be rescored under any thresholds. This
# is what makes threshold calibration cheap and free of extra simulation noise.
bayes_trajectories <- function(p_control, p_treatment, look_n, prior, direction,
                               post_draws, streams, workers) {
  nsim <- length(streams)
  increments <- diff(c(0L, look_n))
  n_chunks <- min(nsim, 16L)
  chunks <- split(seq_len(nsim), ceiling(seq_len(nsim) / ceiling(nsim / n_chunks)))

  run_chunk <- function(idx) {
    out <- matrix(NA_real_, nrow = length(idx), ncol = length(look_n))
    for (j in seq_along(idx)) {
      assign(".Random.seed", streams[[idx[j]]], envir = globalenv())
      e_t <- cumsum(stats::rbinom(length(increments), increments, p_treatment))
      e_c <- cumsum(stats::rbinom(length(increments), increments, p_control))
      for (k in seq_along(look_n)) {
        out[j, k] <- posterior_prob_better(
          events_t = e_t[k], n_t = look_n[k],
          events_c = e_c[k], n_c = look_n[k],
          prior = prior, direction = direction,
          draws = post_draws, seed = NULL
        )
      }
    }
    out
  }

  parts <- if (workers > 1L) {
    parallel::mclapply(chunks, run_chunk, mc.cores = workers)
  } else {
    lapply(chunks, run_chunk)
  }
  bad <- vapply(parts, function(p) !is.matrix(p), logical(1))
  if (any(bad)) {
    stop(
      "Simulation failed in ", sum(bad), " of ", length(parts),
      " chunks. First error: ", as.character(parts[[which(bad)[1]]])[1],
      call. = FALSE
    )
  }
  do.call(rbind, parts)
}

# Applies the stopping rule to an already-simulated trajectory matrix.
bayes_decide <- function(ptreat, look_n, efficacy_threshold, futility_threshold) {
  looks <- ncol(ptreat)
  stop_look <- rep(NA_integer_, nrow(ptreat))
  reason <- rep(NA_character_, nrow(ptreat))
  for (k in seq_len(looks)) {
    open <- is.na(stop_look)
    if (!any(open)) break
    eff <- open & ptreat[, k] > efficacy_threshold
    fut <- open & ptreat[, k] < futility_threshold
    stop_look[eff] <- k
    reason[eff] <- "efficacy"
    stop_look[fut] <- k
    reason[fut] <- "futility"
  }
  open <- is.na(stop_look)
  stop_look[open] <- looks
  reason[open] <- "inconclusive"
  list(
    stop_look = stop_look,
    reason = reason,
    n_total = 2L * look_n[stop_look]
  )
}

bayes_summary <- function(dec, looks) {
  nsim <- length(dec$reason)
  reject <- dec$reason == "efficacy"
  p <- mean(reject)
  list(
    reject = p,
    reject_mcse = sqrt(p * (1 - p) / nsim),
    expected_n = mean(dec$n_total),
    expected_n_mcse = stats::sd(dec$n_total) / sqrt(nsim),
    stop_efficacy = vapply(
      seq_len(looks), function(k) mean(reject & dec$stop_look == k), numeric(1)
    ),
    stop_futility = vapply(
      seq_len(looks), function(k) mean(dec$reason == "futility" & dec$stop_look == k), numeric(1)
    ),
    inconclusive = mean(dec$reason == "inconclusive")
  )
}

resolve_workers <- function(workers) {
  workers <- check_count(workers, "workers", min = 1L)
  if (workers > 1L && .Platform$OS.type == "windows") {
    warning(
      "`parallel::mclapply` cannot fork on Windows; falling back to workers = 1. ",
      "Results are unchanged, only slower.",
      call. = FALSE
    )
    workers <- 1L
  }
  workers
}

#' Bayesian adaptive design with simulated operating characteristics
#'
#' Simulates a two-arm trial that looks at accumulating data `looks` times, at
#' equally spaced fractions of `n_max`. At each look the posterior probability
#' that treatment is better is computed by conjugate updating, the trial stops
#' for efficacy if that probability exceeds `efficacy_threshold`, and stops for
#' futility if it falls below `futility_threshold`. Type I error is estimated
#' by rerunning the identical procedure with both arms generated at the control
#' rate.
#'
#' Allocation is 1:1, so the per-arm size is `floor(n_max / 2)`. An odd `n_max`
#' is therefore rounded down by one participant, and the returned `n_total` is
#' the value actually simulated.
#'
#' The returned `alpha` is the *simulated* type I error of this procedure, not
#' a nominal level supplied by the user. Nothing in a Bayesian stopping rule
#' guarantees frequentist error control, so the honest report is what the
#' simulation found. Use [calibrate_bayesian()] to choose a threshold that
#' delivers a target level.
#'
#' Two separate Monte Carlo errors sit in these numbers and they behave
#' differently. `nsim` controls *variance*: how much a reported rate moves from
#' one seed to the next, falling as `1 / sqrt(nsim)`. `post_draws` controls
#' *bias*: the stopping rule compares an estimated posterior probability, not
#' the exact one, against a fixed threshold, and the estimation noise pushes
#' trials across that threshold asymmetrically. The resulting distortion is a
#' property of `post_draws` alone and does not shrink as `nsim` grows, so a
#' large `nsim` on top of a small `post_draws` reports a biased answer with
#' impressive precision. Set `post_draws` for accuracy first, then `nsim` for
#' precision.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param n_max Maximum total sample size across both arms.
#' @param looks Number of analyses, including the final one. `looks = 1` is a
#'   single final analysis with no interim monitoring.
#' @param prior Length-2 numeric vector `c(a0, b0)`, the Beta prior applied
#'   independently to each arm.
#' @param efficacy_threshold Stop and declare efficacy when the posterior
#'   probability that treatment is better exceeds this value.
#' @param futility_threshold Stop for futility when that probability falls
#'   below this value. Must be smaller than `efficacy_threshold`.
#' @param nsim Number of simulated trials per hypothesis. This controls the
#'   variance of the reported rates and nothing else: the Monte Carlo standard
#'   error of every reported rate is at most `0.5 / sqrt(nsim)`. It does not
#'   control the bias described under `post_draws`, which no value of `nsim`
#'   reduces.
#' @param seed Integer seed. Independent L'Ecuyer-CMRG streams are derived from
#'   it, one per simulated trial, so results do not depend on `workers`.
#' @param post_draws Monte Carlo draws used for each posterior probability, at
#'   every look of every simulated trial. This is the accuracy dial of the
#'   design and too small a value corrupts the operating characteristics rather
#'   than merely blurring them. Each posterior probability is a mean of
#'   `post_draws` Bernoulli draws, so its standard error is at most
#'   `0.5 / sqrt(post_draws)`, attained when the true probability is 0.5, which
#'   is where the null distribution of the decision quantity is centred. If that
#'   standard error is not small next to the distance from `efficacy_threshold`
#'   to 0.5, trials from the bulk of the null cross the threshold on estimation
#'   noise alone: at `efficacy_threshold = 0.975` and `post_draws = 2` the
#'   simulated type I error is about 0.33, more than ten times nominal, and it
#'   stays there however large `nsim` is. The function therefore requires that
#'   distance to cover at least ten standard errors, which is
#'   `ceiling((5 / (efficacy_threshold - 0.5))^2)` draws, 111 at the default
#'   threshold and never fewer than 100, and warns below twenty standard errors,
#'   `ceiling((10 / (efficacy_threshold - 0.5))^2)`, which is 444 at the default
#'   threshold. Meeting the minimum is not the same as being accurate: a
#'   secondary effect, that an estimate taking only the values `0/post_draws` to
#'   `post_draws/post_draws` cannot land exactly on the threshold, still moves
#'   the simulated type I error by roughly `1 / ((1 - efficacy_threshold) *
#'   post_draws)` in relative terms, about 9 percent at 444 draws and 1 percent
#'   at the default of 4000. Treat the minimum as the point below which the
#'   answer is wrong, and the default as the point at which it is usable.
#' @param workers Number of forked processes used via [parallel::mclapply()].
#'   Results are identical for any value; only the run time changes.
#' @return An object of class `gi_design` with `type = "bayesian_adaptive"`,
#'   `power` the simulated power under the scenario rates and `alpha` the
#'   simulated type I error under the null. Its `detail` element holds the look
#'   sizes and information rates, the thresholds and prior, the simulated power
#'   and type I error each with a Monte Carlo standard error, expected sample
#'   size under both hypotheses with standard errors, a `stopping` data frame
#'   giving the probability of stopping at each look for efficacy and for
#'   futility under both hypotheses, and the simulation settings including the
#'   seed and RNG kind.
#' @seealso [calibrate_bayesian()], [posterior_prob_better()]
#' @examples
#' sc <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
#' d <- design_bayesian(sc, n_max = 200, looks = 2, nsim = 200, post_draws = 500, seed = 7)
#' d$power
#' d$detail$stopping
#' @export
design_bayesian <- function(scenario, n_max, looks = 3, prior = c(1, 1),
                            efficacy_threshold = 0.975,
                            futility_threshold = 0.10,
                            nsim = 5000, seed = 1, post_draws = 4000,
                            workers = 1) {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`scenario` must be a gi_scenario, as returned by scenario().", call. = FALSE)
  }
  n_max <- check_count(n_max, "n_max", min = 4L)
  looks <- check_count(looks, "looks", min = 1L)
  prior <- check_prior(prior)
  efficacy_threshold <- check_unit(efficacy_threshold, "efficacy_threshold")
  futility_threshold <- check_unit(futility_threshold, "futility_threshold")
  if (futility_threshold >= efficacy_threshold) {
    stop(
      "`futility_threshold` (", futility_threshold, ") must be below ",
      "`efficacy_threshold` (", efficacy_threshold, ").",
      call. = FALSE
    )
  }
  nsim <- check_count(nsim, "nsim", min = 2L)
  post_draws <- check_post_draws(post_draws, efficacy_threshold)
  seed <- check_count(seed, "seed", min = -.Machine$integer.max)
  workers <- resolve_workers(workers)
  direction <- check_direction(scenario$direction)

  n_per_arm <- n_max %/% 2L
  look_n <- bayes_look_sizes(n_per_arm, looks)
  streams <- bayes_stream_sets(seed, nsim)

  snap <- rng_snapshot()
  on.exit(rng_restore(snap), add = TRUE)

  traj_alt <- bayes_trajectories(
    p_control = scenario$control_rate, p_treatment = scenario$treatment_rate,
    look_n = look_n, prior = prior, direction = direction,
    post_draws = post_draws, streams = streams$alt, workers = workers
  )
  traj_null <- bayes_trajectories(
    p_control = scenario$control_rate, p_treatment = scenario$control_rate,
    look_n = look_n, prior = prior, direction = direction,
    post_draws = post_draws, streams = streams$null, workers = workers
  )

  alt <- bayes_summary(
    bayes_decide(traj_alt, look_n, efficacy_threshold, futility_threshold), looks
  )
  null <- bayes_summary(
    bayes_decide(traj_null, look_n, efficacy_threshold, futility_threshold), looks
  )

  stopping <- data.frame(
    look = seq_len(looks),
    n_total = 2L * look_n,
    information_rate = look_n / n_per_arm,
    efficacy_alt = alt$stop_efficacy,
    futility_alt = alt$stop_futility,
    efficacy_null = null$stop_efficacy,
    futility_null = null$stop_futility,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      type = "bayesian_adaptive",
      scenario = scenario,
      alpha = null$reject,
      power = alt$reject,
      n_total = 2L * n_per_arm,
      n_per_arm = n_per_arm,
      engine = "gitrialsim conjugate beta-binomial (simulation)",
      detail = list(
        looks = looks,
        look_n_per_arm = look_n,
        look_n_total = 2L * look_n,
        information_rates = look_n / n_per_arm,
        prior = prior,
        efficacy_threshold = efficacy_threshold,
        futility_threshold = futility_threshold,
        power = alt$reject,
        power_mcse = alt$reject_mcse,
        alpha_null = null$reject,
        alpha_mcse = null$reject_mcse,
        expected_n_alt = alt$expected_n,
        expected_n_alt_mcse = alt$expected_n_mcse,
        expected_n_null = null$expected_n,
        expected_n_null_mcse = null$expected_n_mcse,
        inconclusive_alt = alt$inconclusive,
        inconclusive_null = null$inconclusive,
        stopping = stopping,
        alt_rates = c(control = scenario$control_rate, treatment = scenario$treatment_rate),
        null_rates = c(control = scenario$control_rate, treatment = scenario$control_rate),
        direction = direction,
        nsim = nsim,
        post_draws = post_draws,
        seed = seed,
        workers = workers,
        rng_kind = gi_rng_kind,
        calibration = NULL
      )
    ),
    class = c("gi_design", "list")
  )
}

#' Calibrate the efficacy threshold to a target type I error
#'
#' Searches for the posterior probability threshold at which the simulated type
#' I error of [design_bayesian()] meets `target_alpha`. The null trajectories
#' are simulated once and rescored at every candidate threshold, so the search
#' compares like with like: using common random numbers makes simulated alpha a
#' monotone non-increasing step function of the threshold, which is what allows
#' plain bisection to be used.
#'
#' The search returns the smallest bracketed threshold whose simulated type I
#' error does not exceed the target, so the achieved alpha is at or just below
#' `target_alpha` rather than straddling it.
#'
#' Honest statement of what this gives you: the calibrated threshold is itself
#' a Monte Carlo estimate, not an exact constant. It is calibrated against one
#' particular set of `nsim` simulated null trials, and a different seed will
#' return a slightly different threshold. The achieved alpha is reported with
#' its Monte Carlo standard error, `sqrt(alpha * (1 - alpha) / nsim)`, which
#' falls only as `1 / sqrt(nsim)`: at `nsim = 5000` and a target of 0.025 the
#' standard error is about 0.0022, and reaching 0.0011 needs `nsim = 20000`.
#' Quadrupling `nsim` halves the uncertainty. The uncertainty in the threshold
#' itself is that error divided by the local slope of alpha against threshold,
#' so a design whose alpha changes slowly with the threshold will have a
#' correspondingly less precisely determined threshold. For a submission,
#' recalibrate at large `nsim` and confirm the threshold on an independent
#' seed.
#'
#' `nsim` is not the only thing to raise. The alpha the search is bisecting is
#' itself biased when `post_draws` is small, because the stopping rule compares
#' an estimated posterior probability against the candidate threshold, so the
#' calibration then targets the wrong quantity and no amount of `nsim` corrects
#' it. The requirement and the warning described under the `post_draws`
#' argument of [design_bayesian()] are applied here as well, against
#' `1 - target_alpha`, which is the least favourable threshold the search can
#' return.
#'
#' @param scenario A `gi_scenario`, as returned by [scenario()].
#' @param n_max Maximum total sample size across both arms.
#' @param target_alpha Target one-sided simulated type I error.
#' @param looks Number of analyses, including the final one.
#' @param prior Length-2 numeric vector `c(a0, b0)`, the Beta prior applied
#'   independently to each arm.
#' @param futility_threshold Futility stopping threshold, held fixed during the
#'   search.
#' @param nsim Number of simulated trials per hypothesis. Controls the variance
#'   of the achieved alpha and hence the precision of the calibrated threshold,
#'   not the bias described under `post_draws`.
#' @param seed Integer seed, used for both the calibration run and the final
#'   design, so the reported alpha of the design equals the achieved alpha of
#'   the search.
#' @param post_draws Monte Carlo draws used for each posterior probability. Too
#'   small a value biases every alpha the search evaluates and therefore the
#'   calibrated threshold itself, a bias `nsim` cannot reduce. The minimum and
#'   the recommended value, derived in the `post_draws` entry of
#'   [design_bayesian()], are enforced here against `1 - target_alpha`.
#' @param workers Number of forked processes used via [parallel::mclapply()].
#' @param tol Width of the threshold bracket at which bisection stops.
#' @param max_iter Maximum number of bisection steps.
#' @return A `gi_design` exactly as returned by [design_bayesian()] at the
#'   calibrated threshold, with `detail$calibration` added: the target and
#'   achieved alpha, the Monte Carlo standard error of the achieved alpha, the
#'   calibrated threshold, the final bracket, the number of iterations, and
#'   whether the target was bracketed at all.
#' @seealso [design_bayesian()]
#' @examples
#' sc <- scenario("ercp_acute_cholangitis", control_rate = 0.30, treatment_rate = 0.15)
#' d <- calibrate_bayesian(sc, n_max = 200, looks = 1, nsim = 300, post_draws = 500, seed = 3)
#' d$detail$calibration$threshold
#' d$detail$calibration$achieved_alpha
#' @export
calibrate_bayesian <- function(scenario, n_max, target_alpha = 0.025, looks = 3,
                               prior = c(1, 1), futility_threshold = 0.10,
                               nsim = 5000, seed = 1, post_draws = 4000,
                               workers = 1, tol = 1e-4, max_iter = 40) {
  if (!inherits(scenario, "gi_scenario")) {
    stop("`scenario` must be a gi_scenario, as returned by scenario().", call. = FALSE)
  }
  n_max <- check_count(n_max, "n_max", min = 4L)
  looks <- check_count(looks, "looks", min = 1L)
  prior <- check_prior(prior)
  target_alpha <- check_unit(target_alpha, "target_alpha")
  futility_threshold <- check_unit(futility_threshold, "futility_threshold")
  nsim <- check_count(nsim, "nsim", min = 2L)
  # The threshold the search will land on is not known yet, but it cannot come
  # in below 1 - target_alpha for a single look and only rises with the number
  # of looks. That is therefore the least favourable threshold `post_draws` has
  # to resolve, and checking against it here fails fast rather than after the
  # null trajectories have been simulated.
  post_draws <- check_post_draws(post_draws, 1 - target_alpha)
  seed <- check_count(seed, "seed", min = -.Machine$integer.max)
  workers <- resolve_workers(workers)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }
  max_iter <- check_count(max_iter, "max_iter", min = 1L)
  direction <- check_direction(scenario$direction)

  lo <- 0.5
  hi <- 0.9999
  if (futility_threshold >= lo) {
    stop(
      "`futility_threshold` (", futility_threshold,
      ") must be below the lower end of the efficacy search bracket (", lo, ").",
      call. = FALSE
    )
  }

  n_per_arm <- n_max %/% 2L
  look_n <- bayes_look_sizes(n_per_arm, looks)
  streams <- bayes_stream_sets(seed, nsim)

  snap <- rng_snapshot()
  on.exit(rng_restore(snap), add = TRUE)

  traj_null <- bayes_trajectories(
    p_control = scenario$control_rate, p_treatment = scenario$control_rate,
    look_n = look_n, prior = prior, direction = direction,
    post_draws = post_draws, streams = streams$null, workers = workers
  )
  alpha_at <- function(thr) {
    mean(bayes_decide(traj_null, look_n, thr, futility_threshold)$reason == "efficacy")
  }

  alpha_lo <- alpha_at(lo)
  alpha_hi <- alpha_at(hi)
  iterations <- 0L
  bracketed <- TRUE

  if (alpha_lo <= target_alpha) {
    bracketed <- FALSE
    threshold <- lo
    warning(
      "Simulated type I error is already ", signif(alpha_lo, 3), " at the lowest ",
      "threshold searched (", lo, "), at or below the target of ", target_alpha,
      ". Returning that threshold; the design is more conservative than requested.",
      call. = FALSE
    )
  } else if (alpha_hi > target_alpha) {
    bracketed <- FALSE
    threshold <- hi
    warning(
      "Simulated type I error is still ", signif(alpha_hi, 3), " at the highest ",
      "threshold searched (", hi, "), above the target of ", target_alpha,
      ". Returning that threshold; reduce the number of looks or raise n_max.",
      call. = FALSE
    )
  } else {
    while (hi - lo > tol && iterations < max_iter) {
      mid <- (lo + hi) / 2
      if (alpha_at(mid) > target_alpha) lo <- mid else hi <- mid
      iterations <- iterations + 1L
    }
    threshold <- hi
  }

  achieved <- alpha_at(threshold)
  design <- design_bayesian(
    scenario = scenario, n_max = n_max, looks = looks, prior = prior,
    efficacy_threshold = threshold, futility_threshold = futility_threshold,
    nsim = nsim, seed = seed, post_draws = post_draws, workers = workers
  )
  design$detail$calibration <- list(
    target_alpha = target_alpha,
    achieved_alpha = achieved,
    achieved_alpha_mcse = sqrt(achieved * (1 - achieved) / nsim),
    threshold = threshold,
    bracket = c(lower = lo, upper = hi),
    iterations = iterations,
    bracketed = bracketed,
    tol = tol,
    nsim = nsim,
    seed = seed,
    note = paste(
      "Calibrated by bisection on nsim simulated null trials under common",
      "random numbers. The threshold is a Monte Carlo estimate; its precision",
      "scales as 1/sqrt(nsim)."
    )
  )
  design
}

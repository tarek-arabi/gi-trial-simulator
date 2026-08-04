test_that("rpact, the closed form and stats::power.prop.test agree", {
  # rpact::getSampleSizeRates implements the pooled-variance normal
  # approximation without continuity correction. Pinning it against two
  # independent implementations, one of them base R, is what lets a discrepancy
  # with a published number be attributed to the publication rather than here.
  grid <- expand.grid(
    p_control = c(0.05, 0.10, 0.30, 0.50),
    effect = c(0.3, 0.5),
    power = c(0.80, 0.90)
  )
  for (i in seq_len(nrow(grid))) {
    pc <- grid$p_control[i]
    pt <- pc * (1 - grid$effect[i])
    args <- list(pc, pt, 0.025, grid$power[i])
    n_rpact <- do.call(gi_n_reference, c(args, method = "rpact"))
    n_pooled <- do.call(gi_n_reference, c(args, method = "pooled"))
    n_base <- do.call(gi_n_reference, c(args, method = "power.prop.test"))
    expect_equal(n_rpact, n_pooled, tolerance = 1e-6)
    expect_equal(n_rpact, n_base, tolerance = 1e-6)
  }
})

test_that("a two-sided pack alpha builds the same design as its one-sided half", {
  # Published trials almost always state a two-sided alpha. For a two-arm
  # superiority test the two conventions coincide, since z(1 - 0.05/2) equals
  # z(1 - 0.025), so a pack must be able to record what its source said without
  # changing any number.
  pack <- load_pack("ercp_acute_cholangitis")
  one_sided <- scenario(pack)

  two_sided_pack <- pack
  two_sided_pack$design_defaults$sided <- 2
  two_sided_pack$design_defaults$alpha <- 0.05
  two_sided <- scenario(two_sided_pack)

  d1 <- design_fixed(one_sided)
  d2 <- design_fixed(two_sided)

  expect_equal(d2$alpha, 0.025)
  expect_equal(d1$n_total, d2$n_total)
  expect_equal(d1$detail$critical_value, d2$detail$critical_value)
  expect_equal(d2$detail$sided_source, 2)
  expect_equal(d1$detail$sided_source, 1)
})

test_that("an explicitly supplied alpha is one-sided and is never halved", {
  pack <- load_pack("ercp_acute_cholangitis")
  pack$design_defaults$sided <- 2
  pack$design_defaults$alpha <- 0.05
  s <- scenario(pack)
  expect_equal(design_fixed(s, alpha = 0.025)$alpha, 0.025)
})

test_that("a pack `sided` other than 1 or 2 is rejected", {
  pack <- load_pack("ercp_acute_cholangitis")
  pack$design_defaults$sided <- 3
  expect_error(design_fixed(scenario(pack)), "must be 1 or 2")
})

test_that("Elmunzer 2012 reconciles to the continuity-corrected formula", {
  # Elmunzer BJ et al. A randomized trial of rectal indomethacin to prevent
  # post-ERCP pancreatitis. N Engl J Med 2012;366(15):1414-22. PMID 22494121.
  # "We estimated that 948 patients (474 per study group) would provide a power
  # of at least 80% to detect a 50% reduction in the incidence of post-ERCP
  # pancreatitis, from 10% in the placebo group to 5% in the indomethacin
  # group, on the basis of Fisher's exact test, with a two-sided significance
  # level of 0.05."
  r <- reconcile_sample_size(
    published_n_per_arm = 474,
    p_control = 0.10, p_treat = 0.05,
    alpha_1sided = 0.025, power = 0.80
  )

  expect_identical(attr(r, "attribution"), "continuity")
  expect_true(attr(r, "reproduced"))

  cc <- r$n_per_arm[r$method == "continuity"]
  expect_lt(abs(cc - 474), 1)

  # The engine itself sits about 8% below, which is the whole point: without
  # attribution this reads as a validation failure.
  rp <- r$pct_vs_published[r$method == "rpact"]
  expect_lt(rp, -0.05)
  expect_gt(rp, -0.12)

  # The published n delivers slightly more than the stated 80%, as a
  # conservative continuity correction implies.
  expect_gt(attr(r, "power_at_published"), 0.80)
})

test_that("Levenick 2016 does not reconcile, and is under-sized as published", {
  # Levenick JM et al. Rectal indomethacin does not prevent post-ERCP
  # pancreatitis in consecutive patients. Gastroenterology 2016;150(4):911-917.
  # PMID 26775631. "We estimated that 1,398 patients (699 per study group)
  # would provide a power of 80% to detect a 50% reduction in the rate of PEP
  # from 5% in the placebo group to 2.5% in the indomethacin group using the
  # two-tailed Fisher's Exact test with a two-sided significance of 0.05."
  r <- reconcile_sample_size(
    published_n_per_arm = 699,
    p_control = 0.05, p_treat = 0.025,
    alpha_1sided = 0.025, power = 0.80
  )

  # No method reproduces 699 at the stated two-sided 0.05. Fisher's exact is
  # conservative and requires more, not fewer, so the stated method cannot
  # explain a number that is too small.
  expect_false(attr(r, "reproduced"))
  expect_true(is.na(attr(r, "attribution")))
  expect_true(all(r$pct_vs_published > 0.20))

  # The published n delivers about 69% power, not the stated 80%.
  expect_lt(attr(r, "power_at_published"), 0.72)
  expect_gt(attr(r, "power_at_published"), 0.66)

  # It is however consistent with a one-sided 0.05 calculation, which is the
  # diagnosis the reconciliation exists to produce. Note the claim is about the
  # alpha convention, not about which formula: at this alpha the candidate
  # methods span 694 to 714, so the published 699 cannot identify one of them.
  one_sided <- reconcile_sample_size(
    published_n_per_arm = 699,
    p_control = 0.05, p_treat = 0.025,
    alpha_1sided = 0.05, power = 0.80
  )
  expect_true(attr(one_sided, "reproduced"))
  expect_true(abs(one_sided$pct_vs_published[1]) < 0.02)
})

test_that("reconcile_report reduces a trial set to one row each", {
  trials <- data.frame(
    label = c("Elmunzer 2012", "Levenick 2016"),
    published_n_per_arm = c(474, 699),
    p_control = c(0.10, 0.05),
    p_treat = c(0.05, 0.025),
    alpha_1sided = c(0.025, 0.025),
    power = c(0.80, 0.80),
    stated_method = c("Fisher's exact", "two-tailed Fisher's exact"),
    stringsAsFactors = FALSE
  )
  out <- reconcile_report(trials)

  expect_equal(nrow(out), 2)
  expect_equal(out$label, c("Elmunzer 2012", "Levenick 2016"))
  expect_equal(out$reproduced, c(TRUE, FALSE))
  expect_identical(out$attribution[1], "continuity")
  expect_true(is.na(out$attribution[2]))
  expect_true(all(out$power_stated == 0.80))
})

test_that("reconcile_report rejects a malformed trial set", {
  expect_error(reconcile_report(data.frame(label = "x")), "missing column")
  expect_error(
    reconcile_report(data.frame(
      label = character(0), published_n_per_arm = numeric(0),
      p_control = numeric(0), p_treat = numeric(0),
      alpha_1sided = numeric(0), power = numeric(0)
    )),
    "no rows"
  )
})

test_that("gi_n_reference validates its inputs", {
  expect_error(gi_n_reference(0, 0.05, 0.025, 0.8), "p_control")
  expect_error(gi_n_reference(0.1, 1, 0.025, 0.8), "p_treat")
  expect_error(gi_n_reference(0.1, 0.05, 0.6, 0.8), "alpha_1sided")
  expect_error(gi_n_reference(0.1, 0.05, 0.025, 1.2), "power")
  expect_error(gi_n_reference(0.1, 0.1, 0.025, 0.8), "must differ")
})

test_that("Schoenfeld equals rpact's required events across hazard ratios", {
  # rpact::getSampleSizeSurvival implements Schoenfeld. Pinning it means a
  # disagreement with a published survival calculation can be attributed to the
  # trial's choice of formula rather than left unexplained, exactly as the
  # continuity correction does in the binary stratum.
  for (hr in c(0.5, 0.6, 0.667, 0.7, 0.75, 0.8)) {
    for (pw in c(0.80, 0.90)) {
      expect_equal(
        gi_events_reference(hr, 0.025, pw, method = "rpact"),
        gi_events_reference(hr, 0.025, pw, method = "schoenfeld"),
        tolerance = 1e-6
      )
    }
  }
})

test_that("Freedman sits above Schoenfeld, by more for larger effects", {
  gap <- function(hr) {
    f <- gi_events_reference(hr, 0.025, 0.8, method = "freedman")
    s <- gi_events_reference(hr, 0.025, 0.8, method = "schoenfeld")
    f / s - 1
  }
  expect_gt(gap(0.5), gap(0.8))
  expect_gt(gap(0.5), 0.05)
  expect_lt(gap(0.8), 0.02)
  expect_true(all(vapply(c(0.5, 0.6, 0.7, 0.8), gap, numeric(1)) > 0))
})

test_that("reconcile_events attributes a Freedman-sized figure to Freedman", {
  # A trial that used Freedman at HR 0.5, 80% power, one-sided 0.025 needs 71
  # events where rpact would say 65. Without attribution that reads as a 8%
  # failure; with it, the discrepancy has a name.
  fre <- gi_events_reference(0.5, 0.025, 0.8, method = "freedman")
  r <- reconcile_events(round(fre), hazard_ratio = 0.5, alpha_1sided = 0.025, power = 0.8)
  expect_identical(attr(r, "attribution"), "freedman")
  expect_true(attr(r, "reproduced"))

  sch <- gi_events_reference(0.5, 0.025, 0.8, method = "schoenfeld")
  r2 <- reconcile_events(round(sch), hazard_ratio = 0.5, alpha_1sided = 0.025, power = 0.8)
  expect_true(attr(r2, "attribution") %in% c("rpact", "schoenfeld"))
})

test_that("gi_events_reference validates its inputs", {
  expect_error(gi_events_reference(0, 0.025, 0.8), "hazard_ratio")
  expect_error(gi_events_reference(1, 0.025, 0.8), "no effect")
  expect_error(gi_events_reference(0.7, 0.6, 0.8), "alpha_1sided")
  expect_error(gi_events_reference(0.7, 0.025, 1.5), "power")
  expect_error(reconcile_events(0, 0.7, 0.025, 0.8), "published_events")
})

test_that("the arcsine reference matches pwr, which implements it independently", {
  # Raised during an adversarial review, where a reviewer claimed this formula
  # was out by a factor of two. It is not, but the claim was only refutable by
  # recomputation, so the check is pinned here rather than argued. pwr
  # solves power = pnorm(h * sqrt(n / 2) - z), with Cohen's h = 2asin(sqrt(p1)) -
  # 2asin(sqrt(p2)); the sqrt(n / 2) is the factor the objection missed.
  skip_if_not_installed("pwr")
  grid <- expand.grid(
    p1 = c(0.05, 0.10, 0.25, 0.50, 0.71),
    p2 = c(0.02, 0.15, 0.35, 0.83),
    power = c(0.80, 0.90),
    alpha1 = c(0.025, 0.05)
  )
  grid <- grid[abs(grid$p1 - grid$p2) > 0.01, ]
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    ours <- gi_n_reference(g$p1, g$p2, g$alpha1, g$power, "arcsine")
    theirs <- pwr::pwr.2p.test(
      h = pwr::ES.h(g$p1, g$p2), sig.level = 2 * g$alpha1, power = g$power
    )$n
    # 1e-4 relative is set by pwr's precision, not by what makes this pass:
    # pwr solves for n by uniroot, so it lands within about 1e-5 relative of the
    # closed form. The error this test exists to catch would be a factor of two.
    expect_equal(ours, theirs, tolerance = 1e-4,
      info = sprintf("p1=%.2f p2=%.2f alpha1=%.3f power=%.2f",
                     g$p1, g$p2, g$alpha1, g$power))
  }
})

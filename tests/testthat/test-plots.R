fake_grid <- function() {
  data.frame(
    n_per_arm = rep(c(500, 1000, 1500, 2000), 2),
    power = c(0.31, 0.55, 0.71, 0.82, 0.24, 0.44, 0.60, 0.72),
    power_mcse = rep(0.015, 8),
    design = rep(c("Fixed", "Group sequential"), each = 4),
    stringsAsFactors = FALSE
  )
}

fake_design <- function(type, power, n_total, detail = list()) {
  structure(
    list(
      type = type,
      scenario = NULL,
      alpha = 0.025,
      power = power,
      n_total = n_total,
      n_per_arm = n_total / 2,
      engine = "test fixture",
      detail = detail
    ),
    class = c("gi_design", "list")
  )
}

fake_boundaries <- function(futility = TRUE) {
  out <- data.frame(
    analysis = 1:3,
    information_rate = c(1 / 3, 2 / 3, 1),
    efficacy = c(3.71, 2.51, 1.99)
  )
  if (futility) out$futility <- c(-0.31, 0.98, 1.99)
  out
}

fake_emulator <- function() {
  list(
    predict = function(newdata) {
      0.9 * (1 - exp(-newdata$n_per_arm / 900)) - 4 * (newdata$control_rate - 0.07)^2
    },
    design = data.frame(
      n_per_arm = c(300, 800, 1200, 1800),
      control_rate = c(0.04, 0.09, 0.055, 0.08)
    )
  )
}

expect_untitled_ggplot <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_null(p$labels$title)
  expect_null(p$labels$subtitle)
  expect_no_error(ggplot2::ggplot_build(p))
}

test_that("plot_power_curve returns an untitled ggplot, grouped and ungrouped", {
  skip_if_not_installed("ggplot2")
  grid <- fake_grid()

  expect_untitled_ggplot(plot_power_curve(grid, x = "n_per_arm", y = "power"))
  p <- plot_power_curve(grid, x = "n_per_arm", y = "power", group = "design")
  expect_untitled_ggplot(p)
  expect_identical(p$labels$x, "n_per_arm")
  expect_identical(p$labels$y, "power")
  expect_identical(p$labels$colour, "design")
})

test_that("plot_power_curve draws a Monte Carlo band only when an MCSE column exists", {
  skip_if_not_installed("ggplot2")
  grid <- fake_grid()
  with_mcse <- plot_power_curve(grid, x = "n_per_arm", y = "power")
  without <- plot_power_curve(grid[setdiff(names(grid), "power_mcse")],
    x = "n_per_arm", y = "power"
  )
  expect_length(with_mcse$layers, 3L)
  expect_length(without$layers, 2L)
})

test_that("plot_power_curve errors informatively on a missing column", {
  skip_if_not_installed("ggplot2")
  grid <- fake_grid()

  expect_error(
    plot_power_curve(grid, x = "sample_size", y = "power"),
    "sample_size"
  )
  expect_error(
    plot_power_curve(grid, x = "n_per_arm", y = "rejection_rate"),
    "rejection_rate"
  )
  expect_error(
    plot_power_curve(grid, x = "n_per_arm", y = "power", group = "arm"),
    "'group'"
  )
  expect_error(
    plot_power_curve(grid, x = "design", y = "power"),
    "not numeric"
  )
  expect_error(plot_power_curve(list(a = 1), x = "a", y = "a"), "data.frame")
  expect_error(plot_power_curve(fake_grid()[0, ], x = "n_per_arm", y = "power"), "no rows")
})

test_that("plot_operating_characteristics handles NA target power on adaptive designs", {
  skip_if_not_installed("ggplot2")
  designs <- list(
    Fixed = fake_design("fixed", 0.9, 3400),
    `Group sequential` = fake_design(
      "group_sequential", 0.9, 3600,
      detail = list(expected_n = 2700)
    ),
    Adaptive = fake_design(
      "bayesian_adaptive", NA_real_, 3400,
      detail = list(simulated_power = 0.88, expected_n = 2450)
    )
  )
  p <- plot_operating_characteristics(designs)
  expect_untitled_ggplot(p)

  values <- p$data
  adaptive_power <- values$gi_value[values$gi_label == "Adaptive" &
    values$gi_measure == "Power"]
  expect_equal(adaptive_power, 0.88)
  expect_false(any(is.na(values$gi_value)))
  expect_true(all(c("Analytic", "Simulated") %in% levels(values$gi_kind)))
})

test_that("plot_operating_characteristics credits a fixed design with its full n", {
  skip_if_not_installed("ggplot2")
  p <- plot_operating_characteristics(list(Fixed = fake_design("fixed", 0.9, 3400)))
  expect_untitled_ggplot(p)
  expected <- p$data$gi_value[p$data$gi_measure == "Expected n"]
  expect_equal(expected, 3400)
})

test_that("plot_operating_characteristics drops an empty measure and validates input", {
  skip_if_not_installed("ggplot2")
  p <- plot_operating_characteristics(list(fake_design("group_sequential", 0.9, 3400)))
  expect_untitled_ggplot(p)
  expect_false("Expected n" %in% levels(p$data$gi_measure))

  expect_error(plot_operating_characteristics(list()), "non-empty list")
  expect_error(plot_operating_characteristics(list(1, 2)), "gi_design")
})

test_that("plot_boundaries plots efficacy with and without a futility column", {
  skip_if_not_installed("ggplot2")
  p <- plot_boundaries(fake_boundaries())
  expect_untitled_ggplot(p)
  expect_identical(levels(p$data$gi_bound), c("Efficacy", "Futility"))
  expect_identical(p$labels$x, "Information fraction")
  expect_identical(p$labels$y, "Z statistic")

  q <- plot_boundaries(fake_boundaries(futility = FALSE))
  expect_untitled_ggplot(q)
  expect_identical(levels(q$data$gi_bound), "Efficacy")
})

test_that("plot_boundaries tolerates alternative column names and non-finite bounds", {
  skip_if_not_installed("ggplot2")
  alt <- data.frame(
    timing = c(0.5, 1),
    upper = c(2.96, 1.97),
    lower = c(-Inf, 1.97)
  )
  p <- plot_boundaries(alt)
  expect_untitled_ggplot(p)
  expect_equal(nrow(p$data), 3L)
})

test_that("plot_boundaries rejects designs it cannot draw", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_boundaries(fake_design("fixed", 0.9, 3400)), "group_sequential")
  expect_error(plot_boundaries("not a design"), "gi_design")
  expect_error(
    plot_boundaries(data.frame(analysis = 1:2, z = c(2.9, 2.0))),
    "information fraction"
  )
})

test_that("plot_evsi draws a ribbon from an SE column and copes without one", {
  skip_if_not_installed("ggplot2")
  curve <- data.frame(
    n_per_arm = seq(200, 2000, by = 200),
    evsi = c(120, 210, 275, 320, 350, 368, 378, 383, 385, 386) * 1000,
    se = rep(9000, 10)
  )
  p <- plot_evsi(list(curve = curve))
  expect_untitled_ggplot(p)
  expect_length(p$layers, 2L)
  expect_identical(p$labels$x, "Sample size per arm")

  bare <- plot_evsi(curve[c("n_per_arm", "evsi")])
  expect_untitled_ggplot(bare)
  expect_length(bare$layers, 1L)

  explicit <- curve[c("n_per_arm", "evsi")]
  explicit$ci_lower <- explicit$evsi - 1e4
  explicit$ci_upper <- explicit$evsi + 1e4
  expect_untitled_ggplot(plot_evsi(explicit))
})

test_that("plot_evsi errors when the curve has no usable columns", {
  skip_if_not_installed("ggplot2")
  expect_error(plot_evsi(data.frame(x = 1:3, y = 1:3)), "sample size column")
  expect_error(plot_evsi(data.frame(n = 1:3, y = 1:3)), "EVSI column")
  expect_error(plot_evsi("nope"), "data.frame")
})

test_that("plot_emulator_surface draws a supported surface with training points", {
  skip_if_not_installed("ggplot2")
  fit <- fake_emulator()
  bounds <- list(n_per_arm = c(200, 2000), control_rate = c(0.03, 0.10))
  p <- plot_emulator_surface(fit, bounds, n_grid = 20)
  expect_untitled_ggplot(p)
  expect_equal(nrow(p$data), 400L)
  expect_identical(p$labels$x, "n_per_arm")
  expect_identical(p$labels$fill, "Posterior mean")
  expect_equal(nrow(p$layers[[3]]$data), 4L)

  plain <- plot_emulator_surface(fit$predict, bounds, n_grid = 12)
  expect_untitled_ggplot(plain)
  expect_length(plain$layers, 2L)
})

test_that("plot_emulator_surface validates bounds and grid size", {
  skip_if_not_installed("ggplot2")
  fit <- fake_emulator()
  good <- list(n_per_arm = c(200, 2000), control_rate = c(0.03, 0.10))

  expect_error(plot_emulator_surface(fit, good[1]), "exactly two")
  expect_error(
    plot_emulator_surface(fit, list(c(1, 2), c(3, 4))),
    "two distinct names"
  )
  expect_error(
    plot_emulator_surface(fit, list(n_per_arm = c(2000, 200), control_rate = c(0.03, 0.1))),
    "lower then upper"
  )
  expect_error(plot_emulator_surface(fit, good, n_grid = 1), "'n_grid'")
  expect_error(plot_emulator_surface(NULL, good), "not NULL")
  expect_error(
    plot_emulator_surface(function(newdata) 1, good, n_grid = 5),
    "predictions for 25 grid points"
  )
})

test_that("plot_boundaries reads a real group-sequential design", {
  skip_if_not_installed("ggplot2")
  skip_if(is.null(get0("design_group_sequential", mode = "function")))
  skip_if(is.null(get0("gs_boundaries", mode = "function")))

  design <- design_group_sequential(
    scenario("ercp_acute_cholangitis"),
    k = 3, futility = "nonbinding_obf"
  )
  p <- plot_boundaries(design)
  expect_untitled_ggplot(p)
  expect_identical(levels(p$data$gi_bound), c("Efficacy", "Futility"))
  # The final analysis has no futility bound, so it contributes no point.
  expect_equal(sum(p$data$gi_bound == "Efficacy"), 3L)
  expect_equal(sum(p$data$gi_bound == "Futility"), 2L)
})

test_that("plot_emulator_surface reads a real fitted emulator", {
  skip_if_not_installed("ggplot2")
  skip_if(is.null(get0("fit_emulator", mode = "function")))

  set.seed(11)
  x <- data.frame(
    n_per_arm = stats::runif(20, 200, 2000),
    control_rate = stats::runif(20, 0.03, 0.10)
  )
  y <- 0.9 * (1 - exp(-x$n_per_arm / 900)) - 4 * (x$control_rate - 0.07)^2
  fit <- fit_emulator(x, y)

  p <- plot_emulator_surface(
    fit,
    bounds = list(n_per_arm = c(200, 2000), control_rate = c(0.03, 0.10)),
    n_grid = 15
  )
  expect_untitled_ggplot(p)
  expect_equal(nrow(p$data), 225L)
  expect_equal(nrow(p$layers[[3]]$data), 20L)
})

test_that("gi_theme returns a theme and no figure carries a title", {
  skip_if_not_installed("ggplot2")
  th <- gi_theme()
  expect_s3_class(th, "theme")
  expect_s3_class(gi_theme(base_size = 9), "theme")
  expect_error(gi_theme(base_size = -1), "base_size")

  plots <- list(
    plot_power_curve(fake_grid(), x = "n_per_arm", y = "power", group = "design"),
    plot_operating_characteristics(list(fake_design("fixed", 0.9, 3400))),
    plot_boundaries(fake_boundaries()),
    plot_evsi(data.frame(n_total = c(500, 1000), evsi = c(1e5, 2e5))),
    plot_emulator_surface(
      fake_emulator(),
      list(n_per_arm = c(200, 2000), control_rate = c(0.03, 0.10)),
      n_grid = 10
    )
  )
  for (p in plots) {
    expect_null(p$labels$title)
    expect_null(p$labels$subtitle)
  }
})

test_that("every figure function refuses to run without ggplot2", {
  exported <- c(
    "plot_power_curve", "plot_operating_characteristics", "plot_boundaries",
    "plot_evsi", "plot_emulator_surface", "gi_theme", "save_figure"
  )
  for (nm in exported) {
    first_call <- paste(deparse(body(get(nm))[[2L]]), collapse = " ")
    expect_match(first_call, "fig_require_ggplot2", info = nm)
  }

  without_ggplot2 <- fig_require_ggplot2
  env <- new.env(parent = environment(fig_require_ggplot2))
  env$requireNamespace <- function(...) FALSE
  environment(without_ggplot2) <- env
  expect_error(without_ggplot2("plot_power_curve()"), "install.packages")
  expect_error(without_ggplot2("plot_power_curve()"), "plot_power_curve\\(\\)")
})

test_that("save_figure writes a file at the requested size", {
  skip_if_not_installed("ggplot2")
  p <- plot_power_curve(fake_grid(), x = "n_per_arm", y = "power")
  path <- file.path(tempdir(), "gitrialsim-tests", "power.png")
  on.exit(unlink(dirname(path), recursive = TRUE), add = TRUE)

  out <- save_figure(p, path, width = 4, height = 3, dpi = 150)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  expect_error(save_figure("not a plot", path), "ggplot object")
  expect_error(save_figure(p, character()), "'path'")
  expect_error(save_figure(p, path, width = 0), "'width'")
})

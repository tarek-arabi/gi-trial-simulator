#' Publication figures
#'
#' Figures for manuscripts reporting designs built with this package. Every
#' function returns a ggplot object and writes nothing; only [save_figure()]
#' touches the file system. No function sets a plot title or subtitle: under the
#' figure conventions this project follows, the descriptive caption belongs in
#' the figure legends text, not baked into the image.
#'
#' ggplot2 is a suggested dependency, so each function checks for it and stops
#' with an install hint if it is absent.
#'
#' @name figures
NULL

utils::globalVariables(c(
  "gi_x", "gi_y", "gi_z", "gi_group", "gi_lower", "gi_upper",
  "gi_label", "gi_value", "gi_measure", "gi_kind", "gi_bound"
))

fig_require_ggplot2 <- function(fun) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      fun, " requires the 'ggplot2' package, which is not installed. ",
      "ggplot2 is a suggested dependency of gitrialsim: ",
      "install it with install.packages(\"ggplot2\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Colourblind-safe qualitative palette (Okabe and Ito), reordered so the
# highest-contrast hues are used first and pale yellow is used last.
fig_palette <- function() {
  c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#56B4E9", "#E69F00", "#000000", "#F0E442"
  )
}

fig_scale_colour <- function(n) {
  pal <- fig_palette()
  if (n <= length(pal)) {
    ggplot2::scale_colour_manual(values = pal[seq_len(n)])
  } else {
    ggplot2::scale_colour_viridis_d(end = 0.9)
  }
}

fig_scale_fill <- function(n) {
  pal <- fig_palette()
  if (n <= length(pal)) {
    ggplot2::scale_fill_manual(values = pal[seq_len(n)])
  } else {
    ggplot2::scale_fill_viridis_d(end = 0.9)
  }
}

#' Colourblind-safe discrete scales
#'
#' The palette this package's own figures use, exposed so that a figure built by
#' hand outside the package looks like the ones built inside it. It is the Okabe
#' and Ito qualitative palette, reordered so the highest-contrast hues come first
#' and pale yellow comes last. Beyond eight levels both scales fall back to
#' viridis, which stays discriminable at any number of levels.
#'
#' A `ggplot2` theme cannot set colour scales, so [gi_theme()] alone leaves a
#' plot on ggplot2's default hue palette, which is not colourblind-safe. Use
#' these alongside it.
#'
#' @param n Number of levels the scale has to cover.
#' @return A `ggplot2` scale, to be added to a plot with `+`.
#' @seealso [gi_theme()]
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   d <- data.frame(x = 1:6, y = c(1, 3, 2, 5, 4, 6), g = rep(c("a", "b"), each = 3))
#'   ggplot(d, aes(x, y, colour = g)) +
#'     geom_line() +
#'     scale_colour_gi(2) +
#'     gi_theme()
#' }
#' @name gi_scales
NULL

#' @rdname gi_scales
#' @export
scale_colour_gi <- function(n) {
  fig_require_ggplot2("scale_colour_gi()")
  fig_scale_colour(gi_check_scale_n(n))
}

#' @rdname gi_scales
#' @export
scale_fill_gi <- function(n) {
  fig_require_ggplot2("scale_fill_gi()")
  fig_scale_fill(gi_check_scale_n(n))
}

gi_check_scale_n <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("`n` must be a single whole number of levels, at least 1.", call. = FALSE)
  }
  as.integer(n)
}

fig_pick <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0L) NULL else hit[[1L]]
}

fig_column <- function(data, col, arg) {
  if (!is.character(col) || length(col) != 1L || is.na(col) || !nzchar(col)) {
    stop("'", arg, "' must be a single column name.", call. = FALSE)
  }
  if (!col %in% names(data)) {
    stop(
      "'", arg, "' names the column '", col, "', which is not in 'data'. ",
      "Columns available: ", paste(names(data), collapse = ", "), ".",
      call. = FALSE
    )
  }
  col
}

fig_positive_scalar <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
    stop("'", arg, "' must be a single positive number.", call. = FALSE)
  }
  as.numeric(value)
}

fig_first_number <- function(x, candidates) {
  for (nm in candidates) {
    v <- x[[nm]]
    if (is.numeric(v) && length(v) >= 1L && !all(is.na(v))) {
      v <- v[!is.na(v)]
      return(as.numeric(v[[1L]]))
    }
  }
  NA_real_
}

fig_as_data_frame <- function(x, candidates, arg) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x)) {
    nm <- fig_pick(x, candidates)
    if (!is.null(nm) && is.data.frame(x[[nm]])) {
      return(x[[nm]])
    }
    frames <- Filter(is.data.frame, x)
    if (length(frames)) {
      return(frames[[1L]])
    }
  }
  stop(
    "'", arg, "' must be a data.frame, or a list containing one (looked for ",
    paste(candidates, collapse = ", "), ").",
    call. = FALSE
  )
}

#' A restrained theme for published figures
#'
#' Plain institutional styling: white panel, hairline border, one faint grid,
#' no decorative colour, legend below the panel, type sized to stay legible when
#' a 7 by 5 inch figure is reduced to a single journal column.
#'
#' The theme carries no colour scale, because a ggplot2 theme cannot. The
#' plotting functions in this package apply a colourblind-safe discrete palette
#' (Okabe and Ito) themselves; add
#' `ggplot2::scale_colour_manual(values = gitrialsim:::fig_palette())` to hand
#' built plots that need the same colours.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family. The empty string keeps the device
#'   default, which is what most journals want.
#' @return A ggplot2 theme object, usable as `p + gi_theme()`.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   p <- ggplot2::ggplot(data.frame(n = 1:5, power = c(.2, .4, .6, .7, .8)))
#'   p + ggplot2::geom_point(ggplot2::aes(n, power)) + gi_theme()
#' }
#' @export
gi_theme <- function(base_size = 11, base_family = "") {
  fig_require_ggplot2("gi_theme()")
  base_size <- fig_positive_scalar(base_size, "base_size")
  if (!is.character(base_family) || length(base_family) != 1L) {
    stop("'base_family' must be a single character string.", call. = FALSE)
  }
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.border = ggplot2::element_rect(colour = "grey25", fill = NA, linewidth = 0.4),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.ticks = ggplot2::element_line(colour = "grey25", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "grey15", size = ggplot2::rel(0.9)),
      axis.title = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = "black", size = ggplot2::rel(0.95), hjust = 0),
      # Panels of a facetted figure carry their own axis labels, which run into
      # each other at the default spacing.
      panel.spacing = ggplot2::unit(12, "pt"),
      legend.position = "bottom",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      # The right margin has to clear half of the last axis label, which sits
      # outside the panel whenever a break lands on the panel edge.
      plot.margin = ggplot2::margin(6, 14, 6, 6)
    )
}

#' Plot a performance measure against a design input
#'
#' Draws power, or any other performance measure produced by `simulate_grid()`,
#' against one design input, optionally split by a second one. If `data` carries
#' a Monte Carlo standard error for `y` in a column named `<y>_mcse`,
#' `mcse_<y>` or `mcse`, a plus or minus 1.96 MCSE band is drawn behind the
#' line.
#'
#' Rows whose `x` or `y` value is missing or non-finite cannot be positioned, so
#' they are dropped before plotting, as in [plot_boundaries()]; if that empties
#' the data the function stops and names the column actually responsible. The
#' band is treated separately: a row whose Monte Carlo standard error is missing
#' or non-finite keeps its point and its line segment and simply carries no
#' band there, and a column of entirely missing standard errors gives the same
#' figure as no standard error column at all. A negative standard error is an
#' error rather than a band drawn upside down.
#'
#' @param data A tidy data frame with one row per evaluated design, such as the
#'   one [simulate_grid()] returns.
#' @param x Name of the numeric column holding the design input, for example
#'   `"n_per_arm_max"`.
#' @param y Name of the numeric column holding the performance measure. In the
#'   output of [simulate_grid()] power is reported as `"rejection_rate"`, with
#'   its Monte Carlo standard error in `"rejection_rate_mcse"`.
#' @param group Optional name of a column to split the curves by. Coerced to a
#'   factor and mapped to colour.
#' @return A ggplot object with one line and point series per group, a Monte
#'   Carlo band over the rows that report a usable standard error, the
#'   `gi_theme()` styling applied, axis labels taken from the column names, and
#'   no plot title.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   grid <- data.frame(
#'     n_per_arm_max = rep(c(500, 1000, 1500, 2000), 2),
#'     rejection_rate = c(0.31, 0.55, 0.71, 0.82, 0.24, 0.44, 0.60, 0.72),
#'     rejection_rate_mcse = rep(0.015, 8),
#'     design_type = rep(c("Fixed", "Group sequential"), each = 4)
#'   )
#'   plot_power_curve(
#'     grid,
#'     x = "n_per_arm_max", y = "rejection_rate", group = "design_type"
#'   )
#' }
#' @export
plot_power_curve <- function(data, x, y, group = NULL) {
  fig_require_ggplot2("plot_power_curve()")
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame; got ", class(data)[1L], ".", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'data' has no rows.", call. = FALSE)
  }
  x <- fig_column(data, x, "x")
  y <- fig_column(data, y, "y")
  if (!is.numeric(data[[x]])) {
    stop("'x' names the column '", x, "', which is not numeric.", call. = FALSE)
  }
  if (!is.numeric(data[[y]])) {
    stop("'y' names the column '", y, "', which is not numeric.", call. = FALSE)
  }

  df <- data.frame(gi_x = as.numeric(data[[x]]), gi_y = as.numeric(data[[y]]))
  grouped <- !is.null(group)
  if (grouped) {
    group <- fig_column(data, group, "group")
    df$gi_group <- factor(data[[group]])
  }

  mcse_col <- fig_pick(data, c(paste0(y, "_mcse"), paste0("mcse_", y), "mcse"))
  band <- !is.null(mcse_col) && is.numeric(data[[mcse_col]])
  if (band) {
    mcse <- as.numeric(data[[mcse_col]])
    negative <- which(is.finite(mcse) & mcse < 0)
    if (length(negative)) {
      stop(
        "'", mcse_col, "' holds a negative Monte Carlo standard error at row(s) ",
        paste(utils::head(negative, 5L), collapse = ", "),
        ". An uncertainty band cannot have a negative half width.",
        call. = FALSE
      )
    }
    half <- 1.96 * mcse
    df$gi_lower <- df$gi_y - half
    df$gi_upper <- df$gi_y + half
  }
  # A row is plottable whenever it has a position. A missing standard error
  # costs that row its band, never its point and its line segment.
  keep <- is.finite(df$gi_x) & is.finite(df$gi_y)
  if (!any(keep)) {
    culprits <- c(x, y)[c(!any(is.finite(df$gi_x)), !any(is.finite(df$gi_y)))]
    if (length(culprits) == 0L) culprits <- c(x, y)
    stop(
      "no finite values left to plot: every row has a missing or non-finite ",
      "value in ", paste0("'", culprits, "'", collapse = " or "), ".",
      call. = FALSE
    )
  }
  df <- df[keep, , drop = FALSE]
  df <- df[order(df$gi_x), , drop = FALSE]

  if (band) {
    band_df <- df[is.finite(df$gi_lower) & is.finite(df$gi_upper), , drop = FALSE]
    band <- nrow(band_df) > 0L
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = gi_x, y = gi_y))
  if (band) {
    p <- p + if (grouped) {
      ggplot2::geom_ribbon(
        data = band_df,
        ggplot2::aes(ymin = gi_lower, ymax = gi_upper, fill = gi_group),
        alpha = 0.15, colour = NA
      )
    } else {
      ggplot2::geom_ribbon(
        data = band_df,
        ggplot2::aes(ymin = gi_lower, ymax = gi_upper),
        fill = "grey60", alpha = 0.25
      )
    }
  }
  if (grouped) {
    n_levels <- nlevels(df$gi_group)
    p <- p +
      ggplot2::geom_line(ggplot2::aes(colour = gi_group), linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(colour = gi_group), size = 1.7) +
      fig_scale_colour(n_levels) +
      fig_scale_fill(n_levels) +
      ggplot2::labs(colour = group, fill = group)
  } else {
    p <- p +
      ggplot2::geom_line(linewidth = 0.6, colour = "#111111") +
      ggplot2::geom_point(size = 1.7, colour = "#111111")
  }
  p + ggplot2::labs(x = x, y = y) + gi_theme()
}

fig_design_labels <- function(designs) {
  labels <- names(designs)
  if (is.null(labels)) labels <- rep("", length(designs))
  labels[is.na(labels)] <- ""
  blank <- !nzchar(labels)
  if (any(blank)) {
    labels[blank] <- vapply(
      designs[blank],
      function(d) as.character(d$type %||% "design")[1L],
      character(1)
    )
  }
  make.unique(labels, sep = " ")
}

#' Compare the operating characteristics of several designs
#'
#' Places power, maximum total sample size and expected total sample size for a
#' set of designs on a common dot chart, one panel per measure so the different
#' units never share an axis.
#'
#' Power means two different things depending on how a design was evaluated, and
#' the two are never put on one axis. For a design solved in closed form the
#' `power` field is the target the sample size was solved for, a design input,
#' and it is drawn in a panel labelled `Power (analytic target)`. For a design
#' evaluated by simulation, which is what [design_bayesian()] returns and what
#' any design carrying a missing `power` is treated as, the number is the power
#' the simulation achieved, an output, and it is drawn in a separate panel
#' labelled `Power (simulated)`. Colour and point shape repeat the distinction
#' so it survives a greyscale reading of the panel strips.
#'
#' Power is read from the design's own `power` field whenever that is a single
#' non-missing number, which covers every design this package builds, including
#' the simulated power [design_bayesian()] reports. A design that leaves `power`
#' missing, as a hand-assembled or third-party design may, falls back to the
#' first number among `detail$simulated_power`, `detail$empirical_power`,
#' `detail$power` and `detail$prob_reject`, and is marked as simulated.
#'
#' Expected sample size is the value under the alternative, read from the first
#' of `detail$expected_n`, `detail$expected_sample_size`,
#' `detail$expected_n_total`, `detail$expected_n_alt`, `detail$expected_n_h1`,
#' `detail$en`, `detail$asn` or `detail$ess` that the design carries. A `fixed`
#' design reporting none is credited with its full sample size, since it has no
#' interim analysis at which to stop early; any other design reporting none
#' simply has no point in that panel, and the panel is dropped if no design
#' reports one.
#'
#' @param designs A list of `gi_design` objects, or a single `gi_design`. List
#'   names are used as the design labels; unnamed elements fall back to their
#'   `type`, deduplicated.
#' @return A ggplot object: designs on the vertical axis, measure value on the
#'   horizontal axis, panelled by measure, with the `gi_theme()` styling and no
#'   plot title. Analytic target power and simulated achieved power occupy
#'   separate panels and are also separated by colour and point shape.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   fixed <- structure(
#'     list(
#'       type = "fixed", alpha = 0.025, power = 0.9, n_total = 3400,
#'       n_per_arm = 1700, engine = "rpact::getSampleSizeRates", detail = list()
#'     ),
#'     class = c("gi_design", "list")
#'   )
#'   adaptive <- structure(
#'     list(
#'       type = "bayesian_adaptive", alpha = 0.025, power = NA_real_,
#'       n_total = 3400, n_per_arm = 1700, engine = "monte carlo",
#'       detail = list(simulated_power = 0.88, expected_n = 2450)
#'     ),
#'     class = c("gi_design", "list")
#'   )
#'   plot_operating_characteristics(list(Fixed = fixed, Adaptive = adaptive))
#' }
#' @export
plot_operating_characteristics <- function(designs) {
  fig_require_ggplot2("plot_operating_characteristics()")
  if (inherits(designs, "gi_design")) designs <- list(designs)
  if (!is.list(designs) || length(designs) == 0L) {
    stop("'designs' must be a non-empty list of gi_design objects.", call. = FALSE)
  }
  ok <- vapply(designs, inherits, logical(1), what = "gi_design")
  if (!all(ok)) {
    stop(
      "'designs' must contain only gi_design objects; element(s) ",
      paste(which(!ok), collapse = ", "), " are not.",
      call. = FALSE
    )
  }

  labels <- fig_design_labels(designs)
  power <- vapply(designs, function(d) {
    p <- d$power
    if (is.numeric(p) && length(p) == 1L && !is.na(p)) {
      return(as.numeric(p))
    }
    fig_first_number(
      d$detail %||% list(),
      c("simulated_power", "empirical_power", "power", "prob_reject")
    )
  }, numeric(1))
  n_max <- vapply(designs, function(d) {
    n <- d$n_total
    if (is.numeric(n) && length(n) == 1L && !is.na(n)) {
      return(as.numeric(n))
    }
    fig_first_number(d$detail %||% list(), c("n_total", "max_n", "n_max"))
  }, numeric(1))
  n_exp <- vapply(seq_along(designs), function(i) {
    d <- designs[[i]]
    value <- fig_first_number(
      d$detail %||% list(),
      c(
        "expected_n", "expected_sample_size", "expected_n_total",
        "expected_n_alt", "expected_n_h1", "en", "asn", "ess"
      )
    )
    # A fixed design has no interim analyses, so it always runs to its full
    # size; reporting that keeps the panel comparable across design types.
    if (is.na(value) && identical(d$type, "fixed")) value <- n_max[[i]]
    value
  }, numeric(1))
  kind <- vapply(designs, function(d) {
    if (identical(d$type, "bayesian_adaptive") ||
      (is.numeric(d$power) && length(d$power) == 1L && is.na(d$power))) {
      "Simulated"
    } else {
      "Analytic"
    }
  }, character(1))

  # Target power is a design input and simulated power is a design output. They
  # are not the same quantity, so they get their own panels rather than one
  # shared axis on which a reader would compare them.
  measures <- c(
    "Power (analytic target)", "Power (simulated)", "Maximum n", "Expected n"
  )
  power_measure <- ifelse(
    kind == "Simulated", "Power (simulated)", "Power (analytic target)"
  )
  rows <- data.frame(
    gi_label = factor(rep(labels, times = 3L), levels = rev(labels)),
    gi_measure = factor(
      c(power_measure, rep(c("Maximum n", "Expected n"), each = length(designs))),
      levels = measures
    ),
    gi_value = c(power, n_max, n_exp),
    gi_kind = factor(rep(kind, times = 3L), levels = c("Analytic", "Simulated")),
    stringsAsFactors = FALSE
  )
  rows <- rows[is.finite(rows$gi_value), , drop = FALSE]
  if (nrow(rows) == 0L) {
    stop(
      "none of the designs reported power, maximum n or expected n.",
      call. = FALSE
    )
  }
  rows$gi_measure <- droplevels(rows$gi_measure)
  rows$gi_kind <- droplevels(rows$gi_kind)

  p <- ggplot2::ggplot(rows, ggplot2::aes(x = gi_value, y = gi_label))
  if (nlevels(rows$gi_kind) > 1L) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = gi_kind, shape = gi_kind), size = 2.4) +
      fig_scale_colour(nlevels(rows$gi_kind)) +
      ggplot2::scale_shape_manual(values = c(Analytic = 16L, Simulated = 17L)) +
      ggplot2::labs(colour = "Value source", shape = "Value source")
  } else {
    p <- p + ggplot2::geom_point(size = 2.4, colour = "#111111")
  }
  p +
    ggplot2::facet_wrap(~gi_measure, scales = "free_x") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(x = NULL, y = NULL) +
    gi_theme()
}

#' Plot group-sequential stopping boundaries
#'
#' Efficacy and, where the design has them, futility boundaries on the z scale
#' against information fraction.
#'
#' @param design A `gi_design` of type `group_sequential`, or the data frame
#'   `gs_boundaries()` returns, which is accepted directly so figures can be
#'   redrawn from a stored table. Column names are resolved leniently: the
#'   information fraction may be `information_rate`, `information_fraction` or
#'   `timing`; the efficacy boundary `efficacy`, `efficacy_z` or `upper`; the
#'   futility boundary `futility`, `futility_z` or `lower`. A design with no
#'   futility column plots the efficacy boundary alone.
#' @return A ggplot object with one line and point series per boundary, the
#'   `gi_theme()` styling applied, and no plot title.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   bounds <- data.frame(
#'     analysis = 1:3,
#'     information_rate = c(1 / 3, 2 / 3, 1),
#'     efficacy = c(3.71, 2.51, 1.99),
#'     futility = c(-0.31, 0.98, 1.99)
#'   )
#'   plot_boundaries(bounds)
#' }
#' @export
plot_boundaries <- function(design) {
  fig_require_ggplot2("plot_boundaries()")
  if (is.data.frame(design)) {
    bnd <- design
  } else {
    if (!inherits(design, "gi_design")) {
      stop(
        "'design' must be a gi_design object, or the data.frame returned by ",
        "gs_boundaries(); got ", class(design)[1L], ".",
        call. = FALSE
      )
    }
    type <- as.character(design$type %||% NA_character_)[1L]
    if (!identical(type, "group_sequential")) {
      stop(
        "'design' must be a group_sequential design; this one has type '",
        type, "'.",
        call. = FALSE
      )
    }
    boundary_fun <- get0("gs_boundaries", mode = "function")
    if (is.null(boundary_fun)) {
      stop(
        "gs_boundaries() is not available, so boundaries cannot be extracted ",
        "from 'design'. Pass the boundary data.frame directly instead.",
        call. = FALSE
      )
    }
    bnd <- boundary_fun(design)
  }
  if (!is.data.frame(bnd) || nrow(bnd) == 0L) {
    stop("the boundary table is not a non-empty data.frame.", call. = FALSE)
  }

  info <- fig_pick(bnd, c(
    "information_rate", "information_fraction", "info_fraction",
    "information_frac", "information", "timing"
  ))
  if (is.null(info)) {
    stop(
      "the boundary table has no information fraction column. ",
      "Columns available: ", paste(names(bnd), collapse = ", "), ".",
      call. = FALSE
    )
  }
  eff <- fig_pick(bnd, c(
    "efficacy", "efficacy_z", "z_efficacy", "upper", "upper_bound", "critical_value"
  ))
  if (is.null(eff)) {
    stop(
      "the boundary table has no efficacy boundary column. ",
      "Columns available: ", paste(names(bnd), collapse = ", "), ".",
      call. = FALSE
    )
  }
  fut <- fig_pick(bnd, c("futility", "futility_z", "z_futility", "lower", "lower_bound"))

  long <- data.frame(
    gi_x = as.numeric(bnd[[info]]),
    gi_y = as.numeric(bnd[[eff]]),
    gi_bound = "Efficacy",
    stringsAsFactors = FALSE
  )
  if (!is.null(fut)) {
    long <- rbind(long, data.frame(
      gi_x = as.numeric(bnd[[info]]),
      gi_y = as.numeric(bnd[[fut]]),
      gi_bound = "Futility",
      stringsAsFactors = FALSE
    ))
  }
  # Non-binding futility bounds are reported as -Inf or NA at analyses where
  # stopping is not allowed; those rows are dropped rather than clipped.
  long <- long[is.finite(long$gi_x) & is.finite(long$gi_y), , drop = FALSE]
  if (nrow(long) == 0L) {
    stop("the boundary table has no finite boundary values to plot.", call. = FALSE)
  }
  long$gi_bound <- factor(long$gi_bound, levels = c("Efficacy", "Futility"))
  long$gi_bound <- droplevels(long$gi_bound)
  long <- long[order(long$gi_bound, long$gi_x), , drop = FALSE]

  ggplot2::ggplot(long, ggplot2::aes(x = gi_x, y = gi_y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_line(
      ggplot2::aes(colour = gi_bound, linetype = gi_bound),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = gi_bound), size = 1.9) +
    fig_scale_colour(nlevels(long$gi_bound)) +
    ggplot2::labs(
      x = "Information fraction",
      y = "Z statistic",
      colour = "Boundary",
      linetype = "Boundary"
    ) +
    gi_theme()
}

#' Plot an EVSI curve
#'
#' Expected value of sample information against the size of the trial being
#' costed, with a Monte Carlo uncertainty ribbon where the curve carries one.
#'
#' Rows whose sample size or EVSI value is missing or non-finite cannot be
#' positioned, so they are dropped before plotting, as in [plot_boundaries()];
#' if that empties the curve the function stops. The ribbon is treated
#' separately: a row whose interval bounds are missing or non-finite keeps its
#' point on the curve and simply carries no ribbon there, and a curve whose
#' bounds are all unusable is drawn as a bare line. A negative standard error is
#' an error rather than a ribbon drawn upside down.
#'
#' @param voi_curve_result The result of the package's EVSI curve routine: a
#'   data frame, or a list containing one under `curve`, `evsi_curve`, `evsi`,
#'   `data` or `results`. The sample size column may be `n`, `n_total`,
#'   `n_per_arm` or `sample_size`; the value column `evsi`, `EVSI` or `value`.
#'   The ribbon uses explicit `lower`/`upper` (or `ci_lower`/`ci_upper`) columns
#'   when present, otherwise plus or minus 1.96 times a `se`, `mcse` or
#'   `evsi_se` column, and is omitted when the curve reports neither.
#' @return A ggplot object with the EVSI curve, its Monte Carlo ribbon where
#'   available, the `gi_theme()` styling, and no plot title.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   voi <- list(curve = data.frame(
#'     n_per_arm = seq(200, 2000, by = 200),
#'     evsi = c(120, 210, 275, 320, 350, 368, 378, 383, 385, 386) * 1000,
#'     se = rep(9000, 10)
#'   ))
#'   plot_evsi(voi)
#' }
#' @export
plot_evsi <- function(voi_curve_result) {
  fig_require_ggplot2("plot_evsi()")
  df <- fig_as_data_frame(
    voi_curve_result,
    c("curve", "evsi_curve", "evsi", "data", "results"),
    "voi_curve_result"
  )
  if (nrow(df) == 0L) {
    stop("'voi_curve_result' holds an empty curve.", call. = FALSE)
  }

  n_col <- fig_pick(df, c(
    "n", "n_total", "n_per_arm", "sample_size", "trial_n", "n_future", "n_new"
  ))
  if (is.null(n_col)) {
    stop(
      "the EVSI curve has no sample size column. Columns available: ",
      paste(names(df), collapse = ", "), ".",
      call. = FALSE
    )
  }
  v_col <- fig_pick(df, c("evsi", "EVSI", "value", "mean", "estimate"))
  if (is.null(v_col)) {
    stop(
      "the EVSI curve has no EVSI column. Columns available: ",
      paste(names(df), collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- data.frame(
    gi_x = as.numeric(df[[n_col]]),
    gi_y = as.numeric(df[[v_col]])
  )
  lo_col <- fig_pick(df, c("lower", "ci_lower", "evsi_lower", "lwr", "conf_low"))
  hi_col <- fig_pick(df, c("upper", "ci_upper", "evsi_upper", "upr", "conf_high"))
  se_col <- fig_pick(df, c("se", "mcse", "evsi_se", "evsi_mcse", "se_evsi", "mc_se"))
  ribbon <- FALSE
  if (!is.null(lo_col) && !is.null(hi_col)) {
    out$gi_lower <- as.numeric(df[[lo_col]])
    out$gi_upper <- as.numeric(df[[hi_col]])
    inverted <- which(is.finite(out$gi_lower) & is.finite(out$gi_upper) &
      out$gi_lower > out$gi_upper)
    if (length(inverted)) {
      stop(
        "'", lo_col, "' exceeds '", hi_col, "' at row(s) ",
        paste(utils::head(inverted, 5L), collapse = ", "),
        ". An uncertainty interval cannot run downwards.",
        call. = FALSE
      )
    }
    ribbon <- TRUE
  } else if (!is.null(se_col) && is.numeric(df[[se_col]])) {
    se <- as.numeric(df[[se_col]])
    negative <- which(is.finite(se) & se < 0)
    if (length(negative)) {
      stop(
        "'", se_col, "' holds a negative standard error at row(s) ",
        paste(utils::head(negative, 5L), collapse = ", "),
        ". An uncertainty ribbon cannot have a negative half width.",
        call. = FALSE
      )
    }
    half <- 1.96 * se
    out$gi_lower <- out$gi_y - half
    out$gi_upper <- out$gi_y + half
    ribbon <- TRUE
  }
  out <- out[is.finite(out$gi_x) & is.finite(out$gi_y), , drop = FALSE]
  if (nrow(out) == 0L) {
    stop("the EVSI curve has no finite points to plot.", call. = FALSE)
  }
  out <- out[order(out$gi_x), , drop = FALSE]
  # A missing bound costs that row its ribbon, never its place on the curve.
  if (ribbon) {
    ribbon_df <- out[is.finite(out$gi_lower) & is.finite(out$gi_upper), , drop = FALSE]
    ribbon <- nrow(ribbon_df) > 0L
  }

  x_label <- switch(n_col,
    n_per_arm = "Sample size per arm",
    n_total = "Total sample size",
    sample_size = "Total sample size",
    n = "Total sample size",
    n_col
  )

  p <- ggplot2::ggplot(out, ggplot2::aes(x = gi_x, y = gi_y))
  if (ribbon) {
    p <- p + ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(ymin = gi_lower, ymax = gi_upper),
      fill = "grey60", alpha = 0.25
    )
  }
  p +
    ggplot2::geom_line(linewidth = 0.6, colour = "#111111") +
    ggplot2::labs(x = x_label, y = "Expected value of sample information") +
    gi_theme()
}

fig_emulator_predict <- function(fit, newdata) {
  out <- if (is.function(fit)) {
    fit(newdata)
  } else if (is.list(fit) && is.function(fit[["predict"]])) {
    fit[["predict"]](newdata)
  } else {
    res <- try(stats::predict(fit, newdata = newdata), silent = TRUE)
    if (inherits(res, "try-error")) {
      res <- try(stats::predict(fit, newdata), silent = TRUE)
    }
    if (inherits(res, "try-error")) {
      stop(
        "could not predict from 'fit': it is not a function, has no predict ",
        "element, and no predict() method accepted a data.frame of new inputs.",
        call. = FALSE
      )
    }
    res
  }
  if (is.numeric(out) && is.null(dim(out))) {
    return(as.numeric(out))
  }
  if (is.matrix(out)) {
    return(as.numeric(out[, 1L]))
  }
  if (is.list(out)) {
    nm <- fig_pick(out, c(
      "mean", "fit", "prediction", "pred", "posterior_mean", "estimate", "y", "value"
    ))
    if (!is.null(nm) && is.numeric(out[[nm]])) {
      return(as.numeric(out[[nm]]))
    }
  }
  stop(
    "could not read a posterior mean out of the emulator prediction; expected ",
    "a numeric vector or a list holding one under mean, fit or prediction.",
    call. = FALSE
  )
}

fig_training_points <- function(fit, x_name, y_name) {
  if (is.function(fit) || !is.list(fit)) {
    return(NULL)
  }
  for (nm in c("design", "X", "x", "train", "training", "training_data", "inputs", "data", "grid", "points")) {
    cand <- fit[[nm]]
    if (is.matrix(cand)) cand <- as.data.frame(cand)
    if (is.data.frame(cand) && all(c(x_name, y_name) %in% names(cand))) {
      out <- data.frame(
        gi_x = as.numeric(cand[[x_name]]),
        gi_y = as.numeric(cand[[y_name]])
      )
      out <- out[is.finite(out$gi_x) & is.finite(out$gi_y), , drop = FALSE]
      if (nrow(out)) {
        return(out)
      }
    }
  }
  NULL
}

#' Plot an emulator posterior mean surface
#'
#' Evaluates a fitted emulator on a regular grid over two of its inputs and
#' draws the posterior mean as a filled surface with contour lines. Where the
#' fit carries its training design, those points are overlaid so a reader can
#' see which parts of the surface are supported by actual simulations and which
#' are extrapolation.
#'
#' @param fit A fitted emulator. Any of these work: a function taking a data
#'   frame of inputs and returning predictions; a list with a `predict` element
#'   that is such a function; or an object with a `predict()` method. A
#'   prediction may be a numeric vector, a matrix whose first column is the
#'   mean, or a list holding the mean under `mean`, `fit` or `prediction`.
#'   Training points are looked for in `design`, `X`, `train`, `inputs`, `data`
#'   or `grid`, and are simply omitted if none of those hold a data frame with
#'   both plotted inputs.
#' @param bounds A named list of exactly two ranges, each a numeric vector of
#'   length two giving the lower and upper plotting limit. Names must match the
#'   input names the emulator expects. The first is the horizontal axis.
#' @param n_grid Number of grid points per axis, so the surface is evaluated at
#'   `n_grid^2` points. Must be a whole number of at least 2; a fractional value
#'   is rejected rather than quietly truncated. Defaults to 60.
#' @return A ggplot object: a raster of the posterior mean with contour lines
#'   over it, training points overlaid where available, the `gi_theme()` styling
#'   and a continuous viridis fill scale, and no plot title.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   fit <- list(
#'     predict = function(newdata) {
#'       0.9 * (1 - exp(-newdata$n_per_arm / 900)) - 4 * (newdata$control_rate - 0.07)^2
#'     },
#'     design = data.frame(
#'       n_per_arm = c(300, 800, 1200, 1800),
#'       control_rate = c(0.04, 0.09, 0.055, 0.08)
#'     )
#'   )
#'   plot_emulator_surface(
#'     fit,
#'     bounds = list(n_per_arm = c(200, 2000), control_rate = c(0.03, 0.10)),
#'     n_grid = 30
#'   )
#' }
#' @export
plot_emulator_surface <- function(fit, bounds, n_grid = 60) {
  fig_require_ggplot2("plot_emulator_surface()")
  if (is.null(fit)) {
    stop("'fit' must be a fitted emulator, not NULL.", call. = FALSE)
  }
  if (!is.list(bounds) || length(bounds) != 2L) {
    stop(
      "'bounds' must be a list of exactly two ranges, one per plotted input.",
      call. = FALSE
    )
  }
  nms <- names(bounds)
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    stop(
      "'bounds' must have two distinct names matching the emulator inputs to plot.",
      call. = FALSE
    )
  }
  for (i in seq_along(bounds)) {
    b <- bounds[[i]]
    if (!is.numeric(b) || length(b) != 2L || any(!is.finite(b)) || b[1L] >= b[2L]) {
      stop(
        "'bounds$", nms[i], "' must be two finite numbers, lower then upper.",
        call. = FALSE
      )
    }
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) ||
    n_grid < 2 || n_grid > .Machine$integer.max || n_grid != trunc(n_grid)) {
    stop(
      "'n_grid' must be a single whole number of at least 2. A fractional ",
      "number of grid points is rejected rather than truncated.",
      call. = FALSE
    )
  }
  n_grid <- as.integer(n_grid)

  x_name <- nms[1L]
  y_name <- nms[2L]
  grid <- expand.grid(
    stats::setNames(
      list(
        seq(bounds[[1L]][1L], bounds[[1L]][2L], length.out = n_grid),
        seq(bounds[[2L]][1L], bounds[[2L]][2L], length.out = n_grid)
      ),
      c(x_name, y_name)
    ),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  z <- fig_emulator_predict(fit, grid)
  if (length(z) != nrow(grid)) {
    stop(
      "the emulator returned ", length(z), " predictions for ", nrow(grid),
      " grid points.",
      call. = FALSE
    )
  }
  surface <- data.frame(gi_x = grid[[x_name]], gi_y = grid[[y_name]], gi_z = as.numeric(z))
  # An infinite prediction is not missing, and it colours nothing: a surface
  # with no finite value at all has no fill scale to build.
  if (!any(is.finite(surface$gi_z))) {
    stop("the emulator returned no finite predictions on this grid.", call. = FALSE)
  }

  p <- ggplot2::ggplot(surface, ggplot2::aes(x = gi_x, y = gi_y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = gi_z), interpolate = TRUE) +
    ggplot2::geom_contour(
      ggplot2::aes(z = gi_z),
      colour = "white", linewidth = 0.25, alpha = 0.6
    ) +
    ggplot2::scale_fill_viridis_c()

  train <- fig_training_points(fit, x_name, y_name)
  if (!is.null(train)) {
    p <- p + ggplot2::geom_point(
      data = train, shape = 21, size = 1.8, stroke = 0.4,
      colour = "white", fill = "black"
    )
  }
  p +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      barwidth = ggplot2::unit(28, "mm"),
      barheight = ggplot2::unit(3, "mm"),
      ticks.colour = "white",
      frame.colour = "grey25"
    )) +
    ggplot2::labs(x = x_name, y = y_name, fill = "Posterior mean") +
    gi_theme()
}

#' Save a figure at print quality
#'
#' Thin wrapper over [ggplot2::ggsave()] with defaults suited to a
#' single-column journal figure: 7 by 5 inches at 300 dpi on a white
#' background. The output format follows the file extension. Missing parent
#' directories are created.
#'
#' @param plot A ggplot object, such as any of the figures in this package.
#' @param path File to write. The extension chooses the device, so use `.pdf`
#'   or `.eps` for vector output and `.png` or `.tiff` for raster.
#' @param width,height Figure size in inches.
#' @param dpi Resolution for raster devices.
#' @return The normalised path, invisibly.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   grid <- data.frame(n_per_arm = c(500, 1000, 1500), power = c(0.3, 0.6, 0.8))
#'   p <- plot_power_curve(grid, x = "n_per_arm", y = "power")
#'   out <- save_figure(p, file.path(tempdir(), "power.png"))
#'   file.exists(out)
#' }
#' @export
save_figure <- function(plot, path, width = 7, height = 5, dpi = 300) {
  fig_require_ggplot2("save_figure()")
  if (!inherits(plot, "ggplot")) {
    stop("'plot' must be a ggplot object; got ", class(plot)[1L], ".", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("'path' must be a single non-empty file path.", call. = FALSE)
  }
  width <- fig_positive_scalar(width, "width")
  height <- fig_positive_scalar(height, "height")
  dpi <- fig_positive_scalar(dpi, "dpi")

  target_dir <- dirname(path)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ggplot2::ggsave(
    filename = path, plot = plot,
    width = width, height = height, units = "in", dpi = dpi, bg = "white"
  )
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

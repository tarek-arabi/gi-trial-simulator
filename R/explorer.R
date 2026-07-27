#' Launch the interactive trial design explorer
#'
#' Opens a Shiny application for exploring how required sample size, stopping
#' boundaries and operating characteristics respond to the assumed event rates
#' and the choice of design. Sample sizes and boundaries come from rpact;
#' Bayesian operating characteristics are simulated, so that panel is slower.
#'
#' The application is a demonstration of the package's public interface. It
#' calls only exported functions, so nothing it does is unavailable from the
#' console.
#'
#' @param ... Passed to [shiny::runApp()], for example `port` or `launch.browser`.
#' @return Called for its side effect. Does not return until the app is closed.
#' @examples
#' if (interactive()) {
#'   run_explorer()
#' }
#' @export
run_explorer <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The explorer needs the 'shiny' package. Install it with install.packages('shiny').",
      call. = FALSE
    )
  }
  app_dir <- system.file("shiny", package = "gitrialsim")
  if (!nzchar(app_dir)) {
    stop("Could not locate the explorer app inside the installed package.", call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}

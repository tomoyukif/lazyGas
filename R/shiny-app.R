#' Launch the lazyGas Shiny explorer
#'
#' @param ... Arguments passed to \code{shiny::runApp()}.
#' @export
runLazyGasShiny <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required.", call. = FALSE)
  }
  app_dir <- system.file("shiny", package = "lazyGas")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    dev_dir <- file.path(getwd(), "inst", "shiny")
    if (dir.exists(dev_dir)) {
      app_dir <- dev_dir
    } else {
      stop("Shiny app directory not found.", call. = FALSE)
    }
  }
  shiny::runApp(app_dir, ...)
}

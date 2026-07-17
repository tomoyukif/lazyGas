#' Launch the lazyGas pipeline runner (Shiny)
#'
#' Opens a lightweight Shiny application to run the core GWAS pipeline
#' (\code{buildLazyGas()} → \code{assignPheno()} → \code{runLazyGas()}) from
#' filesystem paths. Use [runLazyGasExplorer()] afterward to explore and rank
#' candidates.
#'
#' @param gds_path Optional default GDS path pre-filled in the Shiny UI.
#' @param ... Arguments passed to \code{shiny::runApp()}.
#'
#' @details
#' Requires the Suggested package \code{shiny}. Phenotype CSV must contain a
#' sample id column named \code{id} or \code{ID}. The companion store
#' (\code{{stem}.lazygas/}) is created next to the GDS file.
#'
#' @seealso [runLazyGas()], [runLazyGasExplorer()], [runExplorerDemo()]
#'
#' @export
runLazyGasRunner <- function(gds_path = NULL, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required for GFF import.", call. = FALSE)
  }
  if (is.null(gds_path)) {
    demo <- file.path(getwd(), "demo_output", "explorer", "sample.gds")
    if (file.exists(demo)) {
      gds_path <- normalizePath(demo, winslash = "/", mustWork = FALSE)
    }
  }
  if (!is.null(gds_path) && nzchar(gds_path)) {
    options(lazygas.runner.gds_path = gds_path)
  }
  app_dir <- system.file("shiny", "runner", package = "lazyGas")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    dev_dir <- file.path(getwd(), "inst", "shiny", "runner")
    if (dir.exists(dev_dir)) {
      app_dir <- dev_dir
    } else {
      stop("Shiny runner directory not found.", call. = FALSE)
    }
  }
  shiny::runApp(app_dir, ...)
}

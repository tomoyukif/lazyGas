#' Launch the lazyGas phenotype explorer (Shiny)
#'
#' Opens an interactive Shiny application for phenotype-guided candidate-gene
#' exploration with optional local LLM explanations. The UI also embeds the
#' same GWAS overview and per-gene locus views used by
#' [makeInteractiveDashboard()] (phenotype distribution, Manhattan, peaks,
#' haplotype plots, and on-demand variant viewer).
#'
#' @param gds_path Optional default GDS path pre-filled in the Shiny UI.
#' @param ... Arguments passed to \code{shiny::runApp()}.
#'
#' @details
#' Requires the Suggested package \code{shiny}. Interactive plots use Imports
#' \code{plotly} / \code{reactable} / \code{htmltools}. For the variant viewer
#' tab, provide a GFF path in the sidebar (demo default:
#' \code{inst/extdata/demo_annotation.gff}). Stored SnpEff tables in the
#' companion store are preferred; an optional SnpEff GDS path is used when
#' stored annotations are missing.
#'
#' @seealso [runLazyGasRunner()], [makeInteractiveDashboard()], [runExplorerDemo()],
#'   [rankPhenotypeCandidates()]
#'
#' @export
runLazyGasExplorer <- function(gds_path = NULL, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required.", call. = FALSE)
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "Package 'plotly' is required for explorer GWAS / locus plots. ",
      "Install it with install.packages(\"plotly\").",
      call. = FALSE
    )
  }
  if (is.null(gds_path)) {
    demo <- file.path(getwd(), "demo_output", "explorer", "sample.gds")
    if (file.exists(demo)) {
      gds_path <- normalizePath(demo, winslash = "/", mustWork = FALSE)
    }
  }
  if (!is.null(gds_path) && nzchar(gds_path)) {
    options(lazygas.explorer.gds_path = gds_path)
  }
  app_dir <- system.file("shiny", "explorer", package = "lazyGas")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    dev_dir <- file.path(getwd(), "inst", "shiny", "explorer")
    if (dir.exists(dev_dir)) {
      app_dir <- dev_dir
    } else {
      stop("Shiny explorer directory not found.", call. = FALSE)
    }
  }
  shiny::runApp(app_dir, ...)
}

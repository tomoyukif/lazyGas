#' Resolve GDS and companion-store paths for a lazyGas project
#'
#' Given a GDS file path, returns the GDS path and the companion Parquet/SQLite
#' store path (\code{{stem}.lazygas/} by default). Used by Shiny apps and
#' scripts when the GDS file and companion store live in the same folder.
#'
#' @param gds_fn Path to a \code{.gds} file.
#' @param companion_path Optional explicit companion-store path (folder ending
#'   in \code{.lazygas} or a \code{.lazygas.sqlite} file).
#'
#' @return A list with \code{gds} and \code{companion}.
#' @export
resolveLazyGasPaths <- function(gds_fn, companion_path = NULL) {
  if (!is.character(gds_fn) || length(gds_fn) != 1L || !nzchar(trimws(gds_fn))) {
    stop("'gds_fn' must be a non-empty character string.", call. = FALSE)
  }
  gds_fn <- normalizePath(trimws(gds_fn), winslash = "/", mustWork = FALSE)
  if (!file.exists(gds_fn)) {
    stop("GDS file not found: ", gds_fn, call. = FALSE)
  }

  companion <- NULL
  if (!is.null(companion_path) && length(companion_path) == 1L &&
      nzchar(trimws(companion_path))) {
    companion <- normalizePath(trimws(companion_path), winslash = "/", mustWork = FALSE)
    if (!file.exists(companion)) {
      stop("Companion store not found: ", companion, call. = FALSE)
    }
  } else {
    stem <- sub("\\.gds$", "", gds_fn, ignore.case = TRUE)
    candidates <- unique(c(
      paste0(stem, ".lazygas"),
      file.path(dirname(gds_fn), paste0(basename(stem), ".lazygas")),
      paste0(stem, ".lazygas.sqlite"),
      file.path(dirname(gds_fn), paste0(basename(stem), ".lazygas.sqlite"))
    ))
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) {
      companion <- hit[1L]
    }
  }

  list(gds = gds_fn, companion = companion)
}

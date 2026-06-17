################################################################################
#' Run the lazyGas association pipeline
#'
#' Orchestrates scan, peak calling, recalculation, candidate listing, and optional
#' HTML reporting in one call.
#'
#' @param object A \code{LazyGas} object with phenotype data assigned.
#' @param steps Character vector of steps to run. One or more of \code{"scan"},
#'   \code{"peakcall"}, \code{"recalc"}, \code{"candidate"}, \code{"dashboard"},
#'   \code{"summary"}.
#' @param resume If \code{TRUE}, skip steps whose outputs already exist in the
#'   companion store.
#' @param gff,gff_fn \code{GRanges} or path to GFF for \code{listCandidate()} and
#'   dashboard steps.
#' @param snpeff Optional \code{snpeff_gds} object.
#' @param ann Optional annotation \code{data.frame}.
#' @param out_fn Output HTML path for dashboard/summary steps.
#' @param recalc Use recalculated peaks for candidate listing and reports.
#' @param seed Optional random seed recorded in pipeline metadata.
#' @param ... Arguments passed to \code{scanAssoc()}, \code{callPeakBlock()},
#'   \code{recalcAssoc()}, \code{listCandidate()}, or reporting functions.
#'
#' @return The \code{LazyGas} object (invisibly).
#' @export
runLazyGas <- function(object,
                       steps = c("scan", "peakcall", "recalc", "candidate"),
                       resume = TRUE,
                       gff = NULL,
                       gff_fn = NULL,
                       snpeff = NULL,
                       ann = NULL,
                       out_fn = NULL,
                       recalc = TRUE,
                       seed = NULL,
                       ...) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  steps <- match.arg(
    steps,
    c("scan", "peakcall", "recalc", "candidate", "dashboard", "summary"),
    several.ok = TRUE
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(gff) && !is.null(gff_fn)) {
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop("Package 'rtracklayer' is required to read gff_fn.", call. = FALSE)
    }
    gff <- rtracklayer::import.gff(gff_fn)
  }

  dots <- list(...)
  scan_args <- dots[names(dots) %in% c(
    "formula", "null_formula", "conv_fun", "fixed_effect", "geno_format",
    "kruskal", "method"
  )]
  peak_args <- dots[names(dots) %in% c("signif", "threshold", "limit_peakcall", "n_threads")]
  recalc_args <- dots[names(dots) %in% c("n_threads", "refine_position", "grouping_threshold")]
  cand_args <- dots[names(dots) %in% c("recalc")]
  report_args <- dots[names(dots) %in% c("what", "genes", "peak_id")]

  ran <- character()
  skipped <- character()

  if ("scan" %in% steps) {
    if (resume && .store_section_exists(object, "scan")) {
      skipped <- c(skipped, "scan")
    } else {
      message("runLazyGas: scanAssoc()")
      do.call(scanAssoc, c(list(object = object), scan_args))
      ran <- c(ran, "scan")
    }
  }

  if ("peakcall" %in% steps) {
    if (resume && .store_section_exists(object, "peakcall")) {
      skipped <- c(skipped, "peakcall")
    } else {
      if (!.store_section_exists(object, "scan")) {
        stop("No scan data. Run step 'scan' first.", call. = FALSE)
      }
      message("runLazyGas: callPeakBlock()")
      do.call(callPeakBlock, c(list(object = object), peak_args))
      ran <- c(ran, "peakcall")
    }
  }

  if ("recalc" %in% steps) {
    if (resume && .store_section_exists(object, "recalc")) {
      skipped <- c(skipped, "recalc")
    } else {
      if (!.store_section_exists(object, "peakcall")) {
        stop("No peakcall data. Run step 'peakcall' first.", call. = FALSE)
      }
      message("runLazyGas: recalcAssoc()")
      do.call(recalcAssoc, c(list(object = object), recalc_args))
      ran <- c(ran, "recalc")
    }
  }

  if ("candidate" %in% steps) {
    if (is.null(gff)) {
      stop("gff or gff_fn is required for step 'candidate'.", call. = FALSE)
    }
    pheno_names <- getPheno(object)$pheno_names
    all_exist <- all(vapply(
      pheno_names,
      function(pn) .store_dataset_exists(object, "candidate", pn),
      logical(1L)
    ))
    if (resume && all_exist) {
      skipped <- c(skipped, "candidate")
    } else {
      if (recalc && !.store_section_exists(object, "recalc")) {
        stop("No recalc data. Run step 'recalc' first or set recalc = FALSE.",
             call. = FALSE)
      }
      if (!recalc && !.store_section_exists(object, "peakcall")) {
        stop("No peakcall data. Run step 'peakcall' first.", call. = FALSE)
      }
      message("runLazyGas: listCandidate()")
      do.call(
        listCandidate,
        c(
          list(object = object, gff = gff, snpeff = snpeff, ann = ann, recalc = recalc),
          cand_args
        )
      )
      ran <- c(ran, "candidate")
    }
  }

  if ("summary" %in% steps || "dashboard" %in% steps) {
    if (is.null(out_fn)) {
      out_fn <- file.path(
        dirname(.store_gds_fn(object)),
        "lazygas_report.html"
      )
    }
    pheno <- getPheno(object)$pheno_names[1]
    what <- if (!is.null(report_args$what)) {
      report_args$what
    } else if ("dashboard" %in% steps) {
      c("scan_png", "peakcall", "recalc", "candidate")
    } else {
      c("scan_png", "peakcall", "recalc", "candidate")
    }
    if ("dashboard" %in% steps) {
      if (is.null(gff)) {
        stop("gff or gff_fn is required for step 'dashboard'.", call. = FALSE)
      }
      message("runLazyGas: makeInteractiveDashboard() -> ", out_fn)
      do.call(
        makeInteractiveDashboard,
        c(
          list(
            object = object,
            pheno = pheno,
            gff = gff,
            out_fn = out_fn,
            snpeff = snpeff,
            ann = ann,
            recalc = recalc,
            what = what
          ),
          report_args[setdiff(names(report_args), "what")]
        )
      )
    } else {
      message("runLazyGas: makeInteractiveSummary() -> ", out_fn)
      do.call(
        makeInteractiveSummary,
        c(
          list(object = object, pheno = pheno, out_fn = out_fn, what = what),
          report_args[setdiff(names(report_args), "what")]
        )
      )
    }
    ran <- c(ran, intersect(steps, c("dashboard", "summary")))
  }

  hist <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    steps_requested = paste(steps, collapse = ","),
    steps_ran = paste(unique(ran), collapse = ","),
    steps_skipped = paste(skipped, collapse = ","),
    seed = if (is.null(seed)) NA else seed,
    stringsAsFactors = FALSE
  )
  .store_append_pipeline_history(object, hist)

  invisible(object)
}

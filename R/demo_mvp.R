################################################################################
#' Resolve lazyGas demo extdata directory
#'
#' @noRd
.demo_extdata_dir <- function() {
  extdata <- system.file("extdata", package = "lazyGas")
  if (!nzchar(extdata) || !file.exists(file.path(extdata, "demo_annotation.gff"))) {
    dev_inst <- file.path(getwd(), "inst", "extdata")
    if (file.exists(file.path(dev_inst, "demo_annotation.gff"))) {
      extdata <- dev_inst
    }
  }
  if (!nzchar(extdata)) {
    stop(
      "Could not find lazyGas extdata. ",
      "Install the package or run devtools::load_all() from the package root.",
      call. = FALSE
    )
  }
  extdata
}

#' Run MVP feature demo (pipeline, QC, multi-trait, fine-mapping, orthologs)
#'
#' Uses bundled \code{sample.gds}, multi-trait phenotypes
#' (\code{demo_pheno_multitrait.csv}), annotation, SnpEff, and ortholog tables.
#' Writes an HTML report and prints summaries for the new APIs introduced in
#' lazyGas v0.6+.
#'
#' @param out_dir Output directory (created if missing). Default:
#'   \code{demo_output/mvp} under the working directory.
#' @param overwrite If \code{TRUE}, rebuild companion-store results.
#' @param seed Random seed recorded in pipeline metadata.
#'
#' @return Invisibly, a list with \code{lg}, \code{gff}, \code{snpeff},
#'   \code{ann}, \code{ortholog}, \code{report_html}, and \code{out_dir}.
#' @export
#'
#' @seealso [runDashboardDemo()]
#'
runMvpDemo <- function(out_dir = file.path(getwd(), "demo_output", "mvp"),
                       overwrite = TRUE,
                       seed = 42L) {
  if (!requireNamespace("GBScleanR", quietly = TRUE)) {
    stop("Package 'GBScleanR' is required.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required.", call. = FALSE)
  }

  extdata <- .demo_extdata_dir()
  paths <- list(
    gds = file.path(extdata, "sample.gds"),
    pheno = file.path(extdata, "demo_pheno_multitrait.csv"),
    gff = file.path(extdata, "demo_annotation.gff"),
    vcf = file.path(extdata, "demo_snpeff.vcf"),
    ann = file.path(extdata, "demo_ann.csv"),
    ortholog = file.path(extdata, "demo_orthologs.csv")
  )
  for (fn in unlist(paths)) {
    if (!file.exists(fn)) {
      stop("Demo file not found: ", fn, call. = FALSE)
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gds_copy <- file.path(out_dir, "sample.gds")
  file.copy(paths$gds, gds_copy, overwrite = TRUE)

  message("Loading demo annotation ...")
  gff <- rtracklayer::import.gff(paths$gff)
  ann <- read.csv(paths$ann, stringsAsFactors = FALSE)
  ortholog <- read.csv(paths$ortholog, stringsAsFactors = FALSE)

  snpeff_gds_fn <- file.path(out_dir, "demo_snpeff.gds")
  if (!file.exists(snpeff_gds_fn)) {
    message("Converting demo SnpEff VCF to GDS ...")
    snpeff2gds(vcf_fn = paths$vcf, out_fn = snpeff_gds_fn, verbose = FALSE)
  }
  snpeff <- open_snpeff(snpeff_gds_fn)

  message("Building LazyGas object ...")
  lg <- buildLazyGas(
    gds_fn = gds_copy,
    load_filter = TRUE,
    overwrite = overwrite
  )
  on.exit(try(closeGDS(lg, verbose = FALSE), silent = TRUE), add = TRUE)

  pheno <- read.csv(paths$pheno)
  trait_names <- c("fruit_weight", "fruit_length", "flowering_time")
  lg <- assignPheno(object = lg, pheno = pheno, rename = trait_names)

  conv_fun <- makeConvFun(geno_format = "dosage", n_levels = 3L)
  report_html <- file.path(out_dir, "lazygas_mvp_report.html")

  message("Running runLazyGas() pipeline ...")
  lg <- runLazyGas(
    object = lg,
    steps = c("scan", "peakcall", "recalc", "candidate", "summary"),
    resume = !overwrite,
    gff = gff,
    snpeff = snpeff,
    ann = ann,
    out_fn = report_html,
    recalc = TRUE,
    seed = seed,
    formula = "add + dom",
    conv_fun = conv_fun,
    geno_format = "dosage",
    limit_peakcall = 3L,
    n_threads = 1L,
    what = c("scan_png", "qq", "qc", "peakcall", "recalc", "cross_trait", "candidate")
  )

  message("\n--- GWAS QC (per trait) ---")
  qc <- summarizeGWASQC(object = lg, store = TRUE)
  print(qc)

  message("\n--- Cross-trait peak clustering ---")
  mt <- plotMultiTraitOverview(lg, recalc = TRUE)
  print(summarizeCrossTraitPeaks(mt$clusters))
  if (nrow(mt$shared_genes) > 0L) {
    message("Shared candidate genes across traits:")
    print(mt$shared_genes)
  }

  message("\n--- Fine-mapping (first trait, peak 1) ---")
  peak_id <- 1L
  cond <- conditionalAssoc(
    object = lg,
    pheno = trait_names[1],
    peak_id = peak_id,
    recalc = TRUE,
    k_step = 2L
  )
  print(head(cond[order(cond$P_conditional), ], 5))

  cred <- calcCredibleSet(
    object = lg,
    pheno = trait_names[1],
    peak_id = peak_id,
    use_conditional = TRUE,
    recalc = TRUE
  )
  message("Credible set summary:")
  print(attr(cred, "summary"))

  message("\n--- Ortholog annotation (", trait_names[1], ") ---")
  cand <- lazyData(lg, dataset = "candidate", pheno = trait_names[1])
  cand_ortho <- annotateOrthologs(candidate = cand, ortholog = ortholog)
  ortho_show <- merge(
    cand_ortho[, c("Gene_ID", "ortholog_id", "ortholog_score"), drop = FALSE],
    ortholog[, c("gene_id", "species"), drop = FALSE],
    by.x = "Gene_ID",
    by.y = "gene_id",
    all.x = TRUE
  )
  ortho_show <- ortho_show[!is.na(ortho_show$ortholog_id), , drop = FALSE]
  print(ortho_show)
  print(head(summarizeOrthologMatches(cand_ortho), 5))

  message("\nDone.")
  message("  Report: ", normalizePath(report_html, mustWork = FALSE))
  message("  Companion store: ", normalizePath(
    paste0(sub("\\.gds$", "", gds_copy, ignore.case = TRUE), ".lazygas"),
    mustWork = FALSE
  ))
  message("  Re-run in R: runMvpDemo()")
  message("  Shiny: runLazyGasShiny() after loading the GDS from ", out_dir)

  invisible(list(
    lg = lg,
    gff = gff,
    snpeff = snpeff,
    ann = ann,
    ortholog = ortholog,
    qc = qc,
    multitrait = mt,
    conditional = cond,
    credible_set = cred,
    candidates_ortholog = cand_ortho,
    report_html = report_html,
    out_dir = normalizePath(out_dir, mustWork = FALSE)
  ))
}

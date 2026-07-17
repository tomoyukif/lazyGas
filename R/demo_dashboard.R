#' Run makeInteractiveDashboard() demo
#'
#' Builds a full lazyGas workflow on bundled \code{sample.gds}, lists synthetic
#' candidate genes, and writes an HTML dashboard.
#'
#' @param out_dir Output directory (created if missing). Default: \code{"demo_output"}
#'   under the current working directory.
#' @param pheno_name Phenotype column name used in the analysis.
#' @param overwrite_gds If \code{TRUE}, rebuild the companion store on a copy of
#'   \code{sample.gds}.
#'
#' @return Invisibly, the path to the generated HTML file.
#' @export
runDashboardDemo <- function(out_dir = file.path(getwd(), "demo_output"),
                             pheno_name = "Demo_trait",
                             overwrite_gds = TRUE) {
  if (!requireNamespace("GBScleanR", quietly = TRUE)) {
    stop("Package 'GBScleanR' is required.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required to load the demo GFF.", call. = FALSE)
  }

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

  gds_src <- file.path(extdata, "sample.gds")
  pheno_fn <- file.path(extdata, "pheno.csv")
  gff_fn <- file.path(extdata, "demo_annotation.gff")
  vcf_fn <- file.path(extdata, "demo_snpeff.vcf")
  ann_fn <- file.path(extdata, "demo_ann.csv")

  for (fn in c(gds_src, pheno_fn, gff_fn, vcf_fn, ann_fn)) {
    if (!file.exists(fn)) {
      stop("Demo file not found: ", fn, call. = FALSE)
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gds_copy <- file.path(out_dir, "sample.gds")
  file.copy(gds_src, gds_copy, overwrite = TRUE)

  message("Loading demo annotation ...")
  gff <- rtracklayer::import.gff(gff_fn)
  ann <- read.csv(ann_fn, stringsAsFactors = FALSE)

  snpeff_gds_fn <- file.path(out_dir, "demo_snpeff.gds")
  if (!file.exists(snpeff_gds_fn)) {
    message("Converting demo SnpEff VCF to GDS ...")
    snpeff2gds(vcf_fn = vcf_fn, out_fn = snpeff_gds_fn, verbose = FALSE)
  }
  snpeff <- open_snpeff(snpeff_gds_fn)

  message("Running association pipeline on sample.gds ...")
  lg <- buildLazyGas(
    gds_fn = gds_copy,
    load_filter = TRUE,
    overwrite = overwrite_gds
  )
  on.exit(try(closeGDS(lg, verbose = FALSE), silent = TRUE), add = TRUE)

  pheno <- read.csv(pheno_fn)
  lg <- assignPheno(object = lg, pheno = pheno, rename = pheno_name)

  scanAssoc(object = lg, formula = "add", geno_format = "dosage")
  callPeakBlock(
    object = lg,
    signif = 0.05,
    threshold = 0.8,
    limit_peakcall = 3L
  )
  recalcAssoc(object = lg, refine_position = FALSE)

  message("Listing candidate genes ...")
  listCandidate(
    object = lg,
    gff = gff,
    snpeff = snpeff,
    ann = ann,
    recalc = TRUE
  )

  candidate <- lazyData(lg, dataset = "candidate", pheno = pheno_name)
  message("Candidate genes found:")
  print(candidate[, c("peak_ID", "Gene_ID", "Gene_chr", "negLog10P", "Name")])

  out_fn <- file.path(out_dir, "lazygas_dashboard.html")
  message("Writing dashboard to ", out_fn, " ...")
  makeInteractiveDashboard(
    object = lg,
    pheno = pheno_name,
    gff = gff,
    out_fn = out_fn,
    snpeff = snpeff,
    ann = ann,
    recalc = TRUE,
    what = c("scan_png", "peakcall", "recalc", "candidate")
  )

  message("Done. Open in a browser:\n  ", normalizePath(out_fn, mustWork = FALSE))
  invisible(out_fn)
}

.skip_without_sample_gds <- function() {
  gds_fn <- system.file("extdata", "sample.gds", package = "lazyGas")
  if (!nzchar(gds_fn)) {
    gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
  }
  if (!nzchar(gds_fn)) {
    testthat::skip("sample.gds not found (lazyGas or GBScleanR)")
  }
  gds_fn
}

.skip_without_demo_extdata <- function() {
  gff_fn <- system.file("extdata", "demo_annotation.gff", package = "lazyGas")
  vcf_fn <- system.file("extdata", "demo_snpeff.vcf", package = "lazyGas")
  ann_fn <- system.file("extdata", "demo_ann.csv", package = "lazyGas")
  pheno_fn <- system.file("extdata", "pheno.csv", package = "lazyGas")
  if (!nzchar(gff_fn) || !file.exists(gff_fn)) {
    testthat::skip("demo_annotation.gff not in package extdata")
  }
  if (!nzchar(vcf_fn) || !file.exists(vcf_fn)) {
    testthat::skip("demo_snpeff.vcf not in package extdata")
  }
  if (!nzchar(ann_fn) || !file.exists(ann_fn)) {
    testthat::skip("demo_ann.csv not in package extdata")
  }
  if (!nzchar(pheno_fn) || !file.exists(pheno_fn)) {
    testthat::skip("pheno.csv not in package extdata")
  }
  list(gff_fn = gff_fn, vcf_fn = vcf_fn, ann_fn = ann_fn, pheno_fn = pheno_fn)
}

.close_lg <- function(lg) {
  if (!missing(lg) && inherits(lg, "LazyGas")) {
    try(GBScleanR::closeGDS(lg, verbose = FALSE), silent = TRUE)
  }
}

.copy_sample_gds <- function(gds_fn) {
  dest <- tempfile(fileext = ".gds")
  file.copy(gds_fn, dest, overwrite = TRUE)
  dest
}

.run_sample_pipeline <- function(pheno_name = "test_trait",
                                limit_peakcall = 3L,
                                overwrite = TRUE) {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = overwrite,
    lazygas_store = "parquet"
  )

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = pheno_name)

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  lazyGas::scanAssoc(
    object = lg,
    formula = "add + dom",
    conv_fun = conv_fun,
    geno_format = "dosage"
  )
  lazyGas::callPeakBlock(
    object = lg,
    signif = 0.05,
    threshold = 0.8,
    limit_peakcall = limit_peakcall
  )
  lazyGas::recalcAssoc(object = lg, refine_position = FALSE, n_threads = 1L)

  gff <- rtracklayer::import.gff(demo$gff_fn)
  ann <- read.csv(demo$ann_fn, stringsAsFactors = FALSE)
  snpeff_gds_fn <- tempfile(pattern = "lazygas_snpeff_", fileext = ".gds")
  on.exit(unlink(snpeff_gds_fn), add = TRUE)
  lazyGas::snpeff2gds(vcf_fn = demo$vcf_fn, out_fn = snpeff_gds_fn, verbose = FALSE)
  snpeff <- lazyGas::open_snpeff(snpeff_gds_fn)

  lazyGas::listCandidate(
    object = lg,
    gff = gff,
    snpeff = snpeff,
    ann = ann,
    recalc = TRUE
  )

  list(
    lg = lg,
    pheno_name = pheno_name,
    gff = gff,
    ann = ann,
    snpeff = snpeff,
    gds_fn = gds_fn
  )
}

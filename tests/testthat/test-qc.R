test_that("GWAS QC helpers work on scan results", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = TRUE)
  on.exit(.close_lg(lg), add = TRUE)

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "qc_trait")

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  suppressMessages(
    lazyGas::scanAssoc(
      object = lg,
      formula = "add + dom",
      conv_fun = conv_fun,
      geno_format = "dosage"
    )
  )

  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = "qc_trait")
  lambda <- lazyGas::calcGenomicInflation(scan_dat$P.model)
  expect_true(is.finite(lambda))
  expect_gt(lambda, 0)

  p_qq <- lazyGas::plotQQ(lg, pheno = "qc_trait")
  expect_s3_class(p_qq, "ggplot")

  qc <- lazyGas::summarizeGWASQC(lg, pheno = "qc_trait")
  expect_equal(qc$pheno, "qc_trait")
  expect_equal(qc$lambda_gc, lambda)
  expect_true("small_n_warning" %in% names(qc))

  qc2 <- lazyGas::lazyData(lg, dataset = "qc", kind = "summary")
  expect_equal(nrow(qc2), 1L)
})

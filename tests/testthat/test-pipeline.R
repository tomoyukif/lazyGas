test_that("runLazyGas runs core steps with resume", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "parquet"
  )
  on.exit(.close_lg(lg), add = TRUE)

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "pipe_trait")

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  gff <- rtracklayer::import.gff(demo$gff_fn)

  suppressMessages(
    lg <- lazyGas::runLazyGas(
      object = lg,
      steps = c("scan", "peakcall", "recalc"),
      resume = FALSE,
      formula = "add + dom",
      conv_fun = conv_fun,
      geno_format = "dosage",
      limit_peakcall = 2L,
      n_threads = 1L
    )
  )

  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = "pipe_trait")
  expect_true(!is.null(scan_dat))

  hist <- lazyGas::lazyData(lg, dataset = "pipeline")
  expect_true(is.data.frame(hist))
  expect_gt(nrow(hist), 0L)

  lg2 <- lazyGas::runLazyGas(
    object = lg,
    steps = c("scan"),
    resume = TRUE,
    formula = "add + dom",
    conv_fun = conv_fun,
    geno_format = "dosage"
  )
  expect_identical(lg2, lg)
})

test_that("buildLazyGas parquet store and scan round-trip", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "parquet"
  )
  on.exit(.close_lg(lg), add = TRUE)

  pheno_fn <- system.file("extdata", "pheno.csv", package = "lazyGas")
  skip_if_not(nzchar(pheno_fn), "pheno.csv not in package")
  pheno <- read.csv(pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "test_trait")

  lazyGas::scanAssoc(object = lg, formula = "add", geno_format = "dosage")

  scan_dat <- lazyGas::lazyData(object = lg, dataset = "scan", pheno = "test_trait")
  expect_true(is.data.frame(scan_dat))
  expect_true("FDR" %in% names(scan_dat))
  expect_true("negLog10P" %in% names(scan_dat))
  expect_equal(nrow(scan_dat), GBScleanR::nmar(lg))

  store_path <- lg@store$path
  expect_true(dir.exists(store_path))
  meta_path <- file.path(store_path, "meta.json")
  expect_true(file.exists(meta_path))
  meta <- jsonlite::fromJSON(meta_path)
  expect_equal(meta$lazygas_schema_version, 1L)
})

test_that("scanAssoc accepts fixed_effect without null_formula", {
  skip_if_not_installed("GBScleanR")
  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = TRUE)
  on.exit(.close_lg(lg), add = TRUE)

  pheno_fn <- system.file("extdata", "pheno.csv", package = "lazyGas")
  skip_if_not(nzchar(pheno_fn), "pheno.csv not in package")
  pheno <- read.csv(pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "trait1")

  n <- GBScleanR::nsam(lg)
  fe <- data.frame(pc1 = rnorm(n))
  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3)
  suppressMessages(
    expect_error(
      lazyGas::scanAssoc(
        object = lg,
        formula = "add + dom + pc1",
        fixed_effect = fe,
        conv_fun = conv_fun,
        geno_format = "dosage",
        null_formula = NULL
      ),
      NA
    )
  )
})

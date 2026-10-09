test_that("listCandidate always writes wide and simple candidates", {
  skip_if_not_installed("rtracklayer")
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

  demo <- .skip_without_demo_extdata()
  pheno_name <- "outlist_trait"
  pheno <- read.csv(demo$pheno_fn)
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = pheno_name)

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  suppressMessages({
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
      limit_peakcall = 2L
    )
    lazyGas::recalcAssoc(object = lg, refine_position = FALSE, n_threads = 1L)
  })

  gff <- rtracklayer::import.gff(demo$gff_fn)
  ann <- read.csv(demo$ann_fn, stringsAsFactors = FALSE)

  suppressMessages(
    lazyGas::listCandidate(
      object = lg,
      gff = gff,
      ann = ann,
      recalc = TRUE,
      gff_fn = demo$gff_fn,
      ann_fn = demo$ann_fn
    )
  )

  simple <- lazyGas::lazyData(lg, dataset = "simple_candidate", pheno = pheno_name)
  expect_true(is.data.frame(simple))
  expect_gt(nrow(simple), 0L)
  expect_true(all(c("Gene_ID", "Gene_chr", "Gene_start") %in% names(simple)))

  wide <- lazyGas::lazyData(lg, dataset = "candidate", pheno = pheno_name)
  expect_true(is.data.frame(wide))
  expect_gt(nrow(wide), 0L)

  stem <- sub("\\.gds$", "", gds_fn, ignore.case = TRUE)
  man_path <- file.path(paste0(stem, ".lazygas"), "inputs_manifest.json")
  expect_true(file.exists(man_path))
})

test_that("buildLazyGas rejects non-parquet store modes", {
  skip_if_not_installed("GBScleanR")
  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)
  expect_error(
    lazyGas::buildLazyGas(gds_fn = gds_fn, lazygas_store = "sqlite"),
    "parquet"
  )
})

test_that("assignPheno accepts sample_id", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)
  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE
  )
  on.exit(.close_lg(lg), add = TRUE)

  demo <- .skip_without_demo_extdata()
  pheno <- read.csv(demo$pheno_fn, stringsAsFactors = FALSE)
  names(pheno)[names(pheno) == "id"] <- "sample_id"
  expect_silent(
    lg2 <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "sid_trait")
  )
  expect_true("sid_trait" %in% lazyGas::getPheno(lg2)$pheno_names)
})

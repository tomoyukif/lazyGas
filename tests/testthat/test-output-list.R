test_that("listCandidate output_list=FALSE writes simple only", {
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
      output_list = FALSE
    )
  )

  simple <- lazyGas::lazyData(lg, dataset = "simple_candidate", pheno = pheno_name)
  expect_true(is.data.frame(simple))
  expect_gt(nrow(simple), 0L)
  expect_true(all(c("Gene_ID", "Gene_chr", "Gene_start") %in% names(simple)))

  wide <- lazyGas::lazyData(lg, dataset = "candidate", pheno = pheno_name)
  expect_null(wide)

  # Default TRUE restores legacy wide table
  suppressMessages(
    lazyGas::listCandidate(
      object = lg,
      gff = gff,
      ann = ann,
      recalc = TRUE,
      output_list = TRUE
    )
  )
  wide2 <- lazyGas::lazyData(lg, dataset = "candidate", pheno = pheno_name)
  expect_true(is.data.frame(wide2))
  expect_gt(nrow(wide2), 0L)
  simple2 <- lazyGas::lazyData(lg, dataset = "simple_candidate", pheno = pheno_name)
  expect_true(is.data.frame(simple2))
  expect_gt(nrow(simple2), 0L)
})

test_that("runLazyGas skips summary/dashboard when output_list=FALSE", {
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
  lg <- lazyGas::assignPheno(object = lg, pheno = pheno, rename = "pipe_outlist")
  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  gff <- rtracklayer::import.gff(demo$gff_fn)

  msgs <- capture.output(
    suppressMessages(
      lazyGas::runLazyGas(
        object = lg,
        steps = c("scan", "peakcall", "recalc", "candidate", "summary"),
        resume = FALSE,
        formula = "add + dom",
        conv_fun = conv_fun,
        geno_format = "dosage",
        limit_peakcall = 2L,
        n_threads = 1L,
        gff = gff,
        output_list = FALSE
      )
    ),
    type = "message"
  )
  hist <- lazyGas::lazyData(lg, dataset = "pipeline")
  expect_true(is.data.frame(hist))
  expect_true(grepl("summary", hist$steps_skipped[nrow(hist)], fixed = TRUE) ||
                grepl("summary", paste(msgs, collapse = "\n"), fixed = TRUE))

  simple <- lazyGas::lazyData(lg, "simple_candidate", "pipe_outlist")
  expect_true(is.data.frame(simple))
  expect_null(lazyGas::lazyData(lg, "candidate", "pipe_outlist"))
})

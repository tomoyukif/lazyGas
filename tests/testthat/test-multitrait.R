test_that("multitrait overview runs after full pipeline", {
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
  pheno$trait_b <- pheno$pheno + stats::rnorm(nrow(pheno), sd = 0.1)
  lg <- lazyGas::assignPheno(
    object = lg,
    pheno = pheno[, c("id", "pheno", "trait_b")],
    rename = c("trait_a", "trait_b")
  )

  conv_fun <- lazyGas::makeConvFun(geno_format = "dosage", n_levels = 3L)
  suppressMessages(
    lazyGas::scanAssoc(
      object = lg,
      formula = "add + dom",
      conv_fun = conv_fun,
      geno_format = "dosage"
    )
  )
  lazyGas::callPeakBlock(object = lg, signif = 0.05, threshold = 0.8, limit_peakcall = 2L)
  lazyGas::recalcAssoc(object = lg, n_threads = 1L)

  clusters <- lazyGas::clusterCrossTraitPeaks(lg, recalc = TRUE, store = TRUE)
  expect_true(is.data.frame(clusters))
  if (nrow(clusters) > 0L) {
    expect_true(all(c("cluster_id", "merge_reason", "trait") %in% names(clusters)))
  }

  mt <- lazyGas::plotMultiTraitOverview(lg, recalc = TRUE)
  expect_s3_class(mt$heatmap, "ggplot")
  expect_true(is.data.frame(mt$summary))
})

test_that("multitrait clustering merges peaks at the same locus", {
  skip_if_not_installed("GBScleanR")

  gds_fn <- .copy_sample_gds(.skip_without_sample_gds())
  on.exit(unlink(gds_fn), add = TRUE)

  lg <- lazyGas::buildLazyGas(
    gds_fn = gds_fn,
    load_filter = TRUE,
    overwrite = TRUE,
    lazygas_store = "parquet"
  )
  on.exit(.close_lg(lg), add = TRUE)

  peaks <- data.frame(
    trait = c("trait_a", "trait_b", "trait_c"),
    peak_ID = 1L,
    peak_variant_ID = "1378",
    chr = 2L,
    pos = 13920070L,
    negLog10P = c(10, 9, 8),
    peak_key = c("trait_a:1", "trait_b:1", "trait_c:1"),
    stringsAsFactors = FALSE
  )
  clusters <- lazyGas:::.multitrait_cluster_peaks(
    object = lg,
    peaks = peaks,
    dist_threshold = 500000L,
    r2_threshold = 0.6
  )
  expect_equal(length(unique(clusters$cluster_id)), 1L)
  expect_false(any(clusters$merge_reason == "singleton"))
})

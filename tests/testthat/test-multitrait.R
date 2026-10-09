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
    heat_df <- lazyGas:::.multitrait_heatmap_df(lg, clusters)
    expect_gt(length(levels(heat_df$chr)), length(unique(clusters$chr)))
    expect_equal(
      nrow(heat_df),
      length(unique(clusters$trait)) * length(levels(heat_df$chr))
    )
  }

  mt <- lazyGas::plotMultiTraitOverview(lg, recalc = TRUE)
  expect_s3_class(mt$heatmap, "ggplot")
  expect_true(is.data.frame(mt$summary))
  expect_equal(
    names(mt$shared_genes),
    c("cluster_id", "Gene_ID", "trait", "n_traits")
  )
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
  expect_equal(
    names(clusters),
    c(
      "cluster_id", "merge_reason", "trait", "peak_ID", "peak_variant_ID",
      "chr", "pos", "negLog10P"
    )
  )
})

test_that("multitrait empty shared_genes has stable columns", {
  empty <- lazyGas:::.multitrait_empty_shared_genes()
  expect_equal(
    names(empty),
    c("cluster_id", "Gene_ID", "trait", "n_traits")
  )
  expect_equal(nrow(empty), 0L)
})

test_that("CS overlap Jaccard and labels follow plan rules", {
  expect_equal(lazyGas:::.multitrait_cs_jaccard(c("a", "b"), c("b", "c")), 1 / 3)
  expect_equal(lazyGas:::.multitrait_cs_jaccard(character(), character()), NA_real_)

  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 2L,
      lead_a_in_cs_b = TRUE,
      lead_b_in_cs_a = FALSE,
      resolution_a = "high",
      resolution_b = "moderate",
      cs_fraction_a = 0.1,
      cs_fraction_b = 0.2
    ),
    "shared-compatible"
  )
  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 1L,
      lead_a_in_cs_b = FALSE,
      lead_b_in_cs_a = FALSE,
      resolution_a = "high",
      resolution_b = "high",
      cs_fraction_a = 0.1,
      cs_fraction_b = 0.1
    ),
    "possible-shared"
  )
  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 0L,
      lead_a_in_cs_b = FALSE,
      lead_b_in_cs_a = FALSE,
      resolution_a = "high",
      resolution_b = "moderate",
      cs_fraction_a = 0.1,
      cs_fraction_b = 0.1
    ),
    "independent-nearby"
  )
  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 5L,
      lead_a_in_cs_b = TRUE,
      lead_b_in_cs_a = TRUE,
      resolution_a = "low",
      resolution_b = "low",
      cs_fraction_a = 0.9,
      cs_fraction_b = 0.9
    ),
    "inconclusive"
  )
  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 3L,
      lead_a_in_cs_b = TRUE,
      lead_b_in_cs_a = FALSE,
      resolution_a = "low",
      resolution_b = "high",
      cs_fraction_a = 0.85,
      cs_fraction_b = 0.9
    ),
    "inconclusive"
  )
  expect_equal(
    lazyGas:::.multitrait_label_cs_overlap(
      n_overlap = 1L,
      lead_a_in_cs_b = FALSE,
      lead_b_in_cs_a = FALSE,
      resolution_a = "low",
      resolution_b = "high",
      cs_fraction_a = 0.5,
      cs_fraction_b = 0.2
    ),
    "possible-shared"
  )

  empty <- lazyGas:::.multitrait_empty_cs_overlap()
  expect_true(is.data.frame(empty))
  expect_equal(nrow(empty), 0L)
  expect_true(all(c(
    "cluster_id", "n_overlap", "jaccard", "share_label",
    "lead_a_in_cs_b", "lead_b_in_cs_a", "resolution_a", "resolution_b"
  ) %in% names(empty)))
})

test_that("compareCrossTraitCredibleSets labels synthetic CS pairs", {
  cs_a <- list(
    members = c("m1", "m2", "m3"),
    resolution = "high",
    cs_fraction = 0.1,
    n_in_set = 3L
  )
  cs_b <- list(
    members = c("m2", "m4"),
    resolution = "high",
    cs_fraction = 0.1,
    n_in_set = 2L
  )
  row <- lazyGas:::.multitrait_cs_overlap_row(
    cluster_id = 1L,
    trait_a = "trait_a",
    peak_ID_a = 1L,
    lead_a = "m1",
    cs_a = cs_a,
    trait_b = "trait_b",
    peak_ID_b = 2L,
    lead_b = "m2",
    cs_b = cs_b
  )
  expect_equal(row$n_overlap, 1L)
  expect_true(row$lead_b_in_cs_a)
  expect_false(row$lead_a_in_cs_b)
  expect_equal(row$share_label, "shared-compatible")
})

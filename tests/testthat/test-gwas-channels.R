test_that("gwas marker channel scores follow dist and peak formulas", {
  blocks <- data.frame(
    variant_ID = c(1L, 2L, 3L),
    peak_variant_ID = 2L,
    Chr = "1",
    Pos = c(100L, 200L, 400L),
    dist2peak = c(-100, 0, 200),
    negLog10P = c(4, 8, 2),
    stringsAsFactors = FALSE
  )
  ms <- lazyGas:::.gwas_marker_channel_scores(blocks)
  expect_equal(ms$d_max, 200)
  expect_equal(ms$lead_nlp, 8)
  # lead: (200-0)/200 = 1; far: (200-200)/200 = 0; mid: (200-100)/200 = 0.5
  expect_equal(ms$markers$dist_score, c(0.5, 1, 0))
  # nlp/lead: 4/8, 8/8, 2/8
  expect_equal(ms$markers$peak_score, c(0.5, 1, 0.25))
})

test_that("gwas_dist/gwas_peak take max over attributed markers", {
  cand <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("gA", "gB"),
    Gene_chr = "1",
    Gene_start = c(90L, 390L),
    Gene_end = c(110L, 410L),
    dist2peak = c(-110, 190),
    negLog10P = c(4, 2),
    stringsAsFactors = FALSE
  )
  # No LazyGas object: fallback uses candidate rows as synthetic markers
  batch <- lazyGas:::.evidence_gwas_channels_batch(
    gene_ids = c("gA", "gB"),
    candidate = cand,
    object = NULL,
    pheno_name = NULL
  )
  # D_max among |dist2peak| = 190; gA: (190-110)/190; gB: (190-190)/190 = 0
  expect_equal(batch$gA$gwas_dist$score, (190 - 110) / 190, tolerance = 1e-8)
  expect_equal(batch$gB$gwas_dist$score, 0, tolerance = 1e-8)
  # lead = max nlp = 4; gA: 4/4=1; gB: 2/4=0.5
  expect_equal(batch$gA$gwas_peak$score, 1, tolerance = 1e-8)
  expect_equal(batch$gB$gwas_peak$score, 0.5, tolerance = 1e-8)
})

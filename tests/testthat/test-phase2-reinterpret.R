test_that("ann_specificity distinguishes strong vs weak-only", {
  strong <- c("leaf width", "lateral leaf", "amylose")
  weak <- c("leaf", "development", "protein", "domain")
  expect_equal(
    lazyGas:::.ann_specificity("Control of lateral leaf growth", strong, weak),
    1.0
  )
  expect_equal(
    lazyGas:::.ann_specificity("S-domain receptor like kinase in leaf tissue", strong, weak),
    0.2
  )
  expect_equal(
    lazyGas:::.ann_specificity("Granule-bound starch synthase; amylose synthesis", strong, weak),
    0.67
  )
  expect_equal(
    lazyGas:::.ann_specificity("leaf width QTL candidate protein", strong, weak),
    1.0
  )
  expect_equal(lazyGas:::.ann_specificity("hypothetical conserved gene", strong, weak), 0)
})

test_that("function champion requires strong hit; weak-only yields NA", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("g_weak", "g_strong", "g_far"),
    dist2peak = c(100, 500000, 50),
    HIGH = c(8, 0, 2),
    MODERATE = c(1, 1, 1),
    score_annotation = c(1.0, 0.8, 0.1),
    score_snpeff = c(1, 0.3, 1),
    score_gwas = c(0.9, 0.8, 0.9),
    score_expression = c(1, 0, 0.5),
    composite_score = c(0.99, 0.70, 0.60),
    RAP_Note = c(
      "S-domain receptor kinase; leaf drought response",
      "Serine protease, Control of lateral leaf growth",
      "Unknown"
    ),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq_test_f",
      text = "leaf width leaf morphology",
      trait_keywords = c("leaf width", "leaf morphology", "lateral leaf"),
      synonym_phrases = list(),
      related_phrases = list(list(phrase = "photosynthesis", from_user = "leaf")),
      tissues = "leaf",
      developmental_stage = character(),
      conditions = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(
    ranked,
    query = q,
    rules = reinterpretDefaultRules(dist_gate_mode = "off")
  )
  expect_s3_class(ri, "PeakReinterpretation")
  expect_equal(champion(ri, 1, "function"), "g_strong")
  expect_equal(champion(ri, 1, "composite"), "g_weak")
  expect_equal(champion(ri, 1, "proximity"), "g_far")
  expect_gt(nrow(conflictTable(ri)), 0)

  # related-only phrase must not create a strong seat by itself
  ranked2 <- data.frame(
    peak_ID = 1L,
    Gene_ID = "g_rel",
    dist2peak = 10,
    HIGH = 0,
    MODERATE = 0,
    score_annotation = 1,
    score_snpeff = 0,
    score_gwas = 0.5,
    score_expression = 0,
    composite_score = 0.5,
    RAP_Note = "photosynthesis related protein",
    stringsAsFactors = FALSE
  )
  ri2 <- reinterpretPeakCandidates(ranked2, query = q)
  expect_true(is.na(champion(ri2, 1, "function")))
  ch <- ri2$peaks[["1"]]$champions
  fun <- ch[ch$view == "function", , drop = FALSE]
  expect_equal(fun$reason_code[1], "no_eligible")
  expect_match(fun$reason_detail[1], "no function hypothesis", fixed = TRUE)
})

test_that("distance gate off keeps far strong genes in function view", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("near_weak", "far_strong"),
    dist2peak = c(100, 2e6),
    HIGH = c(1, 0),
    MODERATE = c(0, 0),
    score_annotation = c(0.9, 0.85),
    score_snpeff = c(1, 0.3),
    score_gwas = c(0.9, 0.8),
    score_expression = c(0, 0),
    composite_score = c(0.95, 0.7),
    RAP_Note = c(
      "leaf development kinase",
      "Control of lateral leaf growth"
    ),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq_dist",
      text = "leaf width",
      trait_keywords = c("leaf width", "lateral leaf"),
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = character(),
      developmental_stage = character(),
      conditions = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(
    ranked,
    query = q,
    rules = reinterpretDefaultRules(dist_gate_mode = "off")
  )
  expect_equal(champion(ri, 1, "function"), "far_strong")
  gs <- geneSetForReport(ri, 1, mode = "union_champions")
  expect_true("far_strong" %in% gs)
  expect_true("near_weak" %in% gs)
  html <- formatPeakReinterpretation(ri, peak_id = 1, format = "html")
  expect_match(html, "far_strong", fixed = TRUE)
  expect_match(html, "<th>view</th>", fixed = TRUE)
})

test_that("weak-only peak yields empty function champion", {
  ranked <- data.frame(
    peak_ID = 2L,
    Gene_ID = c("a", "b"),
    dist2peak = c(10, 20),
    HIGH = c(1, 0),
    MODERATE = c(0, 1),
    score_annotation = c(1, 0.5),
    score_snpeff = c(1, 0.5),
    score_gwas = c(0.8, 0.7),
    score_expression = c(0.2, 0),
    composite_score = c(0.9, 0.4),
    RAP_Note = c("leaf protein domain family", "plant development"),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq_weak",
      text = "leaf width",
      trait_keywords = c("leaf width"),
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = character(),
      developmental_stage = character(),
      conditions = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(ranked, query = q)
  expect_true(is.na(champion(ri, 2, "function")))
  df <- as.data.frame(ri)
  expect_true(any(df$view == "function" & is.na(df$Gene_ID)))
})

test_that("impact champion uses HIGH_at_var when present", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("g_window_only", "g_at_var"),
    dist2peak = c(10, 100),
    HIGH = c(9, 0),
    MODERATE = c(0, 0),
    HIGH_at_var = c(0, 1),
    MODERATE_at_var = c(0, 2),
    score_annotation = c(0.1, 0.1),
    score_snpeff = c(1, 0.5),
    score_gwas = c(0.5, 0.5),
    score_expression = c(0, 0),
    composite_score = c(0.5, 0.5),
    RAP_Note = c("x", "y"),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq",
      text = "trait",
      trait_keywords = "trait",
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = character(),
      developmental_stage = character(),
      conditions = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(ranked, query = q)
  expect_identical(champion(ri, 1, "impact"), "g_at_var")
  df <- as.data.frame(ri)
  imp <- df[df$view == "impact" & !is.na(df$Gene_ID), , drop = FALSE]
  expect_match(imp$reason_detail[1L], "HIGH=1")
})

test_that("snpeff_for_peak keeps peak-block-linked rows only", {
  snp <- data.frame(
    Gene_ID = c("G1", "G1", "G1", "G1"),
    Chr = "1",
    Pos = c(10, 10, 20, 30),
    Annotation_Impact = c("HIGH", "HIGH", "MODERATE", "LOW"),
    Allele = "A",
    Annotation = "x",
    Feature_ID = c("t1", "t1", "t2", "t3"),
    negLog10P = c(5, "", 4, NA),
    stringsAsFactors = FALSE
  )
  markers <- data.frame(Chr = "1", Pos = c(10, 20), stringsAsFactors = FALSE)
  out <- lazyGas:::.snpeff_for_peak(snpeff = snp, peak_markers = markers)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$Pos %in% c(10, 20)))
  expect_true(all(is.finite(suppressWarnings(as.numeric(out$negLog10P)))))
  expect_identical(sort(out$Annotation_Impact), c("HIGH", "MODERATE"))
})

test_that("function champion uses score_annotation threshold only", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("g_word", "g_action", "g_far"),
    dist2peak = c(100, 500000, 50),
    HIGH = c(8, 0, 2),
    MODERATE = c(1, 1, 1),
    score_annotation = c(0.25, 1.0, 0.1),
    score_snpeff = c(1, 0.3, 1),
    score_gwas_dist = c(0.9, 0.8, 0.9),
    score_gwas_peak = c(0.9, 0.8, 0.9),
    score_expression = c(1, 0, 0.5),
    composite_score = c(0.99, 0.70, 0.60),
    RAP_Note = c("leaf drought response", "regulates flowering", "Unknown"),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq_test_f",
      trait_text = "leaf width leaf morphology",
      trait_keywords = c("leaf width", "leaf morphology"),
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = "leaf",
      developmental_stage = character(),
      conditions = character(),
      action_phrases = character(),
      noun_phrases = character(),
      words = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(
    ranked,
    query = q,
    rules = reinterpretDefaultRules(
      dist_gate_mode = "off",
      min_annotation_score = 0.5
    )
  )
  expect_s3_class(ri, "PeakReinterpretation")
  expect_equal(champion(ri, 1, "function"), "g_action")
  expect_equal(champion(ri, 1, "composite"), "g_word")
  expect_equal(champion(ri, 1, "proximity"), "g_far")
  expect_gt(nrow(conflictTable(ri)), 0)

  ranked2 <- data.frame(
    peak_ID = 1L,
    Gene_ID = "g_low",
    dist2peak = 10,
    HIGH = 0,
    MODERATE = 0,
    score_annotation = 0.25,
    score_snpeff = 0,
    score_gwas_dist = 0.5,
    score_gwas_peak = 0.5,
    score_expression = 0,
    composite_score = 0.5,
    RAP_Note = "photosynthesis related protein",
    stringsAsFactors = FALSE
  )
  ri2 <- reinterpretPeakCandidates(
    ranked2,
    query = q,
    rules = reinterpretDefaultRules(min_annotation_score = 0.5)
  )
  expect_true(is.na(champion(ri2, 1, "function")))
  ch <- ri2$peaks[["1"]]$champions
  fun <- ch[ch$view == "function", , drop = FALSE]
  expect_equal(fun$reason_code[1], "no_eligible")
  expect_match(fun$reason_detail[1], "score_annotation", fixed = TRUE)
})

test_that("distance gate off keeps far high-annotation genes in function view", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("near_low", "far_high"),
    dist2peak = c(100, 2e6),
    HIGH = c(1, 0),
    MODERATE = c(0, 0),
    score_annotation = c(0.25, 1.0),
    score_snpeff = c(1, 0.3),
    score_gwas_dist = c(0.9, 0.8),
    score_gwas_peak = c(0.9, 0.8),
    score_expression = c(0, 0),
    composite_score = c(0.95, 0.7),
    RAP_Note = c("leaf development kinase", "regulates flowering"),
    stringsAsFactors = FALSE
  )
  q <- structure(
    list(
      query_id = "pq_dist",
      trait_text = "leaf width",
      trait_keywords = c("leaf width"),
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = character(),
      developmental_stage = character(),
      conditions = character(),
      action_phrases = character(),
      noun_phrases = character(),
      words = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
  ri <- reinterpretPeakCandidates(
    ranked,
    query = q,
    rules = reinterpretDefaultRules(dist_gate_mode = "off", min_annotation_score = 0.5)
  )
  expect_equal(champion(ri, 1, "function"), "far_high")
})

test_that("geneSetForReport unions champions", {
  ranked <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("a", "b", "c"),
    dist2peak = c(10, 20, 30),
    HIGH = c(1, 0, 0),
    MODERATE = c(0, 1, 0),
    score_annotation = c(1, 0.5, 0),
    score_snpeff = c(1, 0.5, 0),
    score_gwas_dist = c(0.8, 0.7, 0.1),
    score_gwas_peak = c(0.8, 0.7, 0.1),
    score_expression = c(0, 1, 0),
    composite_score = c(0.9, 0.85, 0.1),
    Description = c("regulates flowering", "flowering time", "other"),
    stringsAsFactors = FALSE
  )
  q <- phenotypeQuery("flowering", use_llm = FALSE, noun_phrases = "flowering time")
  ri <- reinterpretPeakCandidates(
    ranked,
    query = q,
    rules = reinterpretDefaultRules(min_annotation_score = 0.5)
  )
  gs <- geneSetForReport(ri, peak_id = 1L)
  expect_true(length(gs) >= 1L)
  expect_true(all(gs %in% ranked$Gene_ID))
})

test_that("legacy min_specificity alias still sets threshold", {
  rules <- reinterpretDefaultRules(min_specificity = 0.75)
  expect_equal(rules$min_annotation_score, 0.75)
  expect_equal(rules$min_specificity, 0.75)
})

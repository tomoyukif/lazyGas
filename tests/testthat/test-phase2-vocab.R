.make_pq <- function(trait_text,
                     trait_keywords,
                     synonym_phrases = list(),
                     related_phrases = list(),
                     tissues = character(),
                     stage = character(),
                     conditions = character(),
                     parsed_by = "llm") {
  structure(
    list(
      trait_text = trait_text,
      trait_keywords = trait_keywords,
      tissues = tissues,
      developmental_stage = stage,
      conditions = conditions,
      species = NULL,
      synonym_phrases = synonym_phrases,
      related_phrases = related_phrases,
      query_id = "pq_test",
      parsed_by = parsed_by
    ),
    class = c("PhenotypeQuery", "list")
  )
}

test_that("A1 TRUE normalize drops species and gene symbols", {
  skip_if_not_installed("lazyGas")
  parsed <- list(
    trait_keywords = c("Heading date", "HESO1", "rice", "flowering time"),
    tissues = "leaf",
    developmental_stage = character(),
    conditions = character(),
    species = NULL
  )
  out <- lazyGas:::.phenotype_query_normalize_true(
    parsed = parsed,
    trait_text = "Heading date / flowering time HESO1 in rice",
    tissues = "leaf",
    species = NULL
  )
  expect_equal(out$parsed_by, "llm")
  expect_true("heading date" %in% out$trait_keywords)
  expect_true("flowering time" %in% out$trait_keywords)
  expect_false("heso1" %in% out$trait_keywords)
  expect_false("rice" %in% out$trait_keywords)
  expect_true("leaf" %in% out$tissues)
  expect_false("leaf" %in% out$trait_keywords)
})

test_that("A1 FALSE path leaves synonym/related empty", {
  skip_if_not_installed("lazyGas")
  q <- phenotypeQuery("amylose content", use_llm = FALSE)
  expect_equal(q$parsed_by, "fallback")
  expect_equal(q$synonym_phrases, list())
  expect_equal(q$related_phrases, list())
})

test_that("A2 linked phrases normalize from_user and drop gene tokens", {
  skip_if_not_installed("lazyGas")
  items <- list(
    list(phrase = "Starch Synthase", from_user = "Amylose content", source = "llm"),
    list(phrase = "Waxy", from_user = "amylose content", source = "llm"),
    list(phrase = "Os06g0133000", from_user = "amylose content")
  )
  out <- lazyGas:::.phenotype_normalize_linked_phrases(
    items,
    trait_text = "amylose content",
    user_phrases = "amylose content"
  )
  phrases <- vapply(out, `[[`, character(1L), "phrase")
  expect_true("starch synthase" %in% phrases)
  expect_false("waxy" %in% phrases)
  expect_false(any(grepl("^os06g", phrases, ignore.case = TRUE)))
  fu <- out[[which(phrases == "starch synthase")]]$from_user
  expect_equal(fu, "amylose content")
})

test_that("A2/A4 family unmatched and related from_user", {
  skip_if_not_installed("lazyGas")
  q <- .make_pq(
    trait_text = "amylose content / leaf width",
    trait_keywords = c("amylose content", "leaf width"),
    related_phrases = list(
      list(
        phrase = "starch synthase",
        from_user = "amylose content",
        source = "llm"
      )
    )
  )
  # Only related hit for amylose; leaf width family empty
  det <- lazyGas:::.annotation_keyword_family_details(
    text = "Granule-bound starch synthase; seed starch metabolism",
    query = q
  )[[1L]]
  expect_equal(det$matched_keywords, character())
  expect_equal(length(det$matched_related), 1L)
  expect_equal(det$matched_related[[1L]]$phrase, "starch synthase")
  expect_equal(det$matched_related[[1L]]$from_user, "amylose content")
  expect_equal(det$unmatched_keywords, "leaf width")
  expect_false("starch synthase" %in% det$unmatched_keywords)

  terms <- lazyGas:::.keyword_terms_from_query(q)
  expect_true("starch synthase" %in% terms)
  expect_true("amylose content" %in% terms)
})

test_that("A5 skeleton: Waxy related hit; HESO1/NAL1 no related; gene scrub", {
  skip_if_not_installed("lazyGas")

  # Gene symbols scrubbed from TRUE keywords
  scrub <- lazyGas:::.phenotype_filter_true_keywords(
    c("heading date", "HESO1", "NAL1", "Waxy", "Os01g0846450"),
    trait_text = "heading date HESO1 NAL1"
  )
  expect_equal(scrub, "heading date")

  waxy_q <- .make_pq(
    trait_text = "amylose content Waxy Os06g0133000",
    trait_keywords = "amylose content",
    related_phrases = list(
      list(
        phrase = "starch synthase",
        from_user = "amylose content",
        source = "llm"
      )
    )
  )
  waxy_ann <- paste(
    "Granule-bound starch synthase; amylose synthesis; Wx protein"
  )
  waxy_det <- lazyGas:::.annotation_keyword_family_details(
    text = waxy_ann,
    query = waxy_q
  )[[1L]]
  expect_equal(waxy_det$matched_related[[1L]]$phrase, "starch synthase")
  expect_equal(waxy_det$matched_related[[1L]]$from_user, "amylose content")
  # Query terms must not include gene self-match tokens
  expect_false("waxy" %in% lazyGas:::.keyword_terms_from_query(waxy_q))

  heso_q <- .make_pq(
    trait_text = "heading date HESO1",
    trait_keywords = "heading date",
    related_phrases = list()
  )
  heso_ann <- paste(
    "Homolog of Arabidopsis thaliana HEN1 suppressor 1, Heading date",
    "expressed protein miRNA uridylyltransferase HESO1"
  )
  heso_det <- lazyGas:::.annotation_keyword_family_details(
    text = heso_ann,
    query = heso_q
  )[[1L]]
  expect_equal(heso_det$matched_related, list())
  expect_true("heading date" %in% heso_det$matched_keywords)
  expect_equal(heso_det$unmatched_keywords, character())

  nal_q <- .make_pq(
    trait_text = "leaf width NAL1",
    trait_keywords = "leaf width",
    related_phrases = list()
  )
  nal_det <- lazyGas:::.annotation_keyword_family_details(
    text = "NAL1 regulates leaf width and photosynthesis",
    query = nal_q
  )[[1L]]
  expect_equal(nal_det$matched_related, list())
  expect_true("leaf width" %in% nal_det$matched_keywords)

  # FALSE path: no related in evidence
  false_q <- phenotypeQuery("amylose content", use_llm = FALSE)
  false_det <- lazyGas:::.annotation_keyword_family_details(
    text = waxy_ann,
    query = false_q
  )[[1L]]
  expect_equal(false_det$matched_related, list())
})

test_that("fallback query forces llm relevance off in annotation batch", {
  skip_if_not_installed("lazyGas")
  cand <- data.frame(
    Gene_ID = "g1",
    Description = "fruit weight development",
    stringsAsFactors = FALSE
  )
  q <- phenotypeQuery("fruit weight", use_llm = FALSE)
  out <- lazyGas:::.evidence_annotation_batch(
    gene_ids = "g1",
    candidate = cand,
    query = q,
    search_text = "fruit weight",
    use_keyword = TRUE,
    use_llm_relevance = TRUE
  )
  expect_equal(out$g1$details$llm_relevance_score, 0)
  expect_true(out$g1$details$keyword_score > 0)
})

test_that("A2 linked phrase cap is 20", {
  skip_if_not_installed("lazyGas")
  expect_equal(lazyGas:::.phenotype_linked_phrase_cap(), 20L)
  many <- lapply(seq_len(30), function(i) {
    list(phrase = paste0("phrase_", i), from_user = "fruit weight", source = "llm")
  })
  out <- lazyGas:::.phenotype_normalize_linked_phrases(
    many,
    trait_text = "fruit weight",
    user_phrases = "fruit weight"
  )
  expect_lte(length(out), 20L)
})

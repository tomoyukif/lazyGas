make_test_candidate <- function() {
  data.frame(
    peak_ID = rep(1L, 4L),
    Gene_ID = c("g1", "g2", "g3", "g4"),
    Gene_chr = "1",
    Gene_start = 1:4,
    dist2peak = c(100, 200, 300, 400),
    negLog10P = 3,
    Description = c(
      "fruit weight development and ripening",
      "root hair elongation",
      "seed weight regulation",
      "photosystem II reaction center"
    ),
    GO = c(
      "fruit development",
      "root development",
      "seed development",
      "photosynthesis"
    ),
    stringsAsFactors = FALSE
  )
}

test_that("searchCandidateGenes is exported", {
  skip_if_not_installed("lazyGas")
  expect_true("searchCandidateGenes" %in% getNamespaceExports("lazyGas"))
})

.load_search_fn <- function() {
  env <- .load_search_env()
  env$searchCandidateGenes
}

.load_search_env <- function() {
  fn_path <- file.path(
    testthat::test_path(),
    "..", "..", "R", "search_candidate_genes.R"
  )
  fn_path <- normalizePath(fn_path, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$lazyData <- function(object, dataset, pheno) NULL
  source(fn_path, local = env)
  env
}

test_that("NULL ann_cols uses all non-core annotation columns", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  out <- searchCandidateGenes(
    candidate = cand,
    query = "photosynthesis",
    mode = "keyword",
    keyword_match = "any"
  )
  expect_true("g4" %in% out$Gene_ID)
  expect_equal(
    attr(out, "searchCandidateGenes")$ann_cols,
    c("Description", "GO")
  )
})

test_that("keyword search requires all tokens when keyword_match is all", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  out <- searchCandidateGenes(
    candidate = cand,
    query = "fruit weight",
    ann_cols = "Description",
    mode = "keyword",
    keyword_match = "all"
  )
  expect_equal(out$Gene_ID, "g1")
  expect_equal(out$keyword_score, 1)
})

test_that("keyword any matches partial token overlap", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  out <- searchCandidateGenes(
    candidate = cand,
    query = "fruit photosynthesis",
    ann_cols = c("Description", "GO"),
    mode = "keyword",
    keyword_match = "any"
  )
  expect_true("g1" %in% out$Gene_ID)
  expect_true("g4" %in% out$Gene_ID)
})

test_that("dedupe_genes keeps best score per Gene_ID", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  cand <- rbind(cand, cand[1, , drop = FALSE])
  cand$peak_ID <- c(cand$peak_ID[1:4], 2L)
  cand$dist2peak[5] <- 50
  out <- searchCandidateGenes(
    candidate = cand,
    query = "fruit",
    ann_cols = "Description",
    mode = "keyword",
    keyword_match = "any",
    dedupe_genes = TRUE
  )
  expect_equal(sum(out$Gene_ID == "g1"), 1L)
})

test_that("missing ann_cols errors", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  expect_error(
    searchCandidateGenes(
      candidate = cand,
      query = "fruit",
      ann_cols = "NotAColumn",
      mode = "keyword"
    ),
    "not found"
  )
})

test_that("semantic mode is rejected after LSA removal", {
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  expect_error(
    searchCandidateGenes(
      candidate = cand,
      query = "fruit",
      ann_cols = "Description",
      mode = "semantic"
    ),
    "keyword"
  )
})

test_that("stopword in does not match containing", {
  env <- .load_search_env()
  ph_ann <- "Similar to pleckstrin homology domain-containing protein 1."
  score <- env$.keyword_match_score(
    text = ph_ann,
    query = "heading date / flowering time HESO1 in rice flowering",
    match = "any"
  )
  expect_equal(score, 0)

  heso_ann <- paste(
    "Homolog of Arabidopsis thaliana HEN1 suppressor 1, Heading date",
    "expressed protein miRNA uridylyltransferase HESO1",
    "TO:0000137 - days to heading, TO:0002616 - flowering time"
  )
  score_h <- env$.keyword_match_score(
    text = heso_ann,
    query = "heading date flowering time",
    match = "any"
  )
  expect_true(score_h > 0.5)

  det <- env$.keyword_match_details(
    text = ph_ann,
    terms = c("heading date", "flowering time", "flowering")
  )[[1L]]
  expect_equal(det$keyword_score, 0)
  expect_equal(det$matched_keywords, character())
  expect_true(all(c("heading date", "flowering time", "flowering") %in% det$unmatched_keywords))
})

test_that("phrase match prefers multi-word trait keywords", {
  env <- .load_search_env()
  q <- list(
    trait_text = "heading date / flowering time HESO1 in rice",
    trait_keywords = c("heading date", "flowering time"),
    tissues = character(),
    developmental_stage = "flowering",
    conditions = character()
  )
  class(q) <- c("PhenotypeQuery", "list")
  terms <- env$.keyword_terms_from_query(q)
  expect_true("heading date" %in% terms)
  expect_true("flowering time" %in% terms)
  expect_false("in" %in% terms)
})

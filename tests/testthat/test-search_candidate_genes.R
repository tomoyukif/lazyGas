.load_keyword_env <- function() {
  fn_path <- file.path(
    testthat::test_path(),
    "..", "..", "R", "search_candidate_genes.R"
  )
  fn_path <- normalizePath(fn_path, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  source(fn_path, local = env)
  env
}

test_that("searchCandidateGenes is not exported", {
  ns_path <- normalizePath(
    file.path(testthat::test_path(), "..", "..", "NAMESPACE"),
    mustWork = TRUE
  )
  ns <- readLines(ns_path)
  expect_false(any(grepl("export\\(searchCandidateGenes\\)", ns)))
})

test_that("resolve_ann_cols errors on missing columns", {
  env <- .load_keyword_env()
  cand <- data.frame(Description = "fruit", stringsAsFactors = FALSE)
  expect_error(
    env$.resolve_ann_cols(candidate = cand, ann_cols = "NotAColumn"),
    "not found"
  )
})

test_that("stopword in does not match containing", {
  env <- .load_keyword_env()
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
  env <- .load_keyword_env()
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

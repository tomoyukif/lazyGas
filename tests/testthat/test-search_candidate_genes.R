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
  if (requireNamespace("lazyGas", quietly = TRUE)) {
    return(get("searchCandidateGenes", envir = asNamespace("lazyGas")))
  }
  fn_path <- file.path(
    testthat::test_path(),
    "..", "..", "R", "search_candidate_genes.R"
  )
  fn_path <- normalizePath(fn_path, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  env$lazyData <- function(object, dataset, pheno) NULL
  source(fn_path, local = env)
  env$searchCandidateGenes
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

test_that("semantic search runs when text2vec is available", {
  skip_if_not_installed("text2vec")
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  out <- searchCandidateGenes(
    candidate = cand,
    query = "grain yield fruit mass",
    ann_cols = c("Description", "GO"),
    mode = "semantic",
    min_score = 0,
    top_n = 3L
  )
  expect_true(nrow(out) >= 1L)
  expect_true(all(c("semantic_score", "match_score") %in% names(out)))
  expect_true("g1" %in% out$Gene_ID || "g3" %in% out$Gene_ID)
})

test_that("semantic mode errors without text2vec", {
  if (requireNamespace("text2vec", quietly = TRUE)) {
    skip("text2vec is installed")
  }
  searchCandidateGenes <- .load_search_fn()
  cand <- make_test_candidate()
  expect_error(
    searchCandidateGenes(
      candidate = cand,
      query = "fruit",
      ann_cols = "Description",
      mode = "semantic"
    ),
    "text2vec"
  )
})

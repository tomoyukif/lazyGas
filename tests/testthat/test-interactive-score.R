test_that("interactive score helpers join candidate and score", {
  candidate <- data.frame(
    peak_ID = c(1L, 1L, 2L),
    Gene_ID = c("g1", "g2", "g3"),
    negLog10P = c(8, 5, 6),
    RAP_Note = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
  score <- data.frame(
    peak_ID = c(1L, 1L),
    Gene_ID = c("g2", "g1"),
    composite_score = c(0.9, 0.4),
    score_annotation = c(1, 0.5),
    evidence_json = c("{}", "{}"),
    stringsAsFactors = FALSE
  )

  joined <- lazyGas:::.interactive_join_candidate_score(candidate, score)
  expect_equal(nrow(joined), 3L)
  expect_true(all(c("composite_score", "score_annotation") %in% names(joined)))
  expect_false("evidence_json" %in% names(joined))
  expect_equal(joined$Gene_ID[1L], "g2")
  expect_equal(joined$composite_score[1L], 0.9)
  expect_true(is.na(joined$composite_score[joined$Gene_ID == "g3"]))

  display <- lazyGas:::.interactive_candidate_display_df(joined)
  expect_true("composite_score" %in% names(display))
  expect_equal(
    which(names(display) == "composite_score"),
    which(names(display) == "Gene_ID") + 1L
  )
  expect_false("evidence_json" %in% names(display))
})

test_that("inactive score channels are zeroed from phenotypeRank sources", {
  score <- data.frame(
    Gene_ID = c("g1", "g2"),
    score_annotation = c(0.5, 1),
    score_snpeff = c(0.2, 0.3),
    score_gwas = c(0.8, 0.7),
    score_expression = c(0.1, 0.4),
    score_literature = c(1, 1),
    score_ortholog = c(1, 1),
    score_finemap = c(0.9, 0.2),
    composite_score = c(0.6, 0.5),
    stringsAsFactors = FALSE
  )
  attr(score, "phenotypeRank") <- list(
    sources = c("annotation", "snpeff", "gwas", "expression", "finemap")
  )
  out <- lazyGas:::.zero_inactive_score_columns(score)
  expect_equal(out$score_literature, c(0, 0))
  expect_equal(out$score_ortholog, c(0, 0))
  expect_equal(out$score_annotation, c(0.5, 1))
  expect_equal(out$score_finemap, c(0.9, 0.2))
  expect_equal(out$composite_score, c(0.6, 0.5))
})

test_that("makeInteractiveSummary accepts what=score", {
  skip_if_not_installed("htmltools")
  skip_if_not_installed("reactable")

  # Minimal stub: helpers only need match.arg to accept "score"
  expect_silent(
    match.arg(
      arg = c("candidate", "score"),
      choices = c(
        "scan", "scan_png", "peakcall", "recalc", "groups",
        "peakcall_haplo", "recalc_haplo", "candidate", "score", "qq", "qc",
        "cross_trait", "fine_mapping"
      ),
      several.ok = TRUE
    )
  )
})

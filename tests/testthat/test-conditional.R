test_that("conditional and credible set run on a peak", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "cond_trait")
  lg <- ctx$lg
  on.exit(.close_lg(lg), add = TRUE)

  recalc <- lazyGas::lazyData(object = lg, dataset = "recalc", pheno = ctx$pheno_name)
  peak_id <- recalc$peak_ID[1]

  cond <- lazyGas::conditionalAssoc(
    object = lg,
    pheno = ctx$pheno_name,
    peak_id = peak_id,
    recalc = TRUE,
    k_step = 1L,
    store = TRUE
  )
  expect_true(is.data.frame(cond))
  expect_true("P_conditional" %in% names(cond))
  expect_gt(nrow(cond), 0L)

  cred <- lazyGas::calcCredibleSet(
    object = lg,
    pheno = ctx$pheno_name,
    peak_id = peak_id,
    use_conditional = TRUE,
    recalc = TRUE,
    store = TRUE
  )
  expect_true("PIP" %in% names(cred))
  expect_true(any(cred$in_credible_set))
  expect_true(!is.null(attr(cred, "summary")))

  p_cs <- lazyGas::plotCredibleSet(cred)
  expect_s3_class(p_cs, "ggplot")

  cred_stored <- lazyGas::lazyData(
    lg,
    dataset = "credible_set",
    pheno = ctx$pheno_name,
    kind = paste0("peak_", peak_id)
  )
  expect_equal(nrow(cred_stored), nrow(cred))
})

test_that("annotateOrthologs merges ortholog table", {
  cand <- data.frame(
    peak_ID = 1L,
    Gene_ID = c("g1", "g2"),
    stringsAsFactors = FALSE
  )
  ortho <- data.frame(
    gene_id = "g1",
    ortholog_id = "og1",
    score = 0.9,
    stringsAsFactors = FALSE
  )
  out <- lazyGas::annotateOrthologs(candidate = cand, ortholog = ortho)
  expect_equal(out$ortholog_id[out$Gene_ID == "g1"], "og1")
  expect_true(is.na(out$ortholog_id[out$Gene_ID == "g2"]))

  summ <- lazyGas::summarizeOrthologMatches(out)
  expect_equal(summ$ortholog_id, "og1")
  expect_equal(summ$n_genes, 1L)

  p <- lazyGas::plotOrthologSummary(out)
  expect_s3_class(p, "ggplot")
})

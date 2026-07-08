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

test_that("summarizeFineMapping reports resolution and independent signals", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "fm_trait")
  lg <- ctx$lg
  on.exit(.close_lg(lg), add = TRUE)

  recalc <- lazyGas::lazyData(object = lg, dataset = "recalc", pheno = ctx$pheno_name)
  peak_id <- recalc$peak_ID[1]

  fm <- lazyGas::runFineMapping(
    object = lg,
    pheno = ctx$pheno_name,
    peak_ids = peak_id,
    recalc = TRUE,
    store = TRUE
  )
  expect_true(is.data.frame(fm$summary))
  expect_true(all(c(
    "n_independent_signals", "n_stepwise_leads", "n_ld_artifacts",
    "n_ambiguous_signals", "n_credible_set", "resolution",
    "credible_set_message", "conditional_message",
    "stepwise_lead_variants", "ld_artifact_variants"
  ) %in% names(fm$summary)))
  expect_true(fm$summary$resolution[1] %in% c("high", "moderate", "low"))
  expect_true(fm$summary$n_independent_signals[1] <= fm$summary$n_stepwise_leads[1])

  summ <- lazyGas::summarizeFineMapping(
    object = lg,
    pheno = ctx$pheno_name,
    peak_ids = peak_id,
    run_if_missing = FALSE
  )
  expect_equal(nrow(summ), 1L)

  cond_sum <- lazyGas::summarizeConditionalSignals(
    object = lg,
    pheno = ctx$pheno_name,
    peak_id = peak_id
  )
  expect_equal(cond_sum$n_independent_signals >= 1L, TRUE)
  expect_true("n_stepwise_leads" %in% names(cond_sum))
})

test_that("assessIndependentSignals classifies identical stepwise leads as LD artifacts", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  # Synthetic conditional table: three stepwise leads with identical step-1 stats.
  cond <- data.frame(
    variant_ID = c("A", "B", "C", "B", "C", "C"),
    beta = c(1.5, 1.5, 1.5, 0.01, 0.01, 0.01),
    se = c(0.2, 0.2, 0.2, 0.3, 0.3, 0.3),
    P_conditional = c(1e-8, 1e-8, 1e-8, 1e-6, 1e-6, 1e-6),
    step = c(1L, 1L, 1L, 2L, 2L, 3L),
    conditioned_on = c("", "", "", "A", "A", "A,B"),
    stringsAsFactors = FALSE
  )
  # Make B and C the stepwise leads at steps 2 and 3 (lowest P in their step).
  # Already true above.

  assessed <- lazyGas:::.assess_independent_leads(
    object = NULL,
    cond = cond,
    p_threshold = 0.05,
    r2_threshold = 0.8
  )
  expect_equal(assessed$stepwise, c("A", "B", "C"))
  expect_equal(assessed$leads$verdict, c("supported", "ld_artifact", "ld_artifact"))
  expect_true(all(assessed$pairs$stats_identical))

  summ <- lazyGas:::.conditional_summary_from_assessment(
    peak_id = 1L,
    n_variants_in_locus = 3L,
    assessment = assessed,
    p_threshold = 0.05,
    r2_threshold = 0.8
  )
  expect_equal(summ$n_stepwise_leads, 3L)
  expect_equal(summ$n_independent_signals, 1L)
  expect_equal(summ$n_ld_artifacts, 2L)
  expect_equal(summ$independent_lead_variants, "A")
  expect_equal(summ$ld_artifact_variants, "B,C")
  expect_match(summ$message, "LD artifact")
})

test_that("assessIndependentSignals runs on a LazyGas peak", {
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("rtracklayer")

  ctx <- .run_sample_pipeline(pheno_name = "assess_trait")
  lg <- ctx$lg
  on.exit(.close_lg(lg), add = TRUE)

  recalc <- lazyGas::lazyData(object = lg, dataset = "recalc", pheno = ctx$pheno_name)
  peak_id <- recalc$peak_ID[1]

  out <- lazyGas::assessIndependentSignals(
    object = lg,
    pheno = ctx$pheno_name,
    peak_id = peak_id,
    recalc = TRUE,
    store = FALSE
  )
  expect_true(is.list(out))
  expect_true(all(c("pairs", "leads", "summary") %in% names(out)))
  expect_true(is.data.frame(out$summary))
  expect_true(out$summary$n_independent_signals <= out$summary$n_stepwise_leads)
  if (nrow(out$leads) > 0L) {
    expect_true(all(out$leads$verdict %in% c("supported", "ld_artifact", "ambiguous")))
  }
})

test_that("E1 sentence split and payload keys", {
  sents <- lazyGas:::.e_split_sentences(
    c("Encodes GBSSI. Mutations reduce amylose.", "No delimiter cell")
  )
  expect_true(length(sents) >= 3L)
  expect_true(any(grepl("GBSSI", sents)))

  p <- lazyGas:::.e1_build_payload(
    gene_id = "G1",
    name = "Waxy",
    trait_text = "amylose",
    free_text = "Encodes granule-bound starch synthase. Unrelated sentence here.",
    keyword_details = list(
      matched_keywords = "amylose",
      matched_synonym = list(),
      matched_related = list(list(phrase = "starch synthase", from_user = "amylose")),
      unmatched_keywords = character(),
      matched_context = character()
    )
  )
  expect_identical(p$Gene_ID, "G1")
  expect_true(length(p$matched_sentences) >= 1L)
  expect_true(any(grepl("starch synthase", p$matched_sentences, ignore.case = TRUE)))
  pr <- lazyGas:::.e1_prompts(p)
  expect_identical(pr$return_key, "keywords")
})

test_that("E2/E3 payloads and return keys; E3 has no resolution", {
  snp <- data.frame(
    Gene_ID = "G1",
    Annotation_Impact = c("HIGH", "LOW"),
    Annotation = c("stop_gained", "synonymous"),
    stringsAsFactors = FALSE
  )
  p2 <- lazyGas:::.e2_build_payload("G1", "n", snp)
  expect_equal(p2$HIGH, 1L)
  expect_identical(p2$worst_impact, "HIGH")
  expect_identical(lazyGas:::.e2_prompts(p2)$return_key, "snpeff")

  cand <- data.frame(
    Gene_ID = "G1", dist2peak = 100, max_PIP = 0.2, nearest_cs_PIP = 0.15,
    stringsAsFactors = FALSE
  )
  cred <- data.frame(
    Chr = "chr1", Pos = 1:3, PIP = c(0.2, 0.1, 0.05),
    in_credible_set = TRUE, stringsAsFactors = FALSE
  )
  p3 <- lazyGas:::.e3_build_payload("G1", "n", "gwas", cand, cred, peak_id = "1")
  expect_false("resolution" %in% names(p3))
  expect_false("resolution_for_e9" %in% names(p3))
  expect_true("max_PIP" %in% names(p3))
  expect_identical(lazyGas:::.e3_prompts(p3)$return_key, "gwas")

  p3q <- lazyGas:::.e3_build_payload("G1", "n", "qtl", cand, cred, peak_id = "1")
  expect_true("nearest_cs_PIP" %in% names(p3q))
  expect_null(p3q$max_PIP)
})

test_that("SnpEff keeps worst ANN per site; no multi-transcript double-count", {
  snp <- data.frame(
    Gene_ID = rep("G1", 6L),
    Chr = rep("chr1", 6L),
    Pos = rep(100L, 6L),
    Allele = rep("A", 6L),
    Annotation_Impact = c("HIGH", "HIGH", "HIGH", "LOW", "MODIFIER", "MODERATE"),
    Annotation = c(
      "frameshift_variant", "frameshift_variant", "frameshift_variant",
      "synonymous_variant", "upstream_gene_variant", "missense_variant"
    ),
    Feature_ID = paste0("TX", 1:6),
    `HGVS.p` = c(
      "p.Arg38fs", "p.Arg38fs", "p.Arg38fs",
      "p.Arg38Arg", "", "p.Arg38His"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  red <- lazyGas:::.snpeff_worst_ann_per_site(snp)
  expect_equal(nrow(red), 1L)
  expect_identical(as.character(red$Annotation_Impact), "HIGH")

  ic <- lazyGas:::.e2_impact_counts(snp)
  expect_equal(unname(ic$counts[["HIGH"]]), 1L)
  expect_equal(unname(ic$counts[["MODERATE"]]), 0L)
  expect_equal(unname(ic$counts[["LOW"]]), 0L)
  expect_equal(unname(ic$counts[["MODIFIER"]]), 0L)
  expect_equal(length(ic$high_effects), 1L)
  expect_match(ic$high_effects[[1L]]$effect, "frameshift")

  snp2 <- rbind(
    snp,
    data.frame(
      Gene_ID = "G1",
      Chr = "chr1",
      Pos = 200L,
      Allele = "T",
      Annotation_Impact = "LOW",
      Annotation = "synonymous_variant",
      Feature_ID = "TX9",
      `HGVS.p` = "p.Leu10Leu",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  ic2 <- lazyGas:::.e2_impact_counts(snp2)
  expect_equal(unname(ic2$counts[["HIGH"]]), 1L)
  expect_equal(unname(ic2$counts[["LOW"]]), 1L)

  collapsed <- lazyGas:::.collapseSnpEff(transform(snp2, negLog10P = 5))
  expect_equal(as.integer(collapsed$HIGH[1L]), 1L)
  expect_equal(as.integer(collapsed$LOW[1L]), 1L)
  expect_equal(as.integer(collapsed$MODERATE[1L]), 0L)
})

test_that("E5–E7 empty skip; E8/E9/E10 return keys", {
  expect_null(lazyGas:::.e5_build_payload("G1", ann_row = list(), go_cols = "GO"))
  p5 <- lazyGas:::.e5_build_payload(
    "G1", ann_row = list(GO = "GO:0001; starch binding"), go_cols = "GO"
  )
  expect_true(length(p5$terms) >= 1L)
  expect_identical(lazyGas:::.e5_prompts(p5)$return_key, "go")
  expect_identical(lazyGas:::.e6_prompts(list())$return_key, "kegg")
  expect_identical(lazyGas:::.e7_prompts(list())$return_key, "domains")

  p8 <- lazyGas:::.e8_build_payload(
    "G1", "n",
    prose = list(keywords = "k"),
    coded = list(HIGH = 0L)
  )
  expect_identical(lazyGas:::.e8_prompts(p8)$return_key, "summary")

  p10 <- lazyGas:::.e10_build_payload(
    gene_id = "G1",
    name = "n",
    trait_text = "amylose",
    prose = list(keywords = "starch"),
    summary_text = "brief",
    coded = list(HIGH = 0L, dist2peak_bp = 10)
  )
  expect_identical(p10$summary, "brief")
  expect_identical(lazyGas:::.e10_prompts(p10)$return_key, "interpretation")
  expect_match(lazyGas:::.e10_prompts(p10)$system, "Plausibility")

  p9 <- lazyGas:::.e9_build_payload(
    "G1", "n",
    coded = list(HIGH = 0L, max_PIP = 0.01, dist2peak_bp = 10),
    resolution = "low",
    n_cs_variants = 20,
    mapping_mode = "gwas"
  )
  expect_identical(p9$peak$resolution, "low")
  expect_identical(p9$peak$resolution_label, "not resolved")
  expect_identical(lazyGas:::.e9_prompts(p9)$return_key, "validity")
})

test_that("Yanai tau and B5 number allowlist", {
  expect_equal(lazyGas:::.yanai_tau(c(10, 0, 0)), 1)
  expect_true(is.na(lazyGas:::.yanai_tau(5)))

  allow <- c(0.2, 10, 3)
  found <- lazyGas:::.e_extract_numbers("PIP 0.2 at 10 bp; n=3")
  expect_true(all(lazyGas:::.e_numbers_allowed(found, allow)))
  found_bad <- lazyGas:::.e_extract_numbers("score 0.99 invented")
  expect_false(all(lazyGas:::.e_numbers_allowed(found_bad, allow)))

  # Label digits from payload remain allowlisted; invented quantities are not.
  allow_lab <- lazyGas:::.e_collect_allowed_numbers(
    list(domains = list(list(name = "Glycosyl transferase family 1", source = "pfam")))
  )
  expect_true(1 %in% allow_lab)
  expect_false(all(lazyGas:::.e_numbers_allowed(c(2.4), allow_lab)))

  rule <- lazyGas:::.e_no_digits_prose_rule()
  expect_match(rule, "quantitative", ignore.case = TRUE)
  expect_match(lazyGas:::.e7_prompts(list(Gene_ID = "G1", domains = list()))$system, "quantitative", ignore.case = TRUE)
  expect_match(lazyGas:::.e5_prompts(list(Gene_ID = "G1", terms = list()))$system, "quantitative", ignore.case = TRUE)
  expect_match(lazyGas:::.e6_prompts(list(Gene_ID = "G1", entries = list()))$system, "quantitative", ignore.case = TRUE)
})

test_that("E4 names high-expression tissues when specific", {
  mat <- matrix(
    c(100, 90, 5, 4, 3, 2),
    nrow = 1L,
    dimnames = list("G1", c("leaf_a", "leaf_b", "root_a", "root_b", "seed_a", "seed_b"))
  )
  smap <- data.frame(
    sample = colnames(mat),
    tissue = c("leaf", "leaf", "root", "root", "seed", "seed"),
    stringsAsFactors = FALSE
  )
  p4 <- lazyGas:::.e4_build_payload("G1", expr_mat = mat, sample_map = smap)
  expect_true(p4$specificity$tissue_summary %in% c("high", "moderate"))
  expect_true(p4$specificity$summary %in% c("high", "moderate"))
  expect_true("leaf" %in% p4$high_tissues)
  expect_identical(p4$tissues_by_level[1L], "leaf")
  expect_identical(p4$matrix_scale, "abundance")
  expect_true(is.finite(p4$relative_level$cohort_size))
  pr <- lazyGas:::.e4_prompts(p4)
  expect_match(pr$system, "high_expression_tissues")
  tmpl <- lazyGas:::.e_template_text("expression", "G1", p4)
  expect_match(tmpl, "leaf")
})

test_that("E4 combo tau drops query-axis not_available; Z nulls tau", {
  mat <- matrix(
    c(100, 90, 5, 4, 50, 40),
    nrow = 1L,
    dimnames = list("G1", paste0("s", 1:6))
  )
  smap <- data.frame(
    sample = colnames(mat),
    tissue = c("leaf", "leaf", "root", "root", "leaf", "leaf"),
    stage = c("veg", "veg", "veg", "veg", "not_available", "not_available"),
    condition = c("control", "control", "control", "control", "drought", "drought"),
    stringsAsFactors = FALSE
  )
  p <- lazyGas:::.e4_build_payload(
    "G1", expr_mat = mat, sample_map = smap,
    cfg = list(query_axes = c("tissue", "stage", "condition"))
  )
  # columns 5–6 dropped from combo because stage is not_available on query axis
  expect_true(is.finite(p$specificity$combination_tau))
  expect_equal(p$n_states_combination, 2L)

  pz <- lazyGas:::.e4_build_payload(
    "G1", expr_mat = mat, sample_map = smap,
    cfg = list(matrix_scale = "row_zscore")
  )
  expect_null(pz$specificity$combination_tau)
  expect_null(pz$specificity$summary)
  expect_null(pz$relative_level$band)
  expect_identical(pz$relative_level$reason, "row_zscore_height_not_applicable")
})

test_that("E4 relative height can restrict to same biotype", {
  mat <- matrix(
    c(
      100, 100,
      10, 10,
      50, 50
    ),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("G1", "G2", "G3"), c("leaf", "root"))
  )
  bt <- c(G1 = "protein_coding", G2 = "lncRNA", G3 = "protein_coding")
  p <- lazyGas:::.e4_build_payload(
    "G1", expr_mat = mat, sample_map = NULL,
    cfg = list(gene_biotype = bt)
  )
  expect_equal(p$relative_level$cohort_size, 2L)
  expect_match(p$relative_level$basis, "protein_coding")
})

test_that("E3/E8/E9 prompts carry QTL ABF and peak-relative cues", {
  p3 <- list(mode = "qtl", peak = list(n_cs_variants = 3L), Gene_ID = "G1",
             dist2peak_bp = 10, negLog10P = 5, nearest_cs_PIP = 0.2)
  expect_match(lazyGas:::.e3_prompts(p3)$system, "ABF")
  expect_match(
    lazyGas:::.e8_prompts(list(Gene_ID = "G1", prose = list(), coded = list()))$system,
    "peak-relative"
  )
  expect_match(
    lazyGas:::.e_template_text("summary", "G1", list(coded = list(dist2peak_bp = 1))),
    "peak-relative"
  )
  p9 <- list(
    peak = list(n_cs_variants = 2L, resolution = "low",
                resolution_label = "not resolved"),
    Gene_ID = "G1", mapping_mode = "qtl",
    dist2peak_bp = 1, negLog10P = 2, max_PIP = NULL, nearest_cs_PIP = 0.1,
    worst_impact = "LOW", HIGH = 0L, MODERATE = 0L, LOW = 1L, MODIFIER = 0L,
    domain_overlap_n = 0L
  )
  expect_match(lazyGas:::.e9_prompts(p9)$system, "ABF")
})

test_that("E gene loop use_llm=FALSE templates and HTML order", {
  skip_if_not_installed("yaml")
  wd <- tempfile("lazygas_e_")
  dir.create(wd)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  writeLazyGasExploreConfig(wd, mapping_mode = "gwas", top_n = Inf)

  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1000, 1100), end = c(1000, 2000)),
    strand = c("+", "+"),
    type = c("gene", "CDS"),
    ID = c("G1", "cds1"),
    Parent = c("", "G1")
  )
  cand <- data.frame(
    Gene_ID = "G1", Name = "geneA", Chr = "chr1", dist2peak = 5,
    stringsAsFactors = FALSE
  )
  cred <- data.frame(
    Chr = "chr1", Pos = 1500, PIP = 0.5, in_credible_set = TRUE,
    stringsAsFactors = FALSE
  )
  ann <- data.frame(
    Gene_ID = "G1",
    Description = "Involved in flowering time regulation.",
    GO = "flowering involvement",
    stringsAsFactors = FALSE
  )
  map <- classifyAnnColumns(
    ann[, c("Description", "GO"), drop = FALSE],
    use_llm = FALSE,
    overrides = c(GO = "go")
  )
  lazyGas:::.write_ann_column_map(wd, map)

  e_run <- lazyGas:::.report_e_run_peak(
    work_dir = wd,
    mapping_mode = "gwas",
    simple_candidates = cand,
    credible_set = cred,
    gff = gr,
    ann = ann,
    snpeff = NULL,
    use_llm = FALSE,
    llm = NULL,
    resolution = "low",
    top_n = Inf,
    cfg = list(id_map_download = FALSE)
  )
  expect_equal(nrow(e_run$genes), 1L)
  gr1 <- e_run$gene_results[["G1"]]
  expect_true(nzchar(gr1$results$keywords$text))
  expect_true(nzchar(gr1$results$summary$text))
  expect_true(nzchar(gr1$results$interpretation$text))
  expect_true(nzchar(gr1$results$validity$text))
  expect_null(gr1$results$snpeff$text) # use_llm=FALSE → no E2/E3 prose
  expect_true(!isTRUE(gr1$results$go$skipped))

  html <- lazyGas:::.report_evidence_html_e(e_run)
  expect_match(html, "GWAS / position")
  expect_match(html, "Variant effect \\(SnpEff\\)")
  expect_match(html, "Keywords / annotation")
  expect_match(html, "Brief summary")
  expect_match(html, "Integrated interpretation")
  expect_match(html, "Brief validity")
  # GO before summary before interpretation
  expect_true(regexpr("GO", html) < regexpr("Brief summary", html))
  expect_true(
    regexpr("Brief summary", html) < regexpr("Integrated interpretation", html)
  )
})

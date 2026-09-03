test_that("peakcall region bounds use min/max Pos per peak_ID", {
  peakcall <- data.frame(
    peak_ID = c(1L, 1L, 1L, 2L, 2L),
    Pos = c(36131042, 36487522, 39056110, 100, 200),
    stringsAsFactors = FALSE
  )
  bounds <- lazyGas:::.peakcall_region_bounds(peakcall)
  expect_equal(nrow(bounds), 2L)
  b1 <- bounds[as.character(bounds$peak_ID) == "1", , drop = FALSE]
  expect_equal(b1$region_start, 36131042)
  expect_equal(b1$region_end, 39056110)
})

test_that("report i18n is English-only (B4)", {
  en <- lazyGas:::.report_i18n("en")
  ja <- lazyGas:::.report_i18n("ja")
  expect_equal(en$candidate_list, "Candidate list")
  expect_equal(ja$candidate_list, "Candidate list")
  expect_equal(ja$basic_info, "Basic information")
  expect_equal(ja$resolvable, "resolvable")
  expect_equal(ja$score_note, en$score_note)
  expect_equal(lazyGas:::.report_fmt_score(c(0.1, NA, 1)), c("0.100", "0.000", "1.000"))
  expect_equal(lazyGas:::.report_name_cell(NA_character_), "")
  expect_equal(lazyGas:::.report_name_cell("  "), "")
  expect_equal(lazyGas:::.report_name_cell("Waxy"), "Waxy")
})

test_that("B4 translation prompt keeps IDs and technical terms", {
  msg <- lazyGas:::.report_translate_prose_ja_system()
  expect_match(msg, "Translate the JSON string values from English to Japanese", fixed = TRUE)
  expect_match(msg, "gene IDs", fixed = TRUE)
  expect_match(msg, "Keep scientific / technical terms in English", fixed = TRUE)
  expect_match(msg, "matched keyword tokens exactly", fixed = TRUE)
  expect_match(msg, "SnpEff impact labels", fixed = TRUE)
})

test_that("ja + use_llm=FALSE keeps English chrome and template prose", {
  ev <- list(
    gwas = list(details = list(dist2peak = 100, negLog10P = 5)),
    snpeff = list(details = list(worst_impact = "MODERATE", MODERATE = 2)),
    annotation = list(details = list(matched_keywords = list("fruit")))
  )
  ranked <- data.frame(
    Gene_ID = "Os01g00001",
    Name = "Waxy",
    evidence_json = jsonlite::toJSON(ev, auto_unbox = TRUE),
    stringsAsFactors = FALSE
  )
  html <- lazyGas:::.report_evidence_html(
    rank_sub = ranked,
    query = list(text = "fruit weight", query_id = "q"),
    language = "ja",
    use_llm = FALSE,
    llm = NULL
  )
  expect_match(html, "GWAS / position", fixed = TRUE)
  expect_match(html, "Keywords / annotation", fixed = TRUE)
  expect_match(html, "Matched keywords: fruit", fixed = TRUE)
  expect_false(grepl("候補一覧|基礎情報|キーワード／|enable use_llm", html))
})

test_that("translate rejects incomplete JSON and keeps English", {
  eng <- list(
    g1 = list(
      keywords = "Matched keywords: fruit.",
      expression = NULL,
      validity = "Check distance."
    )
  )
  # Force llmChat failure via bogus endpoint
  out <- suppressWarnings(
    lazyGas:::.report_translate_prose_ja(
      eng,
      llm = list(
        model = "none",
        base_url = "http://127.0.0.1:9",
        timeout = 1
      )
    )
  )
  expect_identical(out, eng)
})

test_that("candidate table HTML scrolls and uses 3 decimals", {
  ranked <- data.frame(
    Gene_ID = c("g1", "g2"),
    Name = c("Waxy", NA_character_),
    score_annotation = c(0.92, 0.11),
    score_snpeff = c(0.4, 0.85),
    score_gwas = c(0.71, 0.55),
    score_expression = c(0.1, 0),
    score_finemap = c(0.8, 0.2),
    composite_score = c(0.61, 0.48),
    stringsAsFactors = FALSE
  )
  html <- lazyGas:::.report_candidate_table_html(
    ranked,
    sources = c("annotation", "snpeff", "gwas", "expression", "finemap"),
    language = "en"
  )
  expect_match(html, "candidate-table-scroll", fixed = TRUE)
  expect_match(html, "overflow-x: auto", fixed = TRUE)
  expect_match(html, "0.920", fixed = TRUE)
  expect_match(html, ">Waxy<", fixed = TRUE)
  expect_match(html, "Scores are a weighted sum", fixed = TRUE)
  # empty Name cell still present
  expect_match(html, "<td></td>", fixed = TRUE)
})

test_that("rank result csv/rds round-trip preserves rows and attributes", {
  ranked <- data.frame(
    peak_ID = c(1L, 1L, 2L),
    Gene_ID = c("g1", "g2", "g3"),
    composite_score = c(0.9, 0.4, 0.1),
    evidence_json = c("{\"a\":1}", "{\"a\":2}", "{\"a\":3}"),
    stringsAsFactors = FALSE
  )
  attr(ranked, "phenotypeRank") <- list(query_id = "q1", n_genes = 3L)

  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  csv_path <- file.path(tmp, "rank.csv")
  rds_path <- file.path(tmp, "rank.rds")
  lazyGas:::.write_rank_result_files(ranked, csv_path, rds_path)

  from_rds <- lazyGas:::.read_rank_result(rds_path)
  expect_equal(from_rds$Gene_ID, ranked$Gene_ID)
  expect_equal(attr(from_rds, "phenotypeRank")$query_id, "q1")

  from_csv <- lazyGas:::.read_rank_result(csv_path)
  expect_equal(from_csv$Gene_ID, ranked$Gene_ID)
  expect_equal(from_csv$composite_score, ranked$composite_score)
})

test_that("bind_rank_results concatenates peaks and sets n_genes", {
  a <- data.frame(peak_ID = 1L, Gene_ID = "g1", composite_score = 0.8,
                  stringsAsFactors = FALSE)
  b <- data.frame(peak_ID = 2L, Gene_ID = "g2", composite_score = 0.3,
                  stringsAsFactors = FALSE)
  attr(a, "phenotypeRank") <- list(query_id = "q1", pheno = "x", n_genes = 1L)
  out <- lazyGas:::.bind_rank_results(list(a, b))
  expect_equal(nrow(out), 2L)
  expect_equal(out$Gene_ID, c("g1", "g2"))
  expect_equal(attr(out, "phenotypeRank")$n_genes, 2L)
  expect_equal(attr(out, "phenotypeRank")$query_id, "q1")
})

test_that("gwas evidence no longer adds PIP*5", {
  rows <- data.frame(
    Gene_ID = "g1",
    negLog10P = 5,
    dist2peak = 0,
    PIP = 1,
    stringsAsFactors = FALSE
  )
  ev <- lazyGas:::.evidence_gwas("g1", rows)
  # score = 5 + proximity(1) = 6; PIP must not add +5
  expect_equal(ev$score, 6)
  expect_null(ev$details$PIP)
})

test_that("finemap evidence uses max PIP and sources always include finemap", {
  lookup <- data.frame(
    peak_ID = c("1", "1", "2"),
    Gene_ID = c("g1", "g2", "g1"),
    max_PIP = c(0.82, 0.11, 0.99),
    stringsAsFactors = FALSE
  )
  ev1 <- lazyGas:::.evidence_finemap(
    "g1",
    data.frame(peak_ID = "1", Gene_ID = "g1", stringsAsFactors = FALSE),
    lookup
  )
  expect_equal(ev1$score, 0.82)
  expect_equal(ev1$details$max_PIP, 0.82)

  ev_miss <- lazyGas:::.evidence_finemap(
    "g9",
    data.frame(peak_ID = "1", Gene_ID = "g9", stringsAsFactors = FALSE),
    lookup
  )
  expect_equal(ev_miss$score, 0)

  sw <- lazyGas:::.report_sources_with_finemap(
    sources = c("annotation", "gwas"),
    weights = c(annotation = 0.3, gwas = 0.25)
  )
  expect_true("finemap" %in% sw$sources)
  expect_equal(unname(sw$weights[["finemap"]]), 0.5)
})

test_that("evidence HTML writes coded metrics and omits empty expression", {
  ev <- list(
    gwas = list(details = list(dist2peak = 1200, negLog10P = 8.5)),
    snpeff = list(details = list(worst_impact = "HIGH", HIGH = 1, MODERATE = 0)),
    finemap = list(details = list(max_PIP = 0.77)),
    annotation = list(details = list(
      matched_keywords = list("fruit"),
      unmatched_keywords = list("weight")
    ))
  )
  ranked <- data.frame(
    Gene_ID = "Os01g00001",
    Name = "Waxy",
    dist2peak = 1200,
    negLog10P = 8.5,
    composite_score = 0.5,
    evidence_json = jsonlite::toJSON(ev, auto_unbox = TRUE),
    stringsAsFactors = FALSE
  )
  html <- lazyGas:::.report_evidence_html(
    rank_sub = ranked,
    query = list(text = "fruit weight", query_id = "q"),
    language = "en",
    use_llm = FALSE,
    llm = NULL
  )
  expect_match(html, "GWAS / position", fixed = TRUE)
  expect_match(html, "Variant effect (SnpEff)", fixed = TRUE)
  expect_match(html, "max_PIP (95% CS ∩ SnpEff) = 0.770", fixed = TRUE)
  expect_match(html, "dist2peak = 1,200 bp", fixed = TRUE)
  expect_match(html, "HIGH = 1", fixed = TRUE)
  expect_match(html, "Matched keywords: fruit", fixed = TRUE)
  expect_false(grepl("Unmatched", html, fixed = TRUE))
  expect_false(grepl("\\bweight\\b", html))
  expect_false(grepl("Expression", html, fixed = TRUE))
  expect_false(grepl("composite_score", html, fixed = TRUE))
  expect_false(grepl("score_annotation", html, fixed = TRUE))
})

test_that("report template prose lists matched keywords only", {
  expect_equal(
    lazyGas:::.report_gene_template_prose(list(
      ann_matched = c("fruit", "starch"),
      has_expression = FALSE,
      expr_snippets = character()
    ))$keywords,
    "Matched keywords: fruit, starch."
  )
  expect_equal(
    lazyGas:::.report_gene_template_prose(list(
      ann_matched = character(),
      has_expression = FALSE,
      expr_snippets = character()
    ))$keywords,
    "No matched phenotype keywords."
  )
  payload <- lazyGas:::.report_llm_evidence_payload(
    data.frame(
      Gene_ID = "g1",
      Name = "",
      evidence_json = jsonlite::toJSON(
        list(annotation = list(details = list(
          matched_keywords = list("fruit"),
          unmatched_keywords = list("weight")
        ))),
        auto_unbox = TRUE
      ),
      stringsAsFactors = FALSE
    )
  )
  expect_null(payload[[1]]$evidence$annotation$details$unmatched_keywords)
  expect_equal(
    unlist(payload[[1]]$evidence$annotation$details$matched_keywords),
    "fruit"
  )
})

test_that("llm_report writes HTML with table and reuses rank rds", {
  skip_if_not_installed("lazyGas")
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")

  res <- tryCatch(
    .run_sample_pipeline(pheno_name = "llm_rank_reuse"),
    error = function(e) {
      testthat::skip(paste("sample pipeline unavailable:", conditionMessage(e)))
    }
  )
  on.exit(.close_lg(res$lg), add = TRUE)

  cand <- lazyData(res$lg, dataset = "candidate", pheno = res$pheno_name)
  skip_if(is.null(cand) || nrow(cand) == 0L, "no candidate rows")

  peak_ids <- unique(as.character(cand$peak_ID))
  has_cs <- vapply(peak_ids, function(pid) {
    !is.null(tryCatch(
      lazyData(
        res$lg,
        dataset = "credible_set",
        pheno = res$pheno_name,
        kind = paste0("peak_", pid)
      ),
      error = function(e) NULL
    ))
  }, logical(1))
  skip_if(!any(has_cs), "no credible set available for sample pipeline")

  ranked <- cand
  ranked$composite_score <- seq(from = 1, to = 0, length.out = nrow(ranked))
  ranked$score_annotation <- ranked$composite_score
  ranked$score_snpeff <- 0
  ranked$score_gwas <- 0
  ranked$score_expression <- 0
  ranked$score_finemap <- 0
  ranked$evidence_json <- "{}"
  attr(ranked, "phenotypeRank") <- list(
    query_id = "reuse_q",
    pheno = res$pheno_name,
    n_genes = nrow(ranked)
  )

  out_dir <- tempfile("ai_report_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  q <- phenotypeQuery("fruit weight", tissues = "fruit", use_llm = FALSE)
  rep1 <- llm_report(
    object = res$lg,
    pheno = res$pheno_name,
    query = q,
    language = "en",
    top_n = 2L,
    use_llm = FALSE,
    query_use_llm = FALSE,
    out_dir = out_dir,
    rank_result = ranked,
    peak_id = peak_ids[has_cs][1],
    download_tenor = FALSE
  )
  expect_true(grepl("\\.html$", rep1$path))
  expect_true(file.exists(rep1$path))
  expect_match(rep1$html, "<!DOCTYPE html>", fixed = TRUE)
  expect_match(rep1$html, "Candidate list", fixed = TRUE)
  expect_match(rep1$html, "candidate-table-scroll", fixed = TRUE)
  expect_match(rep1$html, "Credible set", fixed = TRUE)
  expect_false(grepl("Evidence verification", rep1$html, fixed = TRUE))
  expect_null(rep1$verification)
  expect_true(file.exists(rep1$rank_csv))
  expect_true(file.exists(rep1$rank_rds))

  rep2 <- llm_report(
    object = res$lg,
    pheno = res$pheno_name,
    query = q,
    language = "en",
    top_n = 2L,
    use_llm = FALSE,
    query_use_llm = FALSE,
    out_dir = out_dir,
    rank_result = rep1$rank_rds,
    peak_id = peak_ids[has_cs][1],
    download_tenor = FALSE
  )
  expect_equal(nrow(rep2$rank_result), nrow(ranked))
  expect_equal(attr(rep2$rank_result, "phenotypeRank")$query_id, "reuse_q")
})

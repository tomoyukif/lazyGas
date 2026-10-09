
.test_ann_column_map <- function(candidate = NULL) {
  cols <- c("Description", "GO", "Name")
  if (!is.null(candidate)) {
    cols <- intersect(cols, names(candidate))
  }
  if (!length(cols)) {
    cols <- "Description"
  }
  data.frame(
    column = cols,
    class = ifelse(cols == "GO", "go", "free_description"),
    use_for_ranking = TRUE,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

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

test_that("phenotypeQuery builds structured query without LLM", {
  skip_if_not_installed("lazyGas")
  q <- phenotypeQuery(
    "fruit weight at maturity",
    tissues = "fruit",
    stage = "maturation",
    use_llm = FALSE
  )
  expect_s3_class(q, "PhenotypeQuery")
  expect_true(nzchar(q$query_id))
  expect_true("fruit" %in% q$tissues)
  expect_equal(q$parsed_by, "fallback")
  expect_null(q$species)
  expect_equal(q$synonym_phrases, list())
  expect_equal(q$related_phrases, list())
  expect_output(print(q), "PhenotypeQuery:")
  expect_output(print(q), "parsed_by: fallback")
})

test_that("phenotypeQuery round-trips through store without print error", {
  skip_if_not_installed("lazyGas")
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")

  res <- .run_sample_pipeline(pheno_name = "explorer_pq")
  on.exit(.close_lg(res$lg), add = TRUE)

  q <- phenotypeQuery(
    "fruit weight at maturity",
    tissues = "fruit",
    stage = "maturation",
    use_llm = FALSE,
    object = res$lg,
    save = TRUE
  )
  expect_null(q$species)

  ids <- lazyData(res$lg, dataset = "phenotype_query")
  expect_true(q$query_id %in% ids)

  loaded <- lazyData(res$lg, dataset = "phenotype_query", kind = q$query_id)
  expect_s3_class(loaded, "PhenotypeQuery")
  expect_null(loaded$species)
  expect_equal(loaded$tissues, "fruit")
  expect_equal(loaded$developmental_stage, "maturation")
  expect_output(print(loaded), "PhenotypeQuery:")
  expect_false(any(grepl("species:", capture.output(print(loaded)))))
})

test_that("rankPhenotypeCandidates ranks demo candidates", {
  skip_if_not_installed("lazyGas")
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  res <- .run_sample_pipeline(pheno_name = "explorer_trait")
  on.exit(.close_lg(res$lg), add = TRUE)

  q <- phenotypeQuery(
    "fruit weight development",
    tissues = "fruit",
    use_llm = FALSE,
    object = res$lg,
    save = FALSE
  )

  extdata <- system.file("extdata", package = "lazyGas")
  expr_df <- read.csv(
    file.path(extdata, "demo_expression.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(expr_df) <- expr_df$Gene_ID
  expr_mat <- as.matrix(expr_df[, setdiff(names(expr_df), "Gene_ID"), drop = FALSE])
  expr_meta <- read.csv(
    file.path(extdata, "demo_expression_meta.csv"),
    stringsAsFactors = FALSE
  )
  rownames(expr_meta) <- expr_meta$sample
  ortho_graph <- file.path(extdata, "demo_ortholog_graph.csv")

  cand <- lazyData(res$lg, dataset = "candidate", pheno = res$pheno_name)
  ranked <- rankPhenotypeCandidates(
    object = res$lg,
    pheno = res$pheno_name,
    query = q,
    sources = c("annotation", "gwas_dist", "gwas_peak", "expression", "ortholog"),
    ann_column_map = .test_ann_column_map(cand),
    expression_matrix = expr_mat,
    expression_sample_map = expr_meta,
    ortholog_graph = ortho_graph,
    top_n = 10L,
    save = TRUE
  )

  expect_true(nrow(ranked) > 0L)
  expect_true("composite_score" %in% names(ranked))
  expect_true("rank_within_peak" %in% names(ranked))
  expect_true(all(ranked$composite_score >= 0 & ranked$composite_score <= 1))
  if ("peak_ID" %in% names(ranked)) {
    for (pid in unique(ranked$peak_ID)) {
      sub <- ranked[ranked$peak_ID == pid | (is.na(pid) & is.na(ranked$peak_ID)), , drop = FALSE]
      expect_true(sub$composite_score[1L] >= sub$composite_score[nrow(sub)])
      expect_equal(sub$rank_within_peak, seq_len(nrow(sub)))
    }
  } else {
    expect_true(ranked$composite_score[1L] >= ranked$composite_score[nrow(ranked)])
  }

  stored <- lazyData(res$lg, dataset = "phenotype_rank", pheno = res$pheno_name)
  expect_false(is.null(stored))
  expect_equal(nrow(stored), nrow(ranked))

  ranked_all <- rankPhenotypeCandidates(
    object = res$lg,
    pheno = res$pheno_name,
    query = q,
    sources = c("annotation", "gwas_dist", "gwas_peak", "expression", "ortholog"),
    ann_column_map = .test_ann_column_map(cand),
    expression_matrix = expr_mat,
    expression_sample_map = expr_meta,
    ortholog_graph = ortho_graph,
    top_n = NULL,
    save = FALSE
  )
  expect_true(nrow(ranked_all) >= nrow(ranked))
  # No global Gene_ID collapse: multi-peak genes may appear more than once
  expect_true(nrow(ranked_all) >= length(unique(ranked_all$Gene_ID)))
})

test_that("collectGeneEvidence returns annotation scores", {
  skip_if_not_installed("lazyGas")
  cand <- make_test_candidate()
  q <- phenotypeQuery(
    "fruit weight",
    tissues = "fruit",
    use_llm = FALSE,
    noun_phrases = "fruit weight",
    words = c("fruit", "weight")
  )
  ev <- collectGeneEvidence(
    gene_ids = cand$Gene_ID,
    query = q,
    candidate = cand,
    sources = c("annotation", "gwas_dist", "gwas_peak"),
    ann_column_map = .test_ann_column_map(cand),
    use_cache = FALSE
  )
  expect_length(ev, nrow(cand))
  expect_equal(ev[["g1"]]$annotation$score, 0.5)
  expect_true(ev[["g1"]]$annotation$score >= ev[["g2"]]$annotation$score)
})

test_that("annotation tier prefers action over noun over word", {
  skip_if_not_installed("lazyGas")
  cand <- data.frame(
    Gene_ID = c("a", "n", "w", "z"),
    Description = c(
      "protein that regulates flowering time in rice",
      "involved in flowering time control",
      "flowering related gene",
      "photosystem subunit"
    ),
    stringsAsFactors = FALSE
  )
  q <- phenotypeQuery(
    "flowering",
    action_phrases = "regulates flowering",
    noun_phrases = "flowering time",
    words = "flowering"
  )
  ev <- collectGeneEvidence(
    gene_ids = cand$Gene_ID,
    query = q,
    candidate = cand,
    sources = "annotation",
    ann_column_map = .test_ann_column_map(cand),
    use_cache = FALSE
  )
  expect_equal(ev[["a"]]$annotation$score, 1)
  expect_equal(ev[["n"]]$annotation$score, 0.5)
  expect_equal(ev[["w"]]$annotation$score, 0.25)
  expect_equal(ev[["z"]]$annotation$score, 0)
})

test_that("network literature/ortholog scores are 1 or 0", {
  skip_if_not_installed("lazyGas")
  cand <- data.frame(
    Gene_ID = c("g1", "g2"),
    Name = c("GENE1", "GENE2"),
    Description = c("x", "y"),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("g1", "fruit"),
    to = c("fruit", "weight"),
    stringsAsFactors = FALSE
  )
  q <- phenotypeQuery("fruit weight", tissues = "fruit", words = "fruit")
  ev <- collectGeneEvidence(
    gene_ids = cand$Gene_ID,
    query = q,
    candidate = cand,
    sources = c("literature", "ortholog"),
    literature_graph = edges,
    ortholog_graph = edges,
    use_cache = FALSE
  )
  expect_equal(ev[["g1"]]$literature$score, 1)
  expect_equal(ev[["g2"]]$literature$score, 0)
  expect_equal(ev[["g1"]]$ortholog$score, 1)
  expect_equal(ev[["g2"]]$ortholog$score, 0)
})


test_that("phenotype query list via lazyData", {
  skip_if_not_installed("lazyGas")
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  res <- .run_sample_pipeline(pheno_name = "pq_store_trait")
  on.exit(.close_lg(res$lg), add = TRUE)

  q <- phenotypeQuery("root development", tissues = "root", object = res$lg, save = TRUE)
  ids <- lazyData(res$lg, dataset = "phenotype_query")
  expect_true(q$query_id %in% ids)
  q2 <- lazyData(res$lg, dataset = "phenotype_query", kind = q$query_id)
  expect_equal(q2$trait_text, q$trait_text)
})

test_that("restorePhenoFromStore recovers phenotypes from companion store", {
  skip_if_not_installed("lazyGas")
  skip_if_not_installed("GBScleanR")
  skip_if_not_installed("arrow")

  res <- .run_sample_pipeline(pheno_name = "restore_trait")
  on.exit(.close_lg(res$lg), add = TRUE)

  expect_equal(getPheno(res$lg)$pheno_names, "restore_trait")

  res$lg@lazydata$pheno <- NULL
  res$lg@lazydata$pheno_names <- NULL
  res$lg@lazydata$pheno_type <- NULL

  lg2 <- lazyGas::restorePhenoFromStore(res$lg, warn_stub = FALSE)
  expect_equal(getPheno(lg2)$pheno_names, "restore_trait")
})

test_that("resolveLazyGasPaths finds companion store", {
  skip_if_not_installed("lazyGas")
  demo_gds <- system.file("extdata", "sample.gds", package = "lazyGas")
  if (!nzchar(demo_gds)) {
    demo_gds <- system.file("extdata", "sample.gds", package = "GBScleanR")
  }
  skip_if_not(nzchar(demo_gds) && file.exists(demo_gds), "sample.gds missing")

  paths <- resolveLazyGasPaths(demo_gds)
  expect_true(file.exists(paths$gds))
})





test_that("runner Shiny app is a thin pipeline launcher", {
  skip_if_not_installed("lazyGas")
  candidates <- c(
    system.file("shiny", "runner", "app.R", package = "lazyGas"),
    normalizePath(
      file.path(testthat::test_path(), "..", "..", "inst", "shiny", "runner", "app.R"),
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      file.path(getwd(), "inst", "shiny", "runner", "app.R"),
      winslash = "/",
      mustWork = FALSE
    )
  )
  app_path <- candidates[nzchar(candidates) & file.exists(candidates)][1L]
  skip_if(is.na(app_path) || !nzchar(app_path), "runner app.R not found")
  src <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("runLazyGas", src, fixed = TRUE))
  expect_true(grepl("assignPheno", src, fixed = TRUE))
  expect_true(grepl("runLazyGasExplorer", src, fixed = TRUE))
  expect_false(grepl("rankPhenotypeCandidates", src, fixed = TRUE))
  expect_error(parse(file = app_path), NA)
})

test_that("explorer Shiny app embeds dashboard plot tabs", {
  skip_if_not_installed("lazyGas")
  candidates <- c(
    system.file("shiny", "explorer", "app.R", package = "lazyGas"),
    normalizePath(
      file.path(testthat::test_path(), "..", "..", "inst", "shiny", "explorer", "app.R"),
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      file.path(getwd(), "inst", "shiny", "explorer", "app.R"),
      winslash = "/",
      mustWork = FALSE
    ),
    normalizePath(
      file.path(getwd(), "..", "..", "inst", "shiny", "explorer", "app.R"),
      winslash = "/",
      mustWork = FALSE
    )
  )
  app_path <- candidates[nzchar(candidates) & file.exists(candidates)][1L]
  skip_if(is.na(app_path) || !nzchar(app_path), "explorer app.R not found")
  src <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("GWAS overview", src, fixed = TRUE))
  expect_true(grepl("Locus / variants", src, fixed = TRUE))
  expect_true(grepl("plotManhattan", src, fixed = TRUE))
  expect_true(grepl("haploPlot", src, fixed = TRUE))
  expect_true(grepl("getVariantViewerData", src, fixed = TRUE))
  expect_true(grepl("plotly::plotlyOutput", src, fixed = TRUE))
  expect_error(parse(file = app_path), NA)
})

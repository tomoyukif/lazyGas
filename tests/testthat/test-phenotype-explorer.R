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
  ortho <- read.csv(file.path(extdata, "demo_orthologs.csv"), stringsAsFactors = FALSE)

  ranked <- rankPhenotypeCandidates(
    object = res$lg,
    pheno = res$pheno_name,
    query = q,
    sources = c("annotation", "gwas", "expression", "ortholog"),
    expression_matrix = expr_mat,
    expression_meta = expr_meta,
    ortholog_table = ortho,
    top_n = 10L,
    save = TRUE
  )

  expect_true(nrow(ranked) > 0L)
  expect_true("composite_score" %in% names(ranked))
  expect_true(all(ranked$composite_score >= 0 & ranked$composite_score <= 1))
  expect_true(ranked$composite_score[1L] >= ranked$composite_score[nrow(ranked)])

  stored <- lazyData(res$lg, dataset = "phenotype_rank", pheno = res$pheno_name)
  expect_false(is.null(stored))
  expect_equal(nrow(stored), nrow(ranked))
})

test_that("explainPhenotypeCandidates returns template without LLM", {
  skip_if_not_installed("lazyGas")
  cand <- make_test_candidate()
  q <- phenotypeQuery("fruit weight at maturity", tissues = "fruit", use_llm = FALSE)
  ranked <- cand
  ranked$composite_score <- c(0.9, 0.4, 0.2, 0.1)
  ranked$score_annotation <- ranked$composite_score
  ranked$score_gwas <- 0
  ranked$score_expression <- 0
  ranked$score_literature <- 0
  ranked$score_ortholog <- 0
  ranked$evidence_json <- vapply(seq_len(nrow(ranked)), function(i) {
    jsonlite::toJSON(
      list(
        annotation = list(
          source = "annotation",
          score = ranked$composite_score[i],
          snippets = ranked$Description[i]
        )
      ),
      auto_unbox = TRUE
    )
  }, character(1))
  ranked$query_id <- q$query_id
  attr(ranked, "phenotypeRank") <- list(query_id = q$query_id, pheno = "test")

  txt <- explainPhenotypeCandidates(
    rank_result = ranked,
    query = q,
    use_llm = FALSE,
    language = "en"
  )
  expect_true(grepl("Phenotype exploration", txt))
  expect_true(grepl("g1", txt))
})

test_that("collectGeneEvidence returns annotation scores", {
  skip_if_not_installed("lazyGas")
  cand <- make_test_candidate()
  q <- phenotypeQuery("fruit weight", tissues = "fruit", use_llm = FALSE)
  ev <- collectGeneEvidence(
    gene_ids = cand$Gene_ID,
    query = q,
    candidate = cand,
    sources = c("annotation", "gwas"),
    use_cache = FALSE
  )
  expect_length(ev, nrow(cand))
  expect_true(ev[["g1"]]$annotation$score >= ev[["g2"]]$annotation$score)
})

test_that("llmHealthCheck fails gracefully without server", {
  skip_if_not_installed("lazyGas")
  ok <- llmHealthCheck(base_url = "http://127.0.0.1:1", timeout = 1)
  expect_false(isTRUE(ok))
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

test_that("answerPhenotypeQuestion compares genes in Japanese without LLM", {
  skip_if_not_installed("lazyGas")
  cand <- make_test_candidate()
  q <- phenotypeQuery("fruit weight at maturity", tissues = "fruit", use_llm = FALSE)
  ranked <- cand[1:2, , drop = FALSE]
  ranked$composite_score <- c(1, 0.15)
  ranked$score_annotation <- c(1, 0)
  ranked$score_gwas <- c(1, 0.5)
  ranked$score_expression <- c(0.6, 0.6)
  ranked$score_literature <- c(1, 1)
  ranked$score_ortholog <- c(0.94, 0.88)
  ranked$evidence_json <- vapply(seq_len(2L), function(i) {
    jsonlite::toJSON(
      list(
        annotation = list(
          source = "annotation",
          score = ranked$score_annotation[i],
          snippets = if (i == 1L) "fruit weight development" else "root development"
        ),
        ortholog = list(
          source = "ortholog",
          score = ranked$score_ortholog[i],
          snippets = if (i == 1L) "Os02g01010" else "Os02g01020"
        )
      ),
      auto_unbox = TRUE
    )
  }, character(1))

  ans <- answerPhenotypeQuestion(
    rank_result = ranked,
    question = "どちらの遺伝子がより候補として信頼できる？",
    query = q,
    use_llm = FALSE,
    language = "ja"
  )
  expect_true(grepl("g1", ans))
  expect_true(grepl("信頼", ans))
  expect_false(grepl("### 1\\.", ans))
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

test_that("describeOllamaStatus reports unreachable server", {
  skip_if_not_installed("lazyGas")
  st <- describeOllamaStatus(base_url = "http://127.0.0.1:1", model = "llama3.2:3b")
  expect_true(grepl("Not running|not configured|CLI not found", st))
})

test_that("lazyGasOllamaHome returns a directory path", {
  skip_if_not_installed("lazyGas")
  path <- lazyGasOllamaHome()
  expect_type(path, "character")
  expect_true(nzchar(path))
})

test_that("Apptainer Ollama helpers are exported", {
  skip_if_not_installed("lazyGas")
  exports <- getNamespaceExports("lazyGas")
  expect_true(all(c(
    "configureLazyGasOllama",
    "buildLazyGasOllamaSif",
    "lazyGasOllamaHome"
  ) %in% exports))
})

test_that("Phase4 exports are registered", {
  skip_if_not_installed("lazyGas")
  exports <- getNamespaceExports("lazyGas")
  expect_true(all(c(
    "phenotypeQuery",
    "collectGeneEvidence",
    "rankPhenotypeCandidates",
    "answerPhenotypeQuestion",
    "llmChat",
    "llmHealthCheck",
    "restorePhenoFromStore",
    "resolveLazyGasPaths",
    "startOllamaServer",
    "startLocalLLM",
    "pullOllamaModel",
    "describeOllamaStatus",
    "configureLazyGasOllama",
    "buildLazyGasOllamaSif",
    "lazyGasOllamaHome",
    "runExplorerDemo",
    "runLazyGasExplorer",
    "runLazyGasRunner"
  ) %in% exports))
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

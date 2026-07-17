################################################################################
#' Run phenotype explorer demo
#'
#' Executes a minimal GWAS-to-candidate workflow on bundled demo data, ranks
#' genes by phenotype-relevant evidence, and optionally writes an explanation
#' report. Does not start the Shiny app; use [runLazyGasExplorer()] for the UI.
#'
#' @param out_dir Output directory (default \code{demo_output/explorer}).
#' @param overwrite Rebuild companion-store results.
#' @param trait_text Phenotype description for ranking.
#' @param tissues Tissue / organ hint.
#' @param stage Developmental stage hint.
#' @param use_llm Use local LLM for query parsing and explanation.
#' @param write_report If \code{TRUE}, save a Markdown explanation to
#'   \code{out_dir/phenotype_explorer_report.md}.
#'
#' @return Invisibly, a list with \code{lg}, \code{query}, \code{ranked},
#'   \code{explanation}, and \code{pheno_name}.
#' @export
#'
#' @seealso [runLazyGasExplorer()], [rankPhenotypeCandidates()]
runExplorerDemo <- function(out_dir = file.path(getwd(), "demo_output", "explorer"),
                            overwrite = TRUE,
                            trait_text = "fruit weight at maturity",
                            tissues = "fruit",
                            stage = "maturation",
                            use_llm = FALSE,
                            write_report = TRUE) {
  if (!requireNamespace("GBScleanR", quietly = TRUE)) {
    stop("Package 'GBScleanR' is required.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required.", call. = FALSE)
  }

  extdata <- .demo_extdata_dir()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gds_copy <- file.path(out_dir, "sample.gds")
  file.copy(file.path(extdata, "sample.gds"), gds_copy, overwrite = TRUE)

  gff <- rtracklayer::import.gff(file.path(extdata, "demo_annotation.gff"))
  ann <- read.csv(file.path(extdata, "demo_ann.csv"), stringsAsFactors = FALSE)
  ortho <- read.csv(file.path(extdata, "demo_orthologs.csv"), stringsAsFactors = FALSE)

  snpeff_gds_fn <- file.path(out_dir, "demo_snpeff.gds")
  if (!file.exists(snpeff_gds_fn) || overwrite) {
    snpeff2gds(
      vcf_fn = file.path(extdata, "demo_snpeff.vcf"),
      out_fn = snpeff_gds_fn,
      verbose = FALSE,
      overwrite = overwrite
    )
  }
  snpeff <- open_snpeff(snpeff_gds_fn)

  lg <- buildLazyGas(
    gds_fn = gds_copy,
    load_filter = TRUE,
    overwrite = overwrite,
    lazygas_store = "parquet"
  )
  on.exit(try(closeGDS(lg, verbose = FALSE), silent = TRUE), add = TRUE)

  pheno <- read.csv(file.path(extdata, "demo_pheno_multitrait.csv"))
  trait_names <- c("fruit_weight", "fruit_length", "flowering_time")
  lg <- assignPheno(object = lg, pheno = pheno, rename = trait_names)
  pheno_name <- trait_names[1L]

  conv_fun <- makeConvFun(geno_format = "dosage", n_levels = 3L)
  lg <- runLazyGas(
    object = lg,
    steps = c("scan", "peakcall", "recalc", "candidate"),
    resume = !overwrite,
    gff = gff,
    snpeff = snpeff,
    ann = ann,
    recalc = TRUE,
    seed = 42L,
    formula = "add + dom",
    conv_fun = conv_fun,
    geno_format = "dosage",
    limit_peakcall = 3L,
    n_threads = 1L
  )

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

  q <- phenotypeQuery(
    text = trait_text,
    tissues = tissues,
    stage = stage,
    use_llm = use_llm,
    object = lg,
    save = TRUE
  )

  ranked <- rankPhenotypeCandidates(
    object = lg,
    pheno = pheno_name,
    query = q,
    expression_matrix = expr_mat,
    expression_meta = expr_meta,
    ortholog_table = ortho,
    sources = c("annotation", "gwas", "expression", "ortholog"),
    top_n = 20L,
    save = TRUE
  )

  explanation <- explainPhenotypeCandidates(
    rank_result = ranked,
    query = q,
    object = lg,
    top_n = 10L,
    use_llm = use_llm,
    language = "en"
  )

  if (write_report) {
    report_fn <- file.path(out_dir, "phenotype_explorer_report.md")
    writeLines(explanation, report_fn)
    message("Wrote ", report_fn)
  }

  message("\n--- Top phenotype-ranked candidates (", pheno_name, ") ---")
  show_cols <- intersect(
    names(ranked),
    c("Gene_ID", "Name", "composite_score", "score_annotation", "score_gwas", "score_expression")
  )
  print(ranked[, show_cols, drop = FALSE])
  message("Launch UI: runLazyGasExplorer()")

  invisible(list(
    lg = lg,
    query = q,
    ranked = ranked,
    explanation = explanation,
    pheno_name = pheno_name,
    out_dir = normalizePath(out_dir, mustWork = FALSE)
  ))
}

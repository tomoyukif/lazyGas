################################################################################
#' Rank candidate genes by phenotype-relevant evidence
#'
#' Integrates annotation, GWAS, expression, literature, and ortholog evidence
#' into a composite score for each candidate gene near GWAS peaks.
#'
#' @param object A \code{LazyGas} object with candidate genes in the companion
#'   store.
#' @param pheno Phenotype name or index.
#' @param query A [PhenotypeQuery] object or character string. Character input
#'   is parsed via [phenotypeQuery()] with \code{use_llm = FALSE}.
#' @param sources Evidence sources to use (see [collectGeneEvidence()]).
#' @param weights Named numeric vector of source weights (must sum to 1 when
#'   all sources are used; missing names default to equal weight among active
#'   sources).
#' @param top_n Return at most this many rows after ranking.
#' @param expression_matrix,expression_meta Passed to [collectGeneEvidence()].
#' @param ortholog_table,ortholog_cols Ortholog data for ortholog evidence.
#' @param use_cache Use the companion-store evidence cache.
#' @param save Save ranked results to the companion store (\code{phenotype_rank}
#'   section).
#'
#' @return A \code{data.frame} sorted by \code{composite_score} with per-source
#'   score columns and \code{evidence_json}. Attributes include
#'   \code{phenotypeRank} metadata.
#' @export
#'
#' @seealso [phenotypeQuery()], [collectGeneEvidence()], [explainPhenotypeCandidates()]
#'
#' @examples
#' \dontrun{
#' q <- phenotypeQuery("fruit weight", tissues = "fruit", use_llm = FALSE)
#' ranked <- rankPhenotypeCandidates(lg, pheno = "Fruit weight", query = q)
#' head(ranked[, c("Gene_ID", "composite_score", "score_annotation", "score_gwas")])
#' }
rankPhenotypeCandidates <- function(object,
                                    pheno,
                                    query,
                                    sources = c("annotation", "gwas", "expression",
                                                "literature", "ortholog"),
                                    weights = c(
                                      annotation = 0.35,
                                      gwas = 0.25,
                                      expression = 0.2,
                                      literature = 0.15,
                                      ortholog = 0.05
                                    ),
                                    top_n = 30L,
                                    expression_matrix = NULL,
                                    expression_meta = NULL,
                                    ortholog_table = NULL,
                                    ortholog_cols = c(
                                      gene_id = "gene_id",
                                      ortholog_id = "ortholog_id",
                                      score = "score"
                                    ),
                                    use_cache = TRUE,
                                    save = TRUE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }

  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  if (is.null(candidate) || nrow(candidate) == 0L) {
    stop("No candidate data found for phenotype '", pheno_name, "'.",
         call. = FALSE)
  }
  if (!"Gene_ID" %in% names(candidate)) {
    stop("Candidate data lacks Gene_ID column.", call. = FALSE)
  }

  if (is.character(query)) {
    query <- phenotypeQuery(
      text = query,
      use_llm = FALSE,
      object = object,
      save = save
    )
  } else if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  } else if (save && !is.null(query$query_id)) {
    .store_write_phenotype_query(object, query)
  }

  sources <- match.arg(sources, several.ok = TRUE)
  active_weights <- weights[names(weights) %in% sources]
  if (!length(active_weights)) {
    active_weights <- rep(1 / length(sources), length(sources))
    names(active_weights) <- sources
  }
  if (any(active_weights < 0)) {
    stop("'weights' must be non-negative.", call. = FALSE)
  }
  w_sum <- sum(active_weights)
  if (w_sum <= 0) {
    stop("'weights' must have a positive sum.", call. = FALSE)
  }
  active_weights <- active_weights / w_sum

  gene_ids <- unique(as.character(candidate$Gene_ID))
  evidence <- collectGeneEvidence(
    gene_ids = gene_ids,
    query = query,
    candidate = candidate,
    object = object,
    pheno = pheno_name,
    sources = sources,
    expression_matrix = expression_matrix,
    expression_meta = expression_meta,
    ortholog_table = ortholog_table,
    ortholog_cols = ortholog_cols,
    use_cache = use_cache
  )

  ann_scores <- vapply(gene_ids, function(g) .evidence_source_score(evidence[[g]], "annotation"), numeric(1))
  gwas_scores <- vapply(gene_ids, function(g) .evidence_source_score(evidence[[g]], "gwas"), numeric(1))
  expr_scores <- vapply(gene_ids, function(g) .evidence_source_score(evidence[[g]], "expression"), numeric(1))
  lit_scores <- vapply(gene_ids, function(g) .evidence_source_score(evidence[[g]], "literature"), numeric(1))
  ortho_scores <- vapply(gene_ids, function(g) .evidence_source_score(evidence[[g]], "ortholog"), numeric(1))

  score_map <- list(
    annotation = ann_scores,
    gwas = gwas_scores,
    expression = expr_scores,
    literature = lit_scores,
    ortholog = ortho_scores
  )

  raw_gwas <- vapply(gene_ids, function(g) {
    if (is.null(evidence[[g]]$gwas)) {
      return(0)
    }
    as.numeric(evidence[[g]]$gwas$score)
  }, numeric(1))
  gwas_norm <- .normalize_score_vec(raw_gwas)
  names(gwas_norm) <- gene_ids

  composite <- rep(0, length(gene_ids))
  names(composite) <- gene_ids
  for (src in sources) {
    sc <- score_map[[src]]
    names(sc) <- gene_ids
    if (src == "gwas") {
      sc <- gwas_norm
    } else {
      sc <- .normalize_score_vec(sc)
      names(sc) <- gene_ids
    }
    composite <- composite + active_weights[src] * sc
  }

  score_df <- data.frame(
    Gene_ID = gene_ids,
    score_annotation = .normalize_score_vec(ann_scores),
    score_gwas = gwas_norm,
    score_expression = .normalize_score_vec(expr_scores),
    score_literature = .normalize_score_vec(lit_scores),
    score_ortholog = .normalize_score_vec(ortho_scores),
    composite_score = as.numeric(composite),
    evidence_json = vapply(gene_ids, function(g) .evidence_to_json(evidence[[g]]), character(1L)),
    query_id = query$query_id,
    stringsAsFactors = FALSE
  )

  out <- merge(candidate, score_df, by = "Gene_ID", all.x = FALSE, sort = FALSE)
  out <- out[order(-out$composite_score, out$Gene_ID), , drop = FALSE]
  rownames(out) <- NULL

  if ("Gene_ID" %in% names(out)) {
    ord <- order(-out$composite_score, out$Gene_ID)
    out <- out[ord, , drop = FALSE]
    out <- out[!duplicated(out$Gene_ID), , drop = FALSE]
  }

  top_n <- as.integer(top_n)[1L]
  if (is.finite(top_n) && top_n > 0L && nrow(out) > top_n) {
    out <- out[seq_len(top_n), , drop = FALSE]
  }

  attr(out, "phenotypeRank") <- list(
    query_id = query$query_id,
    pheno = pheno_name,
    sources = sources,
    weights = active_weights,
    n_genes = nrow(out)
  )

  if (save) {
    .store_write_section_df(
      object,
      "phenotype_rank",
      paste0("rank_", .store_safe_name(query$query_id)),
      out,
      pheno_name
    )
    meta <- .store_read_meta(object)
    prev <- meta$phenotype_rank_latest
    if (is.null(prev)) {
      prev <- list()
    }
    prev[[pheno_name]] <- query$query_id
    meta$phenotype_rank_latest <- prev
    .store_write_meta(object, meta)
  }

  out
}

################################################################################
#' Rank candidate genes by phenotype-relevant evidence
#'
#' Integrates annotation, SnpEff, GWAS, and expression (and optionally
#' literature / ortholog) evidence into a composite score for each candidate
#' gene near GWAS peaks.
#'
#' Phase 1 defaults: sources
#' \code{c("annotation","snpeff","gwas","expression")} with weights
#' annotation 0.30 / snpeff 0.35 / gwas 0.25 / expression 0.10.
#'
#' @param object A \code{LazyGas} object with candidate genes in the companion
#'   store.
#' @param pheno Phenotype name or index.
#' @param query A [PhenotypeQuery] object or character string. Character input
#'   is parsed via [phenotypeQuery()] with \code{use_llm = FALSE} unless
#'   \code{query_use_llm = TRUE}.
#' @param sources Evidence sources to use (see [collectGeneEvidence()]).
#' @param weights Named numeric vector of source weights (renormalized over
#'   active sources). Default: [phase1DefaultWeights()].
#' @param top_n Return at most this many rows after ranking.
#' @param candidate Optional candidate \code{data.frame}; when \code{NULL},
#'   candidates are loaded from the companion store.
#' @param expression_matrix,expression_meta Passed to [collectGeneEvidence()].
#' @param ortholog_table,ortholog_cols Ortholog data for ortholog evidence.
#' @param use_cache Use the companion-store evidence cache.
#' @param save Save ranked results to the companion store (\code{phenotype_rank}
#'   section).
#' @param use_keyword,use_llm_relevance,download_tenor,query_use_llm
#'   Passed through to evidence / query helpers.
#' @param use_semantic Ignored; LSA annotation scoring was removed.
#' @param llm_model,llm_base_url,llm_timeout LLM settings for annotation
#'   relevance (when enabled).
#'
#' @return A \code{data.frame} sorted by \code{composite_score} with per-source
#'   score columns and \code{evidence_json}. Attributes include
#'   \code{phenotypeRank} metadata.
#' @export
#'
#' @seealso [phenotypeQuery()], [collectGeneEvidence()], [llm_report()]
#'
#' @examples
#' \dontrun{
#' q <- phenotypeQuery("heading date", tissues = "leaf", use_llm = FALSE)
#' ranked <- rankPhenotypeCandidates(lg, pheno = "Heading date", query = q)
#' head(ranked[, c("Gene_ID", "composite_score", "score_annotation", "score_snpeff")])
#' }
rankPhenotypeCandidates <- function(object,
                                    pheno,
                                    query,
                                    sources = phase1DefaultSources(),
                                    weights = phase1DefaultWeights(),
                                    top_n = 30L,
                                    candidate = NULL,
                                    expression_matrix = NULL,
                                    expression_meta = NULL,
                                    ortholog_table = NULL,
                                    ortholog_cols = c(
                                      gene_id = "gene_id",
                                      ortholog_id = "ortholog_id",
                                      score = "score"
                                    ),
                                    use_cache = TRUE,
                                    save = TRUE,
                                    use_keyword = NULL,
                                    use_semantic = NULL,
                                    use_llm_relevance = NULL,
                                    download_tenor = TRUE,
                                    query_use_llm = FALSE,
                                    llm_model = NULL,
                                    llm_base_url = NULL,
                                    llm_timeout = NULL) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }

  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  if (is.null(candidate)) {
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  }
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
      use_llm = isTRUE(query_use_llm),
      object = object,
      save = save
    )
  } else if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  } else if (save && !is.null(query$query_id)) {
    .store_write_phenotype_query(object, query)
  }

  sources <- match.arg(sources, choices = .EVIDENCE_SOURCE_CHOICES, several.ok = TRUE)
  if (is.null(weights) || !length(weights)) {
    weights <- phase1DefaultWeights()
  }
  active_weights <- setNames(rep(NA_real_, length(sources)), sources)
  known <- intersect(names(weights), sources)
  if (length(known)) {
    active_weights[known] <- as.numeric(weights[known])
  }
  if (anyNA(active_weights)) {
    fill <- mean(active_weights, na.rm = TRUE)
    if (!is.finite(fill) || fill <= 0) {
      fill <- 1
    }
    active_weights[is.na(active_weights)] <- fill
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
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  if (!length(gene_ids)) {
    stop(
      "No usable Gene_ID values in candidate table for phenotype '",
      pheno_name, "'.",
      call. = FALSE
    )
  }
  # Drop rows without Gene_ID before scoring / merge
  candidate <- candidate[
    !is.na(candidate$Gene_ID) & nzchar(as.character(candidate$Gene_ID)),
    ,
    drop = FALSE
  ]
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
    use_cache = use_cache,
    use_keyword = use_keyword,
    use_semantic = use_semantic,
    use_llm_relevance = use_llm_relevance,
    download_tenor = download_tenor,
    llm_model = llm_model,
    llm_base_url = llm_base_url,
    llm_timeout = llm_timeout
  )

  score_raw <- list()
  for (src in .EVIDENCE_SOURCE_CHOICES) {
    score_raw[[src]] <- vapply(
      gene_ids,
      function(g) .evidence_source_score(evidence[[g]], src),
      numeric(1)
    )
  }

  raw_gwas <- vapply(gene_ids, function(g) {
    if (is.null(evidence[[g]]$gwas)) {
      return(0)
    }
    as.numeric(evidence[[g]]$gwas$score)
  }, numeric(1))
  gwas_norm <- .normalize_score_vec(raw_gwas)
  names(gwas_norm) <- gene_ids

  score_norm <- list()
  for (src in .EVIDENCE_SOURCE_CHOICES) {
    if (src == "gwas") {
      score_norm[[src]] <- gwas_norm
    } else {
      sc <- .normalize_score_vec(score_raw[[src]])
      names(sc) <- gene_ids
      score_norm[[src]] <- sc
    }
  }

  composite <- rep(0, length(gene_ids))
  names(composite) <- gene_ids
  for (src in sources) {
    composite <- composite + active_weights[src] * score_norm[[src]]
  }

  score_df <- data.frame(
    Gene_ID = gene_ids,
    score_annotation = score_norm$annotation,
    score_snpeff = score_norm$snpeff,
    score_gwas = gwas_norm,
    score_expression = score_norm$expression,
    score_literature = score_norm$literature,
    score_ortholog = score_norm$ortholog,
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

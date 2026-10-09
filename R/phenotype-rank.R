################################################################################
#' Rank candidate genes by phenotype-relevant evidence
#'
#' Integrates annotation, SnpEff, GWAS, and expression (and optionally
#' literature / ortholog) evidence into a composite score for each candidate
#' gene **within each GWAS peak**. Channel scores are already on \code{[0, 1]}
#' and are **not** re-normalized within the peak. Inactive channels are stored
#' as 0.
#' The same gene on two peaks keeps two rows. There is no cross-peak
#' \code{Gene_ID} deduplication.
#'
#' Defaults: sources \code{c("snpeff","gwas_dist","gwas_peak")}. Default channel
#' weights come from \code{inst/config/ai-lazygas-models.yaml} (via
#' [phase1DefaultWeights()]) and are renormalized over active sources.
#'
#' @param object A \code{LazyGas} object with candidate genes in the companion
#'   store.
#' @param pheno Phenotype name or index.
#' @param query A [PhenotypeQuery] object or character string. Character input
#'   is parsed via [phenotypeQuery()].
#' @param sources Evidence sources to use (see [collectGeneEvidence()]).
#' @param weights Named numeric vector of source weights (renormalized over
#'   active sources). Default: [phase1DefaultWeights()].
#' @param top_n Return at most this many rows **per peak** after ranking.
#'   \code{NULL} keeps every scored gene in each peak.
#' @param candidate Optional candidate \code{data.frame}; when \code{NULL},
#'   candidates are loaded from the companion store. When \code{peak_ID} is
#'   present, ranking runs separately for each peak.
#' @param ann_column_map Required when \code{annotation} is in \code{sources}
#'   (see [readAnnColumnMap()]).
#' @param expression_matrix,expression_sample_map Required when
#'   \code{expression} is in \code{sources} (see [readExpressionSampleMap()]).
#' @param expression_meta Deprecated alias for \code{expression_sample_map}.
#' @param literature_graph,ortholog_graph Edge tables (\code{from}, \code{to})
#'   required when the matching source is active (see [collectGeneEvidence()]).
#' @param ortholog_table,ortholog_cols Deprecated; use \code{ortholog_graph}.
#' @param use_cache Use the companion-store evidence cache.
#' @param save Save ranked results to the companion store (\code{phenotype_rank}
#'   section).
#' @param use_llm_relevance,download_tenor,query_use_llm,llm_model,llm_base_url,llm_timeout
#'   Ignored; retained for call compatibility.
#' @param gff,gff_windows Passed to [collectGeneEvidence()] for finemap
#'   GFF-window fallback when SnpEff is absent.
#'
#' @return A \code{data.frame} with one row per \code{(peak_ID, Gene_ID)}
#'   (when \code{peak_ID} exists), columns including \code{rank_within_peak},
#'   \code{composite_score}, per-source scores, and \code{evidence_json}.
#'   Rows are ordered by \code{peak_ID}, then \code{rank_within_peak}.
#'   Attributes include \code{phenotypeRank} metadata.
#' @export
#'
#' @seealso [phenotypeQuery()], [collectGeneEvidence()], [readAnnColumnMap()],
#'   [readExpressionSampleMap()]
#'
#' @examples
#' \dontrun{
#' q <- phenotypeQuery("heading date", tissues = "leaf")
#' ranked <- rankPhenotypeCandidates(
#'   lg, pheno = "Heading date", query = q,
#'   ann_column_map = "ann_column_map.csv"
#' )
#' head(ranked[, c("peak_ID", "Gene_ID", "rank_within_peak", "composite_score")])
#' }
rankPhenotypeCandidates <- function(object,
                                    pheno,
                                    query,
                                    sources = phase1DefaultSources(),
                                    weights = phase1DefaultWeights(),
                                    top_n = 30L,
                                    candidate = NULL,
                                    ann_column_map = NULL,
                                    expression_matrix = NULL,
                                    expression_sample_map = NULL,
                                    expression_meta = NULL,
                                    literature_graph = NULL,
                                    ortholog_graph = NULL,
                                    ortholog_table = NULL,
                                    ortholog_cols = c(
                                      gene_id = "gene_id",
                                      ortholog_id = "ortholog_id",
                                      score = "score"
                                    ),
                                    use_cache = TRUE,
                                    save = TRUE,
                                    use_llm_relevance = NULL,
                                    download_tenor = FALSE,
                                    query_use_llm = FALSE,
                                    llm_model = NULL,
                                    llm_base_url = NULL,
                                    llm_timeout = NULL,
                                    gff = NULL,
                                    gff_windows = NULL) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }

  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  if (is.null(candidate)) {
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  }
  if (is.null(candidate) || nrow(candidate) == 0L) {
    has_simple <- .store_dataset_exists(object, "simple_candidate", pheno_name)
    stop(
      "No wide candidate data found for phenotype '", pheno_name, "'. ",
      if (has_simple) {
        "A simple_candidate table exists but ranking requires the wide candidate. "
      } else {
        ""
      },
      "Run listCandidate() (wide tables are always written) before ranking.",
      call. = FALSE
    )
  }
  if (!"Gene_ID" %in% names(candidate)) {
    stop("Candidate data lacks Gene_ID column.", call. = FALSE)
  }

  if (is.character(query)) {
    query <- phenotypeQuery(
      text = query,
      object = object,
      save = save
    )
  } else if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  } else if (save && !is.null(query$query_id)) {
    .store_write_phenotype_query(object, query)
  }

  if (is.null(expression_sample_map) && !is.null(expression_meta)) {
    expression_sample_map <- expression_meta
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

  candidate <- candidate[
    !is.na(candidate$Gene_ID) & nzchar(as.character(candidate$Gene_ID)),
    ,
    drop = FALSE
  ]
  if (!nrow(candidate)) {
    stop(
      "No usable Gene_ID values in candidate table for phenotype '",
      pheno_name, "'.",
      call. = FALSE
    )
  }

  peak_groups <- .rank_peak_group_ids(candidate)
  parts <- vector("list", length(peak_groups$keys))
  for (i in seq_along(peak_groups$keys)) {
    cand_peak <- candidate[peak_groups$index[[i]], , drop = FALSE]
    part <- .rank_phenotype_candidates_block(
      object = object,
      pheno_name = pheno_name,
      query = query,
      candidate = cand_peak,
      sources = sources,
      active_weights = active_weights,
      ann_column_map = ann_column_map,
      expression_matrix = expression_matrix,
      expression_sample_map = expression_sample_map,
      literature_graph = literature_graph,
      ortholog_graph = ortholog_graph,
      ortholog_table = ortholog_table,
      ortholog_cols = ortholog_cols,
      use_cache = use_cache,
      gff = gff,
      gff_windows = gff_windows
    )
    part$rank_within_peak <- seq_len(nrow(part))
    if (!is.null(top_n)) {
      tn <- as.integer(top_n)[1L]
      if (is.finite(tn) && tn > 0L && nrow(part) > tn) {
        part <- part[seq_len(tn), , drop = FALSE]
      }
    }
    parts[[i]] <- part
  }

  out <- .bind_rank_peak_parts(parts)
  rownames(out) <- NULL

  attr(out, "phenotypeRank") <- list(
    query_id = query$query_id,
    pheno = pheno_name,
    sources = sources,
    weights = active_weights,
    n_rows = nrow(out),
    n_genes = length(unique(as.character(out$Gene_ID))),
    n_peaks = length(peak_groups$keys),
    ranking_unit = "peak"
  )

  if (save) {
    .store_write_phenotype_rank(
      object = object,
      rank_result = out,
      pheno_name = pheno_name,
      query_id = query$query_id
    )
  }

  out
}

#' Split candidate rows into peak groups (order of first appearance)
#' @noRd
.rank_peak_group_ids <- function(candidate) {
  if (!"peak_ID" %in% names(candidate)) {
    return(list(
      keys = list(NA),
      index = list(seq_len(nrow(candidate)))
    ))
  }
  keys <- unique(candidate$peak_ID)
  index <- lapply(keys, function(pid) {
    if (is.na(pid)) {
      which(is.na(candidate$peak_ID))
    } else {
      which(!is.na(candidate$peak_ID) & candidate$peak_ID == pid)
    }
  })
  list(keys = as.list(keys), index = index)
}

#' Bind per-peak rank tables; order by peak_ID then rank_within_peak
#' @noRd
.bind_rank_peak_parts <- function(parts) {
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0L, parts)
  if (!length(parts)) {
    stop("Ranking produced no rows for any peak.", call. = FALSE)
  }
  out <- do.call(rbind, parts)
  if ("peak_ID" %in% names(out) && "rank_within_peak" %in% names(out)) {
    out <- out[order(out$peak_ID, out$rank_within_peak, out$Gene_ID), , drop = FALSE]
  } else if ("rank_within_peak" %in% names(out)) {
    out <- out[order(out$rank_within_peak, out$Gene_ID), , drop = FALSE]
  }
  out
}

#' Score and rank one peak's candidate rows (no Gene_ID dedup)
#' @noRd
.rank_phenotype_candidates_block <- function(object,
                                             pheno_name,
                                             query,
                                             candidate,
                                             sources,
                                             active_weights,
                                             ann_column_map = NULL,
                                             expression_matrix = NULL,
                                             expression_sample_map = NULL,
                                             literature_graph = NULL,
                                             ortholog_graph = NULL,
                                             ortholog_table = NULL,
                                             ortholog_cols = NULL,
                                             use_cache = TRUE,
                                             gff = NULL,
                                             gff_windows = NULL,
                                             collect_extra = list()) {
  gene_ids <- unique(as.character(candidate$Gene_ID))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  if (!length(gene_ids)) {
    return(candidate[FALSE, , drop = FALSE])
  }

  evidence_args <- c(
    list(
      gene_ids = gene_ids,
      query = query,
      candidate = candidate,
      object = object,
      pheno = pheno_name,
      sources = sources,
      ann_column_map = ann_column_map,
      expression_matrix = expression_matrix,
      expression_sample_map = expression_sample_map,
      literature_graph = literature_graph,
      ortholog_graph = ortholog_graph,
      ortholog_table = ortholog_table,
      ortholog_cols = ortholog_cols,
      use_cache = use_cache,
      gff = gff,
      gff_windows = gff_windows
    ),
    collect_extra
  )
  evidence <- do.call(collectGeneEvidence, evidence_args)

  # Raw channel scores are already on [0, 1]. No per-peak min–max.
  # Inactive sources are forced to 0 in the output table (not collected/normalized).
  score_out <- list()
  for (src in .EVIDENCE_SOURCE_CHOICES) {
    if (src %in% sources) {
      sc <- vapply(
        gene_ids,
        function(g) .evidence_source_score(evidence[[g]], src),
        numeric(1)
      )
      names(sc) <- gene_ids
      score_out[[src]] <- sc
    } else {
      score_out[[src]] <- setNames(rep(0, length(gene_ids)), gene_ids)
    }
  }

  composite <- rep(0, length(gene_ids))
  names(composite) <- gene_ids
  for (src in sources) {
    composite <- composite + active_weights[src] * score_out[[src]]
  }

  score_df <- data.frame(
    Gene_ID = gene_ids,
    score_annotation = score_out$annotation,
    score_snpeff = score_out$snpeff,
    score_gwas_dist = score_out$gwas_dist,
    score_gwas_peak = score_out$gwas_peak,
    score_expression = score_out$expression,
    score_literature = score_out$literature,
    score_ortholog = score_out$ortholog,
    score_finemap = score_out$finemap,
    composite_score = as.numeric(composite),
    evidence_json = vapply(gene_ids, function(g) .evidence_to_json(evidence[[g]]), character(1L)),
    query_id = query$query_id,
    stringsAsFactors = FALSE
  )

  out <- merge(candidate, score_df, by = "Gene_ID", all.x = FALSE, sort = FALSE)
  out <- out[order(-out$composite_score, out$Gene_ID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Save a phenotype-rank table to the companion store
#' @keywords internal
.store_write_phenotype_rank <- function(object, rank_result, pheno_name, query_id) {
  if (is.null(query_id) || !nzchar(as.character(query_id)[1L])) {
    return(invisible(NULL))
  }
  query_id <- as.character(query_id)[1L]
  .store_write_section_df(
    object,
    "phenotype_rank",
    paste0("rank_", .store_safe_name(query_id)),
    rank_result,
    pheno_name
  )
  meta <- .store_read_meta(object)
  prev <- meta$phenotype_rank_latest
  if (is.null(prev)) {
    prev <- list()
  }
  prev[[pheno_name]] <- query_id
  meta$phenotype_rank_latest <- prev
  .store_write_meta(object, meta)
  invisible(NULL)
}

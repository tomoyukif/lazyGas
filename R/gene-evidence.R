################################################################################
#' Collect multi-source evidence for candidate genes
#'
#' Gathers annotation, SnpEff, GWAS, expression (user matrix), literature,
#' and ortholog evidence for each gene in \code{gene_ids}, optionally using
#' cached results in the companion store.
#'
#' Default sources are \code{snpeff}, \code{gwas_dist}, and \code{gwas_peak}.
#' Other channels (\code{annotation}, \code{expression}, \ldots) remain available
#' when listed in \code{sources}. When \code{annotation} is among \code{sources},
#' scoring uses tiered phrase hits (action / noun / word → 1 / 0.5 / 0.25) and
#' [readAnnColumnMap()] input is required. When \code{expression} is among
#' \code{sources}, both \code{expression_matrix} and
#' [readExpressionSampleMap()] input are required. When \code{literature} or
#' \code{ortholog} is among \code{sources}, the matching edge table
#' (\code{literature_graph} / \code{ortholog_graph} with columns \code{from},
#' \code{to}) is required; score is 1 if a gene node reaches a phenotype node
#' on the undirected graph, else 0.
#'
#' @param gene_ids Character vector of gene identifiers.
#' @param query A [PhenotypeQuery] object or character search string.
#' @param candidate Candidate-gene \code{data.frame} (optional if \code{object}
#'   and \code{pheno} are given).
#' @param object \code{LazyGas} object for cache I/O and candidate lookup.
#' @param pheno Phenotype name or index when reading candidates from \code{object}.
#' @param sources Evidence sources to collect. Any of \code{"annotation"},
#'   \code{"snpeff"}, \code{"gwas_dist"}, \code{"gwas_peak"}, \code{"expression"},
#'   \code{"literature"}, \code{"ortholog"}, \code{"finemap"}.
#' @param ann_column_map Path or [readAnnColumnMap()] result. Required when
#'   \code{annotation} is in \code{sources}.
#' @param expression_matrix Gene x sample expression matrix (\code{data.frame}
#'   or matrix) with row names = gene IDs. Required when \code{expression}
#'   is in \code{sources}.
#' @param expression_sample_map Path or [readExpressionSampleMap()] result
#'   (or legacy meta \code{data.frame} with a \code{sample} column). Required
#'   when \code{expression} is in \code{sources}.
#' @param expression_meta Deprecated alias for \code{expression_sample_map}.
#' @param literature_graph Path or edge \code{data.frame} (\code{from}, \code{to})
#'   for the literature network channel. Required when \code{literature} is in
#'   \code{sources}.
#' @param ortholog_graph Path or edge \code{data.frame} (\code{from}, \code{to})
#'   for the ortholog network channel. Required when \code{ortholog} is in
#'   \code{sources}.
#' @param ortholog_table,ortholog_cols,literature_max Deprecated; ignored.
#'   Use \code{ortholog_graph} / \code{literature_graph}.
#' @param use_cache Read/write evidence cache in the companion store.
#' @param use_llm_relevance,download_tenor,llm_model,llm_base_url,llm_timeout
#'   Ignored; retained for call compatibility.
#' @param gff Optional \code{GRanges} for finemap GFF-window fallback when
#'   SnpEff is missing or does not assign a gene.
#' @param gff_windows Optional precomputed gene windows (from
#'   \code{.gff_gene_windows}); overrides building from \code{gff}.
#'
#' @return A named list keyed by gene ID; each element is a list of evidence
#'   records (\code{source}, \code{score}, \code{snippets}, \code{details}).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [readAnnColumnMap()],
#'   [readExpressionSampleMap()]
collectGeneEvidence <- function(gene_ids,
                                query,
                                candidate = NULL,
                                object = NULL,
                                pheno = NULL,
                                sources = phase1DefaultSources(),
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
                                literature_max = 5L,
                                use_llm_relevance = NULL,
                                download_tenor = FALSE,
                                llm_model = NULL,
                                llm_base_url = NULL,
                                llm_timeout = NULL,
                                gff = NULL,
                                gff_windows = NULL) {
  sources <- match.arg(sources, choices = .EVIDENCE_SOURCE_CHOICES, several.ok = TRUE)
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids)]
  if (!length(gene_ids)) {
    stop("'gene_ids' must contain at least one gene ID.", call. = FALSE)
  }

  search_text <- .phenotype_query_as_search_text(query)
  if (isTRUE(use_llm_relevance)) {
    warning("'use_llm_relevance' is ignored in this build.", call. = FALSE)
  }
  if (isTRUE(download_tenor)) {
    warning(
      "'download_tenor' is ignored; provide expression_matrix and expression_sample_map.",
      call. = FALSE
    )
  }

  if (is.null(expression_sample_map) && !is.null(expression_meta)) {
    expression_sample_map <- expression_meta
  }

  if (is.null(candidate) && !is.null(object)) {
    pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  }
  if (!is.null(candidate) && !"Gene_ID" %in% names(candidate)) {
    stop("'candidate' must contain a Gene_ID column.", call. = FALSE)
  }

  if ("annotation" %in% sources) {
    if (is.null(candidate)) {
      stop(
        "'candidate' (with annotation columns) is required for the annotation channel.",
        call. = FALSE
      )
    }
    ann_column_map <- .validate_ann_column_map(ann_column_map, names(candidate))
  }
  if (!is.null(ortholog_table)) {
    warning(
      "'ortholog_table' is deprecated; pass ortholog_graph (from/to edges) instead.",
      call. = FALSE
    )
  }
  if (!missing(literature_max) && !identical(literature_max, 5L)) {
    warning(
      "'literature_max' is ignored; literature uses literature_graph path reachability.",
      call. = FALSE
    )
  }

  literature_edges <- NULL
  ortholog_edges <- NULL
  if ("literature" %in% sources) {
    literature_edges <- .coerce_evidence_graph(
      literature_graph,
      arg_name = "literature_graph"
    )
  }
  if ("ortholog" %in% sources) {
    ortholog_edges <- .coerce_evidence_graph(
      ortholog_graph,
      arg_name = "ortholog_graph"
    )
  }

  if ("expression" %in% sources) {
    expression_sample_map <- .validate_expression_inputs(
      expression_matrix,
      expression_sample_map
    )
  }

  cache_key <- .evidence_cache_hash(list(
    query = if (inherits(query, "PhenotypeQuery")) {
      unclass(query)
    } else {
      query
    },
    sources = sources,
    use_llm_relevance = FALSE,
    keyword_engine = "annotation_tier_v1",
    ann_columns = if ("annotation" %in% sources) {
      .ann_ranking_columns(ann_column_map)
    } else {
      character()
    },
    literature_graph = if (!is.null(literature_edges)) {
      .evidence_graph_fingerprint(literature_edges)
    } else {
      NULL
    },
    ortholog_graph = if (!is.null(ortholog_edges)) {
      .evidence_graph_fingerprint(ortholog_edges)
    } else {
      NULL
    },
    # Peak-local scoring: same Gene_ID can have different GWAS/SnpEff rows
    peak_ids = if (!is.null(candidate) && "peak_ID" %in% names(candidate)) {
      sort(unique(as.character(candidate$peak_ID)))
    } else {
      character()
    }
  ))

  ann_batch <- NULL
  if ("annotation" %in% sources) {
    ann_batch <- .evidence_annotation_batch(
      gene_ids = gene_ids,
      candidate = candidate,
      search_text = search_text,
      query = query,
      ann_column_map = ann_column_map
    )
  }

  gwas_need <- intersect(c("gwas_dist", "gwas_peak"), sources)
  gwas_batch <- NULL
  if (length(gwas_need)) {
    pheno_name_gwas <- if (!is.null(object)) {
      .determine_phenotype_name(object = object, pheno = pheno)
    } else {
      NULL
    }
    gwas_batch <- .evidence_gwas_channels_batch(
      gene_ids = gene_ids,
      candidate = candidate,
      object = object,
      pheno_name = pheno_name_gwas
    )
  }

  finemap_lookup <- NULL
  if ("finemap" %in% sources) {
    if (is.null(object)) {
      stop("'object' is required when sources include 'finemap'.", call. = FALSE)
    }
    pheno_name_fm <- .determine_phenotype_name(object = object, pheno = pheno)
    finemap_lookup <- .finemap_pip_lookup(
      object = object,
      pheno_name = pheno_name_fm,
      candidate = candidate,
      gff = gff,
      gff_windows = gff_windows
    )
  }

  out <- setNames(vector("list", length(gene_ids)), gene_ids)

  for (gid in gene_ids) {
    cached <- if (use_cache && !is.null(object)) {
      .store_read_evidence_cache(object, gid, cache_key)
    } else {
      NULL
    }
    if (!is.null(cached)) {
      out[[gid]] <- cached
      next
    }

    gene_rows <- if (!is.null(candidate)) {
      candidate[candidate$Gene_ID == gid, , drop = FALSE]
    } else {
      NULL
    }

    ev <- list()
    if ("annotation" %in% sources) {
      ev$annotation <- ann_batch[[gid]]
    }
    if ("snpeff" %in% sources) {
      ev$snpeff <- .evidence_snpeff(gene_id = gid, gene_rows = gene_rows)
    }
    if ("gwas_dist" %in% sources) {
      ev$gwas_dist <- gwas_batch[[gid]]$gwas_dist
    }
    if ("gwas_peak" %in% sources) {
      ev$gwas_peak <- gwas_batch[[gid]]$gwas_peak
    }
    if ("expression" %in% sources) {
      ev$expression <- .evidence_expression(
        gene_id = gid,
        query = query,
        expression_matrix = expression_matrix,
        expression_sample_map = expression_sample_map
      )
    }
    if ("literature" %in% sources) {
      ev$literature <- .evidence_network_channel(
        source = "literature",
        gene_id = gid,
        gene_rows = gene_rows,
        query = query,
        edges = literature_edges
      )
    }
    if ("ortholog" %in% sources) {
      ev$ortholog <- .evidence_network_channel(
        source = "ortholog",
        gene_id = gid,
        gene_rows = gene_rows,
        query = query,
        edges = ortholog_edges
      )
    }
    if ("finemap" %in% sources) {
      ev$finemap <- .evidence_finemap(
        gene_id = gid,
        gene_rows = gene_rows,
        lookup = finemap_lookup
      )
    }

    out[[gid]] <- ev
    if (use_cache && !is.null(object)) {
      .store_write_evidence_cache(object, gid, cache_key, ev)
    }
  }

  attr(out, "collectGeneEvidence") <- list(
    n_genes = length(gene_ids),
    sources = sources,
    cache_key = cache_key,
    search_text = search_text
  )
  out
}

.evidence_cache_hash <- function(x) {
  s <- jsonlite::toJSON(x, auto_unbox = TRUE, null = "null")
  r <- as.integer(charToRaw(s))
  val <- (sum(r) * 31 + length(r) * 17) %% 4294967296
  sprintf("%08x", val)
}

.store_evidence_cache_path <- function(object, gene_id, cache_key) {
  file.path(
    .store_path(object),
    "evidence",
    "cache",
    paste0(.store_safe_name(gene_id), "_", cache_key, ".json")
  )
}

.store_read_evidence_cache <- function(object, gene_id, cache_key) {
  if (.store_is_gds(object)) {
    return(NULL)
  }
  path <- .store_evidence_cache_path(object, gene_id, cache_key)
  if (!file.exists(path)) {
    return(NULL)
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

.store_write_evidence_cache <- function(object, gene_id, cache_key, evidence) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  dir.create(
    file.path(.store_path(object), "evidence", "cache"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  path <- .store_evidence_cache_path(object, gene_id, cache_key)
  jsonlite::write_json(evidence, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

.evidence_annotation_batch <- function(gene_ids,
                                       candidate,
                                       search_text,
                                       query,
                                       ann_column_map) {
  empty <- function(gid) {
    list(
      source = "annotation",
      score = 0,
      snippets = character(),
      details = list(
        keyword_score = 0,
        matched_keywords = character(),
        matched_synonym = list(),
        matched_related = list(),
        unmatched_keywords = character(),
        matched_context = character(),
        match_method = "none"
      )
    )
  }

  out <- setNames(lapply(gene_ids, empty), gene_ids)
  if (is.null(candidate) || !nrow(candidate)) {
    return(out)
  }

  ann_cols <- .ann_ranking_columns(ann_column_map)
  .validate_ann_column_map(ann_column_map, names(candidate))

  texts <- character(length(gene_ids))
  names(texts) <- gene_ids
  snippets <- vector("list", length(gene_ids))
  names(snippets) <- gene_ids
  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[[i]]
    rows <- candidate[candidate$Gene_ID == gid, , drop = FALSE]
    if (!nrow(rows)) {
      texts[[i]] <- ""
      snippets[[i]] <- character()
      next
    }
    ann_text <- .candidate_annotation_text(candidate = rows, ann_cols = ann_cols)
    keep <- nzchar(ann_text)
    texts[[i]] <- if (any(keep)) paste(unique(ann_text[keep]), collapse = " | ") else ""
    snippets[[i]] <- ann_text[keep]
  }

  details_list <- .annotation_tier_score_details(
    text = unname(texts),
    query = query,
    ignore.case = TRUE
  )
  kw_details <- setNames(details_list, gene_ids)

  for (gid in gene_ids) {
    det <- kw_details[[gid]] %||% list(
      keyword_score = 0,
      annotation_tier = "none",
      matched_action_phrases = character(),
      matched_noun_phrases = character(),
      matched_words = character(),
      matched_keywords = character(),
      unmatched_keywords = character()
    )
    kw <- as.numeric(det$keyword_score %||% 0)
    if (!is.finite(kw)) kw <- 0
    out[[gid]] <- list(
      source = "annotation",
      score = kw,
      snippets = snippets[[gid]] %||% character(),
      details = list(
        keyword_score = kw,
        annotation_tier = as.character(det$annotation_tier %||% "none")[1L],
        matched_action_phrases = as.character(
          det$matched_action_phrases %||% character()
        ),
        matched_noun_phrases = as.character(
          det$matched_noun_phrases %||% character()
        ),
        matched_words = as.character(det$matched_words %||% character()),
        matched_keywords = as.character(det$matched_keywords %||% character()),
        unmatched_keywords = as.character(det$unmatched_keywords %||% character()),
        match_method = if (kw > 0) "annotation_tier" else "none"
      )
    )
  }
  out
}


.evidence_snpeff <- function(gene_id, gene_rows) {
  empty <- list(
    source = "snpeff",
    score = 0,
    snippets = character(),
    details = list(worst_impact = "none")
  )
  if (is.null(gene_rows) || nrow(gene_rows) == 0L) {
    return(empty)
  }

  impact_cols <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  counts <- setNames(rep(0, length(impact_cols)), impact_cols)
  for (col in impact_cols) {
    at <- paste0(col, "_at_var")
    if (at %in% names(gene_rows)) {
      counts[[col]] <- sum(as.numeric(gene_rows[[at]]), na.rm = TRUE)
    } else if (col %in% names(gene_rows)) {
      counts[[col]] <- sum(as.numeric(gene_rows[[col]]), na.rm = TRUE)
    }
  }

  if (isTRUE(counts[["HIGH"]] > 0)) {
    score <- 1.0
    worst <- "HIGH"
  } else if (isTRUE(counts[["MODERATE"]] > 0)) {
    score <- 0.5
    worst <- "MODERATE"
  } else if (isTRUE(counts[["LOW"]] > 0) || isTRUE(counts[["MODIFIER"]] > 0)) {
    score <- 0.25
    worst <- if (counts[["LOW"]] > 0) "LOW" else "MODIFIER"
  } else {
    score <- 0.0
    worst <- "none"
  }

  snippets <- character()
  if (score > 0) {
    snippets <- sprintf(
      "worst_impact=%s (HIGH=%s MODERATE=%s LOW=%s MODIFIER=%s)",
      worst,
      counts[["HIGH"]],
      counts[["MODERATE"]],
      counts[["LOW"]],
      counts[["MODIFIER"]]
    )
  }

  list(
    source = "snpeff",
    score = score,
    snippets = snippets,
    details = list(
      worst_impact = worst,
      HIGH = counts[["HIGH"]],
      MODERATE = counts[["MODERATE"]],
      LOW = counts[["LOW"]],
      MODIFIER = counts[["MODIFIER"]]
    )
  )
}

#' Empty GWAS channel evidence record
#' @keywords internal
.evidence_gwas_empty <- function(source) {
  list(source = source, score = 0, snippets = character(), details = list())
}

#' Load peak-block table (recalc, else peakcall)
#' @keywords internal
.gwas_load_peak_blocks <- function(object, pheno_name) {
  if (is.null(object) || is.null(pheno_name) || !nzchar(as.character(pheno_name)[[1L]])) {
    return(NULL)
  }
  pk <- tryCatch(
    lazyData(object = object, dataset = "recalc", pheno = pheno_name),
    error = function(e) NULL
  )
  if (is.null(pk) || !nrow(pk)) {
    pk <- tryCatch(
      lazyData(object = object, dataset = "peakcall", pheno = pheno_name),
      error = function(e) NULL
    )
  }
  if (is.null(pk) || !nrow(pk)) {
    return(NULL)
  }
  pk
}

#' Per-marker GWAS dist/peak scores within one peak block
#' @keywords internal
.gwas_marker_channel_scores <- function(blocks) {
  empty <- list(
    d_max = NA_real_,
    lead_nlp = NA_real_,
    markers = data.frame(
      Chr = character(),
      Pos = numeric(),
      abs_dist = numeric(),
      negLog10P = numeric(),
      dist_score = numeric(),
      peak_score = numeric(),
      stringsAsFactors = FALSE
    )
  )
  if (is.null(blocks) || !nrow(blocks) || !"Pos" %in% names(blocks)) {
    return(empty)
  }
  pos <- as.numeric(blocks$Pos)
  chr <- if ("Chr" %in% names(blocks)) {
    as.character(blocks$Chr)
  } else {
    rep(NA_character_, nrow(blocks))
  }
  if ("dist2peak" %in% names(blocks)) {
    abs_dist <- abs(as.numeric(blocks$dist2peak))
  } else {
    is_lead <- if ("peak_variant_ID" %in% names(blocks) &&
                   "variant_ID" %in% names(blocks)) {
      as.character(blocks$variant_ID) == as.character(blocks$peak_variant_ID)
    } else {
      rep(FALSE, nrow(blocks))
    }
    lead_pos <- if (any(is_lead, na.rm = TRUE)) {
      pos[which(is_lead)[1L]]
    } else {
      NA_real_
    }
    abs_dist <- if (is.finite(lead_pos)) abs(pos - lead_pos) else rep(NA_real_, length(pos))
  }
  nlp <- if ("negLog10P" %in% names(blocks)) {
    as.numeric(blocks$negLog10P)
  } else {
    rep(NA_real_, nrow(blocks))
  }

  d_max <- suppressWarnings(max(abs_dist[is.finite(abs_dist)], na.rm = TRUE))
  if (!is.finite(d_max)) {
    d_max <- NA_real_
  }

  is_lead <- if ("peak_variant_ID" %in% names(blocks) &&
                 "variant_ID" %in% names(blocks)) {
    as.character(blocks$variant_ID) == as.character(blocks$peak_variant_ID)
  } else {
    # Prefer tallest marker when lead IDs are absent.
    if (any(is.finite(nlp))) {
      seq_along(nlp) == which.max(replace(nlp, !is.finite(nlp), -Inf))
    } else {
      rep(FALSE, nrow(blocks))
    }
  }
  lead_nlp <- if (any(is_lead & is.finite(nlp), na.rm = TRUE)) {
    nlp[which(is_lead & is.finite(nlp))[1L]]
  } else {
    suppressWarnings(max(nlp[is.finite(nlp)], na.rm = TRUE))
  }
  if (!is.finite(lead_nlp)) {
    lead_nlp <- NA_real_
  }

  dist_score <- if (!is.finite(d_max) || d_max <= 0) {
    ifelse(is.finite(abs_dist) & abs_dist <= 0, 1, 0)
  } else {
    ifelse(is.finite(abs_dist), pmax(0, pmin(1, (d_max - abs_dist) / d_max)), 0)
  }
  peak_score <- if (!is.finite(lead_nlp) || lead_nlp <= 0) {
    rep(0, length(nlp))
  } else {
    ifelse(is.finite(nlp), pmax(0, pmin(1, nlp / lead_nlp)), 0)
  }

  list(
    d_max = d_max,
    lead_nlp = lead_nlp,
    markers = data.frame(
      Chr = chr,
      Pos = pos,
      abs_dist = abs_dist,
      negLog10P = nlp,
      dist_score = as.numeric(dist_score),
      peak_score = as.numeric(peak_score),
      stringsAsFactors = FALSE
    )
  )
}

#' Marker rows attributed to one gene within a peak
#'
#' Preference: SnpEff Gene_ID ∩ peak markers by Chr/Pos; else peak markers
#' whose Pos falls in the gene interval; else a synthetic marker from the
#' candidate row (\code{dist2peak} / \code{negLog10P}).
#' @keywords internal
.gwas_markers_for_gene <- function(gid,
                                   gene_rows,
                                   marker_scores,
                                   snpeff_peak) {
  markers <- marker_scores$markers
  empty <- markers[0, , drop = FALSE]
  if (!nrow(markers)) {
    return(empty)
  }

  if (!is.null(snpeff_peak) && nrow(snpeff_peak) &&
      "Gene_ID" %in% names(snpeff_peak) && "Pos" %in% names(snpeff_peak)) {
    sub <- snpeff_peak[as.character(snpeff_peak$Gene_ID) == as.character(gid), , drop = FALSE]
    if (nrow(sub)) {
      key_m <- paste(as.character(markers$Chr), as.numeric(markers$Pos), sep = "\t")
      sch <- if ("Chr" %in% names(sub)) as.character(sub$Chr) else rep(NA_character_, nrow(sub))
      key_s <- paste(sch, as.numeric(sub$Pos), sep = "\t")
      hit <- markers[key_m %in% unique(key_s), , drop = FALSE]
      if (nrow(hit)) {
        return(hit)
      }
    }
  }

  if (!is.null(gene_rows) && nrow(gene_rows)) {
    g_start <- if ("Gene_start" %in% names(gene_rows)) {
      suppressWarnings(as.numeric(gene_rows$Gene_start[1L]))
    } else {
      NA_real_
    }
    g_end <- if ("Gene_end" %in% names(gene_rows)) {
      suppressWarnings(as.numeric(gene_rows$Gene_end[1L]))
    } else {
      NA_real_
    }
    g_chr <- if ("Gene_chr" %in% names(gene_rows)) {
      as.character(gene_rows$Gene_chr[1L])
    } else if ("Chr" %in% names(gene_rows)) {
      as.character(gene_rows$Chr[1L])
    } else {
      NA_character_
    }
    if (is.finite(g_start) && is.finite(g_end)) {
      lo <- min(g_start, g_end)
      hi <- max(g_start, g_end)
      same_chr <- if (!is.na(g_chr) && nzchar(g_chr) && "Chr" %in% names(markers)) {
        as.character(markers$Chr) == g_chr
      } else {
        rep(TRUE, nrow(markers))
      }
      hit <- markers[same_chr & is.finite(markers$Pos) &
                       markers$Pos >= lo & markers$Pos <= hi, , drop = FALSE]
      if (nrow(hit)) {
        return(hit)
      }
    }

    # Synthetic marker from candidate gene row (no block marker attribution).
    d <- if ("dist2peak" %in% names(gene_rows)) {
      suppressWarnings(min(abs(as.numeric(gene_rows$dist2peak)), na.rm = TRUE))
    } else {
      NA_real_
    }
    nlp <- if ("negLog10P" %in% names(gene_rows)) {
      suppressWarnings(max(as.numeric(gene_rows$negLog10P), na.rm = TRUE))
    } else {
      NA_real_
    }
    d_max <- marker_scores$d_max
    lead_nlp <- marker_scores$lead_nlp
    dist_score <- if (!is.finite(d)) {
      0
    } else if (!is.finite(d_max) || d_max <= 0) {
      if (d <= 0) 1 else 0
    } else {
      max(0, min(1, (d_max - d) / d_max))
    }
    peak_score <- if (!is.finite(nlp) || !is.finite(lead_nlp) || lead_nlp <= 0) {
      0
    } else {
      max(0, min(1, nlp / lead_nlp))
    }
    return(data.frame(
      Chr = if (!is.na(g_chr)) g_chr else NA_character_,
      Pos = NA_real_,
      abs_dist = if (is.finite(d)) d else NA_real_,
      negLog10P = if (is.finite(nlp)) nlp else NA_real_,
      dist_score = dist_score,
      peak_score = peak_score,
      stringsAsFactors = FALSE
    ))
  }

  empty
}

#' Batch gwas_dist / gwas_peak evidence for gene IDs
#' @keywords internal
.evidence_gwas_channels_batch <- function(gene_ids,
                                          candidate = NULL,
                                          object = NULL,
                                          pheno_name = NULL) {
  empty_pair <- function() {
    list(
      gwas_dist = .evidence_gwas_empty("gwas_dist"),
      gwas_peak = .evidence_gwas_empty("gwas_peak")
    )
  }
  out <- setNames(lapply(gene_ids, function(gid) empty_pair()), gene_ids)
  if (!length(gene_ids)) {
    return(out)
  }

  blocks <- .gwas_load_peak_blocks(object, pheno_name)
  snpeff_all <- NULL
  if (!is.null(object) && !is.null(pheno_name)) {
    snpeff_all <- tryCatch(
      lazyData(object = object, dataset = "snpeff", pheno = pheno_name),
      error = function(e) NULL
    )
  }

  peak_ids <- if (!is.null(candidate) && "peak_ID" %in% names(candidate)) {
    unique(as.character(candidate$peak_ID))
  } else if (!is.null(blocks) && "peak_ID" %in% names(blocks)) {
    unique(as.character(blocks$peak_ID))
  } else {
    NA_character_
  }

  for (pid in peak_ids) {
    gene_in_peak <- gene_ids
    cand_peak <- candidate
    if (!is.null(candidate) && "peak_ID" %in% names(candidate) &&
        !is.na(pid) && nzchar(pid)) {
      cand_peak <- candidate[as.character(candidate$peak_ID) == pid, , drop = FALSE]
      gene_in_peak <- intersect(gene_ids, unique(as.character(cand_peak$Gene_ID)))
    }
    if (!length(gene_in_peak)) {
      next
    }

    blk <- blocks
    if (!is.null(blocks) && "peak_ID" %in% names(blocks) &&
        !is.na(pid) && nzchar(as.character(pid))) {
      blk <- blocks[as.character(blocks$peak_ID) == as.character(pid), , drop = FALSE]
    }

    # Fallback when no peak-block markers: synthesize from candidate rows.
    if (is.null(blk) || !nrow(blk)) {
      if (is.null(cand_peak) || !nrow(cand_peak)) {
        next
      }
      nlp_c <- if ("negLog10P" %in% names(cand_peak)) {
        as.numeric(cand_peak$negLog10P)
      } else {
        rep(NA_real_, nrow(cand_peak))
      }
      lead_i <- which.max(replace(nlp_c, !is.finite(nlp_c), -Inf))
      if (!length(lead_i) || !is.finite(lead_i) || lead_i < 1L) {
        lead_i <- 1L
      }
      blk <- data.frame(
        peak_ID = pid,
        variant_ID = seq_len(nrow(cand_peak)),
        peak_variant_ID = as.integer(lead_i),
        Chr = if ("Gene_chr" %in% names(cand_peak)) {
          as.character(cand_peak$Gene_chr)
        } else {
          NA_character_
        },
        Pos = if ("Gene_start" %in% names(cand_peak)) {
          as.numeric(cand_peak$Gene_start)
        } else {
          NA_real_
        },
        dist2peak = if ("dist2peak" %in% names(cand_peak)) {
          as.numeric(cand_peak$dist2peak)
        } else {
          NA_real_
        },
        negLog10P = nlp_c,
        stringsAsFactors = FALSE
      )
    }

    ms <- .gwas_marker_channel_scores(blk)
    snpeff_peak <- NULL
    if (!is.null(snpeff_all) && is.data.frame(snpeff_all) && nrow(snpeff_all)) {
      snpeff_peak <- .snpeff_for_peak(
        snpeff = snpeff_all,
        object = object,
        pheno_name = pheno_name,
        peak_id = pid,
        peak_markers = ms$markers[, c("Chr", "Pos"), drop = FALSE]
      )
    }

    for (gid in gene_in_peak) {
      grow <- if (!is.null(cand_peak) && nrow(cand_peak)) {
        cand_peak[as.character(cand_peak$Gene_ID) == as.character(gid), , drop = FALSE]
      } else {
        NULL
      }
      attributed <- .gwas_markers_for_gene(
        gid = gid,
        gene_rows = grow,
        marker_scores = ms,
        snpeff_peak = snpeff_peak
      )
      if (!nrow(attributed)) {
        next
      }
      dist_sc <- suppressWarnings(max(as.numeric(attributed$dist_score), na.rm = TRUE))
      peak_sc <- suppressWarnings(max(as.numeric(attributed$peak_score), na.rm = TRUE))
      if (!is.finite(dist_sc)) dist_sc <- 0
      if (!is.finite(peak_sc)) peak_sc <- 0
      best_d <- attributed[which.max(attributed$dist_score), , drop = FALSE]
      best_p <- attributed[which.max(attributed$peak_score), , drop = FALSE]
      out[[gid]] <- list(
        gwas_dist = list(
          source = "gwas_dist",
          score = dist_sc,
          snippets = sprintf(
            "gwas_dist = %.3f (abs_dist=%s; D_max=%s)",
            dist_sc,
            format(best_d$abs_dist[1L], scientific = FALSE),
            format(ms$d_max, scientific = FALSE)
          ),
          details = list(
            abs_dist = best_d$abs_dist[1L],
            d_max = ms$d_max,
            n_markers = nrow(attributed)
          )
        ),
        gwas_peak = list(
          source = "gwas_peak",
          score = peak_sc,
          snippets = sprintf(
            "gwas_peak = %.3f (negLog10P=%.2f; lead=%.2f)",
            peak_sc,
            best_p$negLog10P[1L],
            ms$lead_nlp
          ),
          details = list(
            negLog10P = best_p$negLog10P[1L],
            lead_negLog10P = ms$lead_nlp,
            n_markers = nrow(attributed)
          )
        )
      )
    }
  }
  out
}

#' Max PIP per gene: SnpEff assignment first, else GFF window (Phase 2 H1)
#'
#' Ensure credible-set rows have Chr/Pos for finemap PIP→gene mapping
#' @keywords internal
.e_enrich_cs_coords <- function(cred, object, pheno_name, peak_id) {
  if (is.null(cred) || !nrow(cred)) {
    return(cred)
  }
  if ("chr" %in% names(cred) && !"Chr" %in% names(cred)) {
    cred$Chr <- cred$chr
  }
  if ("pos" %in% names(cred) && !"Pos" %in% names(cred)) {
    cred$Pos <- cred$pos
  }
  have_chr <- "Chr" %in% names(cred) &&
    any(nzchar(as.character(cred$Chr)), na.rm = TRUE)
  have_pos <- "Pos" %in% names(cred) &&
    any(is.finite(as.numeric(cred$Pos)), na.rm = TRUE)
  if (have_chr && have_pos) {
    return(cred)
  }
  # Fall back: merge coordinates from peakcall / recalc blocks
  blk <- NULL
  for (ds in c("recalc", "peakcall")) {
    blk <- tryCatch(
      lazyData(object = object, dataset = ds, pheno = pheno_name),
      error = function(e) NULL
    )
    if (!is.null(blk) && nrow(blk)) {
      break
    }
  }
  if (is.null(blk) || !nrow(blk) || !"variant_ID" %in% names(blk)) {
    return(cred)
  }
  if (!is.null(peak_id) && "peak_ID" %in% names(blk)) {
    blk <- blk[as.character(blk$peak_ID) == as.character(peak_id), , drop = FALSE]
  }
  coord_cols <- intersect(c("variant_ID", "Chr", "chr", "Pos", "pos"), names(blk))
  if (!"variant_ID" %in% coord_cols) {
    return(cred)
  }
  coords <- unique(blk[, coord_cols, drop = FALSE])
  if ("chr" %in% names(coords) && !"Chr" %in% names(coords)) {
    coords$Chr <- coords$chr
  }
  if ("pos" %in% names(coords) && !"Pos" %in% names(coords)) {
    coords$Pos <- coords$pos
  }
  keep <- intersect(c("variant_ID", "Chr", "Pos"), names(coords))
  coords <- coords[, keep, drop = FALSE]
  merge(cred, coords, by = "variant_ID", all.x = TRUE, sort = FALSE)
}

#' For each peak's 95% CS: (1) CS ∩ SnpEff \code{Gene_ID} → \code{max(PIP)};
#' (2) genes still unassigned (no SnpEff or no hit) get
#' \code{max(PIP)} of CS variants whose position falls in the gene's GFF window.
#' @keywords internal
.finemap_pip_lookup <- function(object,
                                pheno_name,
                                candidate,
                                gff = NULL,
                                gff_windows = NULL) {
  empty <- data.frame(
    peak_ID = character(),
    Gene_ID = character(),
    max_PIP = numeric(),
    pip_source = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(candidate) || !"peak_ID" %in% names(candidate)) {
    stop(
      "Finemap scoring requires a candidate table with peak_ID.",
      call. = FALSE
    )
  }
  snpeff <- tryCatch(
    lazyData(object = object, dataset = "snpeff", pheno = pheno_name),
    error = function(e) NULL
  )
  have_snpeff <- !is.null(snpeff) && nrow(snpeff) &&
    "Gene_ID" %in% names(snpeff) && "Pos" %in% names(snpeff)
  if (have_snpeff) {
    snpeff$Gene_ID <- as.character(snpeff$Gene_ID)
    snpeff$Chr <- as.character(snpeff$Chr %||% snpeff$chr)
    snpeff$Pos <- as.numeric(snpeff$Pos)
  }

  win <- gff_windows
  if ((is.null(win) || !is.data.frame(win) || !nrow(win)) && !is.null(gff)) {
    win <- .gff_gene_windows(gff)
  }
  have_gff <- !is.null(win) && is.data.frame(win) && nrow(win) &&
    all(c("Gene_ID", "Chr", "window_start", "window_end") %in% names(win))
  if (have_gff) {
    win$Gene_ID <- as.character(win$Gene_ID)
    win$Chr <- as.character(win$Chr)
    win$window_start <- as.numeric(win$window_start)
    win$window_end <- as.numeric(win$window_end)
  }

  if (!have_snpeff && !have_gff) {
    # Neither assignment path available — score 0 (llm_report requires GFF)
    return(empty)
  }

  peak_ids <- unique(as.character(candidate$peak_ID))
  peak_ids <- peak_ids[!is.na(peak_ids) & nzchar(peak_ids)]
  rows <- list()
  for (pid in peak_ids) {
    cred <- tryCatch(
      lazyData(
        object = object,
        dataset = "credible_set",
        pheno = pheno_name,
        kind = paste0("peak_", pid)
      ),
      error = function(e) NULL
    )
    if (is.null(cred) || !nrow(cred) || !"PIP" %in% names(cred)) {
      stop(
        "Credible set required for peak ", pid,
        " (phenotype '", pheno_name, "'). Run calcCredibleSet() first.",
        call. = FALSE
      )
    }
    if (!"in_credible_set" %in% names(cred)) {
      stop("Credible set for peak ", pid, " lacks in_credible_set.", call. = FALSE)
    }
    cred_in <- as.logical(cred$in_credible_set)
    cred_in[is.na(cred_in)] <- FALSE
    cred <- cred[cred_in, , drop = FALSE]
    if (!nrow(cred)) {
      next
    }
    # Newer calcCredibleSet() stores Chr/Pos on the CS table. Merging those
    # rows again with the peak block yields Pos.x/Pos.y and a false "no Pos"
    # failure — enrich only when coordinates are missing (same as E path).
    cred2 <- .e_enrich_cs_coords(
      cred,
      object = object,
      pheno_name = pheno_name,
      peak_id = pid
    )
    if (!"Pos" %in% names(cred2) ||
        !any(is.finite(as.numeric(cred2$Pos)), na.rm = TRUE)) {
      stop(
        "Could not map credible-set variants to positions for peak ", pid, ".",
        call. = FALSE
      )
    }
    if (!"Chr" %in% names(cred2) && "chr" %in% names(cred2)) {
      cred2$Chr <- cred2$chr
    }
    cred2$Chr <- as.character(cred2$Chr)
    cred2$Pos <- as.numeric(cred2$Pos)
    cred2$PIP <- as.numeric(cred2$PIP)

    pip_by_gene <- list()
    src_by_gene <- list()

    # (1) SnpEff-first
    if (have_snpeff) {
      keep_cols <- intersect(
        c("variant_ID", "PIP", "Chr", "Pos"),
        names(cred2)
      )
      merged <- merge(
        cred2[, keep_cols, drop = FALSE],
        snpeff[, c("Gene_ID", "Chr", "Pos"), drop = FALSE],
        by = c("Chr", "Pos")
      )
      if (nrow(merged)) {
        merged$Gene_ID <- as.character(merged$Gene_ID)
        merged$PIP <- as.numeric(merged$PIP)
        agg <- tapply(merged$PIP, merged$Gene_ID, function(x) {
          m <- max(x, na.rm = TRUE)
          if (!is.finite(m)) NA_real_ else m
        })
        for (nm in names(agg)) {
          if (is.finite(agg[[nm]])) {
            pip_by_gene[[nm]] <- as.numeric(agg[[nm]])
            src_by_gene[[nm]] <- "snpeff"
          }
        }
      }
    }

    # (2) GFF window for genes still without a SnpEff PIP
    cand_peak <- candidate[
      as.character(candidate$peak_ID) == pid,
      ,
      drop = FALSE
    ]
    gene_ids <- unique(as.character(cand_peak$Gene_ID))
    gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
    if (have_gff) {
      need <- gene_ids
      if (length(pip_by_gene)) {
        need <- setdiff(need, names(pip_by_gene))
      }
      # Also score genes that appear only via windows (not in cand) when
      # SnpEff path produced nothing for this peak — prefer candidate set.
      if (!length(need) && !length(pip_by_gene)) {
        need <- unique(win$Gene_ID)
      }
      for (gid in need) {
        wrow <- win[win$Gene_ID == gid, , drop = FALSE]
        if (!nrow(wrow)) next
        # same Chr preferred; take first matching window
        hit <- rep(FALSE, nrow(cred2))
        for (wi in seq_len(nrow(wrow))) {
          hit <- hit |
            (cred2$Chr == wrow$Chr[wi] &
               is.finite(cred2$Pos) &
               cred2$Pos >= wrow$window_start[wi] &
               cred2$Pos <= wrow$window_end[wi])
        }
        if (!any(hit)) next
        mp <- max(cred2$PIP[hit], na.rm = TRUE)
        if (is.finite(mp)) {
          pip_by_gene[[gid]] <- mp
          src_by_gene[[gid]] <- "gff_window"
        }
      }
    }

    if (!length(pip_by_gene)) {
      next
    }
    gids <- names(pip_by_gene)
    rows[[length(rows) + 1L]] <- data.frame(
      peak_ID = pid,
      Gene_ID = gids,
      max_PIP = vapply(gids, function(g) pip_by_gene[[g]], numeric(1)),
      pip_source = vapply(gids, function(g) src_by_gene[[g]], character(1)),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(empty)
  }
  do.call(rbind, rows)
}

.evidence_finemap <- function(gene_id, gene_rows, lookup) {
  empty <- list(
    source = "finemap",
    score = 0,
    snippets = character(),
    details = list(max_PIP = NA_real_)
  )
  if (is.null(lookup) || !nrow(lookup)) {
    return(empty)
  }
  peak_id <- if (!is.null(gene_rows) && "peak_ID" %in% names(gene_rows) &&
                 nrow(gene_rows)) {
    as.character(gene_rows$peak_ID[1L])
  } else {
    NA_character_
  }
  hit <- lookup$Gene_ID == as.character(gene_id)
  if (!is.na(peak_id) && nzchar(peak_id)) {
    hit <- hit & as.character(lookup$peak_ID) == peak_id
  }
  if (!any(hit)) {
    return(empty)
  }
  pip <- max(as.numeric(lookup$max_PIP[hit]), na.rm = TRUE)
  if (!is.finite(pip)) {
    return(empty)
  }
  src <- if ("pip_source" %in% names(lookup)) {
    as.character(lookup$pip_source[hit][which.max(as.numeric(lookup$max_PIP[hit]))])
  } else {
    "snpeff"
  }
  snip <- if (identical(src, "gff_window")) {
    sprintf("max_PIP (95%% CS ∩ GFF window) = %.3f", pip)
  } else {
    sprintf("max_PIP (95%% CS ∩ SnpEff) = %.3f", pip)
  }
  list(
    source = "finemap",
    score = max(0, min(1, pip)),
    snippets = snip,
    details = list(max_PIP = pip, peak_ID = peak_id, pip_source = src)
  )
}

.evidence_expression <- function(gene_id,
                                 query,
                                 expression_matrix,
                                 expression_sample_map) {
  empty <- list(source = "expression", score = 0, snippets = character(), details = list())

  mat <- as.matrix(expression_matrix)
  if (is.null(rownames(mat))) {
    return(empty)
  }
  row_hit <- which(rownames(mat) == gene_id)
  if (!length(row_hit)) {
    return(empty)
  }

  expr <- mat[row_hit[1L], , drop = TRUE]
  expr <- as.numeric(expr)
  names(expr) <- colnames(mat)

  keywords <- if (inherits(query, "PhenotypeQuery")) {
    unique(c(tolower(query$tissues), tolower(query$developmental_stage)))
  } else {
    .phenotype_query_tokenize(as.character(query))
  }
  keywords <- keywords[nzchar(keywords)]

  matched_samples <- character()
  meta <- expression_sample_map
  meta$.__sample__ <- as.character(meta$sample)
  text_cols <- intersect(names(meta), c("tissue", "stage", "condition", "organ"))
  if (length(text_cols)) {
    meta_text <- apply(meta[, text_cols, drop = FALSE], 1L, function(row) {
      paste(stats::na.omit(as.character(row)), collapse = " ")
    })
    meta_text <- tolower(meta_text)
    for (kw in keywords) {
      hit <- grepl(kw, meta_text, fixed = TRUE)
      if (any(hit)) {
        matched_samples <- c(matched_samples, meta$.__sample__[hit])
      }
    }
  }

  matched_samples <- unique(matched_samples)
  score <- 0
  snippets <- character()

  if (length(matched_samples)) {
    vals <- expr[matched_samples]
    vals <- vals[is.finite(vals)]
    if (length(vals)) {
      z <- (mean(vals) - mean(expr, na.rm = TRUE)) / (stats::sd(expr, na.rm = TRUE) + 1e-8)
      score <- max(0, min(1, (z + 3) / 6))
      snippets <- c(
        snippets,
        sprintf(
          "expression in %d matching sample(s); mean = %.3f",
          length(vals),
          mean(vals)
        )
      )
    }
  }
  # No tissue/stage keywords → score stays 0 (do not invent a baseline).

  list(
    source = "expression",
    score = score,
    snippets = snippets,
    details = list(matched_samples = matched_samples)
  )
}

#' Coerce path or data.frame to from/to edge table
#' @keywords internal
.coerce_evidence_graph <- function(graph, arg_name = "graph") {
  if (is.null(graph)) {
    stop(
      "'", arg_name, "' is required when the matching evidence source is active ",
      "(edge table with columns 'from' and 'to').",
      call. = FALSE
    )
  }
  if (is.character(graph) && length(graph) == 1L) {
    if (!nzchar(graph) || !file.exists(graph)) {
      stop("'", arg_name, "' file not found: ", graph, call. = FALSE)
    }
    graph <- utils::read.csv(graph, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(graph)) {
    stop("'", arg_name, "' must be a data.frame or path to a CSV.", call. = FALSE)
  }
  nms <- tolower(names(graph))
  from_i <- match("from", nms)
  to_i <- match("to", nms)
  if (is.na(from_i) || is.na(to_i)) {
    stop(
      "'", arg_name, "' must contain columns 'from' and 'to'.",
      call. = FALSE
    )
  }
  edges <- data.frame(
    from = trimws(as.character(graph[[from_i]])),
    to = trimws(as.character(graph[[to_i]])),
    stringsAsFactors = FALSE
  )
  keep <- nzchar(edges$from) & nzchar(edges$to) &
    !is.na(edges$from) & !is.na(edges$to)
  edges <- edges[keep, , drop = FALSE]
  if (!nrow(edges)) {
    stop("'", arg_name, "' has no usable edges.", call. = FALSE)
  }
  edges
}

#' Compact fingerprint for cache keys
#' @keywords internal
.evidence_graph_fingerprint <- function(edges) {
  list(
    n = nrow(edges),
    hash = .evidence_cache_hash(list(
      from = edges$from,
      to = edges$to
    ))
  )
}

#' Build undirected adjacency list from from/to edges
#' @keywords internal
.evidence_graph_adjacency <- function(edges) {
  adj <- new.env(parent = emptyenv())
  add_edge <- function(a, b) {
    if (!exists(a, envir = adj, inherits = FALSE)) {
      assign(a, character(), envir = adj)
    }
    cur <- get(a, envir = adj, inherits = FALSE)
    if (!b %in% cur) {
      assign(a, c(cur, b), envir = adj)
    }
  }
  for (i in seq_len(nrow(edges))) {
    a <- edges$from[[i]]
    b <- edges$to[[i]]
    if (identical(a, b)) {
      next
    }
    add_edge(a, b)
    add_edge(b, a)
  }
  adj
}

#' BFS: any start reaches any target on undirected adjacency
#' @keywords internal
.evidence_graph_reachable <- function(adj, starts, targets) {
  starts <- unique(as.character(starts))
  starts <- starts[nzchar(starts) & !is.na(starts)]
  targets <- unique(as.character(targets))
  targets <- targets[nzchar(targets) & !is.na(targets)]
  if (!length(starts) || !length(targets)) {
    return(FALSE)
  }
  target_set <- new.env(parent = emptyenv())
  for (t in targets) {
    assign(t, TRUE, envir = target_set)
  }
  for (s in starts) {
    if (exists(s, envir = target_set, inherits = FALSE)) {
      return(TRUE)
    }
  }
  visited <- new.env(parent = emptyenv())
  queue <- starts
  qi <- 1L
  while (qi <= length(queue)) {
    node <- queue[[qi]]
    qi <- qi + 1L
    if (exists(node, envir = visited, inherits = FALSE)) {
      next
    }
    assign(node, TRUE, envir = visited)
    if (exists(node, envir = target_set, inherits = FALSE)) {
      return(TRUE)
    }
    if (!exists(node, envir = adj, inherits = FALSE)) {
      next
    }
    nbrs <- get(node, envir = adj, inherits = FALSE)
    for (nb in nbrs) {
      if (!exists(nb, envir = visited, inherits = FALSE)) {
        queue <- c(queue, nb)
      }
    }
  }
  FALSE
}

#' Gene-side node labels for network evidence
#' @keywords internal
.evidence_gene_graph_nodes <- function(gene_id, gene_rows = NULL) {
  nodes <- as.character(gene_id)[1L]
  if (!is.null(gene_rows) && nrow(gene_rows) > 0L && "Name" %in% names(gene_rows)) {
    nm <- as.character(gene_rows$Name)
    nm <- trimws(nm)
    nm <- nm[nzchar(nm) & !is.na(nm)]
    nodes <- c(nodes, nm)
  }
  unique(nodes)
}

#' Phenotype-side node labels for network evidence
#' @keywords internal
.evidence_phenotype_graph_nodes <- function(query) {
  if (is.character(query)) {
    return(unique(.phenotype_query_tokenize(query)))
  }
  if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  }
  parts <- c(
    query$trait_keywords %||% character(),
    query$tissues %||% character(),
    query$developmental_stage %||% character(),
    query$conditions %||% character(),
    query$action_phrases %||% character(),
    query$noun_phrases %||% character(),
    query$words %||% character()
  )
  parts <- trimws(as.character(unlist(parts, use.names = FALSE)))
  unique(parts[nzchar(parts) & !is.na(parts)])
}

#' Literature / ortholog network channel (path reachability → 1/0)
#' @keywords internal
.evidence_network_channel <- function(source,
                                      gene_id,
                                      gene_rows,
                                      query,
                                      edges) {
  empty <- list(
    source = source,
    score = 0,
    snippets = character(),
    details = list(reachable = FALSE, n_edges = nrow(edges))
  )
  gene_nodes <- .evidence_gene_graph_nodes(gene_id, gene_rows)
  pheno_nodes <- .evidence_phenotype_graph_nodes(query)
  if (!length(pheno_nodes)) {
    empty$details$note <- "no phenotype nodes in query"
    return(empty)
  }
  adj <- .evidence_graph_adjacency(edges)
  hit <- .evidence_graph_reachable(adj, gene_nodes, pheno_nodes)
  if (!isTRUE(hit)) {
    return(empty)
  }
  list(
    source = source,
    score = 1,
    snippets = sprintf(
      "path:%s↔{%s}",
      paste(gene_nodes, collapse = "|"),
      paste(utils::head(pheno_nodes, 5L), collapse = ",")
    ),
    details = list(
      reachable = TRUE,
      n_edges = nrow(edges),
      gene_nodes = gene_nodes,
      phenotype_nodes = pheno_nodes
    )
  )
}

.evidence_source_score <- function(evidence_list, source) {
  if (is.null(evidence_list[[source]])) {
    return(0)
  }
  sc <- evidence_list[[source]]$score
  if (!is.finite(sc)) {
    return(0)
  }
  max(0, min(1, sc))
}

.evidence_to_json <- function(evidence_list) {
  jsonlite::toJSON(evidence_list, auto_unbox = TRUE, null = "null")
}

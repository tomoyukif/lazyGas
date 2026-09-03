################################################################################
#' Collect multi-source evidence for candidate genes
#'
#' Gathers annotation, SnpEff, GWAS, expression (TENOR or matrix), literature,
#' and ortholog evidence for each gene in \code{gene_ids}, optionally using
#' cached results in the companion store.
#'
#' Phase 1 default sources are \code{annotation}, \code{snpeff}, \code{gwas},
#' and \code{expression}. Annotation scores combine phrase-aware keyword
#' matching and optional LLM relevance via \code{max} (LSA / text2vec removed).
#'
#' @param gene_ids Character vector of gene identifiers.
#' @param query A [PhenotypeQuery] object or character search string.
#' @param candidate Candidate-gene \code{data.frame} (optional if \code{object}
#'   and \code{pheno} are given).
#' @param object \code{LazyGas} object for cache I/O, candidate lookup, and
#'   TENOR cache path.
#' @param pheno Phenotype name or index when reading candidates from \code{object}.
#' @param sources Evidence sources to collect. Any of \code{"annotation"},
#'   \code{"snpeff"}, \code{"gwas"}, \code{"expression"}, \code{"literature"},
#'   \code{"ortholog"}, \code{"finemap"}.
#' @param expression_matrix Optional gene x sample expression matrix
#'   (\code{data.frame} or matrix) with row names = gene IDs. When \code{NULL}
#'   and \code{expression} is requested, [evaluateTenorExpression()] is used.
#' @param expression_meta Optional \code{data.frame} with one row per sample
#'   column in \code{expression_matrix}; should include \code{tissue} and/or
#'   \code{stage} columns when available.
#' @param ortholog_table Optional ortholog table with gene/ortholog IDs (and
#'   optionally scores) used as one evidence channel in ranking.
#' @param ortholog_cols Named character vector mapping column names
#'   (\code{gene_id}, \code{ortholog_id}, \code{score}).
#' @param use_cache Read/write evidence cache in the companion store.
#' @param literature_max Maximum PubMed hits per gene (when \pkg{rentrez} is
#'   available).
#' @param use_keyword,use_llm_relevance Annotation channel flags
#'   (defaults from \code{inst/config/ai-lazygas-models.yaml}).
#' @param use_semantic Ignored; LSA / semantic annotation scoring was removed.
#'   Retained for call compatibility.
#' @param download_tenor If \code{TRUE}, download missing TENOR CSVs when using
#'   the expression channel without \code{expression_matrix}.
#' @param llm_model,llm_base_url,llm_timeout LLM settings for annotation
#'   relevance scoring.
#'
#' @return A named list keyed by gene ID; each element is a list of evidence
#'   records (\code{source}, \code{score}, \code{snippets}, \code{details}).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [evaluateTenorExpression()]
collectGeneEvidence <- function(gene_ids,
                                query,
                                candidate = NULL,
                                object = NULL,
                                pheno = NULL,
                                sources = phase1DefaultSources(),
                                expression_matrix = NULL,
                                expression_meta = NULL,
                                ortholog_table = NULL,
                                ortholog_cols = c(
                                  gene_id = "gene_id",
                                  ortholog_id = "ortholog_id",
                                  score = "score"
                                ),
                                use_cache = TRUE,
                                literature_max = 5L,
                                use_keyword = NULL,
                                use_semantic = NULL,
                                use_llm_relevance = NULL,
                                download_tenor = TRUE,
                                llm_model = NULL,
                                llm_base_url = NULL,
                                llm_timeout = NULL) {
  sources <- match.arg(sources, choices = .EVIDENCE_SOURCE_CHOICES, several.ok = TRUE)
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids)]
  if (!length(gene_ids)) {
    stop("'gene_ids' must contain at least one gene ID.", call. = FALSE)
  }

  search_text <- .phenotype_query_as_search_text(query)
  flags <- .ai_config_ranking_flags()
  if (is.null(use_keyword)) use_keyword <- flags$use_keyword
  if (is.null(use_llm_relevance)) use_llm_relevance <- flags$use_llm_relevance
  if (isTRUE(use_semantic)) {
    warning(
      "'use_semantic' is ignored; LSA / text2vec annotation scoring was removed.",
      call. = FALSE
    )
  }
  use_semantic <- FALSE
  llm_settings <- .ai_config_llm_settings(
    model = llm_model %||% flags$llm_model,
    base_url = llm_base_url,
    timeout = llm_timeout
  )

  if (is.null(candidate) && !is.null(object)) {
    pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  }
  if (!is.null(candidate) && !"Gene_ID" %in% names(candidate)) {
    stop("'candidate' must contain a Gene_ID column.", call. = FALSE)
  }

  cache_key <- .evidence_cache_hash(list(
    query = if (inherits(query, "PhenotypeQuery")) {
      unclass(query)
    } else {
      query
    },
    sources = sources,
    use_keyword = use_keyword,
    use_semantic = FALSE,
    use_llm_relevance = use_llm_relevance,
    keyword_engine = "phrase_boundary_v1"
  ))

  ann_batch <- NULL
  if ("annotation" %in% sources) {
    ann_batch <- .evidence_annotation_batch(
      gene_ids = gene_ids,
      candidate = candidate,
      search_text = search_text,
      query = query,
      use_keyword = use_keyword,
      use_llm_relevance = use_llm_relevance,
      llm_model = llm_settings$model,
      llm_base_url = llm_settings$base_url,
      llm_timeout = llm_settings$timeout
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
      candidate = candidate
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
    if ("gwas" %in% sources) {
      ev$gwas <- .evidence_gwas(gene_id = gid, gene_rows = gene_rows)
    }
    if ("expression" %in% sources) {
      ev$expression <- .evidence_expression(
        gene_id = gid,
        query = query,
        expression_matrix = expression_matrix,
        expression_meta = expression_meta,
        object = object,
        download_tenor = download_tenor
      )
    }
    if ("literature" %in% sources) {
      ev$literature <- .evidence_literature(
        gene_id = gid,
        gene_rows = gene_rows,
        query = query,
        max_results = literature_max
      )
    }
    if ("ortholog" %in% sources) {
      ev$ortholog <- .evidence_ortholog(
        gene_id = gid,
        ortholog_table = ortholog_table,
        ortholog_cols = ortholog_cols
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

.normalize_score_vec <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (!any(is.finite(x))) {
    return(rep(0, length(x)))
  }
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) < .Machine$double.eps) {
    return(ifelse(is.finite(x), 1, 0))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

.evidence_annotation_batch <- function(gene_ids,
                                       candidate,
                                       search_text,
                                       query,
                                       use_keyword = TRUE,
                                       use_llm_relevance = FALSE,
                                       llm_model = NULL,
                                       llm_base_url = NULL,
                                       llm_timeout = 600) {
  empty <- function(gid) {
    list(
      source = "annotation",
      score = 0,
      snippets = character(),
      details = list(
        keyword_score = 0,
        matched_keywords = character(),
        unmatched_keywords = character(),
        llm_relevance_score = 0,
        match_method = "none"
      )
    )
  }

  out <- setNames(lapply(gene_ids, empty), gene_ids)
  if (is.null(candidate) || !nrow(candidate)) {
    return(out)
  }

  ann_cols <- tryCatch(
    .resolve_ann_cols(candidate = candidate, ann_cols = NULL),
    error = function(e) character()
  )
  if (!length(ann_cols)) {
    return(out)
  }

  # One annotation text per gene (best row later for snippets)
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

  kw_terms <- .keyword_terms_from_query(query)
  if (!length(kw_terms) && nzchar(search_text %||% "")) {
    kw_terms <- .keyword_terms_from_query(search_text)
  }

  kw_details <- vector("list", length(gene_ids))
  names(kw_details) <- gene_ids
  if (isTRUE(use_keyword)) {
    details_list <- .keyword_match_details(
      text = unname(texts),
      terms = kw_terms,
      ignore.case = TRUE
    )
    kw_details <- setNames(details_list, gene_ids)
  } else {
    kw_details <- setNames(
      lapply(gene_ids, function(gid) {
        list(
          keyword_score = 0,
          matched_keywords = character(),
          unmatched_keywords = kw_terms
        )
      }),
      gene_ids
    )
  }

  llm_score <- rep(0, length(gene_ids))
  names(llm_score) <- gene_ids
  if (isTRUE(use_llm_relevance) && any(nzchar(texts))) {
    llm_score <- .llm_annotation_relevance_scores(
      gene_ids = gene_ids,
      texts = texts,
      query = query,
      search_text = search_text,
      model = llm_model,
      base_url = llm_base_url,
      timeout = llm_timeout
    )
  }

  for (gid in gene_ids) {
    det <- kw_details[[gid]] %||% list(
      keyword_score = 0,
      matched_keywords = character(),
      unmatched_keywords = character()
    )
    kw <- as.numeric(det$keyword_score %||% 0)
    llm <- as.numeric(llm_score[[gid]] %||% 0)
    if (!is.finite(kw)) kw <- 0
    if (!is.finite(llm)) llm <- 0
    score <- max(kw, llm)
    methods <- c(
      if (kw > 0) "keyword",
      if (llm > 0) "llm"
    )
    out[[gid]] <- list(
      source = "annotation",
      score = score,
      snippets = snippets[[gid]] %||% character(),
      details = list(
        keyword_score = kw,
        matched_keywords = as.character(det$matched_keywords %||% character()),
        unmatched_keywords = as.character(det$unmatched_keywords %||% character()),
        llm_relevance_score = llm,
        match_method = if (length(methods)) paste(methods, collapse = "+") else "none"
      )
    )
  }
  out
}

.llm_annotation_relevance_scores <- function(gene_ids,
                                             texts,
                                             query,
                                             search_text,
                                             model = NULL,
                                             base_url = NULL,
                                             timeout = 600) {
  scores <- rep(0, length(gene_ids))
  names(scores) <- gene_ids

  if (!llmHealthCheck(base_url = base_url, timeout = 5)) {
    warning(
      "LLM unreachable; llm_relevance_score set to 0 for annotation channel.",
      call. = FALSE
    )
    return(scores)
  }

  payload <- lapply(seq_along(gene_ids), function(i) {
    list(
      gene_id = gene_ids[[i]],
      annotation = texts[[i]]
    )
  })
  # Cap payload size for long candidate lists
  if (length(payload) > 80L) {
    payload <- payload[seq_len(80L)]
  }

  trait <- if (inherits(query, "PhenotypeQuery")) {
    query$trait_text
  } else {
    search_text
  }

  system_msg <- paste(
    "You score how relevant each gene annotation is to a phenotype query.",
    "Return ONLY JSON: {\"scores\":[{\"gene_id\":\"...\",\"score\":0.0}, ...]}",
    "score must be a number from 0 to 1.",
    "Use 1 for clear causal/functional match, 0.5 for plausible, 0 for unrelated."
  )
  user_msg <- jsonlite::toJSON(
    list(phenotype_query = trait, genes = payload),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  raw <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = system_msg),
        list(role = "user", content = user_msg)
      ),
      model = model,
      base_url = base_url,
      json_mode = TRUE,
      timeout = timeout
    ),
    error = function(e) {
      warning("LLM annotation relevance failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (is.null(raw)) {
    return(scores)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(raw, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed)) {
    return(scores)
  }
  items <- parsed$scores %||% parsed
  if (!is.list(items)) {
    return(scores)
  }
  for (item in items) {
    gid <- as.character(item$gene_id %||% NA_character_)
    sc <- suppressWarnings(as.numeric(item$score %||% NA_real_))
    if (!is.na(gid) && gid %in% names(scores) && is.finite(sc)) {
      scores[[gid]] <- max(0, min(1, sc))
    }
  }
  scores
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
    if (col %in% names(gene_rows)) {
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

.evidence_gwas <- function(gene_id, gene_rows) {
  if (is.null(gene_rows) || nrow(gene_rows) == 0L) {
    return(list(source = "gwas", score = 0, snippets = character(), details = list()))
  }

  score <- 0
  snippets <- character()
  details <- list()

  if ("negLog10P" %in% names(gene_rows)) {
    pvals <- as.numeric(gene_rows$negLog10P)
    score <- max(pvals, na.rm = TRUE)
    if (!is.finite(score)) {
      score <- 0
    }
    details$negLog10P <- max(pvals, na.rm = TRUE)
    snippets <- c(snippets, sprintf("negLog10P = %.2f", details$negLog10P))
  }
  if ("dist2peak" %in% names(gene_rows)) {
    d <- min(as.numeric(gene_rows$dist2peak), na.rm = TRUE)
    if (is.finite(d)) {
      details$dist2peak <- d
      snippets <- c(snippets, sprintf("dist2peak = %s bp", d))
      proximity <- 1 / (1 + d / 1e5)
      score <- score + proximity
    }
  }
  list(
    source = "gwas",
    score = score,
    snippets = unique(snippets),
    details = details
  )
}

#' Max PIP per gene from 95% credible set ∩ SnpEff Gene_ID
#' @keywords internal
.finemap_pip_lookup <- function(object, pheno_name, candidate) {
  if (is.null(candidate) || !"peak_ID" %in% names(candidate)) {
    stop(
      "Finemap scoring requires a candidate table with peak_ID.",
      call. = FALSE
    )
  }
  snpeff <- lazyData(object = object, dataset = "snpeff", pheno = pheno_name)
  if (is.null(snpeff) || !nrow(snpeff) || !"Gene_ID" %in% names(snpeff) ||
      !"Pos" %in% names(snpeff)) {
    stop(
      "SnpEff annotations are required for finemap scoring. ",
      "Run listCandidate(..., snpeff = ...) first.",
      call. = FALSE
    )
  }
  snpeff$Gene_ID <- as.character(snpeff$Gene_ID)
  snpeff$Chr <- as.character(snpeff$Chr)
  snpeff$Pos <- as.numeric(snpeff$Pos)

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
    block <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = TRUE)
    if (is.null(block) || !nrow(block)) {
      block <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = FALSE)
    }
    block <- block[as.character(block$peak_ID) == pid, , drop = FALSE]
    if (!nrow(block) || !"variant_ID" %in% names(block)) {
      stop("Peak block missing for peak ", pid, ".", call. = FALSE)
    }
    pos_map <- block[, intersect(c("variant_ID", "Chr", "Pos"), names(block)),
                     drop = FALSE]
    cred2 <- merge(cred, pos_map, by = "variant_ID", all.x = TRUE)
    if (!"Pos" %in% names(cred2) || all(is.na(cred2$Pos))) {
      stop(
        "Could not map credible-set variants to positions for peak ", pid, ".",
        call. = FALSE
      )
    }
    cred2$Chr <- as.character(cred2$Chr)
    cred2$Pos <- as.numeric(cred2$Pos)
    merged <- merge(
      cred2[, c("variant_ID", "PIP", "Chr", "Pos"), drop = FALSE],
      snpeff[, c("Gene_ID", "Chr", "Pos"), drop = FALSE],
      by = c("Chr", "Pos")
    )
    if (!nrow(merged)) {
      next
    }
    merged$Gene_ID <- as.character(merged$Gene_ID)
    merged$PIP <- as.numeric(merged$PIP)
    agg <- tapply(merged$PIP, merged$Gene_ID, function(x) {
      m <- max(x, na.rm = TRUE)
      if (!is.finite(m)) NA_real_ else m
    })
    rows[[length(rows) + 1L]] <- data.frame(
      peak_ID = pid,
      Gene_ID = names(agg),
      max_PIP = as.numeric(agg),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(
      peak_ID = character(),
      Gene_ID = character(),
      max_PIP = numeric(),
      stringsAsFactors = FALSE
    ))
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
  list(
    source = "finemap",
    score = max(0, min(1, pip)),
    snippets = sprintf("max_PIP (95%% CS ∩ SnpEff) = %.3f", pip),
    details = list(max_PIP = pip, peak_ID = peak_id)
  )
}

.evidence_expression <- function(gene_id,
                                 query,
                                 expression_matrix = NULL,
                                 expression_meta = NULL,
                                 object = NULL,
                                 download_tenor = TRUE) {
  empty <- list(source = "expression", score = 0, snippets = character(), details = list())

  # Prefer explicit matrix when provided; otherwise TENOR (Phase 1)
  if (is.null(expression_matrix)) {
    return(evaluateTenorExpression(
      gene_id = gene_id,
      query = query,
      object = object,
      download_if_missing = isTRUE(download_tenor)
    ))
  }

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
  if (!is.null(expression_meta) && nrow(expression_meta) > 0L) {
    meta <- expression_meta
    if (!is.null(rownames(meta))) {
      meta$.__sample__ <- rownames(meta)
    } else if ("sample" %in% names(meta)) {
      meta$.__sample__ <- meta$sample
    } else {
      meta$.__sample__ <- rownames(meta)
    }
    text_cols <- intersect(names(meta), c("tissue", "stage", "condition", "organ"))
    if (length(text_cols)) {
      meta_text <- apply(meta[, text_cols, drop = FALSE], 1L, function(row) {
        paste(na.omit(as.character(row)), collapse = " ")
      })
      meta_text <- tolower(meta_text)
      for (kw in keywords) {
        hit <- grepl(kw, meta_text, fixed = TRUE)
        if (any(hit)) {
          matched_samples <- c(matched_samples, meta$.__sample__[hit])
        }
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
  } else if (length(keywords) == 0L) {
    score <- 0.25
    snippets <- "expression data available (no tissue/stage filter)"
  }

  list(
    source = "expression",
    score = score,
    snippets = snippets,
    details = list(matched_samples = matched_samples)
  )
}

.evidence_literature <- function(gene_id,
                                 gene_rows,
                                 query,
                                 max_results = 5L) {
  empty <- list(source = "literature", score = 0, snippets = character(), details = list())
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    empty$details$note <- "rentrez not installed"
    return(empty)
  }

  symbol <- gene_id
  if (!is.null(gene_rows) && nrow(gene_rows) > 0L) {
    if ("Name" %in% names(gene_rows) && nzchar(gene_rows$Name[1L])) {
      symbol <- gene_rows$Name[1L]
    }
  }

  terms <- if (inherits(query, "PhenotypeQuery")) {
    unique(c(query$trait_keywords, query$tissues, query$developmental_stage))
  } else {
    .phenotype_query_tokenize(as.character(query))
  }
  terms <- terms[nzchar(terms)]
  if (!length(terms)) {
    terms <- "gene"
  }

  q <- paste0(
    "(", symbol, "[Gene Name] OR ", symbol, "[Title/Abstract]) AND (",
    paste(sprintf("%s[Title/Abstract]", terms), collapse = " OR "),
    ")"
  )

  hits <- tryCatch({
    s <- rentrez::entrez_search(db = "pubmed", term = q, retmax = max_results)
    if (s$count == 0L) {
      return(empty)
    }
    summ <- rentrez::entrez_summary(db = "pubmed", id = s$ids)
    if (length(s$ids) == 1L) {
      summ <- list(summ)
    }
    titles <- vapply(summ, function(x) as.character(x$title), character(1L))
    pmids <- vapply(summ, function(x) as.character(x$uid), character(1L))
    list(
      source = "literature",
      score = min(1, length(pmids) / max_results),
      snippets = titles,
      details = list(pmid = pmids, query = q)
    )
  }, error = function(e) {
    list(
      source = "literature",
      score = 0,
      snippets = character(),
      details = list(error = conditionMessage(e))
    )
  })

  hits
}

.evidence_ortholog <- function(gene_id, ortholog_table, ortholog_cols) {
  empty <- list(source = "ortholog", score = 0, snippets = character(), details = list())
  if (is.null(ortholog_table) || !nrow(ortholog_table)) {
    return(empty)
  }

  gene_col <- ortholog_cols[["gene_id"]]
  ortho_col <- ortholog_cols[["ortholog_id"]]
  score_col <- ortholog_cols[["score"]]
  if (!gene_col %in% names(ortholog_table)) {
    return(empty)
  }

  sub <- ortholog_table[ortholog_table[[gene_col]] == gene_id, , drop = FALSE]
  if (!nrow(sub)) {
    return(empty)
  }

  snippets <- as.character(sub[[ortho_col]])
  snippets <- snippets[nzchar(snippets)]
  score <- 0.5
  if (score_col %in% names(sub)) {
    sc <- suppressWarnings(max(as.numeric(sub[[score_col]]), na.rm = TRUE))
    if (is.finite(sc)) {
      score <- max(score, min(1, sc))
    }
  }

  list(
    source = "ortholog",
    score = score,
    snippets = unique(snippets),
    details = list(n_orthologs = length(snippets))
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

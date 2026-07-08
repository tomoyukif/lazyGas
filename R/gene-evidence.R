################################################################################
#' Collect multi-source evidence for candidate genes
#'
#' Gathers annotation, GWAS, expression, literature, and ortholog evidence for
#' each gene in \code{gene_ids}, optionally using cached results in the
#' companion store.
#'
#' @param gene_ids Character vector of gene identifiers.
#' @param query A [PhenotypeQuery] object or character search string.
#' @param candidate Candidate-gene \code{data.frame} (optional if \code{object}
#'   and \code{pheno} are given).
#' @param object \code{LazyGas} object for cache I/O and candidate lookup.
#' @param pheno Phenotype name or index when reading candidates from \code{object}.
#' @param sources Evidence sources to collect. Any of \code{"annotation"},
#'   \code{"gwas"}, \code{"expression"}, \code{"literature"}, \code{"ortholog"}.
#' @param expression_matrix Optional gene x sample expression matrix
#'   (\code{data.frame} or matrix) with row names = gene IDs.
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
#'
#' @return A named list keyed by gene ID; each element is a list of evidence
#'   records (\code{source}, \code{score}, \code{snippets}, \code{details}).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()]
collectGeneEvidence <- function(gene_ids,
                                query,
                                candidate = NULL,
                                object = NULL,
                                pheno = NULL,
                                sources = c("annotation", "gwas", "expression",
                                            "literature", "ortholog"),
                                expression_matrix = NULL,
                                expression_meta = NULL,
                                ortholog_table = NULL,
                                ortholog_cols = c(
                                  gene_id = "gene_id",
                                  ortholog_id = "ortholog_id",
                                  score = "score"
                                ),
                                use_cache = TRUE,
                                literature_max = 5L) {
  sources <- match.arg(sources, several.ok = TRUE)
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids)]
  if (!length(gene_ids)) {
    stop("'gene_ids' must contain at least one gene ID.", call. = FALSE)
  }

  search_text <- .phenotype_query_as_search_text(query)

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
    sources = sources
  ))

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
      ev$annotation <- .evidence_annotation(
        gene_id = gid,
        gene_rows = gene_rows,
        search_text = search_text
      )
    }
    if ("gwas" %in% sources) {
      ev$gwas <- .evidence_gwas(gene_id = gid, gene_rows = gene_rows)
    }
    if ("expression" %in% sources) {
      ev$expression <- .evidence_expression(
        gene_id = gid,
        query = query,
        expression_matrix = expression_matrix,
        expression_meta = expression_meta
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

.evidence_annotation <- function(gene_id, gene_rows, search_text) {
  if (is.null(gene_rows) || nrow(gene_rows) == 0L) {
    return(list(source = "annotation", score = 0, snippets = character(), details = list()))
  }

  ann_cols <- .resolve_ann_cols(candidate = gene_rows, ann_cols = NULL)
  ann_text <- .candidate_annotation_text(candidate = gene_rows, ann_cols = ann_cols)
  kw <- .keyword_match_score(text = ann_text, query = search_text, match = "any")
  score <- if (length(kw)) max(as.numeric(kw), na.rm = TRUE) else 0
  if (!is.finite(score)) {
    score <- 0
  }
  snippets <- ann_text[nzchar(ann_text)]

  list(
    source = "annotation",
    score = score,
    snippets = snippets,
    details = list(match_method = if (score > 0) "keyword" else "none")
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
  if ("PIP" %in% names(gene_rows)) {
    pip <- max(as.numeric(gene_rows$PIP), na.rm = TRUE)
    if (is.finite(pip)) {
      details$PIP <- pip
      snippets <- c(snippets, sprintf("PIP = %.3f", pip))
      score <- score + pip * 5
    }
  }
  for (imp_col in c("HIGH", "MODERATE")) {
    if (imp_col %in% names(gene_rows)) {
      v <- sum(as.numeric(gene_rows[[imp_col]]), na.rm = TRUE)
      if (is.finite(v) && v > 0) {
        details[[imp_col]] <- v
        snippets <- c(snippets, sprintf("%s impact variants: %s", imp_col, v))
        score <- score + 0.1 * v
      }
    }
  }

  list(
    source = "gwas",
    score = score,
    snippets = unique(snippets),
    details = details
  )
}

.evidence_expression <- function(gene_id,
                                 query,
                                 expression_matrix = NULL,
                                 expression_meta = NULL) {
  empty <- list(source = "expression", score = 0, snippets = character(), details = list())
  if (is.null(expression_matrix)) {
    return(empty)
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

################################################################################
#' Search candidate genes by functional annotation
#'
#' Filter and rank genes in a candidate list using phrase-aware keyword
#' matching (stopwords ignored; word-boundary / phrase match).
#'
#' @param candidate A data.frame of candidate genes (e.g. from
#'   [lazyData()] with `dataset = "candidate"`). Used when `object` is
#'   `NULL`.
#' @param object A \code{LazyGas} object. When provided, the candidate list is
#'   read with [lazyData()] (`dataset = "candidate"`). `candidate` is
#'   ignored in that case.
#' @param pheno Phenotype name or index passed to [lazyData()] when
#'   `object` is given. Required if `object` is not `NULL`.
#' @param query Character string describing the phenotype or biological process
#'   of interest (space-separated keywords are allowed).
#' @param ann_cols Character vector of column names in `candidate` that hold
#'   functional annotation text (e.g. GO terms, gene descriptions). If `NULL`
#'   (default), all columns except peak/gene coordinates and SnpEff impact
#'   counts are used (i.e. columns joined from `ann` in [listCandidate()]).
#' @param mode Search mode. Only \code{"keyword"} is supported (LSA / semantic
#'   modes were removed).
#' @param keyword_match For keyword mode, `"all"` requires every content term to
#'   appear in the annotation text; `"any"` requires at least one term.
#' @param ignore.case Passed to keyword matching.
#' @param top_n If not `NULL`, return at most this many rows after ranking by
#'   `match_score`.
#' @param dedupe_genes If `TRUE` and `Gene_ID` is present, keep one row per
#'   gene with the highest `match_score`.
#' @param min_score,n_topics Ignored; retained for call compatibility with
#'   older scripts (semantic / LSA search was removed).
#'
#' @return A subset of the input candidate `data.frame` with columns
#'   `keyword_score`, `match_score`, and `match_method` added, sorted by
#'   decreasing `match_score`.
#'
#' @details
#' Annotation text for each row is the space-separated concatenation of
#' non-empty values in `ann_cols`. Query terms drop English stopwords
#' (e.g. \code{in}, \code{of}) and match with word boundaries (phrases use
#' boundary-aware multi-word patterns).
#'
#' @seealso [lazyData()], [listCandidate()]
#'
#' @export
#'
#' @examples
#' cand <- data.frame(
#'   peak_ID = 1L,
#'   Gene_ID = c("g1", "g2", "g3"),
#'   Gene_chr = "1",
#'   Gene_start = 1:3,
#'   dist2peak = 0,
#'   negLog10P = 3,
#'   Description = c(
#'     "fruit weight development",
#'     "root hair elongation",
#'     "cell wall biosynthesis"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#' searchCandidateGenes(
#'   candidate = cand,
#'   query = "fruit weight",
#'   mode = "keyword",
#'   keyword_match = "all"
#' )
searchCandidateGenes <- function(candidate = NULL,
                                 object = NULL,
                                 pheno = NULL,
                                 query,
                                 ann_cols = NULL,
                                 mode = "keyword",
                                 keyword_match = c("all", "any"),
                                 ignore.case = TRUE,
                                 min_score = 0.1,
                                 top_n = NULL,
                                 dedupe_genes = TRUE,
                                 n_topics = NULL) {
  keyword_match <- match.arg(keyword_match)
  mode <- as.character(mode)[1L]
  if (!identical(mode, "keyword")) {
    stop(
      "searchCandidateGenes() supports mode = \"keyword\" only. ",
      "LSA / semantic search (text2vec) was removed.",
      call. = FALSE
    )
  }
  if (!missing(n_topics) && !is.null(n_topics)) {
    warning("'n_topics' is ignored; semantic / LSA search was removed.",
            call. = FALSE)
  }

  if (!is.null(object)) {
    if (!inherits(object, "LazyGas")) {
      stop("'object' must be a LazyGas object.", call. = FALSE)
    }
    if (is.null(pheno)) {
      stop("'pheno' is required when 'object' is provided.", call. = FALSE)
    }
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno)
    if (is.null(candidate) || nrow(candidate) == 0L) {
      stop("No candidate data found for the given phenotype.", call. = FALSE)
    }
  } else if (is.null(candidate)) {
    stop("Provide either 'object' (LazyGas) or 'candidate' (data.frame).",
         call. = FALSE)
  }

  if (!is.data.frame(candidate)) {
    stop("'candidate' must be a data.frame.", call. = FALSE)
  }
  if (nrow(candidate) == 0L) {
    stop("'candidate' has no rows.", call. = FALSE)
  }
  if (!is.character(query) || length(query) != 1L || !nzchar(trimws(query))) {
    stop("'query' must be a non-empty character string.", call. = FALSE)
  }
  query <- trimws(query)

  ann_cols <- .resolve_ann_cols(candidate = candidate, ann_cols = ann_cols)
  ann_text <- .candidate_annotation_text(candidate = candidate, ann_cols = ann_cols)

  keyword_score <- .keyword_match_score(
    text = ann_text,
    query = query,
    match = keyword_match,
    ignore.case = ignore.case
  )

  kw_hit <- if (keyword_match == "all") {
    keyword_score >= 1
  } else {
    keyword_score > 0
  }

  if (!any(kw_hit)) {
    out <- candidate[0, , drop = FALSE]
    attr(out, "searchCandidateGenes") <- list(
      query = query,
      mode = mode,
      ann_cols = ann_cols,
      n_matched = 0L
    )
    return(out)
  }

  out <- candidate[kw_hit, , drop = FALSE]
  out$keyword_score <- keyword_score[kw_hit]
  out$match_score <- out$keyword_score
  out$match_method <- ifelse(out$keyword_score > 0, "keyword", "none")

  if (dedupe_genes && "Gene_ID" %in% names(out)) {
    out <- .dedupe_candidate_by_gene(out)
  }

  out <- out[order(-out$match_score, out$Gene_ID), , drop = FALSE]
  rownames(out) <- NULL

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)[1L]
    if (top_n > 0L && nrow(out) > top_n) {
      out <- out[seq_len(top_n), , drop = FALSE]
    }
  }

  attr(out, "searchCandidateGenes") <- list(
    query = query,
    mode = mode,
    ann_cols = ann_cols,
    keyword_match = keyword_match,
    n_matched = nrow(out)
  )
  out
}

.candidate_non_ann_cols <- function() {
  c(
    "peak_ID", "Gene_ID", "Gene_chr", "Gene_start", "dist2peak", "negLog10P",
    "HIGH", "MODERATE", "LOW", "MODIFIER",
    "HIGH_at_var", "MODERATE_at_var", "LOW_at_var", "MODIFIER_at_var"
  )
}

.resolve_ann_cols <- function(candidate, ann_cols) {
  if (!is.null(ann_cols) && length(ann_cols)) {
    if (!is.character(ann_cols)) {
      stop("'ann_cols' must be a character vector.", call. = FALSE)
    }
    miss <- setdiff(ann_cols, names(candidate))
    if (length(miss)) {
      stop(
        "Column(s) not found in candidate: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    return(ann_cols)
  }

  cols <- setdiff(names(candidate), .candidate_non_ann_cols())
  if (!length(cols)) {
    stop(
      "No annotation columns found in candidate. ",
      "Run listCandidate() with ann = ..., or specify ann_cols explicitly.",
      call. = FALSE
    )
  }
  cols
}

.candidate_annotation_text <- function(candidate, ann_cols) {
  mat <- candidate[, ann_cols, drop = FALSE]
  apply(mat, 1L, function(row) {
    x <- as.character(row)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (!length(x)) {
      return("")
    }
    paste(x, collapse = " ")
  })
}

#' English function words ignored for annotation keyword matching
#' @keywords internal
.annotation_keyword_stopwords <- function() {
  c(
    "a", "an", "the", "and", "or", "but", "if", "nor", "so", "as", "at", "by",
    "for", "from", "in", "into", "of", "off", "on", "onto", "over", "to", "up",
    "with", "without", "within", "among", "between", "via", "vs", "per", "than",
    "then", "that", "this", "these", "those", "it", "its", "is", "are", "was",
    "were", "be", "been", "being", "not", "no", "also", "such", "can", "will",
    "just", "only", "same", "other", "some", "any", "all", "each", "few", "more",
    "most", "very", "too", "about", "above", "below", "under", "again",
    "further", "once", "here", "there", "when", "where", "why", "how", "out",
    "down", "own", "should", "now", "do", "does", "did", "have", "has", "had"
  )
}

.is_content_keyword_term <- function(term) {
  term <- tolower(trimws(as.character(term)[1L]))
  if (!nzchar(term) || nchar(term) < 2L) {
    return(FALSE)
  }
  if (!grepl("[a-z0-9]", term, perl = TRUE, ignore.case = TRUE)) {
    return(FALSE)
  }
  !(term %in% .annotation_keyword_stopwords())
}

.annotation_tokenize <- function(x) {
  x <- tolower(trimws(as.character(x)[1L]))
  if (!nzchar(x)) {
    return(character())
  }
  terms <- strsplit(x, "[^a-zA-Z0-9]+", perl = TRUE)[[1L]]
  terms <- terms[nzchar(terms) & nchar(terms) >= 2L]
  unique(terms)
}

.filter_keyword_terms <- function(terms) {
  terms <- unique(trimws(as.character(unlist(terms, use.names = FALSE))))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (!length(terms)) {
    return(character())
  }
  keep <- vapply(
    terms,
    function(term) {
      parts <- .annotation_tokenize(term)
      if (!length(parts)) {
        return(.is_content_keyword_term(term))
      }
      any(vapply(parts, .is_content_keyword_term, logical(1L)))
    },
    logical(1L)
  )
  terms[keep]
}

#' Build content keyword / phrase list from a phenotype query
#' @keywords internal
.keyword_terms_from_query <- function(query) {
  if (inherits(query, "PhenotypeQuery")) {
    raw <- c(
      query$trait_keywords,
      query$tissues,
      query$developmental_stage,
      query$conditions
    )
    terms <- .filter_keyword_terms(raw)
    if (!length(terms)) {
      trait <- if (!is.null(query$trait_text)) query$trait_text else ""
      if (nzchar(trimws(as.character(trait)[1L]))) {
        terms <- .filter_keyword_terms(.annotation_tokenize(trait))
      }
    }
    return(terms)
  }
  q <- trimws(as.character(query)[1L])
  if (!nzchar(q)) {
    return(character())
  }
  .filter_keyword_terms(.annotation_tokenize(q))
}

.escape_perl_regex <- function(x) {
  gsub("([.|()\\[\\]{}+*?^$\\\\])", "\\\\\\1", x, perl = TRUE)
}

#' Word-boundary / phrase match for one keyword term
#' @keywords internal
.keyword_term_matches <- function(text, term, ignore.case = TRUE) {
  term <- trimws(as.character(term)[1L])
  if (!nzchar(term)) {
    return(rep(FALSE, length(text)))
  }
  parts <- strsplit(term, "\\s+", perl = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  if (!length(parts)) {
    return(rep(FALSE, length(text)))
  }
  esc <- .escape_perl_regex(parts)
  pat <- if (length(esc) == 1L) {
    paste0("(?i)\\b", esc, "\\b")
  } else {
    paste0("(?i)\\b", paste(esc, collapse = "\\s+"), "\\b")
  }
  if (!isTRUE(ignore.case)) {
    pat <- if (length(esc) == 1L) {
      paste0("\\b", esc, "\\b")
    } else {
      paste0("\\b", paste(esc, collapse = "\\s+"), "\\b")
    }
  }
  grepl(pat, text, perl = TRUE)
}

.keyword_match_score <- function(text,
                                 query,
                                 match = c("all", "any"),
                                 ignore.case = TRUE,
                                 terms = NULL) {
  match <- match.arg(match)
  if (is.null(terms)) {
    terms <- .keyword_terms_from_query(query)
  } else {
    terms <- .filter_keyword_terms(terms)
  }
  n <- length(text)
  if (!length(terms)) {
    return(rep(0, n))
  }

  hits <- vapply(
    terms,
    function(term) .keyword_term_matches(text, term, ignore.case = ignore.case),
    logical(n)
  )
  if (!is.matrix(hits)) {
    hits <- matrix(hits, nrow = n, ncol = length(terms))
  }

  if (match == "all") {
    as.numeric(apply(hits, 1L, all))
  } else {
    apply(hits, 1L, mean)
  }
}

#' Per-text keyword hit details (score, matched, unmatched)
#' @keywords internal
.keyword_match_details <- function(text, terms, ignore.case = TRUE) {
  terms <- .filter_keyword_terms(terms)
  n <- length(text)
  empty <- list(
    keyword_score = 0,
    matched_keywords = character(),
    unmatched_keywords = terms
  )
  if (!length(terms)) {
    return(lapply(seq_len(n), function(i) empty))
  }

  lapply(seq_len(n), function(i) {
    txt <- text[[i]]
    hit <- vapply(
      terms,
      function(term) .keyword_term_matches(txt, term, ignore.case = ignore.case),
      logical(1L)
    )
    matched <- terms[hit]
    unmatched <- terms[!hit]
    list(
      keyword_score = mean(hit),
      matched_keywords = matched,
      unmatched_keywords = unmatched
    )
  })
}

.dedupe_candidate_by_gene <- function(x) {
  if (!"Gene_ID" %in% names(x)) {
    return(x)
  }
  ord <- order(-x$match_score, x$Gene_ID)
  x <- x[ord, , drop = FALSE]
  x[!duplicated(x$Gene_ID), , drop = FALSE]
}

################################################################################
#' Search candidate genes by functional annotation
#'
#' Filter and rank genes in a candidate list using keyword matching and/or
#' semantic similarity (offline LSA via the \pkg{text2vec} package).
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
#' @param mode Search mode: `"keyword"` (fixed substring match),
#'   `"semantic"` (text2vec LSA cosine similarity), or `"both"`.
#' @param keyword_match For keyword mode, `"all"` requires every query token to
#'   appear in the annotation text; `"any"` requires at least one token.
#' @param ignore.case Passed to [grepl()] for keyword matching.
#' @param min_score Minimum cosine similarity for semantic matches (0--1).
#' @param top_n If not `NULL`, return at most this many rows after ranking by
#'   `match_score`.
#' @param dedupe_genes If `TRUE` and `Gene_ID` is present, keep one row per
#'   gene with the highest `match_score`.
#' @param n_topics Number of LSA topics. Default uses up to 30 topics bounded by
#'   corpus size.
#'
#' @return A subset of the input candidate `data.frame` with columns
#'   `keyword_score`, `semantic_score`, `match_score`, and `match_method`
#'   added, sorted by decreasing `match_score`.
#'
#' @details
#' Semantic search requires the suggested package \pkg{text2vec}:
#' `install.packages("text2vec")`.
#'
#' Annotation text for each row is the space-separated concatenation of
#' non-empty values in `ann_cols`.
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
                                 mode = c("keyword", "semantic", "both"),
                                 keyword_match = c("all", "any"),
                                 ignore.case = TRUE,
                                 min_score = 0.1,
                                 top_n = NULL,
                                 dedupe_genes = TRUE,
                                 n_topics = NULL) {
  mode <- match.arg(mode)
  keyword_match <- match.arg(keyword_match)

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

  if (mode %in% c("semantic", "both")) {
    .require_text2vec()
  }

  ann_text <- .candidate_annotation_text(candidate = candidate, ann_cols = ann_cols)

  keyword_score <- rep(0, nrow(candidate))
  semantic_score <- rep(NA_real_, nrow(candidate))

  if (mode %in% c("keyword", "both")) {
    keyword_score <- .keyword_match_score(
      text = ann_text,
      query = query,
      match = keyword_match,
      ignore.case = ignore.case
    )
  }

  if (mode %in% c("semantic", "both")) {
    semantic_score <- .text2vec_semantic_score(
      text = ann_text,
      query = query,
      n_topics = n_topics
    )
  }
  semantic_score[is.na(semantic_score)] <- 0

  kw_hit <- if (keyword_match == "all") {
    keyword_score >= 1
  } else {
    keyword_score > 0
  }
  sem_hit <- !is.na(semantic_score) & semantic_score >= min_score

  keep <- switch(
    mode,
    keyword = kw_hit,
    semantic = sem_hit,
    both = kw_hit | sem_hit
  )

  if (!any(keep)) {
    out <- candidate[0, , drop = FALSE]
    attr(out, "searchCandidateGenes") <- list(
      query = query,
      mode = mode,
      ann_cols = ann_cols,
      n_matched = 0L
    )
    return(out)
  }

  out <- candidate[keep, , drop = FALSE]
  out$keyword_score <- keyword_score[keep]
  out$semantic_score <- semantic_score[keep]

  out$match_score <- switch(
    mode,
    keyword = out$keyword_score,
    semantic = out$semantic_score,
    both = pmax(out$keyword_score, out$semantic_score, na.rm = TRUE)
  )

  out$match_method <- ifelse(
    out$keyword_score > 0 & out$semantic_score >= min_score,
    "both",
    ifelse(out$keyword_score > 0, "keyword",
           ifelse(out$semantic_score >= min_score, "semantic", "none"))
  )
  if (mode == "keyword") {
    out$match_method <- ifelse(out$keyword_score > 0, "keyword", "none")
  } else if (mode == "semantic") {
    out$match_method <- ifelse(out$semantic_score >= min_score, "semantic", "none")
  }

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
    min_score = min_score,
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

.keyword_match_score <- function(text, query, match = c("all", "any"), ignore.case = TRUE) {
  match <- match.arg(match)
  terms <- strsplit(query, "\\s+", perl = TRUE)[[1L]]
  terms <- terms[nzchar(terms)]
  n <- length(text)
  if (!length(terms)) {
    return(rep(0, n))
  }

  text_cmp <- if (ignore.case) tolower(text) else text
  hits <- vapply(
    terms,
    function(term) {
      term_cmp <- if (ignore.case) tolower(term) else term
      grepl(term_cmp, text_cmp, fixed = TRUE)
    },
    logical(n)
  )
  if (!is.matrix(hits)) {
    hits <- matrix(hits, nrow = n, ncol = 1L)
  }

  if (match == "all") {
    as.numeric(apply(hits, 1L, all))
  } else {
    apply(hits, 1L, mean)
  }
}

.require_text2vec <- function() {
  if (!requireNamespace("text2vec", quietly = TRUE)) {
    stop(
      "Package 'text2vec' is required for semantic search. ",
      "Install with install.packages(\"text2vec\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.text2vec_semantic_score <- function(text, query, n_topics = NULL) {
  n <- length(text)
  empty <- !nzchar(text)
  scores <- rep(0, n)

  usable <- which(!empty)
  if (!length(usable)) {
    return(scores)
  }

  docs <- text[usable]
  prep <- get("tolower", envir = asNamespace("base"))
  tok <- get("word_tokenizer", envir = asNamespace("text2vec"))

  it <- text2vec::itoken(
    docs,
    preprocessor = prep,
    tokenizer = tok,
    progressbar = FALSE
  )
  vocab <- text2vec::create_vocabulary(it)
  if (nrow(vocab) == 0L) {
    return(scores)
  }

  vectorizer <- text2vec::vocab_vectorizer(vocab)
  dtm <- text2vec::create_dtm(it, vectorizer)
  if (nrow(dtm) < 1L || ncol(dtm) < 1L) {
    return(scores)
  }

  tfidf <- text2vec::TfIdf$new()
  dtm_tfidf <- text2vec::fit_transform(dtm, tfidf)

  n_row <- nrow(dtm_tfidf)
  n_col <- ncol(dtm_tfidf)
  k <- if (is.null(n_topics)) {
    min(30L, max(2L, min(n_row - 1L, n_col - 1L)))
  } else {
    as.integer(n_topics)[1L]
  }
  k <- max(2L, min(k, n_row - 1L, n_col - 1L))

  if (k >= 2L && n_row >= 2L && n_col >= 2L) {
    lsa <- text2vec::LatentSemanticAnalysis$new(n_topics = k)
    doc_emb <- text2vec::fit_transform(dtm_tfidf, lsa)
  } else {
    lsa <- NULL
    doc_emb <- dtm_tfidf
  }

  query_it <- text2vec::itoken(
    query,
    preprocessor = prep,
    tokenizer = tok,
    progressbar = FALSE
  )
  query_dtm <- text2vec::create_dtm(query_it, vectorizer)
  if (nrow(query_dtm) < 1L || ncol(query_dtm) < 1L) {
    return(scores)
  }

  query_tfidf <- tfidf$transform(query_dtm)
  query_emb <- if (!is.null(lsa)) {
    lsa$transform(query_tfidf)
  } else {
    query_tfidf
  }

  sim <- text2vec::sim2(
    x = doc_emb,
    y = query_emb,
    method = "cosine",
    norm = "l2"
  )
  sim_v <- as.numeric(sim[, 1L])
  sim_v[!is.finite(sim_v)] <- 0
  scores[usable] <- pmax(0, sim_v)
  scores
}

.dedupe_candidate_by_gene <- function(x) {
  if (!"Gene_ID" %in% names(x)) {
    return(x)
  }
  ord <- order(-x$match_score, x$Gene_ID)
  x <- x[ord, , drop = FALSE]
  x[!duplicated(x$Gene_ID), , drop = FALSE]
}

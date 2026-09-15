################################################################################
# Annotation keyword matching helpers (internal).
# Used by collectGeneEvidence() / ranking.
# The public searchCandidateGenes() API was removed (Phase 2 A3).
################################################################################

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
#'
#' Includes user trait keywords, tissues/stage/conditions, and (when present)
#' synonym / related phrases from Phase 2 A2.
#' @keywords internal
.keyword_terms_from_query <- function(query) {
  if (inherits(query, "PhenotypeQuery")) {
    syn <- if (exists(".phenotype_linked_phrase_strings", mode = "function")) {
      .phenotype_linked_phrase_strings(query$synonym_phrases)
    } else {
      character()
    }
    rel <- if (exists(".phenotype_linked_phrase_strings", mode = "function")) {
      .phenotype_linked_phrase_strings(query$related_phrases)
    } else {
      character()
    }
    raw <- c(
      query$trait_keywords,
      query$tissues,
      query$developmental_stage,
      query$conditions,
      syn,
      rel
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

#' Family-level annotation keyword details for ranking / A4 reports
#' @keywords internal
.annotation_keyword_family_details <- function(text, query, ignore.case = TRUE) {
  text <- as.character(text)
  n <- length(text)
  empty <- list(
    keyword_score = 0,
    matched_keywords = character(),
    matched_synonym = list(),
    matched_related = list(),
    unmatched_keywords = character(),
    matched_context = character()
  )
  if (!inherits(query, "PhenotypeQuery")) {
    terms <- .keyword_terms_from_query(query)
    dets <- .keyword_match_details(text, terms, ignore.case = ignore.case)
    return(lapply(dets, function(d) {
      list(
        keyword_score = d$keyword_score,
        matched_keywords = d$matched_keywords,
        matched_synonym = list(),
        matched_related = list(),
        unmatched_keywords = d$unmatched_keywords,
        matched_context = character()
      )
    }))
  }

  user <- .filter_keyword_terms(query$trait_keywords %||% character())
  ctx <- .filter_keyword_terms(c(
    query$tissues %||% character(),
    query$developmental_stage %||% character(),
    query$conditions %||% character()
  ))
  syn <- query$synonym_phrases %||% list()
  rel <- query$related_phrases %||% list()
  if (is.data.frame(syn)) {
    syn <- lapply(seq_len(nrow(syn)), function(i) as.list(syn[i, , drop = FALSE]))
  }
  if (is.data.frame(rel)) {
    rel <- lapply(seq_len(nrow(rel)), function(i) as.list(rel[i, , drop = FALSE]))
  }

  search_terms <- .keyword_terms_from_query(query)
  if (!length(search_terms) && !length(user) && !length(ctx)) {
    return(lapply(seq_len(n), function(i) empty))
  }

  lapply(seq_len(n), function(i) {
    txt <- text[[i]]
    hit_term <- function(term) {
      .keyword_term_matches(txt, term, ignore.case = ignore.case)
    }
    matched_user <- user[vapply(user, hit_term, logical(1L))]
    matched_ctx <- ctx[vapply(ctx, hit_term, logical(1L))]

    matched_syn <- list()
    for (item in syn) {
      ph <- as.character(item$phrase %||% "")[1L]
      if (!nzchar(ph) || !hit_term(ph)) {
        next
      }
      matched_syn[[length(matched_syn) + 1L]] <- list(
        phrase = ph,
        from_user = as.character(item$from_user %||% NA_character_)[1L]
      )
    }
    matched_rel <- list()
    for (item in rel) {
      ph <- as.character(item$phrase %||% "")[1L]
      if (!nzchar(ph) || !hit_term(ph)) {
        next
      }
      matched_rel[[length(matched_rel) + 1L]] <- list(
        phrase = ph,
        from_user = as.character(item$from_user %||% NA_character_)[1L]
      )
    }

    # Family-level unmatched: user phrase with no self/syn/rel hit
    unmatched_user <- character()
    for (u in user) {
      syn_for <- vapply(syn, function(item) {
        identical(
          tolower(as.character(item$from_user %||% "")[1L]),
          tolower(u)
        ) && hit_term(as.character(item$phrase %||% "")[1L])
      }, logical(1L))
      rel_for <- vapply(rel, function(item) {
        identical(
          tolower(as.character(item$from_user %||% "")[1L]),
          tolower(u)
        ) && hit_term(as.character(item$phrase %||% "")[1L])
      }, logical(1L))
      if (!hit_term(u) && !any(syn_for) && !any(rel_for)) {
        unmatched_user <- c(unmatched_user, u)
      }
    }

    score_terms <- unique(c(
      user,
      vapply(syn, function(z) as.character(z$phrase %||% ""), character(1L)),
      vapply(rel, function(z) as.character(z$phrase %||% ""), character(1L)),
      ctx
    ))
    score_terms <- .filter_keyword_terms(score_terms)
    kw_score <- if (!length(score_terms)) {
      0
    } else {
      mean(vapply(score_terms, hit_term, logical(1L)))
    }

    list(
      keyword_score = kw_score,
      matched_keywords = matched_user,
      matched_synonym = matched_syn,
      matched_related = matched_rel,
      unmatched_keywords = unique(unmatched_user),
      matched_context = matched_ctx
    )
  })
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

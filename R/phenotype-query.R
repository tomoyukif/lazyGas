################################################################################
#' Build a structured phenotype exploration query
#'
#' Parses free text (and optional tissue/stage/condition fields) into a
#' \code{PhenotypeQuery} object used by [rankPhenotypeCandidates()] and
#' [collectGeneEvidence()]. Uses a deterministic token-based parser only
#' (no LLM).
#'
#' @param text Character string describing the phenotype of interest.
#' @param tissues Character vector of tissues or organs (optional).
#' @param stage Character vector of developmental stages (optional).
#' @param conditions Character vector of environmental or experimental
#'   conditions (optional).
#' @param species Target species name (optional).
#' @param use_llm Ignored; kept for call compatibility (always deterministic).
#' @param synonym_phrases,related_phrases Ignored in this build (always empty
#'   unless you assign them on the returned object).
#' @param action_phrases Character vector of verb+noun (or similar) phrases for
#'   the annotation channel tier score \code{1} (e.g. \code{"regulates flowering"}).
#' @param noun_phrases Character vector of noun phrases for tier score \code{0.5}
#'   (e.g. \code{"fruit weight"}, \code{"leaf development"}).
#' @param words Character vector of single-token terms for tier score \code{0.25}.
#'   When \code{action_phrases}, \code{noun_phrases}, and \code{words} are all
#'   empty, multi-word \code{trait_keywords} become noun phrases and single-token
#'   keywords become words (backward-compatible fallback).
#' @param llm_model,llm_base_url,llm_timeout Ignored.
#' @param object Optional \code{LazyGas} object; when \code{save = TRUE} the
#'   query is written to the companion store.
#' @param save Save the query JSON to the companion store when \code{object} is
#'   given.
#'
#' @return A \code{PhenotypeQuery} object (list with class attribute).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()]
#'
#' @examples
#' q <- phenotypeQuery(
#'   "fruit weight at maturity",
#'   tissues = "fruit",
#'   stage = "maturation",
#'   action_phrases = "regulates fruit weight",
#'   noun_phrases = c("fruit weight", "fruit development"),
#'   words = c("fruit", "weight")
#' )
#' q$trait_keywords
phenotypeQuery <- function(text,
                           tissues = NULL,
                           stage = NULL,
                           conditions = NULL,
                           species = NULL,
                           use_llm = FALSE,
                           synonym_phrases = NULL,
                           related_phrases = NULL,
                           action_phrases = NULL,
                           noun_phrases = NULL,
                           words = NULL,
                           llm_model = NULL,
                           llm_base_url = NULL,
                           llm_timeout = 600,
                           object = NULL,
                           save = !is.null(object)) {
  if (!is.character(text) || length(text) != 1L || !nzchar(trimws(text))) {
    stop("'text' must be a non-empty character string.", call. = FALSE)
  }
  text <- trimws(text)
  if (isTRUE(use_llm)) {
    warning(
      "'use_llm = TRUE' is not supported in this build; using deterministic parser.",
      call. = FALSE
    )
  }

  query_id <- .phenotype_query_new_id()
  parsed <- .phenotype_query_fallback(
    text = text,
    tissues = tissues,
    stage = stage,
    conditions = conditions,
    species = species
  )

  out <- structure(
    list(
      trait_text = text,
      trait_keywords = parsed$trait_keywords,
      tissues = parsed$tissues,
      developmental_stage = parsed$developmental_stage,
      conditions = parsed$conditions,
      species = .phenotype_query_normalize_scalar(parsed$species),
      synonym_phrases = list(),
      related_phrases = list(),
      action_phrases = .phenotype_query_normalize_phrase_vec(action_phrases),
      noun_phrases = .phenotype_query_normalize_phrase_vec(noun_phrases),
      words = .phenotype_query_normalize_phrase_vec(words),
      query_id = query_id,
      parsed_by = "fallback"
    ),
    class = c("PhenotypeQuery", "list")
  )
  out <- .phenotype_query_fill_annotation_tiers(out)

  if (save && !is.null(object)) {
    if (!inherits(object, "LazyGas")) {
      stop("'object' must be a LazyGas object when save = TRUE.", call. = FALSE)
    }
    .store_write_phenotype_query(object, out)
  }

  out
}

.phenotype_query_has_value <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  if (is.list(x) && !length(x)) {
    return(FALSE)
  }
  x <- as.character(unlist(x, use.names = FALSE))
  any(!is.na(x) & nzchar(x))
}

#' @export
print.PhenotypeQuery <- function(x, ...) {
  cat("PhenotypeQuery:", x$query_id, "\n")
  cat("  trait:", x$trait_text, "\n")
  if (.phenotype_query_has_value(x$tissues)) {
    cat("  tissues:", paste(as.character(unlist(x$tissues)), collapse = ", "), "\n")
  }
  if (.phenotype_query_has_value(x$developmental_stage)) {
    cat(
      "  stage:",
      paste(as.character(unlist(x$developmental_stage)), collapse = ", "),
      "\n"
    )
  }
  if (.phenotype_query_has_value(x$conditions)) {
    cat(
      "  conditions:",
      paste(as.character(unlist(x$conditions)), collapse = ", "),
      "\n"
    )
  }
  if (.phenotype_query_has_value(x$species)) {
    cat("  species:", as.character(unlist(x$species))[1L], "\n")
  }
  if (.phenotype_query_has_value(x$action_phrases)) {
    cat(
      "  action_phrases:",
      paste(as.character(unlist(x$action_phrases)), collapse = "; "),
      "\n"
    )
  }
  if (.phenotype_query_has_value(x$noun_phrases)) {
    cat(
      "  noun_phrases:",
      paste(as.character(unlist(x$noun_phrases)), collapse = "; "),
      "\n"
    )
  }
  if (.phenotype_query_has_value(x$words)) {
    cat("  words:", paste(as.character(unlist(x$words)), collapse = ", "), "\n")
  }
  n_syn <- length(x$synonym_phrases %||% list())
  n_rel <- length(x$related_phrases %||% list())
  if (n_syn || n_rel) {
    cat("  synonym_phrases:", n_syn, " related_phrases:", n_rel, "\n")
  }
  cat("  parsed_by:", x$parsed_by, "\n")
  invisible(x)
}

#' Normalize phrase vectors without splitting on spaces
#' @keywords internal
.phenotype_query_normalize_phrase_vec <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (length(x) == 1L && is.na(x)) {
    return(character())
  }
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

#' Fill annotation tier fields from trait_keywords when user left them empty
#' @keywords internal
.phenotype_query_fill_annotation_tiers <- function(query) {
  has_explicit <- length(query$action_phrases %||% character()) ||
    length(query$noun_phrases %||% character()) ||
    length(query$words %||% character())
  if (isTRUE(has_explicit)) {
    return(query)
  }

  kws <- as.character(query$trait_keywords %||% character())
  kws <- trimws(kws)
  kws <- kws[nzchar(kws)]
  noun <- character()
  word <- character()
  for (kw in kws) {
    parts <- strsplit(kw, "\\s+", perl = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 2L) {
      noun <- c(noun, kw)
    } else if (length(parts) == 1L) {
      word <- c(word, kw)
    }
  }
  query$noun_phrases <- unique(noun)
  query$words <- unique(word)
  query$action_phrases <- character()
  query
}

#' Annotation tier term lists from a PhenotypeQuery (or character)
#' @keywords internal
.annotation_tier_terms_from_query <- function(query) {
  if (is.character(query)) {
    toks <- .phenotype_query_tokenize(query)
    multi <- character()
    single <- character()
    # Prefer multi-word substrings from the original text when possible
    raw <- trimws(as.character(query)[1L])
    if (nzchar(raw)) {
      # Split on / | ; for phrase-like chunks
      chunks <- unlist(strsplit(raw, "[/|;]+", perl = TRUE), use.names = FALSE)
      chunks <- trimws(chunks)
      for (ch in chunks) {
        if (!nzchar(ch)) next
        parts <- strsplit(ch, "\\s+", perl = TRUE)[[1L]]
        parts <- parts[nzchar(parts)]
        if (length(parts) >= 2L) {
          multi <- c(multi, paste(parts, collapse = " "))
        } else if (length(parts) == 1L) {
          single <- c(single, parts)
        }
      }
    }
    if (!length(multi) && !length(single)) {
      single <- toks
    }
    return(list(
      action = character(),
      noun = unique(.filter_keyword_terms(multi)),
      word = unique(.filter_keyword_terms(single))
    ))
  }
  if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  }
  q <- .phenotype_query_fill_annotation_tiers(query)
  list(
    action = unique(.filter_keyword_terms(q$action_phrases %||% character())),
    noun = unique(.filter_keyword_terms(q$noun_phrases %||% character())),
    word = unique(.filter_keyword_terms(q$words %||% character()))
  )
}

#' Convert a PhenotypeQuery to a single search string
#' @param query A \code{PhenotypeQuery} or character string.
#' @keywords internal
.phenotype_query_as_search_text <- function(query) {
  if (is.character(query)) {
    return(trimws(query))
  }
  if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be a PhenotypeQuery object or character string.",
         call. = FALSE)
  }
  parts <- c(
    query$trait_text,
    query$tissues,
    query$developmental_stage,
    query$conditions
  )
  parts <- parts[!is.na(parts) & nzchar(parts)]
  paste(unique(parts), collapse = " ")
}

.phenotype_query_new_id <- function() {
  paste0(
    "pq_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    sample.int(8999L, 1L) + 1000L
  )
}

.phenotype_query_tokenize <- function(x) {
  x <- tolower(trimws(x))
  terms <- strsplit(x, "[^a-zA-Z0-9]+", perl = TRUE)[[1L]]
  terms <- terms[nzchar(terms) & nchar(terms) >= 2L]
  unique(terms)
}

.phenotype_query_normalize_vec <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  if (length(x) == 1L && is.na(x)) {
    return(character())
  }
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (!length(x)) {
    return(character())
  }
  unique(unlist(strsplit(x, "[,;|]+", perl = TRUE)))
}

.phenotype_query_fallback <- function(text,
                                      tissues = NULL,
                                      stage = NULL,
                                      conditions = NULL,
                                      species = NULL) {
  kw <- .phenotype_query_tokenize(text)
  tissues <- .phenotype_query_normalize_vec(tissues)
  stage <- .phenotype_query_normalize_vec(stage)
  conditions <- .phenotype_query_normalize_vec(conditions)

  if (length(tissues)) {
    kw <- unique(c(kw, tolower(tissues)))
  }
  if (length(stage)) {
    kw <- unique(c(kw, tolower(stage)))
  }
  if (length(conditions)) {
    kw <- unique(c(kw, tolower(conditions)))
  }

  list(
    trait_keywords = kw,
    tissues = tissues,
    developmental_stage = stage,
    conditions = conditions,
    species = .phenotype_query_normalize_scalar(species),
    parsed_by = "fallback"
  )
}

.phenotype_query_normalize_scalar <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.list(x) && !length(x)) {
    return(NULL)
  }
  x <- as.character(unlist(x, use.names = FALSE))
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) {
    return(NULL)
  }
  trimws(x[1L])
}

.phenotype_query_normalize_loaded <- function(out) {
  out$trait_text <- if (is.null(out$trait_text)) {
    NA_character_
  } else {
    as.character(unlist(out$trait_text, use.names = FALSE))[1L]
  }
  out$trait_keywords <- .phenotype_query_normalize_vec(out$trait_keywords)
  out$tissues <- .phenotype_query_normalize_vec(out$tissues)
  out$developmental_stage <- .phenotype_query_normalize_vec(
    out$developmental_stage
  )
  out$conditions <- .phenotype_query_normalize_vec(out$conditions)
  out$species <- .phenotype_query_normalize_scalar(out$species)
  out$synonym_phrases <- .phenotype_normalize_linked_phrases(
    out$synonym_phrases %||% list(),
    trait_text = out$trait_text %||% "",
    user_phrases = out$trait_keywords,
    source_default = "user"
  )
  out$related_phrases <- .phenotype_normalize_linked_phrases(
    out$related_phrases %||% list(),
    trait_text = out$trait_text %||% "",
    user_phrases = out$trait_keywords,
    source_default = "user"
  )
  out$query_id <- if (is.null(out$query_id)) {
    NA_character_
  } else {
    as.character(unlist(out$query_id, use.names = FALSE))[1L]
  }
  out$parsed_by <- if (is.null(out$parsed_by)) {
    NA_character_
  } else {
    as.character(unlist(out$parsed_by, use.names = FALSE))[1L]
  }
  out
}

.store_write_phenotype_query <- function(object, query) {
  if (.store_is_gds(object)) {
    return(invisible(NULL))
  }
  dir.create(
    file.path(.store_path(object), "phenotype_query"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  path <- file.path(
    .store_path(object),
    "phenotype_query",
    paste0(.store_safe_name(query$query_id), ".json")
  )
  jsonlite::write_json(
    query,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(path)
}

.store_read_phenotype_query <- function(object, query_id) {
  if (.store_is_gds(object)) {
    return(NULL)
  }
  path <- file.path(
    .store_path(object),
    "phenotype_query",
    paste0(.store_safe_name(query_id), ".json")
  )
  if (!file.exists(path)) {
    return(NULL)
  }
  out <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  out <- .phenotype_query_normalize_loaded(out)
  class(out) <- c("PhenotypeQuery", "list")
  out
}

.store_list_phenotype_queries <- function(object) {
  if (.store_is_gds(object)) {
    return(character())
  }
  dir <- file.path(.store_path(object), "phenotype_query")
  if (!dir.exists(dir)) {
    return(character())
  }
  files <- list.files(dir, pattern = "\\.json$", full.names = FALSE)
  sub("\\.json$", "", files)
}

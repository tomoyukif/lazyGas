################################################################################
#' Build a structured phenotype exploration query
#'
#' Parses free text (and optional tissue/stage/condition fields) into a
#' \code{PhenotypeQuery} object used by [rankPhenotypeCandidates()] and
#' [collectGeneEvidence()]. When \code{use_llm = TRUE} and a local Ollama
#' server is reachable, the query is structured via JSON output from the LLM;
#' otherwise a deterministic token-based parser is used.
#'
#' @param text Character string describing the phenotype of interest.
#' @param tissues Character vector of tissues or organs (optional).
#' @param stage Character vector of developmental stages (optional).
#' @param conditions Character vector of environmental or experimental
#'   conditions (optional).
#' @param species Target species name (optional).
#' @param use_llm If \code{TRUE}, attempt LLM-based structuring via
#'   [llmChat()].
#' @param llm_model Model name passed to Ollama (default from
#'   \code{LAZYGAS_LLM_MODEL} or \code{"llama3.2:3b"}).
#' @param llm_base_url Ollama base URL (default from \code{LAZYGAS_LLM_URL}).
#' @param object Optional \code{LazyGas} object; when \code{save = TRUE} the
#'   query is written to the companion store.
#' @param save Save the query JSON to the companion store when \code{object} is
#'   given.
#'
#' @return A \code{PhenotypeQuery} object (list with class attribute).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [explainPhenotypeCandidates()]
#'
#' @examples
#' q <- phenotypeQuery(
#'   "fruit weight at maturity",
#'   tissues = "fruit",
#'   stage = "maturation",
#'   use_llm = FALSE
#' )
#' q$trait_keywords
phenotypeQuery <- function(text,
                           tissues = NULL,
                           stage = NULL,
                           conditions = NULL,
                           species = NULL,
                           use_llm = FALSE,
                           llm_model = NULL,
                           llm_base_url = NULL,
                           object = NULL,
                           save = !is.null(object)) {
  if (!is.character(text) || length(text) != 1L || !nzchar(trimws(text))) {
    stop("'text' must be a non-empty character string.", call. = FALSE)
  }
  text <- trimws(text)

  query_id <- .phenotype_query_new_id()

  parsed <- if (isTRUE(use_llm)) {
    .phenotype_query_via_llm(
      text = text,
      tissues = tissues,
      stage = stage,
      conditions = conditions,
      species = species,
      llm_model = llm_model,
      llm_base_url = llm_base_url
    )
  } else {
    NULL
  }

  if (is.null(parsed)) {
    parsed <- .phenotype_query_fallback(
      text = text,
      tissues = tissues,
      stage = stage,
      conditions = conditions,
      species = species
    )
  }

  out <- structure(
    list(
      trait_text = text,
      trait_keywords = parsed$trait_keywords,
      tissues = parsed$tissues,
      developmental_stage = parsed$developmental_stage,
      conditions = parsed$conditions,
      species = .phenotype_query_normalize_scalar(parsed$species),
      query_id = query_id,
      parsed_by = parsed$parsed_by
    ),
    class = c("PhenotypeQuery", "list")
  )

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
  cat("  parsed_by:", x$parsed_by, "\n")
  invisible(x)
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

.phenotype_query_via_llm <- function(text,
                                     tissues = NULL,
                                     stage = NULL,
                                     conditions = NULL,
                                     species = NULL,
                                     llm_model = NULL,
                                     llm_base_url = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(NULL)
  }
  llm_ok <- tryCatch(
    llmHealthCheck(base_url = llm_base_url, timeout = 5),
    error = function(e) FALSE
  )
  if (!isTRUE(llm_ok)) {
    warning(
      "Local LLM is not reachable; using fallback query parser.",
      call. = FALSE
    )
    return(NULL)
  }

  user_hints <- list(
    tissues = .phenotype_query_normalize_vec(tissues),
    developmental_stage = .phenotype_query_normalize_vec(stage),
    conditions = .phenotype_query_normalize_vec(conditions),
    species = species
  )
  user_hints <- user_hints[vapply(user_hints, length, integer(1L)) > 0L]

  system_msg <- paste(
    "You extract structured phenotype metadata from user text.",
    "Reply with JSON only, no markdown.",
    "Schema: {\"trait_keywords\": [strings], \"tissues\": [strings],",
    "\"developmental_stage\": [strings], \"conditions\": [strings],",
    "\"species\": string or null}."
  )
  user_msg <- jsonlite::toJSON(
    list(text = text, hints = user_hints),
    auto_unbox = TRUE
  )

  resp <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = system_msg),
        list(role = "user", content = user_msg)
      ),
      model = llm_model,
      base_url = llm_base_url,
      json_mode = TRUE,
      timeout = 90
    ),
    error = function(e) {
      warning("LLM query parsing failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(resp) || !nzchar(trimws(resp))) {
    return(NULL)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(resp, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || !is.list(parsed)) {
    return(NULL)
  }

  list(
    trait_keywords = unique(as.character(unlist(parsed$trait_keywords))),
    tissues = unique(as.character(unlist(parsed$tissues))),
    developmental_stage = unique(as.character(unlist(parsed$developmental_stage))),
    conditions = unique(as.character(unlist(parsed$conditions))),
    species = .phenotype_query_normalize_scalar(parsed$species),
    parsed_by = "llm"
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
  # Use null="null" so missing scalars serialize as JSON null, not {}.
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

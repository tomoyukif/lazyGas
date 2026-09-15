################################################################################
# Phase 2 A1/A2 vocabulary helpers for phenotypeQuery().
# A1: normalize LLM phrase extraction. A2: synonym / related expand (2 LLM calls).
################################################################################

.phenotype_species_denylist <- function() {
  c(
    "rice", "oryza", "sativa", "arabidopsis", "thaliana", "maize", "zea",
    "mays", "wheat", "barley", "sorghum", "soybean", "tomato", "oryza sativa",
    "arabidopsis thaliana", "zea mays"
  )
}

.phenotype_solo_noise_terms <- function() {
  c("content", "date", "time", "index", "value", "level", "rate")
}

.phenotype_meta_noise_terms <- function() {
  c(
    "trait", "phenotype", "gene", "protein", "putative", "homolog",
    "homologue", "similar", "unknown", "predicted"
  )
}

.phenotype_phrase_clean <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- x[!is.na(x) & nzchar(x)]
  unique(x)
}

.phenotype_looks_like_gene_token <- function(term) {
  term <- trimws(as.character(term)[1L])
  if (!nzchar(term)) {
    return(FALSE)
  }
  if (grepl("^os\\d{2}g\\d+$", term, ignore.case = TRUE, perl = TRUE)) {
    return(TRUE)
  }
  if (grepl("^[A-Za-z]{2,}[0-9]+[A-Za-z0-9]*$", term, perl = TRUE) &&
      grepl("[0-9]", term, perl = TRUE) &&
      !grepl("\\s", term, perl = TRUE)) {
    # HESO1, NAL1, GBSS1-style symbols; keep purely alphabetic trait names
    return(TRUE)
  }
  FALSE
}

#' Alphabetic gene nicknames without digits (A1/A5 scrub; provisional)
#' @keywords internal
.phenotype_gene_symbol_denylist <- function() {
  c("waxy", "wx", "heso1", "nal1", "gbss1", "gbssi")
}

.phenotype_is_gene_only_query <- function(trait_text) {
  t <- trimws(as.character(trait_text)[1L])
  if (!nzchar(t)) {
    return(FALSE)
  }
  parts <- strsplit(t, "[^A-Za-z0-9]+", perl = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  length(parts) == 1L && .phenotype_looks_like_gene_token(parts)
}

.phenotype_drop_species_tokens <- function(phrases, species = NULL) {
  denylist <- .phenotype_species_denylist()
  if (!is.null(species) && nzchar(trimws(as.character(species)[1L]))) {
    denylist <- unique(c(denylist, .phenotype_phrase_clean(species)))
  }
  out_species <- NULL
  keep <- character()
  for (p in phrases) {
    if (p %in% denylist) {
      if (is.null(out_species)) {
        out_species <- p
      }
      next
    }
    toks <- strsplit(p, "\\s+", perl = TRUE)[[1L]]
    if (length(toks) && all(toks %in% denylist)) {
      if (is.null(out_species)) {
        out_species <- p
      }
      next
    }
    keep <- c(keep, p)
  }
  list(phrases = unique(keep), species = out_species)
}

.phenotype_filter_true_keywords <- function(phrases, trait_text) {
  phrases <- .phenotype_phrase_clean(phrases)
  gene_only <- .phenotype_is_gene_only_query(trait_text)
  keep <- vapply(phrases, function(p) {
    if (!nzchar(p)) {
      return(FALSE)
    }
    if (p %in% .phenotype_meta_noise_terms()) {
      return(FALSE)
    }
    if (!gene_only &&
        (.phenotype_looks_like_gene_token(p) ||
         p %in% .phenotype_gene_symbol_denylist())) {
      return(FALSE)
    }
    parts <- strsplit(p, "\\s+", perl = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    if (!length(parts)) {
      return(FALSE)
    }
    # Drop English stopword-only phrases
    if (exists(".is_content_keyword_term", mode = "function")) {
      if (!any(vapply(parts, .is_content_keyword_term, logical(1L)))) {
        return(FALSE)
      }
    }
    if (length(parts) == 1L && parts %in% .phenotype_solo_noise_terms()) {
      return(FALSE)
    }
    TRUE
  }, logical(1L))
  phrases[keep]
}

#' Normalize LLM-parsed PhenotypeQuery fields (A1 TRUE path)
#' @keywords internal
.phenotype_query_normalize_true <- function(parsed,
                                            trait_text,
                                            tissues = NULL,
                                            stage = NULL,
                                            conditions = NULL,
                                            species = NULL) {
  kw <- .phenotype_phrase_clean(parsed$trait_keywords)
  # Prefer caller hints for field roles; do not merge them into trait_keywords
  tissues_v <- .phenotype_phrase_clean(
    c(.phenotype_query_normalize_vec(tissues), parsed$tissues)
  )
  stage_v <- .phenotype_phrase_clean(
    c(
      .phenotype_query_normalize_vec(stage),
      parsed$developmental_stage
    )
  )
  cond_v <- .phenotype_phrase_clean(
    c(.phenotype_query_normalize_vec(conditions), parsed$conditions)
  )

  sp_arg <- .phenotype_query_normalize_scalar(species)
  if (is.null(sp_arg)) {
    sp_arg <- .phenotype_query_normalize_scalar(parsed$species)
  }

  dropped <- .phenotype_drop_species_tokens(kw, species = sp_arg)
  kw <- dropped$phrases
  if (is.null(sp_arg) && !is.null(dropped$species)) {
    sp_arg <- dropped$species
  }

  # Also strip species tokens from context fields
  tissues_v <- .phenotype_drop_species_tokens(tissues_v, species = sp_arg)$phrases
  stage_v <- .phenotype_drop_species_tokens(stage_v, species = sp_arg)$phrases
  cond_v <- .phenotype_drop_species_tokens(cond_v, species = sp_arg)$phrases

  kw <- .phenotype_filter_true_keywords(kw, trait_text = trait_text)
  if (!length(kw) && .phenotype_is_gene_only_query(trait_text)) {
    kw <- .phenotype_phrase_clean(trait_text)
  }
  if (!length(kw) && nzchar(trimws(trait_text))) {
    # Last resort: keep multi-word chunks split on slash/comma only
    chunks <- unlist(strsplit(trait_text, "[/|,;]+", perl = TRUE))
    kw <- .phenotype_filter_true_keywords(chunks, trait_text = trait_text)
  }

  list(
    trait_keywords = kw,
    tissues = tissues_v,
    developmental_stage = stage_v,
    conditions = cond_v,
    species = sp_arg,
    parsed_by = "llm"
  )
}

.phenotype_linked_phrase_cap <- function() {
  20L
}

#' Normalize list of {phrase, from_user, source} records
#' @keywords internal
.phenotype_normalize_linked_phrases <- function(x,
                                                trait_text = "",
                                                user_phrases = character(),
                                                source_default = "llm") {
  if (is.null(x) || !length(x)) {
    return(list())
  }
  # Allow character vector shortcut: phrases without from_user
  if (is.character(x)) {
    x <- lapply(x, function(p) {
      list(phrase = p, from_user = NA_character_, source = source_default)
    })
  }
  if (!is.list(x)) {
    return(list())
  }
  # jsonlite may return data.frame
  if (is.data.frame(x)) {
    x <- lapply(seq_len(nrow(x)), function(i) {
      as.list(x[i, , drop = FALSE])
    })
  }
  user_phrases <- .phenotype_phrase_clean(user_phrases)
  out <- list()
  for (item in x) {
    if (is.null(item)) {
      next
    }
    if (is.character(item) && length(item) == 1L) {
      item <- list(phrase = item)
    }
    if (!is.list(item)) {
      next
    }
    phrase <- .phenotype_phrase_clean(item$phrase %||% item[[1L]])[1L]
    if (!nzchar(phrase %||% "")) {
      next
    }
    if ((.phenotype_looks_like_gene_token(phrase) ||
         phrase %in% .phenotype_gene_symbol_denylist()) &&
        !.phenotype_is_gene_only_query(trait_text)) {
      next
    }
    sp <- .phenotype_drop_species_tokens(phrase)
    if (!length(sp$phrases)) {
      next
    }
    phrase <- sp$phrases[1L]
    from_user <- .phenotype_phrase_clean(
      item$from_user %||% item$from %||% NA_character_
    )[1L]
    if (!nzchar(from_user %||% "")) {
      from_user <- if (length(user_phrases)) user_phrases[1L] else NA_character_
    }
    source <- as.character(item$source %||% source_default)[1L]
    out[[length(out) + 1L]] <- list(
      phrase = phrase,
      from_user = from_user,
      source = source
    )
  }
  if (length(out) > .phenotype_linked_phrase_cap()) {
    out <- out[seq_len(.phenotype_linked_phrase_cap())]
  }
  # Dedupe by phrase + from_user
  keys <- vapply(out, function(z) {
    paste(z$phrase, z$from_user %||% "", sep = "\t")
  }, character(1L))
  out[!duplicated(keys)]
}

.phenotype_parse_linked_llm_json <- function(resp,
                                             trait_text,
                                             user_phrases,
                                             source_default = "llm") {
  if (is.null(resp) || !nzchar(trimws(resp))) {
    return(list())
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(resp, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed)) {
    return(list())
  }
  items <- parsed
  if (is.list(parsed) && !is.null(parsed$phrases)) {
    items <- parsed$phrases
  } else if (is.list(parsed) && !is.null(parsed$synonyms)) {
    items <- parsed$synonyms
  } else if (is.list(parsed) && !is.null(parsed$related)) {
    items <- parsed$related
  }
  .phenotype_normalize_linked_phrases(
    items,
    trait_text = trait_text,
    user_phrases = user_phrases,
    source_default = source_default
  )
}

#' Expand synonym phrases via a dedicated LLM call (A2)
#' @keywords internal
.phenotype_query_expand_synonyms <- function(trait_text,
                                             trait_keywords,
                                             llm_model = NULL,
                                             llm_base_url = NULL,
                                             llm_timeout = 600) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list())
  }
  llm_timeout <- as.numeric(llm_timeout)[1L]
  if (!is.finite(llm_timeout) || llm_timeout <= 0) {
    llm_timeout <- 600
  }
  llm_ok <- tryCatch(
    llmHealthCheck(base_url = llm_base_url, timeout = 5),
    error = function(e) FALSE
  )
  if (!isTRUE(llm_ok)) {
    return(list())
  }
  system_msg <- paste(
    "You list synonym phenotype phrases only.",
    "Reply with JSON only: {\"phrases\":[{\"phrase\":\"...\",\"from_user\":\"...\"},...]}",
    "Each phrase must be an alternative name for the same trait/measurement",
    "as from_user (one of the given trait_keywords).",
    "Do not add enzymes, pathways, gene symbols, gene IDs, or species names.",
    "Do not invent unrelated biology. Cap at 20 phrases."
  )
  user_msg <- jsonlite::toJSON(
    list(text = trait_text, trait_keywords = trait_keywords),
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
      timeout = llm_timeout
    ),
    error = function(e) {
      warning("Synonym expand failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  out <- .phenotype_parse_linked_llm_json(
    resp,
    trait_text = trait_text,
    user_phrases = trait_keywords,
    source_default = "llm"
  )
  if (!length(out) && !is.null(resp) && nzchar(trimws(resp))) {
    raw2 <- .llm_scrutiny_pass(
      llm = list(model = llm_model, base_url = llm_base_url, timeout = llm_timeout),
      rules = system_msg,
      prior_output = resp,
      problems = "could not parse phrases array / linked phrase objects",
      extra_user = user_msg,
      timeout = llm_timeout
    )
    if (!inherits(raw2, "error")) {
      out <- .phenotype_parse_linked_llm_json(
        raw2,
        trait_text = trait_text,
        user_phrases = trait_keywords,
        source_default = "llm"
      )
    }
  }
  out
}

#' Expand related phrases via a dedicated LLM call (A2)
#' @keywords internal
.phenotype_query_expand_related <- function(trait_text,
                                            trait_keywords,
                                            llm_model = NULL,
                                            llm_base_url = NULL,
                                            llm_timeout = 600) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list())
  }
  llm_timeout <- as.numeric(llm_timeout)[1L]
  if (!is.finite(llm_timeout) || llm_timeout <= 0) {
    llm_timeout <- 600
  }
  llm_ok <- tryCatch(
    llmHealthCheck(base_url = llm_base_url, timeout = 5),
    error = function(e) FALSE
  )
  if (!isTRUE(llm_ok)) {
    return(list())
  }
  system_msg <- paste(
    "You list related biological search phrases for a phenotype.",
    "Reply with JSON only: {\"phrases\":[{\"phrase\":\"...\",\"from_user\":\"...\"},...]}",
    "Related means processes, molecular functions, or compounds linked to",
    "from_user (one of the given trait_keywords), not mere synonyms.",
    "Do not include gene symbols or gene IDs (Waxy, HESO1, NAL1, Os..g..).",
    "Do not include species-only terms. Cap at 20 phrases."
  )
  user_msg <- jsonlite::toJSON(
    list(text = trait_text, trait_keywords = trait_keywords),
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
      timeout = llm_timeout
    ),
    error = function(e) {
      warning("Related expand failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  out <- .phenotype_parse_linked_llm_json(
    resp,
    trait_text = trait_text,
    user_phrases = trait_keywords,
    source_default = "llm"
  )
  if (!length(out) && !is.null(resp) && nzchar(trimws(resp))) {
    raw2 <- .llm_scrutiny_pass(
      llm = list(model = llm_model, base_url = llm_base_url, timeout = llm_timeout),
      rules = system_msg,
      prior_output = resp,
      problems = "could not parse phrases array / linked phrase objects",
      extra_user = user_msg,
      timeout = llm_timeout
    )
    if (!inherits(raw2, "error")) {
      out <- .phenotype_parse_linked_llm_json(
        raw2,
        trait_text = trait_text,
        user_phrases = trait_keywords,
        source_default = "llm"
      )
    }
  }
  out
}

#' Linked phrase strings from a PhenotypeQuery field
#' @keywords internal
.phenotype_linked_phrase_strings <- function(x) {
  if (is.null(x) || !length(x)) {
    return(character())
  }
  vapply(x, function(z) {
    if (is.list(z)) {
      as.character(z$phrase %||% "")[1L]
    } else {
      as.character(z)[1L]
    }
  }, character(1L))
}

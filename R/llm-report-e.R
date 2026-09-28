# Phase 2 E1–E9: gene×source payloads, Gemma calls, templates, HTML

.E_SOURCE_ORDER <- c(
  "keywords", "snpeff", "gwas", "expression",
  "go", "kegg", "domains", "summary", "interpretation", "validity"
)

#' @keywords internal
.e_resolution_label <- function(resolution) {
  r <- as.character(resolution %||% "low")[1L]
  if (identical(r, "high")) {
    "resolvable"
  } else if (identical(r, "moderate")) {
    "partially resolvable"
  } else {
    "not resolved"
  }
}

#' @keywords internal
.e_template_text <- function(source, gene_id, payload = NULL) {
  gid <- as.character(gene_id)[1L]
  switch(
    source,
    keywords = {
      m <- payload$matched %||% character()
      if (!length(m)) {
        paste0("Keywords did not match free-description text for ", gid, ".")
      } else {
        paste0(
          "Matched keywords for ", gid, ": ",
          paste(m, collapse = ", "), "."
        )
      }
    },
    snpeff = paste0(
      "Variant-effect narrative unavailable for ", gid, " (template)."
    ),
    gwas = paste0(
      "Mapping-context narrative unavailable for ", gid, " (template)."
    ),
    expression = {
      sp <- payload$specificity$summary %||%
        payload$specificity$combination_summary %||%
        payload$specificity$tissue_summary %||% "NA"
      rel <- payload$relative_level$band %||%
        payload$relative_level$summary %||% "NA"
      ht <- as.character(payload$high_tissues %||% character())
      ht <- ht[nzchar(ht)]
      base <- paste0(
        "Expression for ", gid, ": tissue specificity ", sp,
        "; relative level ", rel
      )
      if (length(ht) && !identical(sp, "low")) {
        paste0(base, "; highest in ", paste(ht, collapse = ", "), ".")
      } else {
        paste0(base, ".")
      }
    },
    go = {
      terms <- payload$terms %||% character()
      if (!length(terms)) {
        paste0("No GO terms for ", gid, ".")
      } else {
        paste0("GO terms for ", gid, ": ", paste(terms, collapse = "; "), ".")
      }
    },
    kegg = {
      entries <- payload$entries %||% list()
      if (!length(entries)) {
        paste0("No KEGG entries for ", gid, ".")
      } else {
        paste0("KEGG entries present for ", gid, " (see coded list).")
      }
    },
    domains = {
      doms <- payload$domains %||% list()
      if (!length(doms)) {
        paste0("No domain annotations for ", gid, ".")
      } else {
        paste0("Domain annotations present for ", gid, " (see coded list).")
      }
    },
    summary = {
      coded <- payload$coded
      if (is.null(coded)) {
        paste0(
          "Brief summary unavailable for ", gid,
          " (template). Review coded metrics. Scores are peak-relative."
        )
      } else {
        cq <- .e_coded_qualitative(coded)
        paste0(
          gid, ": proximity ", cq$proximity,
          "; association ", cq$association_signal,
          "; max PIP ", cq$max_PIP,
          "; worst impact ", cq$worst_impact %||% "NA",
          "; expression specificity ",
          cq$expression_specificity %||% "NA",
          {
            ht <- cq$expression_high_tissues
            if (length(ht)) {
              paste0(" (highest in ", paste(ht, collapse = ", "), ")")
            } else {
              ""
            }
          },
          ". Composite cues are peak-relative, not absolute across peaks."
        )
      }
    },
    interpretation = {
      summ <- as.character(payload$summary %||% "")[1L]
      paste0(
        "Integrated interpretation unavailable for ", gid,
        " (template). Plausibility: insufficient_evidence. ",
        if (nzchar(summ)) paste0("Brief summary: ", summ) else
          "Review source prose and coded labels."
      )
    },
    validity = {
      res <- payload$peak$resolution_label %||%
        payload$peak$resolution %||% "unknown"
      paste0(
        gid, ": credible-set resolution ", res,
        ". Review coded SnpEff vs mapping metrics."
      )
    },
    paste0("[template] ", gid, " / ", source)
  )
}

#' Split free-description cells into sentences (E1)
#' @keywords internal
.e_split_sentences <- function(text) {
  text <- as.character(text)
  text <- text[!is.na(text) & nzchar(trimws(text))]
  if (!length(text)) {
    return(character())
  }
  out <- character()
  for (chunk in text) {
    parts <- unlist(
      strsplit(chunk, "(?<=[.?!。？！])\\s+", perl = TRUE),
      use.names = FALSE
    )
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (!length(parts)) {
      out <- c(out, trimws(chunk))
    } else {
      out <- c(out, parts)
    }
  }
  unique(out)
}

#' Collect numeric allowlist from a payload (B5: numbers only)
#'
#' Includes numbers embedded in character fields (e.g. peak_ID "1") so the
#' model may echo identifiers that appeared in the prompt.
#' @keywords internal
.e_collect_allowed_numbers <- function(x, acc = numeric()) {
  if (is.null(x)) {
    return(acc)
  }
  if (is.numeric(x) || is.integer(x)) {
    v <- as.numeric(x)
    return(c(acc, v[is.finite(v)]))
  }
  if (is.logical(x)) {
    return(acc)
  }
  if (is.character(x)) {
    for (s in x) {
      if (is.na(s) || !nzchar(s)) {
        next
      }
      acc <- c(acc, .e_extract_numbers(s))
    }
    return(acc)
  }
  if (is.list(x) || is.data.frame(x)) {
    for (el in x) {
      acc <- .e_collect_allowed_numbers(el, acc)
    }
  }
  unique(acc[is.finite(acc)])
}

#' @keywords internal
.e_extract_numbers <- function(text) {
  text <- as.character(text)[1L]
  if (!nzchar(text)) {
    return(numeric())
  }
  m <- gregexpr(
    "(?<![A-Za-z0-9_])[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?(?![A-Za-z0-9_])",
    text,
    perl = TRUE
  )
  hits <- regmatches(text, m)[[1L]]
  if (!length(hits)) {
    return(numeric())
  }
  as.numeric(hits)
}

#' @keywords internal
.e_numbers_allowed <- function(found, allow, tol = 1e-8) {
  if (!length(found)) {
    return(logical())
  }
  if (!length(allow)) {
    return(rep(FALSE, length(found)))
  }
  vapply(found, function(v) {
    any(abs(allow - v) <= tol * pmax(1, abs(allow)))
  }, logical(1L))
}

#' B5: check numbers; LLM rewrite; if still bad, one scrutiny pass; else fail
#' @keywords internal
.e_b5_check_and_fix <- function(llm, text, allow, return_key) {
  found <- .e_extract_numbers(text)
  ok_mask <- .e_numbers_allowed(found, allow)
  if (!length(found) || all(ok_mask)) {
    return(list(ok = TRUE, text = text, rewritten = FALSE, reason = NULL))
  }
  bad <- found[!ok_mask]
  message(
    "B5: ungrounded number(s) in '", return_key, "': ",
    paste(unique(bad), collapse = ", "),
    " — rewrite..."
  )
  if (is.null(llm)) {
    return(list(ok = FALSE, text = text, rewritten = FALSE, reason = "b5_ungrounded"))
  }
  allow_u <- unique(allow)
  rules <- paste(
    "You correct scientific gene summaries under numeric grounding rules.",
    "EVERY number in the prose must appear in the allowed-numbers list",
    "(typically digits inside provided annotation labels/IDs).",
    "Do not invent quantitative or mathematical numbers",
    "(counts, scores, percentages, p-values, coordinates, distances,",
    "TPM/PIP/tau, EC numbers, version-like decimals not in the allowlist).",
    "Remove or rephrase claims that need disallowed numbers;",
    "use category words for magnitude.",
    "Return JSON only: {\"", return_key, "\": \"...\"}."
  )
  fix_user_extra <- paste0(
    "Allowed numbers (JSON array): ",
    jsonlite::toJSON(allow_u, auto_unbox = TRUE),
    "\nDisallowed numbers found: ",
    jsonlite::toJSON(unique(bad), auto_unbox = TRUE)
  )
  fix_sys <- paste(
    "You correct scientific gene summaries.",
    "Rewrite so EVERY remaining number appears in the allowed list",
    "(label/ID digits only).",
    "Do not invent quantitative numbers. Prefer category words. Return JSON only."
  )
  fix_user <- paste0(
    fix_user_extra,
    "\nOriginal text:\n", text,
    "\nReturn JSON: {\"", return_key, "\": \"...\"}"
  )
  message("B5: LLM rewrite for '", return_key, "'...")
  raw2 <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = fix_sys),
        list(role = "user", content = fix_user)
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = llm$timeout %||% 120
    ),
    error = function(e) e
  )
  if (inherits(raw2, "error")) {
    message("B5: rewrite LLM failed for '", return_key, "'")
    return(list(ok = FALSE, text = text, rewritten = TRUE, reason = "b5_rewrite_error"))
  }
  extract_key <- function(raw) {
    p <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(p) || is.null(p[[return_key]])) {
      return(NULL)
    }
    as.character(p[[return_key]])[1L]
  }
  text2 <- extract_key(raw2)
  if (is.null(text2) || !nzchar(text2)) {
    # Scrutiny: prior JSON + numeric rules
    message("B5: scrutiny (bad JSON after rewrite) for '", return_key, "'...")
    raw_s <- .llm_scrutiny_pass(
      llm = llm,
      rules = rules,
      prior_output = raw2,
      problems = paste0("missing JSON key '", return_key, "' or empty text"),
      extra_user = fix_user_extra,
      timeout = llm$timeout %||% 120
    )
    if (inherits(raw_s, "error")) {
      message("B5: rewrite-JSON scrutiny failed for '", return_key, "'")
      return(list(ok = FALSE, text = text, rewritten = TRUE, reason = "b5_rewrite_bad_json"))
    }
    text2 <- extract_key(raw_s)
    if (is.null(text2) || !nzchar(text2)) {
      message("B5: rewrite-JSON scrutiny still empty for '", return_key, "'")
      return(list(ok = FALSE, text = text, rewritten = TRUE, reason = "b5_rewrite_bad_json"))
    }
    raw2 <- raw_s
  }
  found2 <- .e_extract_numbers(text2)
  if (!length(found2) || all(.e_numbers_allowed(found2, allow))) {
    message("B5: rewrite OK for '", return_key, "'")
    return(list(ok = TRUE, text = text2, rewritten = TRUE, reason = NULL))
  }
  bad2 <- found2[!.e_numbers_allowed(found2, allow)]
  # One scrutiny pass with failed rewrite + allowlist rules
  message(
    "B5: scrutiny (numbers still bad) for '", return_key, "': ",
    paste(unique(bad2), collapse = ", "), "..."
  )
  raw3 <- .llm_scrutiny_pass(
    llm = llm,
    rules = rules,
    prior_output = jsonlite::toJSON(
      setNames(list(text2), return_key),
      auto_unbox = TRUE
    ),
    problems = paste0(
      "disallowed numbers remain: ",
      paste(unique(bad2), collapse = ", ")
    ),
    extra_user = paste0(
      "Allowed numbers (JSON array): ",
      jsonlite::toJSON(allow_u, auto_unbox = TRUE)
    ),
    timeout = llm$timeout %||% 120
  )
  if (inherits(raw3, "error")) {
    message("B5: final scrutiny LLM failed for '", return_key, "'")
    return(list(
      ok = FALSE, text = text, rewritten = TRUE,
      reason = "b5_recheck_failed", bad_after = bad2
    ))
  }
  text3 <- extract_key(raw3)
  if (is.null(text3) || !nzchar(text3)) {
    message("B5: final scrutiny empty for '", return_key, "'")
    return(list(
      ok = FALSE, text = text, rewritten = TRUE,
      reason = "b5_recheck_failed", bad_after = bad2
    ))
  }
  found3 <- .e_extract_numbers(text3)
  if (!length(found3) || all(.e_numbers_allowed(found3, allow))) {
    message("B5: final scrutiny OK for '", return_key, "'")
    return(list(ok = TRUE, text = text3, rewritten = TRUE, reason = NULL))
  }
  message(
    "B5: failed for '", return_key, "' after scrutiny; bad=",
    paste(unique(found3[!.e_numbers_allowed(found3, allow)]), collapse = ", ")
  )
  list(
    ok = FALSE,
    text = text,
    rewritten = TRUE,
    reason = "b5_recheck_failed",
    bad_after = found3[!.e_numbers_allowed(found3, allow)]
  )
}

#' One gene × one source Gemma call (JSON + B5; one scrutiny on bad schema)
#' @keywords internal
.e_call_source <- function(source,
                           gene_id,
                           payload,
                           system_prompt,
                           user_prompt,
                           return_key,
                           use_llm,
                           llm) {
  out <- list(
    gene_id = gene_id,
    source = source,
    ok = FALSE,
    skipped = FALSE,
    text = NULL,
    raw = NULL,
    error = NULL,
    b5 = NULL
  )
  if (!isTRUE(use_llm) || is.null(llm)) {
    out$text <- .e_template_text(source, gene_id, payload)
    out$error <- "llm_disabled"
    return(out)
  }
  message("E: LLM ", source, " for ", gene_id, "...")
  allow <- .e_collect_allowed_numbers(payload)
  raw <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = system_prompt),
        list(role = "user", content = user_prompt)
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = llm$timeout %||% 600
    ),
    error = function(e) e
  )
  if (inherits(raw, "error")) {
    message(
      "E: LLM ", source, " failed for ", gene_id, ": ",
      conditionMessage(raw), " — template fallback"
    )
    out$text <- .e_template_text(source, gene_id, payload)
    out$error <- conditionMessage(raw)
    return(out)
  }
  parse_text <- function(raw_s) {
    parsed <- tryCatch(
      jsonlite::fromJSON(raw_s, simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (is.null(parsed) || !is.list(parsed) || is.null(parsed[[return_key]]) ||
        !nzchar(as.character(parsed[[return_key]])[1L])) {
      return(NULL)
    }
    as.character(parsed[[return_key]])[1L]
  }
  text <- parse_text(raw)
  if (is.null(text)) {
    message("E: schema scrutiny for ", source, " / ", gene_id, "...")
    schema_rules <- paste(
      as.character(system_prompt)[1L],
      "Corrected output must be JSON with non-empty string key",
      paste0("\"", return_key, "\"."),
      "No other top-level keys required."
    )
    raw2 <- .llm_scrutiny_pass(
      llm = llm,
      rules = schema_rules,
      prior_output = raw,
      problems = paste0(
        "invalid JSON or missing/empty key '", return_key, "'"
      ),
      extra_user = as.character(user_prompt)[1L],
      timeout = llm$timeout %||% 600
    )
    if (!inherits(raw2, "error")) {
      text <- parse_text(raw2)
      if (!is.null(text)) {
        raw <- raw2
        message("E: schema scrutiny OK for ", source, " / ", gene_id)
      }
    }
  }
  if (is.null(text)) {
    out$text <- .e_template_text(source, gene_id, payload)
    out$error <- paste0("bad_json_or_missing_", return_key)
    out$raw <- raw
    message("E: bad JSON for ", source, " / ", gene_id, " — template fallback")
    return(out)
  }
  message("B5: check '", return_key, "' for ", gene_id, "...")
  b5 <- .e_b5_check_and_fix(llm, text, allow, return_key)
  out$b5 <- b5
  if (!isTRUE(b5$ok)) {
    out$text <- .e_template_text(source, gene_id, payload)
    out$error <- b5$reason %||% "b5_failed"
    out$raw <- raw
    message("B5: fail '", return_key, "' for ", gene_id, " — template fallback")
    return(out)
  }
  if (isTRUE(b5$rewritten)) {
    message("B5: pass '", return_key, "' for ", gene_id, " (rewritten)")
  } else {
    message("B5: pass '", return_key, "' for ", gene_id)
  }
  out$ok <- TRUE
  out$text <- b5$text
  out$raw <- raw
  out
}

# ---- E1 ----

#' @keywords internal
.e1_build_payload <- function(gene_id,
                              name = NULL,
                              trait_text = NULL,
                              free_text = character(),
                              keyword_details = NULL) {
  sentences <- .e_split_sentences(free_text)
  matched <- as.character(keyword_details$matched_keywords %||% character())
  syn <- keyword_details$matched_synonym %||% list()
  rel <- keyword_details$matched_related %||% list()
  unmatched <- as.character(keyword_details$unmatched_keywords %||% character())
  ctx <- as.character(keyword_details$matched_context %||% character())
  # sentence hit: any phrase from matched / synonym / related / context
  phrases <- unique(c(
    matched,
    vapply(syn, function(x) as.character(x$phrase %||% x[[1L]] %||% ""), character(1L)),
    vapply(rel, function(x) as.character(x$phrase %||% x[[1L]] %||% ""), character(1L)),
    ctx
  ))
  phrases <- phrases[nzchar(phrases)]
  hit_sent <- character()
  if (length(sentences) && length(phrases)) {
    for (s in sentences) {
      sl <- tolower(s)
      if (any(vapply(phrases, function(p) grepl(tolower(p), sl, fixed = TRUE), logical(1L)))) {
        hit_sent <- c(hit_sent, s)
      }
    }
  }
  payload <- list(
    trait_text = trait_text %||% "",
    Gene_ID = gene_id,
    Name = name %||% "",
    matched = matched,
    matched_synonym = syn,
    matched_related = rel,
    unmatched = unmatched,
    matched_context = ctx,
    matched_sentences = unique(hit_sent)
  )
  # drop empty arrays optionally kept for simplicity
  payload
}

#' @keywords internal
.e1_prompts <- function(payload) {
  sys <- paste(
    "You write concise English keyword/annotation evidence for one gene.",
    "linked/matched = direct hits; matched_related = indirect lexical links only.",
    "Use matched_sentences only; do not invent biology.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"keywords\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
    "\nWrite keywords prose; echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "keywords")
}

# ---- E2 ----

#' @keywords internal
.e2_impact_counts <- function(snpeff_gene) {
  counts <- c(HIGH = 0L, MODERATE = 0L, LOW = 0L, MODIFIER = 0L)
  if (is.null(snpeff_gene) || !nrow(snpeff_gene)) {
    return(list(counts = counts, worst = NA_character_, high_effects = list()))
  }
  # One ANN per genomic site (worst impact); no multi-transcript double-count.
  snpeff_gene <- .snpeff_worst_ann_per_site(snpeff_gene)
  impact <- toupper(as.character(
    snpeff_gene$Annotation_Impact %||% snpeff_gene$Impact %||% ""
  ))
  for (k in names(counts)) {
    counts[[k]] <- sum(impact == k, na.rm = TRUE)
  }
  order_imp <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  present <- order_imp[order_imp %in% impact]
  worst <- if (length(present)) present[1L] else NA_character_
  hi <- impact == "HIGH"
  high_effects <- list()
  if (any(hi)) {
    eff <- as.character(snpeff_gene$Annotation %||% snpeff_gene$Effect %||% "")
    hgvs <- character(nrow(snpeff_gene))
    if ("HGVS.p" %in% names(snpeff_gene)) {
      hgvs <- as.character(snpeff_gene[["HGVS.p"]])
    } else if ("HGVS_P" %in% names(snpeff_gene)) {
      hgvs <- as.character(snpeff_gene$HGVS_P)
    } else if ("hgvs_p" %in% names(snpeff_gene)) {
      hgvs <- as.character(snpeff_gene$hgvs_p)
    }
    for (i in which(hi)) {
      high_effects[[length(high_effects) + 1L]] <- list(
        effect = if (length(eff) >= i) eff[i] else "",
        hgvs_p = if (length(hgvs) >= i) hgvs[i] else ""
      )
    }
  }
  list(counts = counts, worst = worst, high_effects = high_effects)
}

#' @keywords internal
.e2_build_payload <- function(gene_id, name = NULL, snpeff_gene = NULL,
                              domain_overlaps = list(),
                              domain_id_match_rate = NA_real_) {
  ic <- .e2_impact_counts(snpeff_gene)
  list(
    Gene_ID = gene_id,
    Name = name %||% "",
    worst_impact = ic$worst,
    HIGH = unname(ic$counts[["HIGH"]]),
    MODERATE = unname(ic$counts[["MODERATE"]]),
    LOW = unname(ic$counts[["LOW"]]),
    MODIFIER = unname(ic$counts[["MODIFIER"]]),
    domain_overlaps = domain_overlaps,
    high_effects = ic$high_effects,
    domain_id_match_rate = domain_id_match_rate
  )
}

#' @keywords internal
.e2_prompts <- function(payload) {
  qual <- list(
    Gene_ID = payload$Gene_ID,
    Name = payload$Name %||% "",
    worst_impact = payload$worst_impact,
    HIGH = .e_count_band(payload$HIGH),
    MODERATE = .e_count_band(payload$MODERATE),
    LOW = .e_count_band(payload$LOW),
    MODIFIER = .e_count_band(payload$MODIFIER),
    high_effects = payload$high_effects,
    domain_overlaps = payload$domain_overlaps,
    domain_match = .e4_band(payload$domain_id_match_rate)
  )
  sys <- paste(
    "You summarize SnpEff variant effects for one gene.",
    "Use only categorical labels and effect names provided.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"snpeff\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite snpeff prose with category words only;",
    " echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "snpeff")
}

#' InterProScan rows + SnpEff AA overlap (same protein/transcript ID only)
#' @keywords internal
.e2_domain_overlaps <- function(snpeff_gene,
                                interpro_df,
                                protein_map,
                                gene_id) {
  empty <- list(
    domain_overlaps = list(),
    high_effects = list(),
    domain_id_match_rate = NA_real_,
    warning = NULL
  )
  if (!is.null(snpeff_gene) && nrow(snpeff_gene)) {
    snpeff_gene <- .snpeff_worst_ann_per_site(snpeff_gene)
  }
  ic_high <- .e2_impact_counts(snpeff_gene)$high_effects
  empty$high_effects <- ic_high

  if (is.null(interpro_df) || !nrow(interpro_df) || is.null(protein_map) ||
      !nrow(protein_map)) {
    return(empty)
  }
  gmap <- protein_map[as.character(protein_map$Gene_ID) == gene_id, , drop = FALSE]
  if (!nrow(gmap)) {
    return(empty)
  }
  pids <- unique(c(
    as.character(gmap$protein_id),
    as.character(gmap$transcript_id)
  ))
  pids <- pids[nzchar(pids)]

  nm <- names(interpro_df)
  low <- tolower(gsub("[ .]+", "_", nm))
  pid_col <- nm[match(c("protein_id", "proteinid", "seqid", "id"), low)][1L]
  st_col <- nm[match(c("start", "start_location", "locstart"), low)][1L]
  en_col <- nm[match(c("stop", "end", "stop_location", "locend"), low)][1L]
  desc_i <- which(low %in% c("description", "signature_description"))
  desc_col <- if (length(desc_i)) nm[desc_i[1L]] else NA_character_
  if (is.na(pid_col) || is.na(st_col) || is.na(en_col)) {
    empty$warning <- "interpro_missing_required_columns"
    return(empty)
  }
  ip <- interpro_df
  ip$._pid <- as.character(ip[[pid_col]])
  ip$._start <- as.numeric(ip[[st_col]])
  ip$._stop <- as.numeric(ip[[en_col]])
  ip$._name <- if (!is.na(desc_col)) as.character(ip[[desc_col]]) else ""
  ip <- ip[ip$._pid %in% pids, , drop = FALSE]
  if (!nrow(ip)) {
    return(empty)
  }

  aa_col <- NULL
  feat_col <- NULL
  if (!is.null(snpeff_gene) && nrow(snpeff_gene)) {
    for (cn in c("Pos.in.AA", "Pos.in.aa", "Protein_position", "aa_pos")) {
      if (cn %in% names(snpeff_gene)) {
        aa_col <- cn
        break
      }
    }
    for (cn in c("Feature_ID", "Feature", "transcript_id", "protein_id")) {
      if (cn %in% names(snpeff_gene)) {
        feat_col <- cn
        break
      }
    }
  }

  impact_rank <- function(x) {
    switch(
      toupper(as.character(x)[1L]),
      HIGH = 4L, MODERATE = 3L, LOW = 2L, MODIFIER = 1L, 0L
    )
  }
  rank_to_lab <- c("none", "MODIFIER", "LOW", "MODERATE", "HIGH")

  .parse_aa <- function(raw) {
    raw <- trimws(as.character(raw))
    if (!nzchar(raw) || raw %in% c(".", "-", "NA")) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(sub("/.*$", "", raw)))
  }

  rate <- NA_real_
  if (!is.null(snpeff_gene) && nrow(snpeff_gene) && !is.null(aa_col)) {
    elig <- integer()
    matched <- integer()
    for (vi in seq_len(nrow(snpeff_gene))) {
      aa <- .parse_aa(snpeff_gene[[aa_col]][vi])
      if (!is.finite(aa)) next
      elig <- c(elig, vi)
      feat <- if (!is.null(feat_col)) {
        as.character(snpeff_gene[[feat_col]][vi])
      } else {
        ""
      }
      if (nzchar(feat) && feat %in% ip$._pid) {
        matched <- c(matched, vi)
      }
    }
    if (length(elig)) {
      rate <- length(unique(matched)) / length(unique(elig))
    }
  }

  overlaps <- list()
  for (di in seq_len(nrow(ip))) {
    d_pid <- ip$._pid[di]
    d_start <- ip$._start[di]
    d_stop <- ip$._stop[di]
    best <- 0L
    if (!is.null(snpeff_gene) && nrow(snpeff_gene) && !is.null(aa_col)) {
      for (vi in seq_len(nrow(snpeff_gene))) {
        aa <- .parse_aa(snpeff_gene[[aa_col]][vi])
        if (!is.finite(aa)) next
        feat <- if (!is.null(feat_col)) {
          as.character(snpeff_gene[[feat_col]][vi])
        } else {
          ""
        }
        if (!identical(feat, d_pid)) next
        if (aa >= d_start && aa <= d_stop) {
          imp <- as.character(
            snpeff_gene$Annotation_Impact[vi] %||% snpeff_gene$Impact[vi] %||% ""
          )
          best <- max(best, impact_rank(imp))
        }
      }
    }
    if (best > 0L) {
      overlaps[[length(overlaps) + 1L]] <- list(
        name = ip$._name[di],
        aa_start = d_start,
        aa_end = d_stop,
        worst_overlapping_impact = rank_to_lab[best + 1L]
      )
    }
  }

  warn <- NULL
  if (is.finite(rate) && rate < 0.8) {
    warn <- paste0(
      "InterProScan ID match rate ", signif(rate, 3),
      " < 80% for gene ", gene_id
    )
    warning(warn, call. = FALSE)
  }
  list(
    domain_overlaps = overlaps,
    high_effects = ic_high,
    domain_id_match_rate = rate,
    warning = warn
  )
}

#' Backward-compatible wrapper
#' @keywords internal
.e2_interpro_overlap <- function(interpro_df, protein_ids) {
  if (is.null(interpro_df) || !nrow(interpro_df) || !length(protein_ids)) {
    return(list(rows = list(), match_rate = NA_real_))
  }
  ov <- .e2_domain_overlaps(
    snpeff_gene = NULL,
    interpro_df = interpro_df,
    protein_map = data.frame(
      Gene_ID = rep("x", length(protein_ids)),
      transcript_id = as.character(protein_ids),
      protein_id = as.character(protein_ids),
      stringsAsFactors = FALSE
    ),
    gene_id = "x"
  )
  list(rows = ov$domain_overlaps, match_rate = ov$domain_id_match_rate)
}

# ---- E3 ----

#' @keywords internal
.e3_build_payload <- function(gene_id,
                              name = NULL,
                              mapping_mode,
                              cand_row,
                              credible_set,
                              peak_id = NULL) {
  csinfo <- .e_cs_chr_minmax(credible_set)
  mode <- as.character(mapping_mode)[1L]
  pip_field <- if (identical(mode, "qtl")) {
    list(nearest_cs_PIP = cand_row$nearest_cs_PIP %||% NA_real_)
  } else {
    list(max_PIP = cand_row$max_PIP %||% NA_real_)
  }
  c(
    list(
      mode = mode,
      peak = list(
        peak_ID = peak_id,
        n_cs_variants = if (!is.null(csinfo)) csinfo$n else NA_integer_,
        cs_region = if (!is.null(csinfo) && identical(mode, "qtl")) {
          list(Chr = csinfo$Chr, start = csinfo$start, end = csinfo$end)
        } else {
          NULL
        }
      ),
      Gene_ID = gene_id,
      Name = name %||% "",
      dist2peak_bp = cand_row$dist2peak %||% NA_real_,
      negLog10P = cand_row$negLog10P %||% cand_row$neglog10p %||% NA_real_
    ),
    pip_field
  )
}

#' @keywords internal
.e3_prompts <- function(payload) {
  qual <- .e3_qualitative_payload(payload)
  mode <- qual$mode %||% "gwas"
  pip_note <- if (identical(mode, "qtl")) {
    paste(
      "Use nearest_cs_PIP band for QTL.",
      "PIP/CS here are Wakefield ABF approximations, not true causal",
      "probabilities (QTL-specific caveat)."
    )
  } else {
    "Use max_PIP band for GWAS when present."
  }
  sys <- paste(
    "You summarize mapping / credible-set context for one gene.",
    pip_note,
    "Use ONLY categorical labels in the payload.",
    .e_no_digits_prose_rule(),
    "Do not narrate CS resolution (reserved for validity).",
    "Return JSON only: {\"gwas\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite gwas prose with category words only;",
    " echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "gwas")
}

# ---- E4 ----

#' Yanai tissue-specificity tau (abundance)
#' @keywords internal
.yanai_tau <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  x <- pmax(x, 0)
  xmax <- max(x)
  if (xmax <= 0) {
    return(0)
  }
  sum(1 - x / xmax) / (length(x) - 1L)
}

#' @keywords internal
.e4_band <- function(x, hi = 0.8, mid = 0.4) {
  x <- as.numeric(x)[1L]
  if (length(x) != 1L || !is.finite(x)) {
    return("not_available")
  }
  if (x >= hi) "high" else if (x >= mid) "moderate" else "low"
}

#' Shared instruction: label digits OK; quantitative numbers forbidden
#'
#' Digits that appear inside provided annotation labels/IDs may be echoed.
#' Invented counts, scores, percentages, coordinates, EC-like decimals, etc.
#' are forbidden (avoids B5 rewrite loops).
#' @keywords internal
.e_no_digits_prose_rule <- function() {
  paste(
    "Numbers that appear only as part of provided annotation labels or IDs",
    "(term names, domain names, pathway names, gene symbols as given) may be",
    "echoed exactly as in the payload.",
    "Do NOT invent or write quantitative/mathematical numbers:",
    "counts, scores, percentages, p-values, coordinates, distances in bp,",
    "TPM, PIP, tau, EC numbers not present in the payload, or other",
    "version-like decimals used as quantities.",
    "For magnitude claims use category words only",
    "(e.g. high / moderate / low / near / moderate_distance / far / one / few / several / many)."
  )
}

#' @keywords internal
.e_count_band <- function(n) {
  n <- as.integer(n)[1L]
  if (!is.finite(n) || n <= 0L) {
    "none"
  } else if (n == 1L) {
    "one"
  } else if (n <= 3L) {
    "few"
  } else if (n <= 10L) {
    "several"
  } else {
    "many"
  }
}

#' @keywords internal
.e3_proximity_band <- function(dist_bp) {
  d <- as.numeric(dist_bp)[1L]
  if (!is.finite(d)) {
    return("not_available")
  }
  if (d <= 0) {
    "at_peak"
  } else if (d < 5000) {
    "very_near_peak"
  } else if (d < 50000) {
    "near_peak"
  } else if (d < 300000) {
    "moderate_distance"
  } else {
    "far_from_peak"
  }
}

#' @keywords internal
.e3_signal_band <- function(neglog10p) {
  x <- as.numeric(neglog10p)[1L]
  if (!is.finite(x)) {
    return("not_available")
  }
  if (x >= 8) {
    "very_strong"
  } else if (x >= 5) {
    "strong"
  } else if (x >= 3) {
    "moderate"
  } else {
    "weak"
  }
}

#' Qualitative GWAS/QTL payload for LLM (no raw metrics)
#' @keywords internal
.e3_qualitative_payload <- function(payload) {
  mode <- payload$mode %||% "gwas"
  pip_qual <- if (identical(mode, "qtl")) {
    list(nearest_cs_PIP = .e4_band(payload$nearest_cs_PIP))
  } else {
    list(max_PIP = .e4_band(payload$max_PIP))
  }
  n_cs <- payload$peak$n_cs_variants %||% NA_integer_
  list(
    mode = mode,
    peak = list(
      peak_ID = payload$peak$peak_ID,
      cs_variant_load = .e_count_band(n_cs),
      cs_region_present = isTRUE(!is.null(payload$peak$cs_region))
    ),
    Gene_ID = payload$Gene_ID,
    Name = payload$Name %||% "",
    proximity = .e3_proximity_band(payload$dist2peak_bp),
    association_signal = .e3_signal_band(payload$negLog10P),
    pip_qual
  )
}

#' Qualitative coded metrics for E8/E9 LLM payloads
#' @keywords internal
.e_coded_qualitative <- function(coded) {
  list(
    proximity = .e3_proximity_band(coded$dist2peak_bp),
    association_signal = .e3_signal_band(coded$negLog10P),
    max_PIP = .e4_band(coded$max_PIP),
    nearest_cs_PIP = .e4_band(coded$nearest_cs_PIP),
    worst_impact = coded$worst_impact,
    HIGH = .e_count_band(coded$HIGH),
    MODERATE = .e_count_band(coded$MODERATE),
    LOW = .e_count_band(coded$LOW),
    MODIFIER = .e_count_band(coded$MODIFIER),
    expression_specificity = coded$expression_specificity,
    expression_relative = coded$expression_relative,
    expression_high_tissues = {
      ht <- as.character(coded$expression_high_tissues %||% character())
      ht <- ht[nzchar(ht)]
      if (length(ht)) ht else NULL
    },
    domain_overlap = .e_count_band(coded$domain_overlap_n),
    domain_match = .e4_band(coded$domain_id_match_rate),
    domain_match_warning = coded$domain_match_warning
  )
}

#' @keywords internal
.e4_load_vocab <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "expression_coarse_vocab.yaml", package = "lazyGas")
  }
  if (!nzchar(path) || !file.exists(path)) {
    return(list(
      tissue = c("leaf", "root", "shoot", "seed", "flower", "endosperm"),
      stage = c("seedling", "vegetative", "reproductive", "ripening"),
      condition = c("control", "drought", "heat", "cold", "salt"),
      other = "other",
      not_available = c("not_available", "NA", "na", "")
    ))
  }
  yaml::read_yaml(path)
}

#' @keywords internal
.e4_map_sample_labels <- function(labels, vocab) {
  labels <- as.character(labels)
  out <- data.frame(
    sample = labels,
    tissue = "not_available",
    stage = "not_available",
    condition = "not_available",
    stringsAsFactors = FALSE
  )
  for (i in seq_along(labels)) {
    lab <- tolower(labels[i])
    if (!nzchar(lab) || lab %in% tolower(vocab$not_available %||% character())) {
      next
    }
    for (ax in c("tissue", "stage", "condition")) {
      terms <- as.character(vocab[[ax]] %||% character())
      hit <- terms[vapply(terms, function(t) {
        grepl(tolower(t), lab, fixed = TRUE)
      }, logical(1L))]
      if (length(hit)) {
        out[[ax]][i] <- hit[1L]
      } else if (identical(ax, "tissue")) {
        out[[ax]][i] <- "other"
      }
    }
  }
  out
}

#' @keywords internal
.e4_warn_row_z <- function(mat, lo = 0.8, hi = 1.2) {
  if (is.null(mat)) {
    return(NULL)
  }
  m <- as.matrix(mat)
  storage.mode(m) <- "numeric"
  rsd <- apply(m, 1L, stats::sd, na.rm = TRUE)
  med <- stats::median(rsd[is.finite(rsd)], na.rm = TRUE)
  if (is.finite(med) && med >= lo && med <= hi) {
    return(paste0(
      "Median row SD=", signif(med, 3),
      " is in [", lo, ",", hi, "] — input may be row-Z; tau requires abundance."
    ))
  }
  NULL
}

#' Tissue labels with highest mean abundance (names only; no values)
#' @keywords internal
.e4_high_tissues <- function(vals, frac = 0.7, max_n = 3L) {
  if (is.null(vals) || !length(vals)) {
    return(character())
  }
  v <- as.numeric(vals)
  names(v) <- names(vals)
  ok <- is.finite(v) & !is.na(names(v)) & nzchar(names(v))
  v <- v[ok]
  if (!length(v)) {
    return(character())
  }
  mx <- max(v, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) {
    return(character())
  }
  hit <- names(v)[v >= mx * frac]
  hit <- hit[order(v[hit], decreasing = TRUE)]
  utils::head(as.character(hit), as.integer(max_n)[1L])
}

#' @keywords internal
.e4_is_unavailable <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | !nzchar(x) | x %in% c("not_available", "na", "null", "nan")
}

#' @keywords internal
.e4_agg_means <- function(vals, labels) {
  vals <- as.numeric(vals)
  labels <- as.character(labels)
  ok <- is.finite(vals) & !is.na(labels) & nzchar(labels)
  if (!any(ok)) {
    return(numeric())
  }
  tapply(vals[ok], labels[ok], mean, na.rm = TRUE)
}

#' Median row SD (finite rows with >= 3 non-NA values and SD > 0)
#' @keywords internal
.e4_row_sd_median <- function(mat) {
  if (is.null(mat)) {
    return(NA_real_)
  }
  m <- as.matrix(mat)
  storage.mode(m) <- "numeric"
  rsd <- apply(m, 1L, function(r) {
    r <- r[is.finite(r)]
    if (length(r) < 3L) {
      return(NA_real_)
    }
    s <- stats::sd(r)
    if (!is.finite(s) || s <= 0) NA_real_ else s
  })
  stats::median(rsd[is.finite(rsd)], na.rm = TRUE)
}

#' Resolve gene biotype from cfg$gene_biotype
#' @keywords internal
.e4_gene_biotype <- function(gene_id, gene_biotype) {
  if (is.null(gene_biotype)) {
    return(NA_character_)
  }
  gid <- as.character(gene_id)[1L]
  if (is.data.frame(gene_biotype) &&
      all(c("Gene_ID", "biotype") %in% names(gene_biotype))) {
    hit <- as.character(
      gene_biotype$biotype[as.character(gene_biotype$Gene_ID) == gid]
    )
    return(if (length(hit)) hit[1L] else NA_character_)
  }
  if (!is.null(names(gene_biotype))) {
    v <- unname(as.character(gene_biotype[gid]))
    return(if (length(v) && nzchar(v[1L])) v[1L] else NA_character_)
  }
  NA_character_
}

#' @keywords internal
.e4_build_payload <- function(gene_id,
                              expr_mat = NULL,
                              sample_map = NULL,
                              cfg = list()) {
  if (is.null(expr_mat) || !gene_id %in% rownames(expr_mat)) {
    return(NULL)
  }
  tau_hi <- cfg$tau_high %||% 0.8
  tau_mid <- cfg$tau_moderate %||% 0.4
  rel_hi <- cfg$relative_high %||% 0.8
  rel_mid <- cfg$relative_moderate %||% 0.4
  matrix_scale <- as.character(cfg$matrix_scale %||% "abundance")[1L]
  warn <- .e4_warn_row_z(
    expr_mat,
    cfg$row_sd_z_warn_lo %||% 0.8,
    cfg$row_sd_z_warn_hi %||% 1.2
  )
  row_sd_med <- .e4_row_sd_median(expr_mat)

  null_spec <- list(
    combination_tau = NULL,
    combination_summary = NULL,
    tissue_tau = NULL,
    tissue_summary = NULL,
    stage_tau = NULL,
    stage_summary = NULL,
    condition_tau = NULL,
    condition_summary = NULL,
    summary = NULL,
    metric = "yanai_tau_combination"
  )
  null_rel <- list(
    percentile = NULL,
    value = NULL,
    band = NULL,
    summary = NULL,
    cohort_n = NULL,
    cohort_size = NULL,
    basis = NULL,
    metric = NULL,
    reason = "row_zscore_height_not_applicable"
  )

  if (identical(matrix_scale, "row_zscore") || identical(matrix_scale, "z")) {
    return(list(
      Gene_ID = gene_id,
      warning = warn %||%
        "matrix_scale=row_zscore; tau and relative height suppressed",
      matrix_scale = "row_zscore",
      row_sd_median = row_sd_med,
      specificity = null_spec,
      relative_level = null_rel,
      high_tissues = character(),
      tissues_by_level = character(),
      n_states_combination = 0L
    ))
  }

  row <- as.numeric(expr_mat[gene_id, , drop = TRUE])
  names(row) <- colnames(expr_mat)

  sm <- sample_map
  if (!is.null(sm) && nrow(sm)) {
    for (ax in c("tissue", "stage", "condition")) {
      if (!ax %in% names(sm)) {
        sm[[ax]] <- "not_available"
      }
    }
    idx <- match(names(row), sm$sample)
    tissue <- as.character(sm$tissue[idx])
    stage <- as.character(sm$stage[idx])
    condition <- as.character(sm$condition[idx])
    tissue[is.na(tissue) | !nzchar(tissue)] <- "not_available"
    stage[is.na(stage) | !nzchar(stage)] <- "not_available"
    condition[is.na(condition) | !nzchar(condition)] <- "not_available"

    query_axes <- cfg$query_axes
    if (is.null(query_axes) || !length(query_axes)) {
      query_axes <- "tissue"
      if (any(!.e4_is_unavailable(stage))) {
        query_axes <- c(query_axes, "stage")
      }
      if (any(!.e4_is_unavailable(condition))) {
        query_axes <- c(query_axes, "condition")
      }
    } else {
      query_axes <- as.character(query_axes)
    }
    keep_combo <- rep(TRUE, length(row))
    if ("tissue" %in% query_axes) {
      keep_combo <- keep_combo & !.e4_is_unavailable(tissue)
    }
    if ("stage" %in% query_axes) {
      keep_combo <- keep_combo & !.e4_is_unavailable(stage)
    }
    if ("condition" %in% query_axes) {
      keep_combo <- keep_combo & !.e4_is_unavailable(condition)
    }
    combo_lab <- paste(tissue, stage, condition, sep = "|")
    combo_means <- .e4_agg_means(row[keep_combo], combo_lab[keep_combo])
    tau_combo <- .yanai_tau(as.numeric(combo_means))

    # tissue single-axis: keep not_available as a level if present
    tissue_means <- .e4_agg_means(row, tissue)
    tau_tissue <- .yanai_tau(as.numeric(tissue_means))

    stage_lab <- ifelse(.e4_is_unavailable(stage), NA_character_, stage)
    stage_means <- .e4_agg_means(row, stage_lab)
    tau_stage <- .yanai_tau(as.numeric(stage_means))

    cond_lab <- ifelse(.e4_is_unavailable(condition), NA_character_, condition)
    cond_means <- .e4_agg_means(row, cond_lab)
    tau_condition <- .yanai_tau(as.numeric(cond_means))
  } else {
    combo_means <- numeric()
    tissue_means <- setNames(row, names(row))
    tau_combo <- .yanai_tau(row)
    tau_tissue <- tau_combo
    tau_stage <- NA_real_
    tau_condition <- NA_real_
  }

  combo_summary <- .e4_band(tau_combo, tau_hi, tau_mid)
  tissue_summary <- .e4_band(tau_tissue, tau_hi, tau_mid)
  stage_summary <- .e4_band(tau_stage, tau_hi, tau_mid)
  condition_summary <- .e4_band(tau_condition, tau_hi, tau_mid)
  # Primary specificity cue: combination tau when states exist, else tissue
  primary_summary <- if (length(combo_means) >= 2L) {
    combo_summary
  } else {
    tissue_summary
  }

  gene_mean <- mean(row, na.rm = TRUE)
  cohort <- rowMeans(expr_mat, na.rm = TRUE)
  names(cohort) <- rownames(expr_mat)
  expressed <- is.finite(cohort) & cohort > 0
  gene_bt <- .e4_gene_biotype(gene_id, cfg$gene_biotype)
  basis <- "expressed_genes"
  if (!is.na(gene_bt) && nzchar(gene_bt)) {
    all_bt <- vapply(names(cohort), function(g) {
      .e4_gene_biotype(g, cfg$gene_biotype)
    }, character(1L))
    same_bt <- !is.na(all_bt) & all_bt == gene_bt
    if (any(expressed & same_bt)) {
      expressed <- expressed & same_bt
      basis <- paste0("expressed_same_biotype:", gene_bt)
    }
  }
  cohort_vals <- cohort[expressed]
  pct <- if (length(cohort_vals) && is.finite(gene_mean)) {
    mean(cohort_vals <= gene_mean, na.rm = TRUE)
  } else {
    NA_real_
  }
  rel_band <- .e4_band(pct, rel_hi, rel_mid)

  high_tissues <- character()
  tissues_by_level <- character()
  if (length(tissue_means)) {
    # Prefer excluding not_available for display names
    disp <- tissue_means
    disp <- disp[!.e4_is_unavailable(names(disp))]
    if (!length(disp)) {
      disp <- tissue_means
    }
    ord <- order(as.numeric(disp), decreasing = TRUE, na.last = TRUE)
    tissues_by_level <- as.character(names(disp)[ord])
    if (!identical(primary_summary, "low") && !is.null(primary_summary)) {
      high_tissues <- .e4_high_tissues(
        disp,
        frac = cfg$high_tissue_frac %||% 0.7,
        max_n = cfg$high_tissue_n %||% 3L
      )
    }
  }

  list(
    Gene_ID = gene_id,
    warning = warn,
    matrix_scale = "abundance",
    row_sd_median = row_sd_med,
    specificity = list(
      combination_tau = tau_combo,
      combination_summary = combo_summary,
      tissue_tau = tau_tissue,
      tissue_summary = tissue_summary,
      stage_tau = tau_stage,
      stage_summary = stage_summary,
      condition_tau = tau_condition,
      condition_summary = condition_summary,
      summary = primary_summary,
      metric = "yanai_tau_combination"
    ),
    relative_level = list(
      percentile = pct,
      value = pct,
      band = rel_band,
      summary = rel_band,
      cohort_n = length(cohort_vals),
      cohort_size = length(cohort_vals),
      basis = basis,
      metric = "percentile_among_expressed_genes",
      reason = NULL
    ),
    high_tissues = high_tissues,
    tissues_by_level = tissues_by_level,
    n_states_combination = length(combo_means)
  )
}

#' @keywords internal
.e4_prompts <- function(payload) {
  # Qualitative-only prompt: omit raw tau / percentile so the model cannot
  # echo digits (avoids B5 rewrite loops on expression prose).
  sp <- payload$specificity$summary %||%
    payload$specificity$combination_summary %||%
    payload$specificity$tissue_summary %||% NA_character_
  ht <- as.character(payload$high_tissues %||% character())
  ht <- ht[nzchar(ht)]
  qual <- list(
    Gene_ID = payload$Gene_ID,
    warning = payload$warning,
    matrix_scale = payload$matrix_scale %||% "abundance",
    tissue_specificity = sp,
    combination_specificity = payload$specificity$combination_summary,
    tissue_axis_specificity = payload$specificity$tissue_summary,
    stage_specificity = payload$specificity$stage_summary,
    condition_specificity = payload$specificity$condition_summary,
    relative_expression = payload$relative_level$band %||%
      payload$relative_level$summary %||% NA_character_,
    high_expression_tissues = if (length(ht)) ht else NULL,
    tissues_by_level = payload$tissues_by_level %||% NULL
  )
  sys <- paste(
    "You summarize RNA-seq tissue specificity and relative expression for one gene.",
    "Use ONLY the categorical labels and tissue names in the payload",
    "(e.g. high / moderate / low).",
    "Prefer combination_specificity when present.",
    "If tissue_specificity is high or moderate AND high_expression_tissues is",
    "non-empty, you MUST name those tissues as where expression is highest.",
    "Do not invent tissue names absent from the payload.",
    "If relative_expression is null, do not compare height to other genes.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"expression\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite expression prose with category words and tissue names only;",
    " echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "expression")
}

# ---- E5 / E6 / E7 parsers ----

#' @keywords internal
.e_ann_cols_by_class <- function(ann_map, class) {
  if (is.null(ann_map)) {
    return(character())
  }
  if (is.data.frame(ann_map)) {
    return(as.character(ann_map$column[ann_map$class == class]))
  }
  cols <- character()
  for (row in ann_map) {
    if (is.list(row) && identical(as.character(row$class), class)) {
      cols <- c(cols, as.character(row$column))
    }
  }
  cols
}

#' @keywords internal
.e_cell_tokens <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) {
    return(character())
  }
  toks <- unlist(strsplit(paste(x, collapse = ";"), "[;|,]"), use.names = FALSE)
  trimws(toks[nzchar(trimws(toks))])
}

#' @keywords internal
.e5_build_payload <- function(gene_id,
                              name = NULL,
                              ann_row = NULL,
                              go_cols = character(),
                              cache_dir = NULL,
                              obo_path = NULL,
                              download = TRUE) {
  tokens <- character()
  if (!is.null(ann_row) && length(go_cols)) {
    for (cn in go_cols) {
      if (cn %in% names(ann_row)) {
        tokens <- c(tokens, .e_cell_tokens(ann_row[[cn]]))
      }
    }
  }
  terms <- .e5_resolve_term_names(
    tokens,
    cache_dir = cache_dir,
    obo_path = obo_path,
    download = download
  )
  if (!length(terms)) {
    return(NULL)
  }
  list(Gene_ID = gene_id, Name = name %||% "", terms = terms)
}

#' @keywords internal
.e5_prompts <- function(payload) {
  sys <- paste(
    "You summarize GO terms for one gene. Use only provided term names/IDs.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"go\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
    "\nWrite go prose; echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "go")
}

#' @keywords internal
.e6_build_payload <- function(gene_id,
                              name = NULL,
                              ann_row = NULL,
                              kegg_cols = character(),
                              cache_dir = NULL,
                              download = TRUE) {
  tokens <- character()
  if (!is.null(ann_row) && length(kegg_cols)) {
    for (cn in kegg_cols) {
      if (cn %in% names(ann_row)) {
        tokens <- c(tokens, .e_cell_tokens(ann_row[[cn]]))
      }
    }
  }
  entries <- .e6_resolve_entries(
    tokens,
    cache_dir = cache_dir,
    download = download
  )
  if (!length(entries)) {
    return(NULL)
  }
  list(Gene_ID = gene_id, Name = name %||% "", entries = entries)
}

#' @keywords internal
.e6_prompts <- function(payload) {
  sys <- paste(
    "You summarize KEGG pathway/KO/EC annotations for one gene.",
    "Use only provided entry names/IDs.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"kegg\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
    "\nWrite kegg prose; echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "kegg")
}

#' @keywords internal
.e7_build_payload <- function(gene_id,
                              name = NULL,
                              ann_row = NULL,
                              dom_cols = character(),
                              cache_dir = NULL,
                              pfam_path = NULL,
                              interpro_path = NULL,
                              download = TRUE) {
  tokens <- character()
  if (!is.null(ann_row) && length(dom_cols)) {
    for (cn in dom_cols) {
      if (cn %in% names(ann_row)) {
        tokens <- c(tokens, .e_cell_tokens(ann_row[[cn]]))
      }
    }
  }
  domains <- .e7_resolve_domains(
    tokens,
    cache_dir = cache_dir,
    pfam_path = pfam_path,
    interpro_path = interpro_path,
    download = download
  )
  if (!length(domains)) {
    return(NULL)
  }
  list(Gene_ID = gene_id, Name = name %||% "", domains = domains)
}

#' @keywords internal
.e7_prompts <- function(payload) {
  sys <- paste(
    "You summarize protein domain annotations for one gene.",
    "Use only provided domain names/sources.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"domains\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
    "\nWrite domains prose; echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "domains")
}

# ---- E8 / E9 ----

#' @keywords internal
.e8_build_payload <- function(gene_id, name = NULL, prose, coded) {
  list(
    Gene_ID = gene_id,
    Name = name %||% "",
    prose = prose,
    coded = coded
  )
}

#' @keywords internal
.e8_prompts <- function(payload) {
  qual <- list(
    Gene_ID = payload$Gene_ID,
    Name = payload$Name %||% "",
    prose = payload$prose,
    coded = .e_coded_qualitative(payload$coded)
  )
  sys <- paste(
    "You write a brief English summary of prior source prose for one gene.",
    "Prefer coded category labels over prose if they conflict. Do not invent facts.",
    "You MUST include one short clause that composite / ranking cues are",
    "peak-relative (within this peak), not absolute strength across peaks.",
    .e_no_digits_prose_rule(),
    "Do not write validity/caveats (that is a separate step).",
    "Return JSON only: {\"summary\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite summary with category words only;",
    " echo label digits only — no invented quantities.",
    " Mention that scores are peak-relative."
  )
  list(system = sys, user = user, return_key = "summary")
}

#' @keywords internal
.e9_build_payload <- function(gene_id,
                              name = NULL,
                              coded,
                              resolution,
                              n_cs_variants,
                              mapping_mode) {
  list(
    peak = list(
      n_cs_variants = n_cs_variants,
      resolution = as.character(resolution %||% "low")[1L],
      resolution_label = .e_resolution_label(resolution)
    ),
    Gene_ID = gene_id,
    Name = name %||% "",
    dist2peak_bp = coded$dist2peak_bp %||% NA_real_,
    negLog10P = coded$negLog10P %||% NA_real_,
    max_PIP = coded$max_PIP,
    nearest_cs_PIP = coded$nearest_cs_PIP,
    worst_impact = coded$worst_impact,
    HIGH = coded$HIGH %||% 0L,
    MODERATE = coded$MODERATE %||% 0L,
    LOW = coded$LOW %||% 0L,
    MODIFIER = coded$MODIFIER %||% 0L,
    domain_overlap_n = coded$domain_overlap_n %||% 0L,
    mapping_mode = mapping_mode
  )
}

#' @keywords internal
.e9_prompts <- function(payload) {
  qual <- list(
    peak = list(
      cs_variant_load = .e_count_band(payload$peak$n_cs_variants),
      resolution = payload$peak$resolution,
      resolution_label = payload$peak$resolution_label
    ),
    Gene_ID = payload$Gene_ID,
    Name = payload$Name %||% "",
    coded = .e_coded_qualitative(list(
      dist2peak_bp = payload$dist2peak_bp,
      negLog10P = payload$negLog10P,
      max_PIP = payload$max_PIP,
      nearest_cs_PIP = payload$nearest_cs_PIP,
      worst_impact = payload$worst_impact,
      HIGH = payload$HIGH,
      MODERATE = payload$MODERATE,
      LOW = payload$LOW,
      MODIFIER = payload$MODIFIER,
      expression_specificity = NULL,
      expression_relative = NULL,
      domain_overlap_n = payload$domain_overlap_n,
      domain_id_match_rate = NA_real_,
      domain_match_warning = NULL
    )),
    mapping_mode = payload$mapping_mode
  )
  sys <- paste(
    "You assess coded contradictions / gaps between SnpEff and mapping metrics",
    "for one gene. Do not use keyword or other prose. Narrate CS resolution here",
    "(resolution + resolution_label). Use ONLY categorical labels.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"validity\": \"...\"}."
  )
  if (identical(as.character(payload$mapping_mode %||% "gwas")[1L], "qtl")) {
    sys <- paste(
      sys,
      "If PIP/CS appear under QTL, note briefly they are Wakefield ABF",
      "approximations, not true causal probabilities."
    )
  }
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite validity with category words only;",
    " echo label digits only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "validity")
}

# ---- E10 interpretation (after Brief summary) ----

#' @keywords internal
.e10_build_payload <- function(gene_id,
                               name = NULL,
                               trait_text = NULL,
                               prose,
                               summary_text = NULL,
                               coded) {
  list(
    Gene_ID = gene_id,
    Name = name %||% "",
    trait_text = as.character(trait_text %||% "")[1L],
    summary = as.character(summary_text %||% "")[1L],
    prose = prose %||% list(),
    coded = coded
  )
}

#' Integrated free-form assessment grounded on E1–E8
#' @keywords internal
.e10_prompts <- function(payload) {
  qual <- list(
    Gene_ID = payload$Gene_ID,
    Name = payload$Name %||% "",
    trait_text = payload$trait_text %||% "",
    summary = payload$summary %||% "",
    prose = payload$prose %||% list(),
    coded = .e_coded_qualitative(payload$coded)
  )
  sys <- paste(
    "You are a plant genetics analyst. Using ONLY the provided brief summary,",
    "prior source prose, and coded category labels for one gene, write a free",
    "but careful interpretation of whether this gene looks related to the",
    "trait and whether it is a plausible causal / prioritization candidate.",
    "Integrate: keyword/function fit, expression specificity and profile",
    "(name high-expression tissues when provided),",
    "domain/structure evidence, and variant-effect impact when present.",
    "Base every claim on the supplied summary/prose/labels; do not invent",
    "pathways, orthologs, phenotypes, or numbers absent from the payload.",
    "A credible set does NOT prove the causal SNP or gene — say so when",
    "resolution is weak or mapping is multi-gene.",
    "End with exactly one line: Plausibility: <label>",
    "where <label> is one of: strong_candidate, plausible, weak,",
    "insufficient_evidence.",
    .e_no_digits_prose_rule(),
    "Return JSON only: {\"interpretation\": \"...\"}."
  )
  user <- paste0(
    "Payload JSON:\n",
    jsonlite::toJSON(qual, auto_unbox = TRUE, null = "null"),
    "\nWrite an integrated interpretation (several sentences OK);",
    " finish with Plausibility: <label>;",
    " use category words only — no invented quantities."
  )
  list(system = sys, user = user, return_key = "interpretation")
}

#' @keywords internal
.e_ann_row <- function(ann, gene_id) {
  if (is.null(ann) || !is.data.frame(ann) || !"Gene_ID" %in% names(ann)) {
    return(NULL)
  }
  hit <- ann[as.character(ann$Gene_ID) == gene_id, , drop = FALSE]
  if (!nrow(hit)) {
    return(NULL)
  }
  as.list(hit[1L, , drop = FALSE])
}

#' @keywords internal
.e_free_text_from_ann <- function(ann_row, free_cols) {
  if (is.null(ann_row) || !length(free_cols)) {
    return(character())
  }
  out <- character()
  for (cn in free_cols) {
    if (cn %in% names(ann_row)) {
      out <- c(out, as.character(ann_row[[cn]]))
    }
  }
  out
}

#' Coerce JSON-decoded keyword fields to character / linked-phrase lists
#' @keywords internal
.e_normalize_keyword_details <- function(det) {
  if (is.null(det) || !is.list(det)) {
    return(NULL)
  }
  as_chr <- function(x) {
    if (is.null(x)) {
      return(character())
    }
    if (is.character(x)) {
      return(x[nzchar(x) & !is.na(x)])
    }
    if (is.list(x)) {
      v <- vapply(x, function(el) {
        if (is.list(el) && !is.null(el$phrase)) {
          as.character(el$phrase)[1L]
        } else {
          as.character(el)[1L]
        }
      }, character(1L))
      return(v[nzchar(v) & !is.na(v)])
    }
    as.character(x)
  }
  list(
    keyword_score = as.numeric(det$keyword_score %||% 0)[1L],
    matched_keywords = as_chr(det$matched_keywords),
    matched_synonym = det$matched_synonym %||% list(),
    matched_related = det$matched_related %||% list(),
    unmatched_keywords = as_chr(det$unmatched_keywords),
    matched_context = as_chr(det$matched_context),
    llm_relevance_score = as.numeric(det$llm_relevance_score %||% 0)[1L],
    match_method = as.character(det$match_method %||% "none")[1L]
  )
}

#' Pull per-gene annotation keyword details from a rank table
#' @keywords internal
.keyword_hits_from_rank_result <- function(rank_df) {
  if (is.null(rank_df) || !nrow(rank_df) || !"Gene_ID" %in% names(rank_df) ||
      !"evidence_json" %in% names(rank_df)) {
    return(NULL)
  }
  out <- list()
  for (i in seq_len(nrow(rank_df))) {
    gid <- as.character(rank_df$Gene_ID[i])
    if (!nzchar(gid) || is.na(gid)) {
      next
    }
    ej <- tryCatch(
      jsonlite::fromJSON(
        as.character(rank_df$evidence_json[i]),
        simplifyVector = FALSE
      ),
      error = function(e) NULL
    )
    det <- .e_normalize_keyword_details(ej$annotation$details)
    if (!is.null(det)) {
      out[[gid]] <- det
    }
  }
  if (!length(out)) {
    return(NULL)
  }
  out
}

#' Run E1–E10 for one gene (skip empty optional sources)
#' @keywords internal
.report_e_gene_sources <- function(gene_id,
                                   name = NULL,
                                   mapping_mode,
                                   cand_row,
                                   ann,
                                   ann_map,
                                   credible_set,
                                   resolution,
                                   snpeff,
                                   protein_map,
                                   interpro = NULL,
                                   expr_mat = NULL,
                                   sample_map = NULL,
                                   trait_text = NULL,
                                   keyword_details = NULL,
                                   peak_id = NULL,
                                   use_llm = FALSE,
                                   llm = NULL,
                                   cfg = list()) {
  ann_row <- .e_ann_row(ann, gene_id)
  free_cols <- .e_ann_cols_by_class(ann_map, "free_description")
  go_cols <- .e_ann_cols_by_class(ann_map, "go")
  kegg_cols <- .e_ann_cols_by_class(ann_map, "kegg")
  dom_cols <- .e_ann_cols_by_class(ann_map, "pfam_interpro")

  results <- list()
  prose <- list(
    keywords = NULL, snpeff = NULL, gwas = NULL, expression = NULL,
    go = NULL, kegg = NULL, domains = NULL
  )

  call_ok <- function(source, payload, prompts, always = TRUE) {
    if (is.null(payload)) {
      results[[source]] <<- list(
        gene_id = gene_id, source = source, ok = FALSE,
        skipped = TRUE, text = NULL, error = "empty"
      )
      return(invisible(NULL))
    }
    # E2/E3: when use_llm=FALSE, plan says coded <ul> only (no prose)
    if (!isTRUE(use_llm) && source %in% c("snpeff", "gwas")) {
      results[[source]] <<- list(
        gene_id = gene_id, source = source, ok = FALSE,
        skipped = FALSE, text = NULL, error = "llm_disabled_no_prose"
      )
      return(invisible(NULL))
    }
    res <- .e_call_source(
      source = source,
      gene_id = gene_id,
      payload = payload,
      system_prompt = prompts$system,
      user_prompt = prompts$user,
      return_key = prompts$return_key,
      use_llm = use_llm,
      llm = llm
    )
    results[[source]] <<- res
    if (!is.null(res$text) && source %in% names(prose)) {
      prose[[source]] <<- res$text
    }
    invisible(res)
  }

  # E1
  free_text <- .e_free_text_from_ann(ann_row, free_cols)
  if (is.null(keyword_details) && length(free_text) &&
      any(nzchar(free_text)) && nzchar(as.character(trait_text %||% "")[1L])) {
    # Fallback when caller did not pass ranking keyword hits
    q_fb <- tryCatch(
      .phenotype_query_fallback(text = as.character(trait_text)[1L]),
      error = function(e) NULL
    )
    if (!is.null(q_fb)) {
      dets <- .annotation_keyword_family_details(
        text = paste(free_text[nzchar(free_text)], collapse = " | "),
        query = structure(
          list(
            trait_text = q_fb$trait_text %||% as.character(trait_text)[1L],
            trait_keywords = q_fb$trait_keywords,
            synonym_phrases = list(),
            related_phrases = list()
          ),
          class = c("PhenotypeQuery", "list")
        ),
        ignore.case = TRUE
      )
      if (length(dets)) {
        keyword_details <- dets[[1L]]
      }
    }
  }
  p1 <- .e1_build_payload(
    gene_id = gene_id,
    name = name,
    trait_text = trait_text,
    free_text = free_text,
    keyword_details = keyword_details
  )
  call_ok("keywords", p1, .e1_prompts(p1))

  # E2
  snp_g <- NULL
  if (!is.null(snpeff) && is.data.frame(snpeff) && nrow(snpeff) &&
      "Gene_ID" %in% names(snpeff)) {
    snp_g <- snpeff[as.character(snpeff$Gene_ID) == gene_id, , drop = FALSE]
  }
  id_dl <- isTRUE(cfg$id_map_download %||% TRUE)
  id_cache <- cfg$id_cache_dir %||% NULL
  ov <- .e2_domain_overlaps(
    snpeff_gene = snp_g,
    interpro_df = interpro,
    protein_map = protein_map,
    gene_id = gene_id
  )
  p2 <- .e2_build_payload(
    gene_id, name, snp_g,
    domain_overlaps = ov$domain_overlaps %||% list(),
    domain_id_match_rate = ov$domain_id_match_rate
  )
  # prefer HIGH effects from overlap helper (same as plan)
  if (length(ov$high_effects)) {
    p2$high_effects <- ov$high_effects
  }
  call_ok("snpeff", p2, .e2_prompts(p2))

  # E3
  p3 <- .e3_build_payload(
    gene_id, name, mapping_mode, cand_row, credible_set, peak_id
  )
  call_ok("gwas", p3, .e3_prompts(p3))

  # E4 (skip if no row)
  p4 <- .e4_build_payload(gene_id, expr_mat, sample_map, cfg = cfg)
  if (!is.null(p4)) {
    call_ok("expression", p4, .e4_prompts(p4))
  } else {
    results$expression <- list(
      gene_id = gene_id, source = "expression", skipped = TRUE, text = NULL
    )
  }

  # E5–E7 (skip empty)
  p5 <- .e5_build_payload(
    gene_id, name, ann_row, go_cols,
    cache_dir = id_cache,
    obo_path = cfg$go_obo_path %||% NULL,
    download = id_dl
  )
  if (!is.null(p5)) call_ok("go", p5, .e5_prompts(p5)) else {
    results$go <- list(skipped = TRUE, text = NULL, source = "go")
  }
  p6 <- .e6_build_payload(
    gene_id, name, ann_row, kegg_cols,
    cache_dir = id_cache,
    download = id_dl
  )
  if (!is.null(p6)) call_ok("kegg", p6, .e6_prompts(p6)) else {
    results$kegg <- list(skipped = TRUE, text = NULL, source = "kegg")
  }
  p7 <- .e7_build_payload(
    gene_id, name, ann_row, dom_cols,
    cache_dir = id_cache,
    pfam_path = cfg$pfam_path %||% NULL,
    interpro_path = cfg$interpro_path %||% NULL,
    download = id_dl
  )
  if (!is.null(p7)) call_ok("domains", p7, .e7_prompts(p7)) else {
    results$domains <- list(skipped = TRUE, text = NULL, source = "domains")
  }

  coded <- list(
    dist2peak_bp = p3$dist2peak_bp,
    negLog10P = p3$negLog10P,
    max_PIP = p3$max_PIP,
    nearest_cs_PIP = p3$nearest_cs_PIP,
    worst_impact = p2$worst_impact,
    HIGH = p2$HIGH,
    MODERATE = p2$MODERATE,
    LOW = p2$LOW,
    MODIFIER = p2$MODIFIER,
    expression_specificity = if (!is.null(p4)) {
      p4$specificity$summary %||%
        p4$specificity$combination_summary %||%
        p4$specificity$tissue_summary
    } else {
      NULL
    },
    expression_relative = if (!is.null(p4)) {
      p4$relative_level$band %||% p4$relative_level$summary
    } else {
      NULL
    },
    expression_high_tissues = if (!is.null(p4)) p4$high_tissues else character(),
    domain_overlap_n = length(ov$domain_overlaps %||% list()),
    domain_overlaps = ov$domain_overlaps %||% list(),
    high_effects = p2$high_effects %||% list(),
    domain_id_match_rate = ov$domain_id_match_rate,
    domain_match_warning = ov$warning
  )

  # E8
  p8 <- .e8_build_payload(gene_id, name, prose, coded)
  res8 <- .e_call_source(
    "summary", gene_id, p8,
    .e8_prompts(p8)$system, .e8_prompts(p8)$user, "summary",
    use_llm, llm
  )
  results$summary <- res8

  # E10 — free interpretation grounded on brief summary + prior prose
  p10 <- .e10_build_payload(
    gene_id = gene_id,
    name = name,
    trait_text = trait_text,
    prose = prose,
    summary_text = res8$text,
    coded = coded
  )
  res10 <- .e_call_source(
    "interpretation", gene_id, p10,
    .e10_prompts(p10)$system, .e10_prompts(p10)$user, "interpretation",
    use_llm, llm
  )
  results$interpretation <- res10

  # E9 — coded only (no E1–E8 prose)
  csinfo <- .e_cs_chr_minmax(credible_set)
  p9 <- .e9_build_payload(
    gene_id, name, coded, resolution,
    n_cs_variants = if (!is.null(csinfo)) csinfo$n else NA_integer_,
    mapping_mode = mapping_mode
  )
  res9 <- .e_call_source(
    "validity", gene_id, p9,
    .e9_prompts(p9)$system, .e9_prompts(p9)$user, "validity",
    use_llm, llm
  )
  results$validity <- res9

  list(results = results, prose = prose, coded = coded, payloads = list(
    e1 = p1, e2 = p2, e3 = p3, e4 = p4, e5 = p5, e6 = p6, e7 = p7,
    e8 = p8, e10 = p10, e9 = p9
  ))
}

#' Orchestrate E for one peak
#' @keywords internal
.report_e_run_peak <- function(work_dir,
                               mapping_mode,
                               simple_candidates,
                               credible_set,
                               gff,
                               ann = NULL,
                               snpeff = NULL,
                               interpro = NULL,
                               expr_mat = NULL,
                               trait_text = NULL,
                               keyword_hits_by_gene = NULL,
                               peak_id = NULL,
                               use_llm = FALSE,
                               llm = NULL,
                               resolution = NULL,
                               top_n = Inf,
                               cfg = list(),
                               gff_windows = NULL,
                               protein_map = NULL) {
  if (is.null(gff_windows) || is.null(protein_map)) {
    tables <- .gff_e_tables(gff, protein_map = protein_map)
    if (is.null(gff_windows)) {
      gff_windows <- tables$windows
    }
    if (is.null(protein_map) || !is.data.frame(protein_map) || !nrow(protein_map)) {
      protein_map <- tables$protein_map
    }
  }
  windows <- gff_windows
  genes <- .report_e_select_genes(
    mapping_mode = mapping_mode,
    simple_candidates = simple_candidates,
    credible_set = credible_set,
    gff_windows = windows,
    snpeff = snpeff,
    top_n = top_n
  )
  ann_map <- .read_ann_column_map(work_dir)
  if (is.null(ann_map) && !is.null(ann)) {
    ann_map <- classifyAnnColumns(ann, use_llm = FALSE)
    .write_ann_column_map(work_dir, ann_map)
  }
  sample_map <- NULL
  if (!is.null(expr_mat) && !is.null(colnames(expr_mat))) {
    sample_map <- .e4_map_sample_labels(colnames(expr_mat), .e4_load_vocab())
  }
  gene_results <- list()
  b5_failures <- list()
  n_genes <- nrow(genes)
  message("E: ", n_genes, " gene(s) for peak ", peak_id %||% "?")
  for (i in seq_len(n_genes)) {
    gid <- genes$Gene_ID[i]
    message("E: gene ", gid, " (", i, "/", n_genes, ")...")
    nm <- if ("Name" %in% names(genes)) as.character(genes$Name[i]) else NULL
    kh <- NULL
    if (!is.null(keyword_hits_by_gene) && gid %in% names(keyword_hits_by_gene)) {
      kh <- keyword_hits_by_gene[[gid]]
    }
    gr <- .report_e_gene_sources(
      gene_id = gid,
      name = nm,
      mapping_mode = mapping_mode,
      cand_row = genes[i, , drop = FALSE],
      ann = ann,
      ann_map = ann_map,
      credible_set = credible_set,
      resolution = resolution,
      snpeff = snpeff,
      protein_map = protein_map,
      interpro = interpro,
      expr_mat = expr_mat,
      sample_map = sample_map,
      trait_text = trait_text,
      keyword_details = kh,
      peak_id = peak_id,
      use_llm = use_llm,
      llm = llm,
      cfg = cfg
    )
    gene_results[[gid]] <- gr
    for (src in names(gr$results)) {
      b5 <- gr$results[[src]]$b5
      if (!is.null(b5) && isFALSE(b5$ok)) {
        b5_failures[[paste(gid, src, sep = "/")]] <- b5
      }
    }
  }
  message("E: peak ", peak_id %||% "?", " gene sources done")
  list(
    genes = genes,
    gene_results = gene_results,
    windows = windows,
    protein_map = protein_map,
    b5_failures = b5_failures
  )
}

#' HTML for E evidence (B2 section order)
#' @keywords internal
.report_evidence_html_e <- function(e_run, language = "en") {
  ui <- .report_i18n(language)
  if (is.null(e_run) || !nrow(e_run$genes)) {
    return(paste0("<p>", .report_escape_html("(no genes)"), "</p>\n"))
  }
  label_map <- c(
    keywords = ui$keywords,
    snpeff = ui$snpeff,
    gwas = ui$gwas_pos,
    expression = ui$expression,
    go = ui$go %||% "GO",
    kegg = ui$kegg %||% "KEGG",
    domains = ui$domains %||% "Domains",
    summary = ui$summary %||% "Brief summary",
    interpretation = ui$interpretation %||% "Integrated interpretation",
    validity = ui$validity
  )
  parts <- character()
  for (i in seq_len(nrow(e_run$genes))) {
    gid <- e_run$genes$Gene_ID[i]
    nm <- if ("Name" %in% names(e_run$genes)) {
      .report_name_cell(e_run$genes$Name[i])
    } else {
      ""
    }
    title <- if (nzchar(nm)) paste0(gid, " (", nm, ")") else gid
    gr <- e_run$gene_results[[gid]]
    if (is.null(gr)) next
    blocks <- character()
    # coded GWAS / SnpEff first (always)
    coded <- gr$coded
    gwas_items <- c(
      sprintf("dist2peak_bp: %s", coded$dist2peak_bp %||% "NA"),
      if (!is.null(coded$max_PIP)) sprintf("max_PIP: %s", coded$max_PIP),
      if (!is.null(coded$nearest_cs_PIP)) {
        sprintf("nearest_cs_PIP: %s", coded$nearest_cs_PIP)
      }
    )
    blocks <- c(
      blocks,
      paste0(
        "<h5>", .report_escape_html(ui$gwas_pos), "</h5>\n<ul>\n",
        paste0("<li>", .report_escape_html(gwas_items), "</li>\n", collapse = ""),
        "</ul>\n"
      )
    )
    if (!is.null(gr$results$gwas$text) && nzchar(gr$results$gwas$text)) {
      blocks <- c(
        blocks,
        paste0("<p>", .report_escape_html(gr$results$gwas$text), "</p>\n")
      )
    }
    blocks <- c(
      blocks,
      paste0(
        "<h5>", .report_escape_html(ui$snpeff), "</h5>\n<ul>\n",
        sprintf(
          "<li>worst_impact: %s; HIGH=%s MODERATE=%s LOW=%s MODIFIER=%s</li>\n",
          .report_escape_html(coded$worst_impact %||% "NA"),
          coded$HIGH %||% 0L, coded$MODERATE %||% 0L,
          coded$LOW %||% 0L, coded$MODIFIER %||% 0L
        ),
        "</ul>\n"
      )
    )
    if (length(coded$domain_overlaps) || length(coded$high_effects)) {
      dom_lis <- character()
      if (length(coded$domain_overlaps)) {
        for (d in coded$domain_overlaps) {
          imp <- as.character(d$worst_overlapping_impact %||% "")[1L]
          if (!nzchar(imp) || identical(tolower(imp), "none")) {
            next
          }
          dom_lis <- c(
            dom_lis,
            sprintf(
              "<li>%s [%s–%s]: %s</li>\n",
              .report_escape_html(d$name %||% ""),
              d$aa_start %||% "?",
              d$aa_end %||% "?",
              .report_escape_html(imp)
            )
          )
        }
      }
      hi_lis <- character()
      if (length(coded$high_effects)) {
        hi_lis <- vapply(coded$high_effects, function(h) {
          sprintf(
            "<li>HIGH: %s %s</li>\n",
            .report_escape_html(h$effect %||% ""),
            .report_escape_html(h$hgvs_p %||% "")
          )
        }, character(1L))
      }
      if (length(dom_lis) || length(hi_lis)) {
        blocks <- c(
          blocks,
          paste0(
            "<ul>\n",
            paste(dom_lis, collapse = ""),
            paste(hi_lis, collapse = ""),
            "</ul>\n"
          )
        )
      }
    }
    if (!is.null(gr$results$snpeff$text) && nzchar(gr$results$snpeff$text)) {
      blocks <- c(
        blocks,
        paste0("<p>", .report_escape_html(gr$results$snpeff$text), "</p>\n")
      )
    }
    for (src in c(
      "keywords", "expression", "go", "kegg", "domains",
      "summary", "interpretation", "validity"
    )) {
      res <- gr$results[[src]]
      if (is.null(res) || isTRUE(res$skipped) || is.null(res$text) || !nzchar(res$text)) {
        next
      }
      lab <- unname(label_map[src]) %||% src
      blocks <- c(
        blocks,
        paste0(
          "<h5>", .report_escape_html(lab), "</h5>\n",
          "<p>", .report_escape_html(res$text), "</p>\n"
        )
      )
    }
    parts <- c(
      parts,
      paste0("<h4>", .report_escape_html(title), "</h4>\n", paste(blocks, collapse = ""))
    )
  }
  if (length(e_run$b5_failures)) {
    parts <- c(
      parts,
      paste0(
        "<p><strong>Numeric grounding failures (B5)</strong>: ",
        .report_escape_html(paste(names(e_run$b5_failures), collapse = ", ")),
        " — B5 checks digits only; non-numeric claims are not verified.",
        "</p>\n"
      )
    )
  }
  paste(parts, collapse = "\n")
}

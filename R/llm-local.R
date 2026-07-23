################################################################################
#' Check connectivity to a local Ollama server
#'
#' @param base_url Ollama base URL.
#' @param timeout Connection timeout in seconds.
#' @return Logical scalar.
#' @export
llmHealthCheck <- function(base_url = NULL, timeout = 5) {
  base_url <- .llm_base_url(base_url)
  resp <- .llm_http_get(
    path = "/api/tags",
    base_url = base_url,
    timeout = timeout
  )
  isTRUE(!is.null(resp))
}

#' Send a chat completion request to a local Ollama server
#'
#' @param messages List of lists with \code{role} and \code{content}.
#' @param model Model name (default from \code{LAZYGAS_LLM_MODEL}).
#' @param base_url Ollama base URL.
#' @param json_mode If \code{TRUE}, request JSON-formatted output.
#' @param timeout Request timeout in seconds.
#'
#' @return Assistant message content as a character string.
#' @export
llmChat <- function(messages,
                    model = NULL,
                    base_url = NULL,
                    json_mode = FALSE,
                    timeout = 120) {
  if (!is.list(messages) || !length(messages)) {
    stop("'messages' must be a non-empty list.", call. = FALSE)
  }
  model <- .llm_model(model)
  base_url <- .llm_base_url(base_url)

  body <- list(
    model = model,
    messages = messages,
    stream = FALSE
  )
  if (isTRUE(json_mode)) {
    body$format <- "json"
  }

  resp <- .llm_http_post(
    path = "/api/chat",
    body = body,
    base_url = base_url,
    timeout = timeout
  )
  if (is.null(resp) || is.null(resp$message$content)) {
    stop("Empty response from local LLM.", call. = FALSE)
  }
  as.character(resp$message$content)
}

#' Explain ranked phenotype candidates using a local LLM
#'
#' Generates a Markdown summary of the top-ranked genes using only evidence
#' collected during ranking. When \code{use_llm = FALSE}, a template-based
#' report is returned.
#'
#' For Phase 1, prefer the public entry point [llm_report()], which ranks
#' candidates, verifies claims, and saves under \code{lazygas/ai_report/}.
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param query A \code{PhenotypeQuery} object (optional; inferred from
#'   \code{rank_result$query_id} when \code{object} is given).
#' @param object \code{LazyGas} object to load saved queries.
#' @param top_n Number of top genes to include in the explanation.
#' @param use_llm Use local LLM when available.
#' @param language Response language: \code{"ja"} or \code{"en"}.
#' @param model,base_url LLM settings passed to [llmChat()].
#' @param timeout Request timeout in seconds (default 600).
#'
#' @return A character string (Markdown).
#' @export
explainPhenotypeCandidates <- function(rank_result,
                                       query = NULL,
                                       object = NULL,
                                       top_n = 10L,
                                       use_llm = TRUE,
                                       language = c("ja", "en"),
                                       model = NULL,
                                       base_url = NULL,
                                       timeout = NULL) {
  language <- match.arg(language)
  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }

  meta <- attr(rank_result, "phenotypeRank")
  if (is.null(query) && !is.null(meta$query_id) && !is.null(object)) {
    query <- .store_read_phenotype_query(object, meta$query_id)
  }
  if (is.null(query)) {
    query <- list(
      trait_text = "phenotype query",
      query_id = if (!is.null(meta)) meta$query_id else "unknown"
    )
    class(query) <- c("PhenotypeQuery", "list")
  }

  llm <- .ai_config_llm_settings(model = model, base_url = base_url, timeout = timeout)
  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)

  if (!isTRUE(use_llm)) {
    return(.explain_phenotype_template(
      query = query,
      evidence_bundle = evidence_bundle,
      language = language
    ))
  }

  if (!llmHealthCheck(base_url = llm$base_url, timeout = 5)) {
    stop(
      "Local LLM is unreachable at ", llm$base_url,
      ". Start the SSH tunnel / Ollama, or set use_llm = FALSE for the template path.",
      call. = FALSE
    )
  }

  lang_note <- if (language == "ja") {
    "回答は日本語で書いてください。"
  } else {
    "Write the response in English."
  }

  system_msg <- paste(
    "You are a plant genetics assistant.",
    "Use ONLY the JSON evidence provided.",
    "Do not recommend genes not in the list.",
    "Cite evidence snippets when explaining each gene.",
    "Always include: (1) keyword/annotation match notes,",
    "(2) SnpEff impact notes when present,",
    "(3) a short validity comment for each gene.",
    "For keyword/annotation notes:",
    "- State whether keyword match is high, low, or none (qualitative only; never print numeric scores).",
    "- List which phenotype keywords matched and which did not, using annotation.details.matched_keywords / unmatched_keywords when present.",
    "- Do NOT say that a 'semantic match' was confirmed or cite LSA/embedding similarity.",
    "- Ignore matches of prepositions, conjunctions, and other non-biological function words (e.g. in, of, to, and); do not treat them as evidence.",
    lang_note
  )

  user_msg <- jsonlite::toJSON(
    list(
      phenotype_query = unclass(query),
      ranked_genes = evidence_bundle
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )

  llmChat(
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    model = llm$model,
    base_url = llm$base_url,
    timeout = llm$timeout
  )
}

#' Generate a Phase 1 AI candidate-gene report
#'
#' Public Phase 1 entry point: run phenotype query structuring (optional LLM),
#' rank candidates with Phase 1 weights/sources, render a Markdown report
#' (keyword + SnpEff + validity comments), verify claims against evidence JSON
#' (unsupported claims are **marked**, not removed), and save under
#' \code{lazygas/ai_report/}.
#'
#' When \code{use_llm = TRUE}, an unreachable LLM raises an error. Use
#' \code{use_llm = FALSE} for the keyword/template path.
#'
#' @param object A \code{LazyGas} object with candidate genes.
#' @param pheno Phenotype name or index.
#' @param query Character phenotype description or [PhenotypeQuery]. If
#'   \code{NULL}, uses \code{pheno} name as the query text.
#' @param language \code{"ja"} or \code{"en"}.
#' @param top_n Number of top genes (default 10).
#' @param sources,weights Ranking channels (defaults: Phase 1).
#' @param use_llm Use local LLM for the report body (and optionally query parse).
#' @param query_use_llm Pass \code{use_llm} to [phenotypeQuery()] when
#'   \code{query} is character.
#' @param model,base_url,timeout LLM connection (default timeout 600s).
#' @param save Write Markdown + verification JSON under [aiReportDir()].
#' @param out_dir Optional report directory override.
#' @param rank_result Optional precomputed [rankPhenotypeCandidates()] output;
#'   skips re-ranking when provided (single-section report).
#' @param peak_id Optional peak ID(s) to report. When \code{rank_result} is
#'   \code{NULL}, candidates are ranked within each peak; default \code{NULL}
#'   uses the top five peaks by lead \code{negLog10P}.
#' @param use_llm_relevance Annotation LLM relevance during ranking.
#' @param download_tenor Download missing TENOR expression CSVs.
#' @param candidate Optional candidate \code{data.frame} override (must include
#'   \code{Gene_ID} / peak columns as produced by \code{listCandidate}). When
#'   \code{NULL}, candidates are loaded from the companion store.
#' @param ... Passed to [rankPhenotypeCandidates()] when ranking.
#'
#' @return A list with \code{markdown}, \code{path}, \code{meta_path},
#'   \code{verification}, \code{rank_result}, and \code{query}.
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [explainPhenotypeCandidates()]
llm_report <- function(object,
                       pheno,
                       query = NULL,
                       language = c("ja", "en"),
                       top_n = 10L,
                       sources = phase1DefaultSources(),
                       weights = phase1DefaultWeights(),
                       use_llm = TRUE,
                       query_use_llm = TRUE,
                       model = NULL,
                       base_url = NULL,
                       timeout = NULL,
                       save = TRUE,
                       out_dir = NULL,
                       rank_result = NULL,
                       peak_id = NULL,
                       use_llm_relevance = NULL,
                       download_tenor = TRUE,
                       candidate = NULL,
                       ...) {
  language <- match.arg(language)
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }

  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  llm <- .ai_config_llm_settings(model = model, base_url = base_url, timeout = timeout)

  if (is.null(query)) {
    query <- pheno_name
  }
  if (is.character(query)) {
    query <- phenotypeQuery(
      text = query,
      use_llm = isTRUE(query_use_llm) && isTRUE(use_llm),
      llm_model = llm$model,
      llm_base_url = llm$base_url,
      object = object,
      save = TRUE
    )
  } else if (!inherits(query, "PhenotypeQuery")) {
    stop("'query' must be NULL, character, or PhenotypeQuery.", call. = FALSE)
  }

  if (is.null(rank_result)) {
    peaks_plan <- .resolve_report_peaks(
      object = object,
      pheno_name = pheno_name,
      peak_id = peak_id
    )
    candidate_all <- if (!is.null(candidate)) {
      candidate
    } else {
      lazyData(
        object = object,
        dataset = "candidate",
        pheno = pheno_name
      )
    }
    if (is.null(candidate_all) || nrow(candidate_all) == 0L) {
      stop("No candidate data found for phenotype '", pheno_name, "'.",
           call. = FALSE)
    }
    rank_top <- max(as.integer(top_n)[1L], 30L)
    rank_args <- list(
      object = object,
      pheno = pheno_name,
      query = query,
      sources = sources,
      weights = weights,
      save = FALSE,
      use_llm_relevance = use_llm_relevance,
      download_tenor = download_tenor,
      query_use_llm = FALSE,
      llm_model = llm$model,
      llm_base_url = llm$base_url,
      llm_timeout = llm$timeout
    )
    rank_args <- c(rank_args, list(...))
    sections <- character()
    rank_parts <- list()
    for (i in seq_len(nrow(peaks_plan))) {
      pid <- peaks_plan$peak_ID[i]
      chr_lab <- peaks_plan$Chr[i]
      pos_lab <- peaks_plan$Pos[i]
      pos_fmt <- format(as.numeric(pos_lab), big.mark = ",", scientific = FALSE,
                        trim = TRUE)
      header <- paste0("## Peak ", pid, " — ", chr_lab, ":", pos_fmt)
      cand_peak <- candidate_all[
        candidate_all$peak_ID == pid,
        ,
        drop = FALSE
      ]
      if (nrow(cand_peak) == 0L) {
        sections <- c(
          sections,
          paste0(header, "\n\n", "No candidate genes listed for this peak.")
        )
        next
      }
      if ("Gene_ID" %in% names(cand_peak)) {
        gid_ok <- !is.na(cand_peak$Gene_ID) & nzchar(as.character(cand_peak$Gene_ID))
        if (!any(gid_ok)) {
          sections <- c(
            sections,
            paste0(
              header, "\n\n",
              "No candidate genes with a usable Gene_ID for this peak."
            )
          )
          next
        }
        cand_peak <- cand_peak[gid_ok, , drop = FALSE]
      }
      rank_peak <- do.call(
        rankPhenotypeCandidates,
        c(
          rank_args,
          list(
            candidate = cand_peak,
            top_n = rank_top
          )
        )
      )
      rank_parts[[length(rank_parts) + 1L]] <- rank_peak
      body <- explainPhenotypeCandidates(
        rank_result = rank_peak,
        query = query,
        object = object,
        top_n = top_n,
        use_llm = use_llm,
        language = language,
        model = llm$model,
        base_url = llm$base_url,
        timeout = llm$timeout
      )
      sections <- c(sections, paste0(header, "\n\n", body))
    }
    markdown_body <- paste(sections, collapse = "\n\n")
    rank_result <- if (length(rank_parts)) {
      do.call(rbind, rank_parts)
    } else {
      NULL
    }
  } else {
    markdown_body <- explainPhenotypeCandidates(
      rank_result = rank_result,
      query = query,
      object = object,
      top_n = top_n,
      use_llm = use_llm,
      language = language,
      model = llm$model,
      base_url = llm$base_url,
      timeout = llm$timeout
    )
  }

  verification <- if (!is.null(rank_result) && nrow(rank_result) > 0L) {
    .verify_report_against_evidence(
      markdown = markdown_body,
      rank_result = rank_result,
      top_n = top_n,
      language = language
    )
  } else {
    list(
      n_checked = 0L,
      n_supported = 0L,
      n_marked = 0L,
      claims = list()
    )
  }

  markdown <- .append_verification_section(
    markdown = markdown_body,
    verification = verification,
    language = language
  )

  report_path <- NULL
  meta_path <- NULL
  if (isTRUE(save)) {
    report_dir <- aiReportDir(object = object, out_dir = out_dir)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    safe_pheno <- .store_safe_name(pheno_name)
    base <- paste0(safe_pheno, "_", stamp)
    report_path <- file.path(report_dir, paste0(base, ".md"))
    meta_path <- file.path(report_dir, paste0(base, ".meta.json"))
    writeLines(markdown, report_path, useBytes = TRUE)
    jsonlite::write_json(
      list(
        pheno = pheno_name,
        query = unclass(query),
        peak_id = peak_id,
        language = language,
        top_n = as.integer(top_n)[1L],
        sources = sources,
        weights = as.list(weights),
        use_llm = isTRUE(use_llm),
        model = llm$model,
        base_url = llm$base_url,
        created_at = stamp,
        report_file = basename(report_path),
        verification = verification,
        rank_meta = attr(rank_result, "phenotypeRank")
      ),
      meta_path,
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null"
    )
  }

  list(
    markdown = markdown,
    path = report_path,
    meta_path = meta_path,
    verification = verification,
    rank_result = rank_result,
    query = query
  )
}

#' Resolve peak table rows for per-peak LLM reports
#' @keywords internal
.resolve_report_peaks <- function(object, pheno_name, peak_id = NULL) {
  peakcall <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = TRUE)
  if (is.null(peakcall) || nrow(peakcall) == 0L) {
    peakcall <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = FALSE)
  }
  if (is.null(peakcall) || nrow(peakcall) == 0L) {
    stop("No peakcall data for phenotype '", pheno_name, "'.", call. = FALSE)
  }

  leads <- peakcall[
    peakcall$peak_variant_ID == peakcall$variant_ID,
    ,
    drop = FALSE
  ]
  if (nrow(leads) == 0L) {
    leads <- do.call(rbind, lapply(split(peakcall, peakcall$peak_ID), function(block) {
      block[which.max(block$negLog10P), , drop = FALSE]
    }))
  }
  ord <- order(-leads$negLog10P, leads$peak_ID)
  leads <- leads[ord, , drop = FALSE]
  leads <- leads[!duplicated(leads$peak_ID), , drop = FALSE]
  out <- data.frame(
    peak_ID = leads$peak_ID,
    Chr = leads$Chr,
    Pos = leads$Pos,
    stringsAsFactors = FALSE
  )

  if (is.null(peak_id)) {
    out <- out[seq_len(min(5L, nrow(out))), , drop = FALSE]
    return(out)
  }

  peak_id <- as.character(peak_id)
  avail <- as.character(out$peak_ID)
  miss <- setdiff(peak_id, avail)
  if (length(miss)) {
    stop(
      "peak_id not found for phenotype '", pheno_name, "': ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  out <- out[match(peak_id, avail), , drop = FALSE]
  out
}

.verify_report_against_evidence <- function(markdown,
                                            rank_result,
                                            top_n = 10L,
                                            language = "ja") {
  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  bundle <- .phenotype_evidence_bundle(sub)

  # Collect searchable evidence tokens per gene
  evidence_text <- setNames(character(nrow(sub)), sub$Gene_ID)
  for (i in seq_along(bundle)) {
    item <- bundle[[i]]
    bits <- character()
    for (src in names(item$evidence %||% list())) {
      block <- item$evidence[[src]]
      if (!is.null(block$snippets)) {
        bits <- c(bits, as.character(unlist(block$snippets)))
      }
      if (!is.null(block$details)) {
        bits <- c(bits, as.character(unlist(block$details)))
      }
    }
    evidence_text[[item$Gene_ID]] <- tolower(paste(bits, collapse = " "))
  }

  # Sentence-like claims
  text <- gsub("\r", "", markdown)
  parts <- unlist(strsplit(text, "(?<=[。．\\.\\!\\?\\n])\\s*", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & nchar(parts) > 12L]

  gene_ids <- as.character(sub$Gene_ID)
  claims <- list()
  n_checked <- 0L
  n_supported <- 0L
  n_marked <- 0L

  for (sent in parts) {
    hit_genes <- gene_ids[vapply(gene_ids, function(g) {
      grepl(g, sent, fixed = TRUE)
    }, logical(1))]
    if (!length(hit_genes)) {
      next
    }
    n_checked <- n_checked + 1L
    supported <- FALSE
    for (g in hit_genes) {
      ev <- evidence_text[[g]]
      # Token overlap: any alphanumeric token length >= 4 from claim in evidence
      toks <- unique(unlist(regmatches(
        tolower(sent),
        gregexpr("[a-z0-9_]{4,}", tolower(sent), perl = TRUE)
      )))
      toks <- setdiff(toks, tolower(gene_ids))
      if (!length(toks)) {
        # Gene mentioned with score-like context counts as supported if gene in rank
        supported <- TRUE
        break
      }
      if (any(vapply(toks, function(t) grepl(t, ev, fixed = TRUE), logical(1)))) {
        supported <- TRUE
        break
      }
    }
    if (supported) {
      n_supported <- n_supported + 1L
      claims[[length(claims) + 1L]] <- list(
        text = sent,
        genes = hit_genes,
        status = "supported"
      )
    } else {
      n_marked <- n_marked + 1L
      claims[[length(claims) + 1L]] <- list(
        text = sent,
        genes = hit_genes,
        status = "unsupported_marked"
      )
    }
  }

  congruence <- if (n_checked > 0L) n_supported / n_checked else NA_real_
  list(
    n_claims_checked = n_checked,
    n_supported = n_supported,
    n_unsupported_marked = n_marked,
    congruence_rate = congruence,
    claims = claims
  )
}

.append_verification_section <- function(markdown, verification, language = "ja") {
  rate <- verification$congruence_rate
  rate_txt <- if (is.finite(rate)) sprintf("%.1f%%", 100 * rate) else "NA"
  marked <- Filter(function(x) identical(x$status, "unsupported_marked"),
                   verification$claims %||% list())

  if (language == "ja") {
    hdr <- c(
      "",
      "---",
      "",
      "## 根拠検証",
      "",
      sprintf("- 照合クレーム数: %s", verification$n_claims_checked),
      sprintf("- 支持されたクレーム: %s", verification$n_supported),
      sprintf("- 根拠なし（マークのみ）: %s", verification$n_unsupported_marked),
      sprintf("- 整合率: %s", rate_txt),
      ""
    )
    if (length(marked)) {
      hdr <- c(hdr, "### 根拠なしとしてマークした記述", "")
      for (m in marked) {
        hdr <- c(hdr, paste0("- ⚠️ ", m$text))
      }
      hdr <- c(hdr, "")
    }
  } else {
    hdr <- c(
      "",
      "---",
      "",
      "## Evidence verification",
      "",
      sprintf("- Claims checked: %s", verification$n_claims_checked),
      sprintf("- Supported: %s", verification$n_supported),
      sprintf("- Unsupported (marked only): %s", verification$n_unsupported_marked),
      sprintf("- Congruence rate: %s", rate_txt),
      ""
    )
    if (length(marked)) {
      hdr <- c(hdr, "### Marked unsupported claims", "")
      for (m in marked) {
        hdr <- c(hdr, paste0("- ⚠️ ", m$text))
      }
      hdr <- c(hdr, "")
    }
  }
  paste0(markdown, "\n", paste(hdr, collapse = "\n"))
}

#' Answer a follow-up question about ranked phenotype candidates
#'
#' Unlike [explainPhenotypeCandidates()], this function responds to a specific
#' user question (e.g. gene comparison) instead of dumping a full evidence
#' report. Uses a local LLM when available; otherwise a rule-based answer.
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param question User question (character).
#' @param query Optional \code{PhenotypeQuery} object.
#' @param object Optional \code{LazyGas} object to load saved queries.
#' @param top_n Number of ranked genes to include as context.
#' @param use_llm Use local LLM when available.
#' @param language \code{"ja"} or \code{"en"}.
#' @param model,base_url LLM settings.
#' @param chat_history Optional prior chat messages (list of \code{role}/\code{content}).
#'
#' @return A character string answering the question.
#' @export
answerPhenotypeQuestion <- function(rank_result,
                                    question,
                                    query = NULL,
                                    object = NULL,
                                    top_n = 10L,
                                    use_llm = TRUE,
                                    language = c("ja", "en"),
                                    model = NULL,
                                    base_url = NULL,
                                    chat_history = NULL) {
  language <- match.arg(language)
  if (!is.character(question) || length(question) != 1L || !nzchar(trimws(question))) {
    stop("'question' must be a non-empty character string.", call. = FALSE)
  }
  question <- trimws(question)

  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }

  meta <- attr(rank_result, "phenotypeRank")
  if (is.null(query) && !is.null(meta$query_id) && !is.null(object)) {
    query <- .store_read_phenotype_query(object, meta$query_id)
  }
  if (is.null(query)) {
    query <- list(
      trait_text = "phenotype query",
      query_id = if (!is.null(meta)) meta$query_id else "unknown"
    )
    class(query) <- c("PhenotypeQuery", "list")
  }

  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)

  if (isTRUE(use_llm) && llmHealthCheck(base_url = base_url, timeout = 5)) {
    lang_note <- if (language == "ja") {
      paste(
        "ユーザーの質問に直接答えてください。",
        "全文レポートや箇条書きの証拠ダンプは不要です。",
        "比較質問なら結論を最初に述べ、根拠を短く説明してください。",
        "日本語で回答してください。"
      )
    } else {
      paste(
        "Answer the user's question directly.",
        "Do not dump a full evidence report.",
        "For comparison questions, state the conclusion first, then brief reasons.",
        "Write in English."
      )
    }

    messages <- list(
      list(
        role = "system",
        content = paste(
          "You are a plant genetics assistant.",
          "Use ONLY the JSON evidence provided.",
          "Do not recommend genes not in the list.",
          lang_note
        )
      )
    )
    if (!is.null(chat_history) && length(chat_history)) {
      for (msg in chat_history) {
        if (!is.null(msg$role) && !is.null(msg$content) && nzchar(msg$content)) {
          messages[[length(messages) + 1L]] <- list(
            role = msg$role,
            content = as.character(msg$content)
          )
        }
      }
    }
    messages[[length(messages) + 1L]] <- list(
      role = "user",
      content = jsonlite::toJSON(
        list(
          question = question,
          phenotype_query = unclass(query),
          ranked_genes = evidence_bundle
        ),
        auto_unbox = TRUE,
        pretty = TRUE
      )
    )

    return(tryCatch(
      llmChat(
        messages = messages,
        model = model,
        base_url = base_url,
        timeout = 180
      ),
      error = function(e) {
        warning("LLM answer failed; using rule-based reply.", call. = FALSE)
        .answer_phenotype_question_rules(
          question = question,
          query = query,
          evidence_bundle = evidence_bundle,
          language = language
        )
      }
    ))
  }

  .answer_phenotype_question_rules(
    question = question,
    query = query,
    evidence_bundle = evidence_bundle,
    language = language
  )
}

.phenotype_evidence_bundle <- function(rank_result) {
  lapply(seq_len(nrow(rank_result)), function(i) {
    row <- rank_result[i, , drop = FALSE]
    ev <- tryCatch(
      jsonlite::fromJSON(row$evidence_json[1L], simplifyVector = FALSE),
      error = function(e) list()
    )
    list(
      rank = i,
      Gene_ID = row$Gene_ID[1L],
      Name = if ("Name" %in% names(row)) row$Name[1L] else NA_character_,
      composite_score = row$composite_score[1L],
      negLog10P = if ("negLog10P" %in% names(row)) row$negLog10P[1L] else NA_real_,
      dist2peak = if ("dist2peak" %in% names(row)) row$dist2peak[1L] else NA_real_,
      scores = list(
        annotation = if ("score_annotation" %in% names(row)) row$score_annotation[1L] else NA_real_,
        snpeff = if ("score_snpeff" %in% names(row)) row$score_snpeff[1L] else NA_real_,
        gwas = if ("score_gwas" %in% names(row)) row$score_gwas[1L] else NA_real_,
        expression = if ("score_expression" %in% names(row)) row$score_expression[1L] else NA_real_,
        literature = if ("score_literature" %in% names(row)) row$score_literature[1L] else NA_real_,
        ortholog = if ("score_ortholog" %in% names(row)) row$score_ortholog[1L] else NA_real_
      ),
      evidence = ev
    )
  })
}

.is_compare_question <- function(question) {
  q <- tolower(question)
  grepl(
    "どちら|どっち|より|比較|信頼|おすすめ|上位|強い|弱い|which|more reliable|better candidate|compare|versus|\\bvs\\b",
    q,
    perl = TRUE
  )
}

.is_ranking_question <- function(question) {
  q <- tolower(question)
  grepl(
    "なぜ|理由|why|how come|rank|ランク|スコア|根拠",
    q,
    perl = TRUE
  )
}

.evidence_snippet_text <- function(item, source) {
  block <- item$evidence[[source]]
  if (is.null(block) || !length(block$snippets)) {
    return(character())
  }
  as.character(unlist(block$snippets))
}

.answer_phenotype_question_rules <- function(question,
                                            query,
                                            evidence_bundle,
                                            language) {
  if (!length(evidence_bundle)) {
    return(if (language == "ja") {
      "候補遺伝子データがありません。先に Rank candidates を実行してください。"
    } else {
      "No ranked genes available. Run Rank candidates first."
    })
  }

  if (.is_compare_question(question)) {
    return(.answer_compare_genes(
      evidence_bundle = evidence_bundle,
      query = query,
      language = language
    ))
  }

  if (.is_ranking_question(question)) {
    return(.answer_ranking_rationale(
      query = query,
      evidence_bundle = evidence_bundle,
      language = language
    ))
  }

  top <- evidence_bundle[[1L]]
  if (language == "ja") {
    paste0(
      "ご質問「", question, "」について:\n",
      "現在のランキング1位は **", top$Gene_ID, "**（総合スコア ",
      sprintf("%.3f", top$composite_score), "）です。",
      "比較の質問には「どちらがより信頼できる？」のように聞いてください。"
    )
  } else {
    paste0(
      "Regarding your question \"", question, "\":\n",
      "The top-ranked gene is **", top$Gene_ID, "** (composite score ",
      sprintf("%.3f", top$composite_score), ").",
      " For comparisons, try asking which gene is more reliable."
    )
  }
}

.answer_compare_genes <- function(evidence_bundle, language, query = NULL) {
  if (length(evidence_bundle) < 2L) {
    g <- evidence_bundle[[1L]]
    if (language == "ja") {
      return(paste0(
        "比較対象は **", g$Gene_ID, "** のみです（総合スコア ",
        sprintf("%.3f", g$composite_score), "）。他候補がランキングに含まれていません。"
      ))
    }
    return(paste0(
      "Only **", g$Gene_ID, "** is available (composite ",
      sprintf("%.3f", g$composite_score), ")."
    ))
  }

  a <- evidence_bundle[[1L]]
  b <- evidence_bundle[[2L]]
  winner <- if (a$composite_score >= b$composite_score) a else b
  loser <- if (identical(winner, a)) b else a

  reasons <- character()
  if (language == "ja") {
    if (winner$scores$annotation > loser$scores$annotation + 0.05) {
      ann_w <- paste(.evidence_snippet_text(winner, "annotation"), collapse = " ")
      trait <- if (!is.null(query$trait_text)) query$trait_text else "表現型"
      reasons <- c(
        reasons,
        paste0(
          "- **機能アノテーション**: ", winner$Gene_ID, " の方が「", trait,
          "」に一致",
          if (nzchar(ann_w)) paste0("（", ann_w, "）") else ""
        )
      )
    }
    if (!is.na(winner$dist2peak) && !is.na(loser$dist2peak)) {
      reasons <- c(
        reasons,
        paste0(
          "- **GWAS 位置**: ", winner$Gene_ID, " dist2peak=",
          format(winner$dist2peak, big.mark = ","), " bp vs ",
          loser$Gene_ID, " ", format(loser$dist2peak, big.mark = ","), " bp"
        )
      )
    } else if (!is.na(winner$negLog10P)) {
      reasons <- c(
        reasons,
        paste0(
          "- **GWAS**: negLog10P ", sprintf("%.2f", winner$negLog10P),
          " vs ", sprintf("%.2f", loser$negLog10P)
        )
      )
    }
    if (winner$scores$expression > loser$scores$expression + 0.01) {
      expr_w <- paste(.evidence_snippet_text(winner, "expression"), collapse = "; ")
      expr_l <- paste(.evidence_snippet_text(loser, "expression"), collapse = "; ")
      reasons <- c(
        reasons,
        paste0("- **発現**: ", winner$Gene_ID, " — ", expr_w, " / ", loser$Gene_ID, " — ", expr_l)
      )
    }
    if (winner$scores$ortholog > loser$scores$ortholog + 0.01) {
      ortho_w <- paste(.evidence_snippet_text(winner, "ortholog"), collapse = ", ")
      ortho_l <- paste(.evidence_snippet_text(loser, "ortholog"), collapse = ", ")
      reasons <- c(
        reasons,
        paste0("- **オルソログ**: ", winner$Gene_ID, " (", ortho_w, ") vs ",
               loser$Gene_ID, " (", ortho_l, ")")
      )
    }
    if (winner$scores$literature > 0 && loser$scores$literature > 0) {
      reasons <- c(
        reasons,
        "- **文献**: 両遺伝子とも汎用的な PubMed ヒットが付いており、候補判断への寄与は限定的です"
      )
    }

    hdr <- paste0(
      "**", winner$Gene_ID, "** の方が候補としてより信頼できます（総合スコア ",
      sprintf("%.3f", winner$composite_score), " vs ",
      sprintf("%.3f", loser$composite_score), "）。\n\n",
      if (length(reasons)) "主な理由:\n" else "",
      paste(reasons, collapse = "\n")
    )
    return(hdr)
  }

  # English
  if (winner$scores$annotation > loser$scores$annotation + 0.05) {
    reasons <- c(reasons, paste0(
      "- **Annotation**: ", winner$Gene_ID, " matches the phenotype query better."
    ))
  }
  if (!is.na(winner$negLog10P)) {
    reasons <- c(reasons, paste0(
      "- **GWAS**: negLog10P ", sprintf("%.2f", winner$negLog10P),
      " vs ", sprintf("%.2f", loser$negLog10P)
    ))
  }
  if (winner$scores$ortholog > loser$scores$ortholog + 0.01) {
    reasons <- c(reasons, paste0(
      "- **Ortholog**: ", winner$Gene_ID, " has a stronger ortholog score."
    ))
  }
  paste0(
    "**", winner$Gene_ID, "** is the more reliable candidate (composite ",
    sprintf("%.3f", winner$composite_score), " vs ",
    sprintf("%.3f", loser$composite_score), ").\n\n",
    if (length(reasons)) paste(paste(reasons, collapse = "\n")) else ""
  )
}

.answer_ranking_rationale <- function(query, evidence_bundle, language) {
  top <- evidence_bundle[[1L]]
  lines <- character()
  if (language == "ja") {
    lines <- c(
      lines,
      paste0("**", top$Gene_ID, "** が1位なのは、総合スコア（",
             sprintf("%.3f", top$composite_score), "）が最も高いためです。")
    )
    if (top$scores$annotation > 0) {
      ann <- paste(.evidence_snippet_text(top, "annotation"), collapse = " ")
      lines <- c(lines, paste0("- アノテーション: ", ann))
    }
    if (!is.na(top$negLog10P)) {
      lines <- c(lines, paste0("- GWAS 有意性: negLog10P = ", sprintf("%.2f", top$negLog10P)))
    }
    return(paste(lines, collapse = "\n"))
  }
  paste0(
    "**", top$Gene_ID, "** ranks first with composite score ",
    sprintf("%.3f", top$composite_score), "."
  )
}

.explain_phenotype_template <- function(query, evidence_bundle, language) {
  hdr <- if (language == "ja") {
    paste0("## 表現型探索結果: ", query$trait_text, "\n\n")
  } else {
    paste0("## Phenotype exploration: ", query$trait_text, "\n\n")
  }

  lines <- character()
  for (item in evidence_bundle) {
    title <- if (language == "ja") {
      sprintf("### %d. %s (総合スコア %.3f)", item$rank, item$Gene_ID, item$composite_score)
    } else {
      sprintf("### %d. %s (composite %.3f)", item$rank, item$Gene_ID, item$composite_score)
    }
    lines <- c(lines, title)

    ev <- item$evidence
    for (src in names(ev)) {
      block <- ev[[src]]
      if (is.null(block) || !length(block$snippets)) {
        next
      }
      snip <- block$snippets
      if (length(snip) > 3L) {
        snip <- snip[seq_len(3L)]
      }
      lines <- c(
        lines,
        sprintf("- **%s** (score %.2f): %s", src, block$score %||% 0, paste(snip, collapse = "; "))
      )
    }
    lines <- c(lines, "")
  }

  paste0(hdr, paste(lines, collapse = "\n"))
}

.llm_model <- function(model) {
  settings <- .ai_config_llm_settings(model = model, base_url = NULL, timeout = NULL)
  settings$model
}

.llm_base_url <- function(base_url) {
  settings <- .ai_config_llm_settings(model = NULL, base_url = base_url, timeout = NULL)
  settings$base_url
}

.llm_http_get <- function(path, base_url, timeout = 10) {
  url <- paste0(base_url, path)
  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    req <- httr2::req_timeout(req, timeout)
    resp <- tryCatch(httr2::req_perform(req), error = function(e) NULL)
    if (is.null(resp)) {
      return(NULL)
    }
    return(httr2::resp_body_json(resp))
  }

  .llm_http_get_base(url, timeout)
}

.llm_http_post <- function(path, body, base_url, timeout = 120) {
  url <- paste0(base_url, path)
  payload <- jsonlite::toJSON(body, auto_unbox = TRUE)

  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    req <- httr2::req_timeout(req, timeout)
    req <- httr2::req_body_raw(req, payload, "application/json")
    resp <- httr2::req_perform(req)
    return(httr2::resp_body_json(resp))
  }

  .llm_http_post_base(url, payload, timeout)
}

.llm_http_get_base <- function(url, timeout) {
  con <- utils::url(url, open = "rt")
  on.exit(close(con), add = TRUE)
  if (timeout > 0) {
    setTimeLimit(elapsed = timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE), add = TRUE)
  }
  txt <- readLines(con, warn = FALSE)
  if (!length(txt)) {
    return(NULL)
  }
  jsonlite::fromJSON(paste(txt, collapse = "\n"), simplifyVector = FALSE)
}

.llm_http_post_base <- function(url, payload, timeout) {
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    curl::handle_setheaders(h, "Content-Type" = "application/json")
    curl::handle_setopt(h, postfields = payload)
    if (timeout > 0) {
      curl::handle_setopt(h, timeout = timeout)
    }
    txt <- curl::curl_fetch_memory(url, handle = h)$content
    return(jsonlite::fromJSON(rawToChar(txt), simplifyVector = FALSE))
  }

  stop(
    "HTTP requests require 'httr2' or 'curl'. Install with install.packages(\"httr2\").",
    call. = FALSE
  )
}

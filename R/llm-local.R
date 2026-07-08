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
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param query A \code{PhenotypeQuery} object (optional; inferred from
#'   \code{rank_result$query_id} when \code{object} is given).
#' @param object \code{LazyGas} object to load saved queries.
#' @param top_n Number of top genes to include in the explanation.
#' @param use_llm Use local LLM when available.
#' @param language Response language: \code{"ja"} or \code{"en"}.
#' @param model,base_url LLM settings passed to [llmChat()].
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
                                       base_url = NULL) {
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

  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)

  if (!isTRUE(use_llm) || !llmHealthCheck(base_url = base_url, timeout = 5)) {
    return(.explain_phenotype_template(
      query = query,
      evidence_bundle = evidence_bundle,
      language = language
    ))
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

  tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = system_msg),
        list(role = "user", content = user_msg)
      ),
      model = model,
      base_url = base_url,
      timeout = 180
    ),
    error = function(e) {
      warning("LLM explanation failed; using template report.", call. = FALSE)
      .explain_phenotype_template(
        query = query,
        evidence_bundle = evidence_bundle,
        language = language
      )
    }
  )
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
        annotation = row$score_annotation[1L],
        gwas = row$score_gwas[1L],
        expression = row$score_expression[1L],
        literature = row$score_literature[1L],
        ortholog = row$score_ortholog[1L]
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

`%||%` <- function(x, y) if (is.null(x)) y else x

.llm_model <- function(model) {
  if (!is.null(model) && nzchar(model)) {
    return(model)
  }
  env <- Sys.getenv("LAZYGAS_LLM_MODEL", unset = "llama3.2:3b")
  if (nzchar(env)) {
    return(env)
  }
  "llama3.2:3b"
}

.llm_base_url <- function(base_url) {
  if (!is.null(base_url) && nzchar(base_url)) {
    return(sub("/+$", "", base_url))
  }
  env <- Sys.getenv("LAZYGAS_LLM_URL", unset = "http://127.0.0.1:11434")
  sub("/+$", "", env)
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

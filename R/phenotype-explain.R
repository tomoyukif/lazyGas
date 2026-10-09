################################################################################
#' Explain ranked phenotype candidates (template report)
#'
#' Generates a Markdown summary of the top-ranked genes using only evidence
#' collected during ranking (deterministic template; no LLM).
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param query A \code{PhenotypeQuery} object (optional; inferred from
#'   \code{rank_result$query_id} when \code{object} is given).
#' @param object \code{LazyGas} object to load saved queries.
#' @param top_n Number of top genes to include in the explanation.
#' @param use_llm Ignored; kept for call compatibility (always template path).
#' @param language Response language: \code{"ja"} or \code{"en"}.
#' @param model,base_url,timeout Ignored; kept for call compatibility.
#' @param save If \code{TRUE}, write Markdown under the companion
#'   \code{reports/} folder (requires \code{object}). Default \code{FALSE}.
#' @param out_fn Optional explicit Markdown path when \code{save = TRUE}.
#'
#' @return A character string (Markdown).
#' @export
#'
#' @seealso [rankPhenotypeCandidates()], [phenotypeQuery()]
explainPhenotypeCandidates <- function(rank_result,
                                       query = NULL,
                                       object = NULL,
                                       top_n = 10L,
                                       use_llm = FALSE,
                                       language = c("ja", "en"),
                                       model = NULL,
                                       base_url = NULL,
                                       timeout = NULL,
                                       save = FALSE,
                                       out_fn = NULL) {
  language <- match.arg(language)
  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }
  if (isTRUE(use_llm)) {
    warning(
      "'use_llm = TRUE' is not supported in this build; using template report.",
      call. = FALSE
    )
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
  md <- .explain_phenotype_template(
    query = query,
    evidence_bundle = evidence_bundle,
    language = language
  )

  if (isTRUE(save)) {
    if (is.null(object)) {
      stop("'object' is required when save = TRUE.", call. = FALSE)
    }
    if (is.null(out_fn) || !nzchar(as.character(out_fn)[1L])) {
      out_fn <- .default_explain_path(
        object,
        query_id = query$query_id %||% meta$query_id,
        pheno = meta$pheno %||% NULL
      )
    } else {
      dir.create(dirname(out_fn), recursive = TRUE, showWarnings = FALSE)
    }
    writeLines(md, out_fn, useBytes = TRUE)
  }
  md
}

#' Answer a free-text question about ranked candidates (rule-based)
#'
#' @param rank_result Output of [rankPhenotypeCandidates()].
#' @param question User question string.
#' @param query Optional [PhenotypeQuery].
#' @param object Optional LazyGas object.
#' @param top_n Genes considered.
#' @param use_llm Ignored.
#' @param language \code{"ja"} or \code{"en"}.
#' @param model,base_url,timeout Ignored.
#' @return Character answer.
#' @export
answerPhenotypeQuestion <- function(rank_result,
                                    question,
                                    query = NULL,
                                    object = NULL,
                                    top_n = 10L,
                                    use_llm = FALSE,
                                    language = c("ja", "en"),
                                    model = NULL,
                                    base_url = NULL,
                                    timeout = NULL) {
  language <- match.arg(language)
  if (!is.character(question) || length(question) != 1L || !nzchar(trimws(question))) {
    stop("'question' must be a non-empty character string.", call. = FALSE)
  }
  if (!is.data.frame(rank_result) || nrow(rank_result) == 0L) {
    stop("'rank_result' must be a non-empty data.frame.", call. = FALSE)
  }
  meta <- attr(rank_result, "phenotypeRank")
  if (is.null(query) && !is.null(meta$query_id) && !is.null(object)) {
    query <- .store_read_phenotype_query(object, meta$query_id)
  }
  if (is.null(query)) {
    query <- list(trait_text = "phenotype query")
    class(query) <- c("PhenotypeQuery", "list")
  }
  top_n <- as.integer(top_n)[1L]
  sub <- rank_result[seq_len(min(top_n, nrow(rank_result))), , drop = FALSE]
  evidence_bundle <- .phenotype_evidence_bundle(rank_result = sub)
  .answer_phenotype_question_rules(
    question = trimws(question),
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


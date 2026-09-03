################################################################################
# Phase 2 B2: HTML report helpers (English canonical; ja = dictionary + translate)

#' Report UI strings (English canonical)
#' @keywords internal
.report_i18n <- function(language = c("en", "ja")) {
  language <- match.arg(language)
  en <- list(
    basic_info = "Basic information",
    input_files = "Input files",
    samples = "Samples",
    markers = "Markers",
    phenotype = "Phenotype",
    non_missing = "Non-missing phenotype",
    candidate_list = "Candidate list",
    evidence = "Evidence",
    score_note = paste(
      "Scores are a weighted sum of values normalized to 0–1",
      "among candidates within this peak"
    ),
    composite = "composite",
    gwas_pos = "GWAS / position",
    snpeff = "Variant effect (SnpEff)",
    keywords = "Keywords / annotation",
    expression = "Expression",
    validity = "Brief validity / caveats",
    credible_set = "Credible set",
    variants = "variants",
    resolvable = "resolvable",
    partially_resolvable = "partially resolvable",
    not_resolved = "not resolved",
    verification = "Evidence verification",
    claims_checked = "Claims checked",
    supported = "Supported",
    unsupported = "Unsupported (marked only)",
    congruence = "Congruence rate",
    marked_claims = "Marked unsupported claims"
  )
  if (identical(language, "en")) {
    return(en)
  }
  ja <- en
  ja$basic_info <- "基礎情報"
  ja$input_files <- "入力ファイル"
  ja$samples <- "サンプル数"
  ja$markers <- "マーカー数"
  ja$phenotype <- "表現型"
  ja$non_missing <- "表現型非欠損"
  ja$candidate_list <- "候補一覧"
  ja$evidence <- "根拠"
  ja$score_note <- "スコアは当該ピーク内の候補間で 0–1 に正規化した値の重み付き和"
  ja$composite <- "総合"
  ja$gwas_pos <- "GWAS・位置"
  ja$snpeff <- "変異影響（SnpEff）"
  ja$keywords <- "キーワード／アノテーション"
  ja$expression <- "発現"
  ja$validity <- "短い妥当性・注意"
  ja$credible_set <- "Credible set"
  ja$variants <- "バリアント"
  ja$resolvable <- "分解可能"
  ja$partially_resolvable <- "部分的に分解可能"
  ja$not_resolved <- "分解不能"
  ja$verification <- "根拠検証"
  ja$claims_checked <- "照合クレーム数"
  ja$supported <- "支持されたクレーム"
  ja$unsupported <- "根拠なし（マークのみ）"
  ja$congruence <- "整合率"
  ja$marked_claims <- "根拠なしとしてマークした記述"
  ja
}

#' @keywords internal
.report_escape_html <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

#' @keywords internal
.report_fmt_score <- function(x) {
  x <- as.numeric(x)
  out <- rep("0.000", length(x))
  ok <- is.finite(x)
  out[ok] <- sprintf("%.3f", x[ok])
  out
}

#' @keywords internal
.report_name_cell <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L) {
    x <- x[1L]
  }
  if (is.na(x) || !nzchar(trimws(x))) {
    return("")
  }
  trimws(x)
}

#' Require SnpEff + credible sets for report peaks
#' @keywords internal
.report_require_snpeff_cs <- function(object, pheno_name, peak_ids) {
  snpeff <- tryCatch(
    lazyData(object = object, dataset = "snpeff", pheno = pheno_name),
    error = function(e) NULL
  )
  if (is.null(snpeff) || !nrow(snpeff) || !"Gene_ID" %in% names(snpeff)) {
    stop(
      "SnpEff annotations are required for llm_report(). ",
      "Run listCandidate(..., snpeff = ...) first.",
      call. = FALSE
    )
  }
  peak_ids <- unique(as.character(peak_ids))
  peak_ids <- peak_ids[!is.na(peak_ids) & nzchar(peak_ids)]
  for (pid in peak_ids) {
    cred <- tryCatch(
      lazyData(
        object = object,
        dataset = "credible_set",
        pheno = pheno_name,
        kind = paste0("peak_", pid)
      ),
      error = function(e) NULL
    )
    if (is.null(cred) || !nrow(cred) || !"PIP" %in% names(cred)) {
      stop(
        "Credible set required for peak ", pid,
        " (phenotype '", pheno_name, "'). Run calcCredibleSet() first.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' Ensure finemap is among ranking sources/weights for llm_report
#' @keywords internal
.report_sources_with_finemap <- function(sources, weights) {
  sources <- unique(c(as.character(sources), "finemap"))
  sources <- match.arg(sources, choices = .EVIDENCE_SOURCE_CHOICES, several.ok = TRUE)
  if (is.null(weights) || !length(weights)) {
    weights <- phase1DefaultWeights()
  }
  weights <- unlist(weights)
  if (!"finemap" %in% names(weights)) {
    weights <- c(weights, finemap = 0.5)
  }
  list(sources = sources, weights = weights)
}

#' Credible-set summary line under a peak heading
#' @keywords internal
.report_cs_comment_html <- function(object, pheno_name, peak_id, language = "en") {
  ui <- .report_i18n(language)
  cred <- lazyData(
    object = object,
    dataset = "credible_set",
    pheno = pheno_name,
    kind = paste0("peak_", peak_id)
  )
  n_in <- if ("in_credible_set" %in% names(cred)) {
    sum(as.logical(cred$in_credible_set), na.rm = TRUE)
  } else {
    NA_integer_
  }
  summ <- attr(cred, "summary")
  if (!is.null(summ$n_variants_in_set)) {
    n_in <- as.integer(summ$n_variants_in_set)[1L]
  }
  if (!is.null(summ$n_variants_in_locus)) {
    n_locus <- as.integer(summ$n_variants_in_locus)[1L]
  } else {
    n_locus <- nrow(cred)
  }
  resolution <- if (!is.null(summ$resolution)) {
    as.character(summ$resolution)[1L]
  } else if (is.finite(n_in) && n_in <= 3L) {
    "high"
  } else if (is.finite(n_in) && n_in <= 10L) {
    "moderate"
  } else {
    "low"
  }
  label <- switch(
    resolution,
    high = ui$resolvable,
    moderate = ui$partially_resolvable,
    ui$not_resolved
  )
  detail <- if (!is.null(summ$message) && nzchar(as.character(summ$message)[1L])) {
    paste0(" (", .report_escape_html(summ$message[1L]), ")")
  } else {
    ""
  }
  sprintf(
    "<p><strong>%s:</strong> %s %s — %s%s</p>",
    .report_escape_html(ui$credible_set),
    .report_escape_html(as.character(n_in)),
    .report_escape_html(ui$variants),
    .report_escape_html(label),
    detail
  )
}

#' Candidate score table (scrollable)
#' @keywords internal
.report_candidate_table_html <- function(rank_sub, sources, language = "en") {
  ui <- .report_i18n(language)
  sources <- intersect(as.character(sources), .EVIDENCE_SOURCE_CHOICES)
  if (!length(sources)) {
    sources <- c("annotation", "snpeff", "gwas", "expression", "finemap")
  }
  score_cols <- paste0("score_", sources)
  missing <- setdiff(score_cols, names(rank_sub))
  for (mc in missing) {
    rank_sub[[mc]] <- 0
  }
  if (!"composite_score" %in% names(rank_sub)) {
    rank_sub$composite_score <- 0
  }
  if (!"Name" %in% names(rank_sub)) {
    rank_sub$Name <- ""
  }

  header <- c("Gene_ID", "Name", sources, ui$composite)
  th <- paste0(
    "<th>", vapply(header, .report_escape_html, character(1L)), "</th>",
    collapse = ""
  )
  body_rows <- character(nrow(rank_sub))
  for (i in seq_len(nrow(rank_sub))) {
    cells <- c(
      .report_escape_html(rank_sub$Gene_ID[i]),
      .report_escape_html(.report_name_cell(rank_sub$Name[i])),
      .report_fmt_score(unlist(rank_sub[i, score_cols, drop = TRUE])),
      .report_fmt_score(rank_sub$composite_score[i])
    )
    body_rows[i] <- paste0(
      "<tr>",
      paste0("<td>", cells, "</td>", collapse = ""),
      "</tr>"
    )
  }
  paste0(
    "<div class=\"candidate-table-scroll\" style=\"overflow-x: auto;\">\n",
    "<table border=\"1\" cellpadding=\"4\" cellspacing=\"0\">\n",
    "<thead><tr>", th, "</tr></thead>\n",
    "<tbody>\n", paste(body_rows, collapse = "\n"), "\n</tbody>\n",
    "</table>\n</div>\n",
    "<p>", .report_escape_html(ui$score_note), "</p>\n"
  )
}

#' Parse evidence list from a rank row
#' @keywords internal
.report_row_evidence <- function(row) {
  if (!"evidence_json" %in% names(row) || is.na(row$evidence_json[1L]) ||
      !nzchar(as.character(row$evidence_json[1L]))) {
    return(list())
  }
  tryCatch(
    jsonlite::fromJSON(as.character(row$evidence_json[1L]), simplifyVector = FALSE),
    error = function(e) list()
  )
}

#' Code-written fact blocks for one gene
#' @keywords internal
.report_gene_fact_blocks <- function(row, language = "en") {
  ui <- .report_i18n(language)
  ev <- .report_row_evidence(row)
  gwas <- ev$gwas$details %||% list()
  snpeff <- ev$snpeff$details %||% list()
  fm <- ev$finemap$details %||% list()

  dist <- if (!is.null(gwas$dist2peak)) {
    gwas$dist2peak
  } else if ("dist2peak" %in% names(row)) {
    row$dist2peak[1L]
  } else {
    NA_real_
  }
  nlp <- if (!is.null(gwas$negLog10P)) {
    gwas$negLog10P
  } else if ("negLog10P" %in% names(row)) {
    row$negLog10P[1L]
  } else {
    NA_real_
  }
  pip <- fm$max_PIP %||% NA_real_

  gwas_lines <- c(
    if (is.finite(as.numeric(dist)[1L])) {
      sprintf("dist2peak = %s bp", format(as.numeric(dist)[1L], big.mark = ",", scientific = FALSE))
    },
    if (is.finite(as.numeric(nlp)[1L])) {
      sprintf("negLog10P = %.2f", as.numeric(nlp)[1L])
    },
    if (is.finite(as.numeric(pip)[1L])) {
      sprintf("max_PIP (95%% CS ∩ SnpEff) = %.3f", as.numeric(pip)[1L])
    }
  )
  if (!length(gwas_lines)) {
    gwas_lines <- "No GWAS / position metrics available."
  }

  snp_lines <- character()
  for (nm in c("HIGH", "MODERATE", "LOW", "MODIFIER")) {
    if (!is.null(snpeff[[nm]])) {
      snp_lines <- c(snp_lines, sprintf("%s = %s", nm, snpeff[[nm]]))
    }
  }
  if (!is.null(snpeff$worst_impact)) {
    snp_lines <- c(sprintf("worst_impact = %s", snpeff$worst_impact), snp_lines)
  }
  if (!length(snp_lines)) {
    snp_lines <- "No SnpEff impact counts available."
  }

  has_expr <- !is.null(ev$expression) &&
    (isTRUE(as.numeric(ev$expression$score %||% 0) > 0) ||
       length(ev$expression$snippets %||% NULL) > 0L)

  list(
    gwas_html = paste0(
      "<h5>", .report_escape_html(ui$gwas_pos), "</h5>\n",
      "<ul>",
      paste0("<li>", .report_escape_html(gwas_lines), "</li>", collapse = ""),
      "</ul>\n"
    ),
    snpeff_html = paste0(
      "<h5>", .report_escape_html(ui$snpeff), "</h5>\n",
      "<ul>",
      paste0("<li>", .report_escape_html(snp_lines), "</li>", collapse = ""),
      "</ul>\n"
    ),
    has_expression = has_expr,
    ann_matched = unlist(ev$annotation$details$matched_keywords %||% list()),
    ann_unmatched = unlist(ev$annotation$details$unmatched_keywords %||% list()),
    expr_snippets = as.character(unlist(ev$expression$snippets %||% list()))
  )
}

#' LLM-free English prose for keyword / expression / validity
#' @keywords internal
.report_gene_template_prose <- function(facts) {
  matched <- facts$ann_matched
  unmatched <- facts$ann_unmatched
  kw <- if (!length(matched) && !length(unmatched)) {
    "No keyword match details were available."
  } else {
    paste0(
      "Keyword match is ",
      if (length(matched)) "present" else "none",
      ". Matched: ",
      if (length(matched)) paste(matched, collapse = ", ") else "(none)",
      ". Unmatched: ",
      if (length(unmatched)) paste(unmatched, collapse = ", ") else "(none)",
      "."
    )
  }
  expr <- if (isTRUE(facts$has_expression)) {
    if (length(facts$expr_snippets)) {
      paste(facts$expr_snippets, collapse = "; ")
    } else {
      "Expression evidence is present."
    }
  } else {
    NULL
  }
  list(
    keywords = kw,
    expression = expr,
    validity = "Automated template narrative; review the coded metrics above."
  )
}

#' Strip ranking scores from evidence bundle for LLM
#' @keywords internal
.report_llm_evidence_payload <- function(rank_sub) {
  lapply(seq_len(nrow(rank_sub)), function(i) {
    row <- rank_sub[i, , drop = FALSE]
    ev <- .report_row_evidence(row)
    list(
      Gene_ID = as.character(row$Gene_ID[1L]),
      Name = .report_name_cell(row$Name[1L]),
      dist2peak = if ("dist2peak" %in% names(row)) row$dist2peak[1L] else NA_real_,
      negLog10P = if ("negLog10P" %in% names(row)) row$negLog10P[1L] else NA_real_,
      evidence = ev
    )
  })
}

#' Ask LLM for English keyword / expression / validity JSON per gene
#' @keywords internal
.report_llm_evidence_prose <- function(rank_sub, query, llm) {
  payload <- .report_llm_evidence_payload(rank_sub)
  system_msg <- paste(
    "You are a plant genetics assistant.",
    "Return ONLY a JSON object keyed by Gene_ID.",
    "Each value must be an object with keys: keywords, expression, validity.",
    "Write all text in English.",
    "Do not invent genes. Do not print composite_score or score_* values.",
    "Do not claim a 'semantic match' or embedding similarity.",
    "For keywords: state high/low/none qualitatively; list matched and unmatched",
    "biological keywords from annotation.details when present.",
    "For expression: summarize expression snippets if present, else use empty string.",
    "For validity: note contradictions or missing information briefly.",
    "Preserve gene IDs and numeric facts exactly when you mention them."
  )
  user_msg <- jsonlite::toJSON(
    list(phenotype_query = unclass(query), genes = payload),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  raw <- llmChat(
    messages = list(
      list(role = "system", content = system_msg),
      list(role = "user", content = user_msg)
    ),
    model = llm$model,
    base_url = llm$base_url,
    json_mode = TRUE,
    timeout = llm$timeout
  )
  parsed <- tryCatch(
    jsonlite::fromJSON(raw, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed) || !is.list(parsed)) {
    warning("LLM evidence JSON parse failed; using templates.", call. = FALSE)
    return(NULL)
  }
  parsed
}

#' Translate English prose blocks to Japanese without changing IDs/numbers
#' @keywords internal
.report_translate_prose_ja <- function(prose_by_gene, llm) {
  system_msg <- paste(
    "Translate the JSON string values from English to Japanese.",
    "Return ONLY JSON with the same keys and structure.",
    "Do not change gene IDs, numbers, punctuation used in IDs, or PIP/distances.",
    "Keep scientific terms when a standard Japanese term is unclear."
  )
  raw <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = system_msg),
        list(
          role = "user",
          content = jsonlite::toJSON(prose_by_gene, auto_unbox = TRUE, pretty = TRUE)
        )
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = llm$timeout
    ),
    error = function(e) {
      warning("Japanese translation failed; keeping English prose. ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(raw)) {
    return(prose_by_gene)
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(raw, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed) || !is.list(parsed)) {
    warning("Japanese translation JSON parse failed; keeping English.", call. = FALSE)
    return(prose_by_gene)
  }
  parsed
}

#' Evidence HTML for genes in one peak
#' @keywords internal
.report_evidence_html <- function(rank_sub,
                                  query,
                                  language = "en",
                                  use_llm = FALSE,
                                  llm = NULL) {
  ui <- .report_i18n(language)
  if (!nrow(rank_sub)) {
    return(paste0("<p>", .report_escape_html("(no genes)"), "</p>\n"))
  }

  prose_map <- list()
  llm_ok <- isTRUE(use_llm) && !is.null(llm) &&
    llmHealthCheck(base_url = llm$base_url, timeout = 5)
  if (llm_ok) {
    prose_map <- .report_llm_evidence_prose(rank_sub, query = query, llm = llm)
  }
  if (is.null(prose_map)) {
    prose_map <- list()
  }

  # Fill missing genes with templates (English)
  for (i in seq_len(nrow(rank_sub))) {
    gid <- as.character(rank_sub$Gene_ID[i])
    facts <- .report_gene_fact_blocks(rank_sub[i, , drop = FALSE], language = "en")
    if (is.null(prose_map[[gid]])) {
      prose_map[[gid]] <- .report_gene_template_prose(facts)
    } else {
      tmpl <- .report_gene_template_prose(facts)
      if (is.null(prose_map[[gid]]$keywords) || !nzchar(prose_map[[gid]]$keywords)) {
        prose_map[[gid]]$keywords <- tmpl$keywords
      }
      if (!isTRUE(facts$has_expression)) {
        prose_map[[gid]]$expression <- NULL
      } else if (is.null(prose_map[[gid]]$expression) ||
                 !nzchar(as.character(prose_map[[gid]]$expression))) {
        prose_map[[gid]]$expression <- tmpl$expression
      }
      if (is.null(prose_map[[gid]]$validity) || !nzchar(prose_map[[gid]]$validity)) {
        prose_map[[gid]]$validity <- tmpl$validity
      }
    }
  }

  if (identical(language, "ja") && llm_ok) {
    prose_map <- .report_translate_prose_ja(prose_map, llm = llm)
  } else if (identical(language, "ja") && !llm_ok) {
    # Dictionary labels only; keep English template prose with a short note
    for (gid in names(prose_map)) {
      prose_map[[gid]]$validity <- paste(
        prose_map[[gid]]$validity,
        "(English template; enable use_llm for Japanese translation.)"
      )
    }
  }

  parts <- character()
  for (i in seq_len(nrow(rank_sub))) {
    row <- rank_sub[i, , drop = FALSE]
    gid <- as.character(row$Gene_ID[1L])
    nm <- .report_name_cell(row$Name[1L])
    title <- if (nzchar(nm)) {
      paste0(gid, " (", nm, ")")
    } else {
      gid
    }
    facts <- .report_gene_fact_blocks(row, language = language)
    prose <- prose_map[[gid]] %||% .report_gene_template_prose(facts)
    blocks <- c(
      facts$gwas_html,
      facts$snpeff_html,
      paste0(
        "<h5>", .report_escape_html(ui$keywords), "</h5>\n",
        "<p>", .report_escape_html(prose$keywords %||% ""), "</p>\n"
      )
    )
    if (isTRUE(facts$has_expression) ||
        (is.character(prose$expression) && nzchar(prose$expression))) {
      blocks <- c(
        blocks,
        paste0(
          "<h5>", .report_escape_html(ui$expression), "</h5>\n",
          "<p>", .report_escape_html(prose$expression %||% ""), "</p>\n"
        )
      )
    }
    blocks <- c(
      blocks,
      paste0(
        "<h5>", .report_escape_html(ui$validity), "</h5>\n",
        "<p>", .report_escape_html(prose$validity %||% ""), "</p>\n"
      )
    )
    parts <- c(
      parts,
      paste0("<h4>", .report_escape_html(title), "</h4>\n", paste(blocks, collapse = ""))
    )
  }
  paste(parts, collapse = "\n")
}

#' Basic information as HTML
#' @keywords internal
.report_basic_info_html <- function(object, pheno_name, language = "en") {
  ui <- .report_i18n(language)
  gds <- .store_gds_fn(object)
  companion <- .store_path(object)
  if (is.null(companion) || !nzchar(as.character(companion)[1L])) {
    companion <- "(none)"
  } else {
    companion <- normalizePath(as.character(companion)[1L], winslash = "/",
                               mustWork = FALSE)
  }
  n_sam <- as.integer(nsam(object))[1L]
  n_mar <- as.integer(nmar(object))[1L]
  pheno_df <- getPheno(object)$pheno
  pheno_vec <- if (!is.null(pheno_df) && pheno_name %in% names(pheno_df)) {
    pheno_df[[pheno_name]]
  } else {
    rep(NA, n_sam)
  }
  n_obs <- sum(!is.na(pheno_vec))
  pct <- if (is.finite(n_sam) && n_sam > 0L) 100 * n_obs / n_sam else NA_real_
  paste0(
    "<h2>", .report_escape_html(ui$basic_info), "</h2>\n",
    "<ul>\n",
    "<li>", .report_escape_html(ui$input_files), ":\n",
    "<ul>\n",
    "<li>GDS: <code>", .report_escape_html(gds), "</code></li>\n",
    "<li>companion: <code>", .report_escape_html(companion), "</code></li>\n",
    "</ul></li>\n",
    "<li>", .report_escape_html(ui$samples), ": ", n_sam, "</li>\n",
    "<li>", .report_escape_html(ui$markers), ": ", n_mar, "</li>\n",
    "<li>", .report_escape_html(ui$phenotype), ": ",
    .report_escape_html(pheno_name), "</li>\n",
    "<li>", .report_escape_html(ui$non_missing), ": ",
    n_obs, " / ", n_sam,
    if (is.finite(pct)) sprintf(" (%.1f%%)", pct) else "",
    "</li>\n",
    "</ul>\n"
  )
}

#' One peak section HTML
#' @keywords internal
.report_peak_section_html <- function(object,
                                      pheno_name,
                                      peak_id,
                                      chr,
                                      pos,
                                      region_start,
                                      region_end,
                                      rank_peak,
                                      sources,
                                      query,
                                      top_n,
                                      language = "en",
                                      use_llm = FALSE,
                                      llm = NULL) {
  ui <- .report_i18n(language)
  header <- sprintf(
    "Peak %s — %s:%s (%s-%s)",
    as.character(peak_id)[1L],
    as.character(chr)[1L],
    .report_fmt_pos(pos),
    .report_fmt_pos(region_start),
    .report_fmt_pos(region_end)
  )
  cs_html <- .report_cs_comment_html(
    object = object,
    pheno_name = pheno_name,
    peak_id = peak_id,
    language = language
  )
  sub <- rank_peak
  if (!is.null(sub) && nrow(sub) > 0L) {
    top_n <- as.integer(top_n)[1L]
    if (is.finite(top_n) && top_n > 0L && nrow(sub) > top_n) {
      sub <- sub[seq_len(top_n), , drop = FALSE]
    }
  } else {
    sub <- data.frame()
  }
  paste0(
    "<h2>", .report_escape_html(header), "</h2>\n",
    cs_html,
    "<h3>", .report_escape_html(ui$candidate_list), "</h3>\n",
    if (nrow(sub)) {
      .report_candidate_table_html(sub, sources = sources, language = language)
    } else {
      "<p>(no candidates)</p>\n"
    },
    "<h3>", .report_escape_html(ui$evidence), "</h3>\n",
    .report_evidence_html(
      rank_sub = sub,
      query = query,
      language = language,
      use_llm = use_llm,
      llm = llm
    )
  )
}

#' Verification section as HTML
#' @keywords internal
.report_verification_html <- function(verification, language = "en") {
  ui <- .report_i18n(language)
  rate <- verification$congruence_rate
  rate_txt <- if (is.finite(rate)) sprintf("%.1f%%", 100 * rate) else "NA"
  marked <- Filter(
    function(x) identical(x$status, "unsupported_marked"),
    verification$claims %||% list()
  )
  lines <- c(
    "<hr/>\n",
    paste0("<h2>", .report_escape_html(ui$verification), "</h2>\n"),
    "<ul>\n",
    sprintf(
      "<li>%s: %s</li>\n",
      .report_escape_html(ui$claims_checked),
      verification$n_claims_checked %||% verification$n_checked %||% 0L
    ),
    sprintf(
      "<li>%s: %s</li>\n",
      .report_escape_html(ui$supported),
      verification$n_supported %||% 0L
    ),
    sprintf(
      "<li>%s: %s</li>\n",
      .report_escape_html(ui$unsupported),
      verification$n_unsupported_marked %||% verification$n_marked %||% 0L
    ),
    sprintf(
      "<li>%s: %s</li>\n",
      .report_escape_html(ui$congruence),
      rate_txt
    ),
    "</ul>\n"
  )
  if (length(marked)) {
    lines <- c(
      lines,
      paste0("<h3>", .report_escape_html(ui$marked_claims), "</h3>\n"),
      "<ul>\n",
      vapply(marked, function(m) {
        paste0("<li>⚠️ ", .report_escape_html(m$text), "</li>\n")
      }, character(1L)),
      "</ul>\n"
    )
  }
  paste(lines, collapse = "")
}

#' Wrap body fragments in a minimal HTML document
#' @keywords internal
.report_html_document <- function(body_html, title = "lazyGas AI report") {
  paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n",
    "<meta charset=\"utf-8\"/>\n",
    "<title>", .report_escape_html(title), "</title>\n",
    "<style>\n",
    "body{font-family:system-ui,sans-serif;margin:1.5rem;line-height:1.45;}\n",
    "table{border-collapse:collapse;}\n",
    "th,td{border:1px solid #ccc;padding:0.35rem 0.5rem;white-space:nowrap;}\n",
    ".candidate-table-scroll{max-width:100%;}\n",
    "code{font-size:0.9em;}\n",
    "</style>\n</head>\n<body>\n",
    body_html,
    "\n</body>\n</html>\n"
  )
}

################################################################################
# Phase 2 F: peak reinterpretation (hypothesis cards)
# Deterministic multi-view champions + conflicts; no LLM reordering.

#' Default reinterpretation views
#'
#' @return Character vector of view names.
#' @export
reinterpretDefaultViews <- function() {
  c("proximity", "function", "impact", "expression", "composite")
}

#' Default rules for [reinterpretPeakCandidates()]
#'
#' Function-view eligibility follows plan alpha: strong terms from the query
#' trait keywords plus synonym phrases only (no manual dictionary); weak-only
#' hits are ineligible; empty function seat is \code{NA}.
#' Distance gating defaults to \code{"off"} (credible-set / peak membership
#' is the spatial filter).
#'
#' @param strong_terms Optional character vector overriding query-derived
#'   strong terms. \code{NULL} = trait keywords + synonym phrases.
#' @param weak_terms Generic terms that never qualify as strong alone.
#' @param min_specificity Minimum \code{ann_specificity} for function view.
#' @param require_strong_hit If \code{TRUE}, function champions need a strong hit.
#' @param dist_gate_mode \code{"off"} (default), \code{"absolute"}, or
#'   \code{"relative"}.
#' @param dist_abs_max Absolute distance cap when gating is on.
#' @param dist_iqr_mult Relative multiplier when \code{dist_gate_mode = "relative"}.
#' @param conflict_if_champions_differ Record conflicts when champions differ.
#' @param soft_flag_composite_only_weak_ann Soft-flag composite winners with
#'   weak-only annotation (stored on support rows; does not change champions).
#'
#' @return A named list of rules.
#' @export
reinterpretDefaultRules <- function(
  strong_terms = NULL,
  weak_terms = c(
    "leaf", "shoot", "plant", "development", "growth",
    "protein", "domain", "family", "putative", "similar"
  ),
  min_specificity = 0.34,
  require_strong_hit = TRUE,
  dist_gate_mode = c("off", "absolute", "relative"),
  dist_abs_max = 5e5,
  dist_iqr_mult = 5,
  conflict_if_champions_differ = TRUE,
  soft_flag_composite_only_weak_ann = TRUE
) {
  dist_gate_mode <- match.arg(dist_gate_mode)
  list(
    strong_terms = strong_terms,
    weak_terms = as.character(weak_terms %||% character()),
    min_specificity = as.numeric(min_specificity)[1L],
    require_strong_hit = isTRUE(require_strong_hit),
    dist_gate_mode = dist_gate_mode,
    dist_abs_max = as.numeric(dist_abs_max)[1L],
    dist_iqr_mult = as.numeric(dist_iqr_mult)[1L],
    conflict_if_champions_differ = isTRUE(conflict_if_champions_differ),
    soft_flag_composite_only_weak_ann = isTRUE(soft_flag_composite_only_weak_ann)
  )
}

#' Reinterpret ranked candidates into per-peak hypothesis champions
#'
#' Builds proximity / function / impact / expression / composite champions
#' without trusting \code{composite_score} alone. Conflicts are listed side by
#' side (no primary-hypothesis label).
#'
#' @param ranked A \code{data.frame} from [rankPhenotypeCandidates()].
#' @param query A [PhenotypeQuery] or \code{NULL} (restored from
#'   \code{attr(ranked, "phenotypeRank")} / \code{query_id} when possible).
#' @param peak_id Optional peak id(s); \code{NULL} = all peaks in \code{ranked}.
#' @param views Character vector; default [reinterpretDefaultViews()].
#' @param rules Rule list; default [reinterpretDefaultRules()].
#' @param top_n_support Rows kept per view in the support table.
#' @param save If \code{TRUE}, write to the companion store (requires
#'   \code{object} and \code{pheno}).
#' @param object Optional \code{LazyGas} for saving / query restore.
#' @param pheno Phenotype name when saving.
#'
#' @return An object of class \code{PeakReinterpretation}.
#' @export
#'
#' @seealso [conflictTable()], [geneSetForReport()], [formatPeakReinterpretation()]
reinterpretPeakCandidates <- function(
  ranked,
  query = NULL,
  peak_id = NULL,
  views = reinterpretDefaultViews(),
  rules = reinterpretDefaultRules(),
  top_n_support = 5L,
  save = FALSE,
  object = NULL,
  pheno = NULL
) {
  if (is.null(ranked) || !is.data.frame(ranked) || nrow(ranked) == 0L) {
    stop("'ranked' must be a non-empty data.frame.", call. = FALSE)
  }
  if (!"Gene_ID" %in% names(ranked)) {
    stop("'ranked' must contain Gene_ID.", call. = FALSE)
  }
  if (is.list(rules) && !is.null(rules$dist_gate_mode)) {
    rules$dist_gate_mode <- match.arg(
      rules$dist_gate_mode,
      c("off", "absolute", "relative")
    )
  } else if (!is.list(rules)) {
    rules <- reinterpretDefaultRules()
  }

  views <- as.character(views)
  views <- views[nzchar(views)]
  if (!length(views)) {
    views <- reinterpretDefaultViews()
  }
  unknown <- setdiff(views, reinterpretDefaultViews())
  if (length(unknown)) {
    stop(
      "Unknown view(s): ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  query <- .reinterpret_resolve_query(ranked, query, object)
  top_n_support <- max(1L, as.integer(top_n_support)[1L])

  df <- ranked
  if (!"peak_ID" %in% names(df)) {
    df$peak_ID <- 1L
  }
  df$Gene_ID <- as.character(df$Gene_ID)
  df$dist2peak <- suppressWarnings(as.numeric(df$dist2peak %||% NA_real_))
  # Prefer peak-block-marker-linked impact counts when present (existing stores).
  if ("HIGH_at_var" %in% names(df)) {
    df$HIGH <- suppressWarnings(as.numeric(df$HIGH_at_var %||% 0))
  } else {
    df$HIGH <- suppressWarnings(as.numeric(df$HIGH %||% 0))
  }
  if ("MODERATE_at_var" %in% names(df)) {
    df$MODERATE <- suppressWarnings(as.numeric(df$MODERATE_at_var %||% 0))
  } else {
    df$MODERATE <- suppressWarnings(as.numeric(df$MODERATE %||% 0))
  }
  if ("LOW_at_var" %in% names(df)) {
    df$LOW <- suppressWarnings(as.numeric(df$LOW_at_var %||% 0))
  }
  if ("MODIFIER_at_var" %in% names(df)) {
    df$MODIFIER <- suppressWarnings(as.numeric(df$MODIFIER_at_var %||% 0))
  }
  df$score_annotation <- suppressWarnings(as.numeric(df$score_annotation %||% 0))
  df$score_snpeff <- suppressWarnings(as.numeric(df$score_snpeff %||% 0))
  df$score_gwas <- suppressWarnings(as.numeric(df$score_gwas %||% 0))
  df$score_expression <- suppressWarnings(as.numeric(df$score_expression %||% 0))
  df$composite_score <- suppressWarnings(as.numeric(df$composite_score %||% 0))
  df$HIGH[is.na(df$HIGH)] <- 0
  df$MODERATE[is.na(df$MODERATE)] <- 0
  for (cn in c(
    "score_annotation", "score_snpeff", "score_gwas", "score_expression",
    "composite_score"
  )) {
    df[[cn]][is.na(df[[cn]])] <- 0
  }

  df$ann_text <- .reinterpret_ann_text(df)
  strong <- .reinterpret_strong_terms(query, rules)
  weak <- unique(tolower(trimws(as.character(rules$weak_terms %||% character()))))
  # Exact weak tokens never count as strong (phrases like "leaf width" stay strong)
  strong <- strong[!tolower(strong) %in% weak]
  spec <- vapply(
    df$ann_text,
    function(txt) .ann_specificity(txt, strong = strong, weak = weak),
    numeric(1L)
  )
  df$ann_specificity <- spec
  df$has_strong_hit <- spec >= 0.67

  peak_ids <- unique(df$peak_ID)
  if (!is.null(peak_id)) {
    peak_ids <- peak_ids[as.character(peak_ids) %in% as.character(peak_id)]
    if (!length(peak_ids)) {
      stop("Requested peak_id not found in ranked table.", call. = FALSE)
    }
  }

  peaks_out <- list()
  for (pid in peak_ids) {
    sub <- df[as.character(df$peak_ID) == as.character(pid), , drop = FALSE]
    sub$dist_gate_pass <- .dist_gate(sub$dist2peak, rules = rules)
    if (isTRUE(rules$soft_flag_composite_only_weak_ann)) {
      sub$soft_flag_weak_ann <- sub$ann_specificity > 0 & sub$ann_specificity < 0.67
    } else {
      sub$soft_flag_weak_ann <- FALSE
    }

    champions <- do.call(rbind, lapply(views, function(v) {
      .pick_champion(sub, view = v, rules = rules)
    }))
    rownames(champions) <- NULL

    support_parts <- lapply(views, function(v) {
      .reinterpret_support_rows(sub, view = v, rules = rules, top_n = top_n_support)
    })
    support <- do.call(rbind, support_parts)
    rownames(support) <- NULL

    conflicts <- .conflict_table(champions, rules = rules)
    gene_set <- .gene_set_from_reinterpret(
      champions = champions,
      support = support,
      mode = "union_champions"
    )

    peaks_out[[as.character(pid)]] <- list(
      peak_ID = pid,
      n_candidates = nrow(sub),
      champions = champions,
      support = support,
      conflicts = conflicts,
      gene_set = gene_set
    )
  }

  pheno_name <- pheno
  if (is.null(pheno_name)) {
    pr <- attr(ranked, "phenotypeRank")
    if (is.list(pr) && !is.null(pr$pheno)) {
      pheno_name <- pr$pheno
    }
  }
  query_id <- if (inherits(query, "PhenotypeQuery")) {
    query$query_id
  } else {
    attr(ranked, "phenotypeRank")$query_id %||% NA_character_
  }

  out <- list(
    meta = list(
      pheno = pheno_name,
      query_id = query_id,
      views = views,
      rules = rules,
      top_n_support = top_n_support,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    peaks = peaks_out
  )
  class(out) <- c("PeakReinterpretation", "list")

  if (isTRUE(save)) {
    if (is.null(object) || !inherits(object, "LazyGas")) {
      stop("'object' LazyGas is required when save = TRUE.", call. = FALSE)
    }
    if (is.null(pheno_name) || !nzchar(as.character(pheno_name)[1L])) {
      stop("'pheno' is required when save = TRUE.", call. = FALSE)
    }
    .store_write_phenotype_reinterpret(
      object = object,
      reinterpret = out,
      pheno_name = as.character(pheno_name)[1L],
      query_id = as.character(query_id)[1L]
    )
  }

  out
}

#' Flatten champions from a PeakReinterpretation
#' @param x A \code{PeakReinterpretation}.
#' @param row.names,optional Unused; required for \code{as.data.frame} S3.
#' @param ... Unused.
#' @return A data.frame.
#' @export
as.data.frame.PeakReinterpretation <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!inherits(x, "PeakReinterpretation")) {
    stop("'x' must be a PeakReinterpretation.", call. = FALSE)
  }
  rows <- lapply(x$peaks, function(pk) {
    ch <- pk$champions
    if (is.null(ch) || !nrow(ch)) {
      return(NULL)
    }
    ch$peak_ID <- pk$peak_ID
    ch
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      peak_ID = integer(),
      view = character(),
      Gene_ID = character(),
      reason_code = character(),
      reason_detail = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("peak_ID", "view", "Gene_ID", "reason_code", "reason_detail")]
}

#' Conflict table from a PeakReinterpretation
#' @param x A \code{PeakReinterpretation}.
#' @return A data.frame of conflicts across peaks.
#' @export
conflictTable <- function(x) {
  if (!inherits(x, "PeakReinterpretation")) {
    stop("'x' must be a PeakReinterpretation.", call. = FALSE)
  }
  rows <- lapply(x$peaks, function(pk) {
    cf <- pk$conflicts
    if (is.null(cf) || !nrow(cf)) {
      return(NULL)
    }
    cf$peak_ID <- pk$peak_ID
    cf
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      peak_ID = integer(),
      type = character(),
      views = character(),
      Gene_IDs = character(),
      note = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c("peak_ID", "type", "views", "Gene_IDs", "note")]
}

#' Gene set for report / E evidence from reinterpretation
#'
#' @param x A \code{PeakReinterpretation}.
#' @param peak_id Peak id.
#' @param mode One of \code{"union_champions"}, \code{"composite_top_n"},
#'   \code{"union_all"}.
#' @param composite_ids Optional Gene_ID vector for composite top_n / union_all.
#' @return Character vector of Gene_ID (no \code{NA}).
#' @export
geneSetForReport <- function(x,
                             peak_id,
                             mode = c("union_champions", "composite_top_n", "union_all"),
                             composite_ids = NULL) {
  mode <- match.arg(mode)
  if (!inherits(x, "PeakReinterpretation")) {
    stop("'x' must be a PeakReinterpretation.", call. = FALSE)
  }
  pk <- x$peaks[[as.character(peak_id)]]
  if (is.null(pk)) {
    stop("peak_id not found in PeakReinterpretation.", call. = FALSE)
  }
  champ_ids <- unique(stats::na.omit(as.character(pk$champions$Gene_ID)))
  support_ids <- if (!is.null(pk$support) && nrow(pk$support)) {
    unique(stats::na.omit(as.character(pk$support$Gene_ID)))
  } else {
    character()
  }
  comp <- unique(stats::na.omit(as.character(composite_ids %||% character())))
  if (identical(mode, "union_champions")) {
    return(unique(c(champ_ids, support_ids)))
  }
  if (identical(mode, "composite_top_n")) {
    return(comp)
  }
  unique(c(champ_ids, support_ids, comp))
}

#' Champion Gene_ID helper
#' @param x PeakReinterpretation
#' @param peak Peak id
#' @param view View name
#' @return Character scalar or \code{NA_character_}
#' @export
champion <- function(x, peak, view) {
  if (!inherits(x, "PeakReinterpretation")) {
    stop("'x' must be a PeakReinterpretation.", call. = FALSE)
  }
  pk <- x$peaks[[as.character(peak)]]
  if (is.null(pk)) {
    return(NA_character_)
  }
  ch <- pk$champions
  hit <- ch[as.character(ch$view) == as.character(view), , drop = FALSE]
  if (!nrow(hit)) {
    return(NA_character_)
  }
  gid <- as.character(hit$Gene_ID[1L])
  if (!nzchar(gid) || is.na(gid)) {
    return(NA_character_)
  }
  gid
}

#' Format hypothesis cards as HTML or Markdown
#'
#' @param x A \code{PeakReinterpretation} or one peak list element.
#' @param peak_id Required when \code{x} is a full PeakReinterpretation.
#' @param language Ignored for chrome (English); accepted for API parity.
#' @param format \code{"html"} or \code{"markdown"}.
#' @return Character scalar.
#' @export
formatPeakReinterpretation <- function(x,
                                       peak_id = NULL,
                                       language = c("en", "ja"),
                                       format = c("html", "markdown")) {
  language <- match.arg(language)
  format <- match.arg(format)
  pk <- if (inherits(x, "PeakReinterpretation")) {
    if (is.null(peak_id)) {
      stop("'peak_id' is required for PeakReinterpretation input.", call. = FALSE)
    }
    x$peaks[[as.character(peak_id)]]
  } else {
    x
  }
  if (is.null(pk) || is.null(pk$champions)) {
    return(if (identical(format, "html")) {
      "<p>(no reinterpretation)</p>\n"
    } else {
      "(no reinterpretation)\n"
    })
  }
  ch <- pk$champions
  if (identical(format, "markdown")) {
    lines <- c(
      "| view | Gene_ID | reason |",
      "|---|---|---|"
    )
    for (i in seq_len(nrow(ch))) {
      gid <- ch$Gene_ID[i]
      if (is.na(gid) || !nzchar(as.character(gid))) {
        gid <- "(none)"
      }
      detail <- ch$reason_detail[i] %||% ch$reason_code[i]
      lines <- c(
        lines,
        sprintf(
          "| %s | %s | %s |",
          ch$view[i],
          gid,
          gsub("\\|", "/", as.character(detail)[1L])
        )
      )
    }
    cf <- pk$conflicts
    if (!is.null(cf) && nrow(cf)) {
      lines <- c(lines, "", "Conflicts:", "")
      for (i in seq_len(nrow(cf))) {
        lines <- c(lines, paste0("- ", cf$note[i]))
      }
    }
    return(paste(c(lines, ""), collapse = "\n"))
  }

  # HTML
  esc <- if (exists(".report_escape_html", mode = "function")) {
    .report_escape_html
  } else {
    function(z) {
      z <- as.character(z)
      z <- gsub("&", "&amp;", z, fixed = TRUE)
      z <- gsub("<", "&lt;", z, fixed = TRUE)
      z <- gsub(">", "&gt;", z, fixed = TRUE)
      z
    }
  }
  rows <- character()
  for (i in seq_len(nrow(ch))) {
    gid <- ch$Gene_ID[i]
    if (is.na(gid) || !nzchar(as.character(gid))) {
      gid <- "(none)"
    }
    detail <- ch$reason_detail[i] %||% ch$reason_code[i]
    rows <- c(
      rows,
      paste0(
        "<tr><td>", esc(ch$view[i]), "</td><td>", esc(gid),
        "</td><td>", esc(detail), "</td></tr>"
      )
    )
  }
  cf_html <- ""
  cf <- pk$conflicts
  if (!is.null(cf) && nrow(cf)) {
    items <- paste0(
      "<li>", esc(cf$note), "</li>",
      collapse = ""
    )
    cf_html <- paste0("<ul>\n", items, "\n</ul>\n")
  }
  paste0(
    "<table class=\"lazygas-hypothesis-cards\">\n",
    "<thead><tr><th>view</th><th>Gene_ID</th><th>reason</th></tr></thead>\n",
    "<tbody>\n", paste(rows, collapse = "\n"), "\n</tbody>\n</table>\n",
    cf_html
  )
}

#' @export
print.PeakReinterpretation <- function(x, ...) {
  cat("PeakReinterpretation\n")
  cat("  pheno:", x$meta$pheno %||% NA, "\n")
  cat("  query_id:", x$meta$query_id %||% NA, "\n")
  cat("  peaks:", length(x$peaks), "\n")
  cat("  views:", paste(x$meta$views %||% character(), collapse = ", "), "\n")
  invisible(x)
}

################################################################################
# Internal helpers

.reinterpret_resolve_query <- function(ranked, query, object) {
  if (inherits(query, "PhenotypeQuery")) {
    return(query)
  }
  if (is.character(query) && nzchar(trimws(query[1L]))) {
    return(phenotypeQuery(text = query[1L], use_llm = FALSE, save = FALSE))
  }
  pr <- attr(ranked, "phenotypeRank")
  qid <- pr$query_id %||% NULL
  if (!is.null(object) && inherits(object, "LazyGas") && !is.null(qid)) {
    q <- tryCatch(
      lazyData(object, dataset = "phenotype_query", kind = qid),
      error = function(e) NULL
    )
    if (inherits(q, "PhenotypeQuery")) {
      return(q)
    }
  }
  # Minimal synthetic query from empty keywords
  structure(
    list(
      query_id = qid %||% "pq_unknown",
      text = "",
      trait_keywords = character(),
      synonym_phrases = list(),
      related_phrases = list(),
      tissues = character(),
      developmental_stage = character(),
      conditions = character()
    ),
    class = c("PhenotypeQuery", "list")
  )
}

.reinterpret_strong_terms <- function(query, rules) {
  if (!is.null(rules$strong_terms) && length(rules$strong_terms)) {
    return(.filter_keyword_terms(as.character(rules$strong_terms)))
  }
  if (!inherits(query, "PhenotypeQuery")) {
    return(character())
  }
  syn <- .phenotype_linked_phrase_strings(query$synonym_phrases %||% list())
  .filter_keyword_terms(c(query$trait_keywords %||% character(), syn))
}

.reinterpret_ann_text <- function(df) {
  skip <- c(
    "peak_ID", "Gene_ID", "Gene_chr", "Gene_start", "dist2peak", "negLog10P",
    "HIGH", "MODERATE", "LOW", "MODIFIER",
    "HIGH_at_var", "MODERATE_at_var", "LOW_at_var", "MODIFIER_at_var",
    "score_annotation", "score_snpeff", "score_gwas", "score_expression",
    "score_literature", "score_ortholog", "score_finemap", "composite_score",
    "evidence_json", "query_id", "ann_text", "ann_specificity", "has_strong_hit",
    "dist_gate_pass", "soft_flag_weak_ann"
  )
  cols <- setdiff(names(df), skip)
  # Prefer known annotation-ish columns first
  prefer <- intersect(
    c(
      "RAP_Note", "MSU_Note", "EGGNog_Description", "MapMan_DESCRIPTION",
      "InterPro_Description", "Oryzabase_Trait_Ontology", "Oryzabase_geneSymbol",
      "RAPDB_geneSymbol", "CGSNL_geneSymbol", "Name", "Description", "description"
    ),
    cols
  )
  cols <- unique(c(prefer, cols))
  if (!length(cols)) {
    return(rep("", nrow(df)))
  }
  apply(df[, cols, drop = FALSE], 1L, function(row) {
    x <- as.character(row)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (!length(x)) {
      return("")
    }
    paste(x, collapse = " ")
  })
}

#' Annotation specificity for function-view gating (plan alpha)
#' @keywords internal
.ann_specificity <- function(text, strong, weak) {
  text <- as.character(text)[1L]
  if (!nzchar(trimws(text %||% ""))) {
    return(0)
  }
  strong <- as.character(strong %||% character())
  weak <- as.character(weak %||% character())
  hit_strong_phrase <- FALSE
  hit_strong <- FALSE
  for (term in strong) {
    if (!nzchar(term)) {
      next
    }
    if (isTRUE(.keyword_term_matches(text, term, ignore.case = TRUE))) {
      hit_strong <- TRUE
      parts <- strsplit(trimws(term), "\\s+", perl = TRUE)[[1L]]
      if (length(parts) >= 2L) {
        hit_strong_phrase <- TRUE
      }
    }
  }
  if (hit_strong_phrase) {
    return(1.0)
  }
  if (hit_strong) {
    return(0.67)
  }
  hit_weak <- FALSE
  for (term in weak) {
    if (!nzchar(term)) {
      next
    }
    if (isTRUE(.keyword_term_matches(text, term, ignore.case = TRUE))) {
      hit_weak <- TRUE
      break
    }
  }
  if (hit_weak) {
    return(0.2)
  }
  0
}

#' Distance gate vector
#' @keywords internal
.dist_gate <- function(dist2peak, rules) {
  n <- length(dist2peak)
  mode <- rules$dist_gate_mode %||% "off"
  if (identical(mode, "off")) {
    return(rep(TRUE, n))
  }
  d <- abs(suppressWarnings(as.numeric(dist2peak)))
  d[is.na(d)] <- Inf
  if (identical(mode, "absolute")) {
    lim <- as.numeric(rules$dist_abs_max)[1L]
    return(d <= lim)
  }
  # relative
  finite <- d[is.finite(d)]
  med <- if (length(finite)) stats::median(finite) else 0
  lim <- max(
    as.numeric(rules$dist_abs_max)[1L],
    as.numeric(rules$dist_iqr_mult)[1L] * med
  )
  d <= lim
}

.pick_champion <- function(df, view, rules) {
  empty <- data.frame(
    view = view,
    Gene_ID = NA_character_,
    reason_code = "no_eligible",
    reason_detail = .reinterpret_empty_detail(view),
    stringsAsFactors = FALSE
  )
  if (is.null(df) || !nrow(df)) {
    return(empty)
  }

  if (identical(view, "function")) {
    thr <- max(
      as.numeric(rules$min_specificity)[1L],
      if (isTRUE(rules$require_strong_hit)) 0.67 else 0
    )
    ok <- df$ann_specificity >= thr
    cand <- df[ok, , drop = FALSE]
    if (!nrow(cand)) {
      empty$reason_detail <- "no function hypothesis (no strong annotation hit)"
      return(empty)
    }
    o <- order(
      -cand$score_annotation,
      abs(cand$dist2peak),
      -cand$HIGH,
      cand$Gene_ID,
      na.last = TRUE
    )
    win <- cand[o[1L], , drop = FALSE]
    return(data.frame(
      view = view,
      Gene_ID = as.character(win$Gene_ID[1L]),
      reason_code = "max_ann_among_strong",
      reason_detail = sprintf(
        "ann_specificity=%.2f; score_annotation=%.3f",
        win$ann_specificity[1L],
        win$score_annotation[1L]
      ),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(view, "proximity")) {
    o <- order(
      abs(df$dist2peak),
      -df$HIGH,
      -df$score_annotation,
      df$Gene_ID,
      na.last = TRUE
    )
    win <- df[o[1L], , drop = FALSE]
    return(data.frame(
      view = view,
      Gene_ID = as.character(win$Gene_ID[1L]),
      reason_code = "min_abs_dist",
      reason_detail = sprintf("|dist|=%s", format(abs(win$dist2peak[1L]), scientific = FALSE)),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(view, "impact")) {
    hi <- df[df$HIGH > 0, , drop = FALSE]
    pool <- if (nrow(hi)) hi else df[df$MODERATE > 0, , drop = FALSE]
    if (!nrow(pool)) {
      empty$reason_detail <- "no HIGH/MODERATE variants"
      return(empty)
    }
    o <- order(
      -pool$HIGH,
      -pool$MODERATE,
      abs(pool$dist2peak),
      -pool$score_annotation,
      pool$Gene_ID,
      na.last = TRUE
    )
    win <- pool[o[1L], , drop = FALSE]
    return(data.frame(
      view = view,
      Gene_ID = as.character(win$Gene_ID[1L]),
      reason_code = if (win$HIGH[1L] > 0) "max_high_then_prox" else "max_moderate_then_prox",
      reason_detail = sprintf("HIGH=%s; MODERATE=%s", win$HIGH[1L], win$MODERATE[1L]),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(view, "expression")) {
    pool <- df[df$score_expression > 0, , drop = FALSE]
    if (!nrow(pool)) {
      empty$reason_detail <- "no positive expression score"
      return(empty)
    }
    o <- order(
      -pool$score_expression,
      abs(pool$dist2peak),
      pool$Gene_ID,
      na.last = TRUE
    )
    win <- pool[o[1L], , drop = FALSE]
    return(data.frame(
      view = view,
      Gene_ID = as.character(win$Gene_ID[1L]),
      reason_code = "max_expr",
      reason_detail = sprintf("score_expression=%.3f", win$score_expression[1L]),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(view, "composite")) {
    o <- order(-df$composite_score, df$Gene_ID, na.last = TRUE)
    win <- df[o[1L], , drop = FALSE]
    return(data.frame(
      view = view,
      Gene_ID = as.character(win$Gene_ID[1L]),
      reason_code = "max_composite",
      reason_detail = sprintf("composite_score=%.3f", win$composite_score[1L]),
      stringsAsFactors = FALSE
    ))
  }

  empty
}

.reinterpret_empty_detail <- function(view) {
  if (identical(view, "function")) {
    return("no function hypothesis (no strong annotation hit)")
  }
  paste0("no eligible gene for view ", view)
}

.reinterpret_support_rows <- function(df, view, rules, top_n) {
  ranked_view <- .reinterpret_order_for_view(df, view = view, rules = rules)
  if (!nrow(ranked_view)) {
    return(data.frame())
  }
  n <- min(as.integer(top_n)[1L], nrow(ranked_view))
  out <- ranked_view[seq_len(n), , drop = FALSE]
  out$view <- view
  out$rank_in_view <- seq_len(n)
  keep <- c(
    "view", "rank_in_view", "Gene_ID", "dist2peak", "HIGH", "MODERATE",
    "score_annotation", "score_snpeff", "score_gwas", "score_expression",
    "composite_score", "ann_specificity", "dist_gate_pass", "soft_flag_weak_ann"
  )
  keep <- intersect(keep, names(out))
  out[, keep, drop = FALSE]
}

.reinterpret_order_for_view <- function(df, view, rules) {
  if (identical(view, "function")) {
    thr <- max(
      as.numeric(rules$min_specificity)[1L],
      if (isTRUE(rules$require_strong_hit)) 0.67 else 0
    )
    ok <- df$ann_specificity >= thr
    cand <- df[ok, , drop = FALSE]
    if (!nrow(cand)) {
      return(cand)
    }
    cand[order(
      -cand$score_annotation,
      abs(cand$dist2peak),
      -cand$HIGH,
      cand$Gene_ID,
      na.last = TRUE
    ), , drop = FALSE]
  } else if (identical(view, "proximity")) {
    df[order(
      abs(df$dist2peak), -df$HIGH, -df$score_annotation, df$Gene_ID,
      na.last = TRUE
    ), , drop = FALSE]
  } else if (identical(view, "impact")) {
    hi <- df[df$HIGH > 0 | df$MODERATE > 0, , drop = FALSE]
    if (!nrow(hi)) {
      return(hi)
    }
    hi[order(
      -hi$HIGH, -hi$MODERATE, abs(hi$dist2peak), -hi$score_annotation, hi$Gene_ID,
      na.last = TRUE
    ), , drop = FALSE]
  } else if (identical(view, "expression")) {
    ex <- df[df$score_expression > 0, , drop = FALSE]
    if (!nrow(ex)) {
      return(ex)
    }
    ex[order(-ex$score_expression, abs(ex$dist2peak), ex$Gene_ID, na.last = TRUE), , drop = FALSE]
  } else {
    df[order(-df$composite_score, df$Gene_ID, na.last = TRUE), , drop = FALSE]
  }
}

.conflict_table <- function(champions, rules) {
  empty <- data.frame(
    type = character(),
    views = character(),
    Gene_IDs = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(rules$conflict_if_champions_differ)) {
    return(empty)
  }
  if (is.null(champions) || nrow(champions) < 2L) {
    return(empty)
  }
  ids <- as.character(champions$Gene_ID)
  # Treat NA as distinct sentinel for mismatch detection
  key <- ifelse(is.na(ids) | !nzchar(ids), "__NA__", ids)
  uniq <- unique(key)
  if (length(uniq) <= 1L) {
    return(empty)
  }
  data.frame(
    type = "champion_mismatch",
    views = paste(champions$view, collapse = "|"),
    Gene_IDs = paste(ifelse(is.na(ids) | !nzchar(ids), "NA", ids), collapse = "|"),
    note = "Hypothesis champions differ across views (listed side by side; no primary label).",
    stringsAsFactors = FALSE
  )
}

.gene_set_from_reinterpret <- function(champions, support, mode = "union_champions") {
  champ_ids <- unique(stats::na.omit(as.character(champions$Gene_ID)))
  champ_ids <- champ_ids[nzchar(champ_ids)]
  support_ids <- if (!is.null(support) && nrow(support) && "Gene_ID" %in% names(support)) {
    unique(stats::na.omit(as.character(support$Gene_ID)))
  } else {
    character()
  }
  support_ids <- support_ids[nzchar(support_ids)]
  unique(c(champ_ids, support_ids))
}

.store_write_phenotype_reinterpret <- function(object, reinterpret, pheno_name, query_id) {
  if (is.null(query_id) || !nzchar(as.character(query_id)[1L]) ||
      identical(as.character(query_id)[1L], "NA")) {
    query_id <- paste0("pq_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  query_id <- as.character(query_id)[1L]
  root <- .store_path(object)
  dir.create(file.path(root, "phenotype_reinterpret"), recursive = TRUE, showWarnings = FALSE)
  fname <- paste0(
    "reinterpret_",
    .store_safe_name(query_id),
    "_",
    .store_safe_name(pheno_name),
    ".rds"
  )
  path <- file.path(root, "phenotype_reinterpret", fname)
  saveRDS(reinterpret, path)
  meta <- .store_read_meta(object)
  prev <- meta$phenotype_reinterpret_latest
  if (is.null(prev)) {
    prev <- list()
  }
  prev[[pheno_name]] <- query_id
  meta$phenotype_reinterpret_latest <- prev
  .store_write_meta(object, meta)
  invisible(path)
}

.store_read_phenotype_reinterpret <- function(object, pheno_name, query_id = NULL) {
  root <- .store_path(object)
  dir <- file.path(root, "phenotype_reinterpret")
  if (!dir.exists(dir)) {
    return(NULL)
  }
  if (is.null(query_id) || !nzchar(as.character(query_id)[1L])) {
    meta <- .store_read_meta(object)
    query_id <- meta$phenotype_reinterpret_latest[[pheno_name]]
  }
  if (is.null(query_id) || !nzchar(as.character(query_id)[1L])) {
    return(NULL)
  }
  path <- file.path(
    dir,
    paste0(
      "reinterpret_",
      .store_safe_name(query_id),
      "_",
      .store_safe_name(pheno_name),
      ".rds"
    )
  )
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

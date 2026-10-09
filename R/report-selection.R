################################################################################
# Evidence-selection / HTML helpers used by reinterpret reports

#' Escape text for HTML attributes and body
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

#' HTML id for an Evidence gene heading
#' @keywords internal
.report_evidence_anchor_id <- function(gene_id) {
  gid <- as.character(gene_id)[1L]
  if (!nzchar(gid) || is.na(gid)) {
    return("evidence-gene-unknown")
  }
  paste0("evidence-gene-", gsub("[^A-Za-z0-9_-]", "_", gid))
}

#' Gene_ID cell: link to Evidence heading when that gene has Evidence
#' @keywords internal
.report_link_gene_id_html <- function(gene_id, evidence_gene_ids = NULL) {
  gid <- as.character(gene_id)[1L]
  if (!nzchar(gid) || is.na(gid)) {
    return(.report_escape_html("(none)"))
  }
  esc <- .report_escape_html(gid)
  ev <- unique(stats::na.omit(as.character(evidence_gene_ids %||% character())))
  if (length(ev) && gid %in% ev) {
    return(sprintf(
      "<a href=\"#%s\">%s</a>",
      .report_escape_html(.report_evidence_anchor_id(gid)),
      esc
    ))
  }
  esc
}

#' Map gene IDs to how they were selected for a report section
#' @keywords internal
.report_evidence_selection_map <- function(composite_ids = character(),
                                           card_ids = character()) {
  composite_ids <- unique(as.character(composite_ids))
  card_ids <- unique(as.character(card_ids))
  composite_ids <- composite_ids[nzchar(composite_ids)]
  card_ids <- card_ids[nzchar(card_ids)]
  both <- intersect(composite_ids, card_ids)
  out <- character()
  out[both] <- "both"
  out[setdiff(composite_ids, both)] <- "composite"
  out[setdiff(card_ids, both)] <- "card"
  out
}

#' Short HTML note describing why a gene was selected
#' @keywords internal
.report_evidence_selection_note_html <- function(gene_id,
                                                 selection_map = NULL,
                                                 language = "en") {
  if (is.null(selection_map) || !length(selection_map)) {
    return("")
  }
  gid <- as.character(gene_id)[1L]
  key <- unname(selection_map[gid])
  if (is.null(key) || !length(key) || is.na(key) || !nzchar(as.character(key)[1L])) {
    return("")
  }
  labels <- if (identical(language, "ja")) {
    list(
      selected_by = "選定理由",
      both = "総合スコア上位かつ仮説カード",
      composite = "総合スコア上位",
      card = "仮説カード"
    )
  } else {
    list(
      selected_by = "Selected by",
      both = "top composite score and hypothesis card",
      composite = "top composite score",
      card = "hypothesis card"
    )
  }
  detail <- switch(
    as.character(key)[1L],
    both = labels$both,
    composite = labels$composite,
    card = labels$card,
    NULL
  )
  if (is.null(detail)) {
    return("")
  }
  paste0(
    "<p><em>",
    labels$selected_by,
    ": ",
    detail,
    "</em></p>\n"
  )
}

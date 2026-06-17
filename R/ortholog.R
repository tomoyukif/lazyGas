################################################################################
#' Annotate candidate genes with ortholog information
#'
#' @param candidate A candidate-gene data.frame, or a \code{LazyGas} object.
#' @param object When \code{candidate} is a \code{LazyGas} object, it is used as
#'   \code{object}; otherwise ignored.
#' @param pheno Phenotype for [lazyData()] when reading candidates from
#'   \code{object}.
#' @param ortholog A data.frame with columns \code{gene_id}, \code{ortholog_id},
#'   and optionally \code{score}.
#' @param species_col,ortholog_col,score_col Column names in \code{ortholog}.
#' @param min_score Minimum ortholog score to retain.
#' @export
annotateOrthologs <- function(candidate = NULL,
                              object = NULL,
                              pheno = NULL,
                              ortholog,
                              species_col = "gene_id",
                              ortholog_col = "ortholog_id",
                              score_col = "score",
                              min_score = 0) {
  if (inherits(candidate, "LazyGas")) {
    object <- candidate
    candidate <- NULL
  }
  if (is.null(candidate)) {
    if (is.null(object)) {
      stop("Provide candidate or object.", call. = FALSE)
    }
    pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
    if (is.null(candidate)) {
      stop("No candidate data found.", call. = FALSE)
    }
  }
  req <- c(species_col, ortholog_col)
  if (!all(req %in% names(ortholog))) {
    stop("ortholog must contain columns: ", paste(req, collapse = ", "),
         call. = FALSE)
  }
  if (!"Gene_ID" %in% names(candidate)) {
    stop("candidate must contain Gene_ID.", call. = FALSE)
  }

  ortho <- ortholog
  if (score_col %in% names(ortho)) {
    ortho <- ortho[is.finite(ortho[[score_col]]) & ortho[[score_col]] >= min_score, , drop = FALSE]
  }
  names(ortho)[names(ortho) == species_col] <- "Gene_ID"
  keep <- c("Gene_ID", ortholog_col)
  if (score_col %in% names(ortho)) {
    keep <- c(keep, score_col)
  }
  ortho <- ortho[, keep, drop = FALSE]
  names(ortho)[names(ortho) == ortholog_col] <- "ortholog_id"
  if (score_col %in% names(ortho)) {
    names(ortho)[names(ortho) == score_col] <- "ortholog_score"
  }

  merged <- merge(candidate, ortho, by = "Gene_ID", all.x = TRUE)
  merged
}

#' Summarize ortholog conservation across candidate lists
#'
#' @param annotated Output of [annotateOrthologs()].
#' @export
summarizeOrthologMatches <- function(annotated) {
  if (!"ortholog_id" %in% names(annotated)) {
    stop("annotated data lacks ortholog_id column.", call. = FALSE)
  }
  sub <- annotated[!is.na(annotated$ortholog_id) & nzchar(annotated$ortholog_id), , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(data.frame())
  }
  if ("trait" %in% names(sub)) {
    out <- aggregate(
      ortholog_id ~ ortholog_id,
      data = sub,
      FUN = function(x) length(unique(x))
    )
    names(out)[2] <- "n_genes"
    trait_col <- aggregate(
      trait ~ ortholog_id,
      data = sub,
      FUN = function(x) paste(sort(unique(x)), collapse = ",")
    )
    out <- merge(out, trait_col, by = "ortholog_id")
  } else {
    out <- as.data.frame(table(sub$ortholog_id), stringsAsFactors = FALSE)
    names(out) <- c("ortholog_id", "n_genes")
  }
  out[order(out$n_genes, decreasing = TRUE), , drop = FALSE]
}

#' Plot ortholog overlap among annotated candidates
#'
#' @param annotated Output of [annotateOrthologs()].
#' @export
plotOrthologSummary <- function(annotated) {
  summ <- summarizeOrthologMatches(annotated)
  if (nrow(summ) == 0L) {
    stop("No ortholog matches to plot.", call. = FALSE)
  }
  top <- head(summ[order(summ$n_genes, decreasing = TRUE), , drop = FALSE], 20L)
  ggplot2::ggplot(top, ggplot2::aes(x = reorder(ortholog_id, n_genes), y = n_genes)) +
    ggplot2::geom_col(fill = "#4C78A8") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Shared orthologs among candidates",
      x = "Ortholog ID",
      y = "Number of candidate genes"
    ) +
    ggplot2::theme_bw()
}

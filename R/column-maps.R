################################################################################
#' Allowed annotation column classes for ranking / evidence
#' @keywords internal
.ANN_COLUMN_CLASSES <- c(
  "free_description",
  "go",
  "kegg",
  "pfam_interpro",
  "numeric",
  "other"
)

#' Annotation classes used as text for keyword ranking by default
#' @keywords internal
.ANN_RANKING_CLASSES <- c(
  "free_description",
  "go",
  "kegg",
  "pfam_interpro"
)

#' Read a manual annotation column map
#'
#' Loads a YAML or CSV defining how annotation table columns are interpreted.
#' Auto-detection and LLM classification are not used.
#'
#' Required columns: \code{column}, \code{class}.
#' Optional: \code{use_for_ranking} (logical; default TRUE for ranking classes),
#' \code{note}.
#'
#' \code{class} must be one of: \code{free_description}, \code{go}, \code{kegg},
#' \code{pfam_interpro}, \code{numeric}, \code{other}.
#'
#' @param path Path to \code{.yaml}/\code{.yml} or \code{.csv}.
#'
#' @return A \code{data.frame} with columns \code{column}, \code{class},
#'   \code{use_for_ranking}, \code{note}.
#' @export
#'
#' @examples
#' \dontrun{
#' map <- readAnnColumnMap(system.file("extdata", "demo_ann_column_map.csv",
#'                                     package = "lazyGas"))
#' }
readAnnColumnMap <- function(path) {
  path <- .column_map_require_path(path, "ann_column_map")
  raw <- .column_map_read_table(path)
  need <- c("column", "class")
  miss <- setdiff(need, names(raw))
  if (length(miss)) {
    stop(
      "ann_column_map must contain columns: ",
      paste(need, collapse = ", "),
      ". Missing: ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  df <- data.frame(
    column = trimws(as.character(raw$column)),
    class = tolower(trimws(as.character(raw$class))),
    stringsAsFactors = FALSE
  )
  if ("use_for_ranking" %in% names(raw)) {
    df$use_for_ranking <- .column_map_as_logical(raw$use_for_ranking)
  } else {
    df$use_for_ranking <- df$class %in% .ANN_RANKING_CLASSES
  }
  df$note <- if ("note" %in% names(raw)) {
    as.character(raw$note)
  } else {
    NA_character_
  }
  df <- df[nzchar(df$column), , drop = FALSE]
  if (!nrow(df)) {
    stop("ann_column_map has no rows.", call. = FALSE)
  }
  bad <- setdiff(unique(df$class), .ANN_COLUMN_CLASSES)
  if (length(bad)) {
    stop(
      "Unknown ann column class(es): ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(.ANN_COLUMN_CLASSES, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(df$column)) {
    stop(
      "Duplicate column entries in ann_column_map: ",
      paste(unique(df$column[duplicated(df$column)]), collapse = ", "),
      call. = FALSE
    )
  }
  df
}

#' Annotation column names used for ranking text
#' @param ann_column_map Output of [readAnnColumnMap()] or a path.
#' @keywords internal
.ann_ranking_columns <- function(ann_column_map) {
  map <- .coerce_ann_column_map(ann_column_map)
  cols <- map$column[as.logical(map$use_for_ranking)]
  cols <- cols[!is.na(cols) & nzchar(as.character(cols))]
  if (!length(cols)) {
    stop(
      "ann_column_map has no columns with use_for_ranking = TRUE.",
      call. = FALSE
    )
  }
  as.character(cols)
}

#' Validate ann map against a candidate / annotation table
#' @keywords internal
.validate_ann_column_map <- function(ann_column_map, table_cols) {
  map <- .coerce_ann_column_map(ann_column_map)
  ranking_cols <- .ann_ranking_columns(map)
  miss <- setdiff(ranking_cols, table_cols)
  if (length(miss)) {
    stop(
      "ann_column_map ranking column(s) not found in table: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(map)
}

.coerce_ann_column_map <- function(ann_column_map) {
  if (is.null(ann_column_map)) {
    stop(
      "'ann_column_map' is required when the annotation evidence channel is used. ",
      "Pass a path or the result of readAnnColumnMap().",
      call. = FALSE
    )
  }
  if (is.character(ann_column_map) && length(ann_column_map) == 1L) {
    return(readAnnColumnMap(ann_column_map))
  }
  if (is.data.frame(ann_column_map)) {
    if (!all(c("column", "class") %in% names(ann_column_map))) {
      stop(
        "ann_column_map data.frame must have columns 'column' and 'class'.",
        call. = FALSE
      )
    }
    if (!"use_for_ranking" %in% names(ann_column_map)) {
      ann_column_map$use_for_ranking <-
        ann_column_map$class %in% .ANN_RANKING_CLASSES
    }
    return(ann_column_map)
  }
  stop(
    "'ann_column_map' must be a file path or data.frame from readAnnColumnMap().",
    call. = FALSE
  )
}

#' Read a manual expression sample map
#'
#' Loads a YAML or CSV mapping expression-matrix sample columns to tissue /
#' stage / condition labels. Required when the expression evidence channel is
#' used. Sample IDs must match \code{colnames(expression_matrix)} exactly.
#'
#' Required column: \code{sample}. At least one of \code{tissue}, \code{stage},
#' \code{condition}, \code{organ} should be present (\code{tissue} recommended).
#' Optional file-level YAML key \code{matrix_scale}: \code{abundance} or
#' \code{row_zscore}.
#'
#' @param path Path to \code{.yaml}/\code{.yml} or \code{.csv}.
#'
#' @return A \code{data.frame} of sample metadata. Attribute
#'   \code{matrix_scale} may be set.
#' @export
readExpressionSampleMap <- function(path) {
  path <- .column_map_require_path(path, "expression_sample_map")
  ext <- tolower(tools::file_ext(path))
  matrix_scale <- "abundance"
  if (ext %in% c("yml", "yaml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required to read YAML column maps.", call. = FALSE)
    }
    obj <- yaml::read_yaml(path)
    if (is.list(obj) && !is.null(obj$matrix_scale)) {
      matrix_scale <- as.character(obj$matrix_scale)[1L]
    }
    if (is.list(obj) && !is.null(obj$samples)) {
      raw <- .column_map_list_to_df(obj$samples)
    } else if (is.list(obj) && !is.null(obj$columns)) {
      # allow same shape as ann map by mistake? no — require samples
      stop(
        "expression_sample_map YAML must contain a top-level 'samples' list.",
        call. = FALSE
      )
    } else {
      raw <- .column_map_list_to_df(obj)
    }
  } else {
    raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!"sample" %in% names(raw)) {
    stop("expression_sample_map must contain a 'sample' column.", call. = FALSE)
  }
  raw$sample <- trimws(as.character(raw$sample))
  raw <- raw[nzchar(raw$sample), , drop = FALSE]
  if (!nrow(raw)) {
    stop("expression_sample_map has no rows.", call. = FALSE)
  }
  if (anyDuplicated(raw$sample)) {
    stop(
      "Duplicate sample entries in expression_sample_map: ",
      paste(unique(raw$sample[duplicated(raw$sample)]), collapse = ", "),
      call. = FALSE
    )
  }
  axis_cols <- intersect(names(raw), c("tissue", "stage", "condition", "organ"))
  if (!length(axis_cols)) {
    stop(
      "expression_sample_map must include at least one of: tissue, stage, condition, organ.",
      call. = FALSE
    )
  }
  matrix_scale <- match.arg(matrix_scale, c("abundance", "row_zscore"))
  attr(raw, "matrix_scale") <- matrix_scale
  raw
}

#' Validate expression matrix against sample map
#' @keywords internal
.validate_expression_inputs <- function(expression_matrix, expression_sample_map) {
  if (is.null(expression_matrix)) {
    stop(
      "'expression_matrix' is required when the expression evidence channel is used.",
      call. = FALSE
    )
  }
  map <- .coerce_expression_sample_map(expression_sample_map)
  mat <- as.matrix(expression_matrix)
  if (is.null(colnames(mat)) || !length(colnames(mat))) {
    stop("'expression_matrix' must have column names (sample IDs).", call. = FALSE)
  }
  miss <- setdiff(colnames(mat), map$sample)
  extra <- setdiff(map$sample, colnames(mat))
  if (length(miss)) {
    stop(
      "expression_matrix column(s) missing from expression_sample_map: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(extra)) {
    warning(
      "expression_sample_map sample(s) not in expression_matrix (ignored): ",
      paste(extra, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(map)
}

.coerce_expression_sample_map <- function(expression_sample_map) {
  if (is.null(expression_sample_map)) {
    stop(
      "'expression_sample_map' is required when the expression evidence channel is used. ",
      "Pass a path or the result of readExpressionSampleMap().",
      call. = FALSE
    )
  }
  if (is.character(expression_sample_map) && length(expression_sample_map) == 1L) {
    return(readExpressionSampleMap(expression_sample_map))
  }
  if (is.data.frame(expression_sample_map)) {
    if (!"sample" %in% names(expression_sample_map)) {
      # accept legacy expression_meta with rownames or sample column already
      if (!is.null(rownames(expression_sample_map)) &&
          !"sample" %in% names(expression_sample_map)) {
        expression_sample_map$sample <- rownames(expression_sample_map)
      } else {
        stop(
          "expression_sample_map data.frame must have a 'sample' column.",
          call. = FALSE
        )
      }
    }
    return(expression_sample_map)
  }
  stop(
    "'expression_sample_map' must be a file path or data.frame from readExpressionSampleMap().",
    call. = FALSE
  )
}

.column_map_require_path <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || !nzchar(trimws(path))) {
    stop("'", label, "' path must be a non-empty string.", call. = FALSE)
  }
  path <- path.expand(trimws(path))
  if (!file.exists(path)) {
    stop("'", label, "' file not found: ", path, call. = FALSE)
  }
  path
}

.column_map_read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yml", "yaml")) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required to read YAML column maps.", call. = FALSE)
    }
    obj <- yaml::read_yaml(path)
    if (is.list(obj) && !is.null(obj$columns)) {
      return(.column_map_list_to_df(obj$columns))
    }
    return(.column_map_list_to_df(obj))
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.column_map_list_to_df <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.list(x) || !length(x)) {
    stop("Column map YAML is empty or invalid.", call. = FALSE)
  }
  # named list of class per column
  if (!is.null(names(x)) && all(vapply(x, function(z) {
    is.character(z) || is.list(z)
  }, logical(1)))) {
    # detect list-of-rows vs named map
    if (is.list(x[[1L]]) && !is.null(x[[1L]]$column)) {
      return(do.call(rbind, lapply(x, function(row) {
        as.data.frame(row, stringsAsFactors = FALSE)
      })))
    }
    if (is.list(x[[1L]]) && !is.null(x[[1L]]$sample)) {
      return(do.call(rbind, lapply(x, function(row) {
        as.data.frame(row, stringsAsFactors = FALSE)
      })))
    }
    # flat column -> class
    if (all(vapply(x, function(z) is.character(z) || is.numeric(z), logical(1)))) {
      return(data.frame(
        column = names(x),
        class = as.character(unlist(x, use.names = FALSE)),
        stringsAsFactors = FALSE
      ))
    }
  }
  do.call(rbind, lapply(x, function(row) {
    as.data.frame(row, stringsAsFactors = FALSE)
  }))
}

.column_map_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x))
  out[x %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[x %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (anyNA(out)) {
    stop("Invalid logical values in use_for_ranking.", call. = FALSE)
  }
  out
}

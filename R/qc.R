################################################################################
#' Genomic inflation factor (lambda GC)
#'
#' @param p_values Numeric vector of p-values (NA omitted).
#' @export
calcGenomicInflation <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  if (length(p_values) < 2L) {
    return(NA_real_)
  }
  chisq <- qchisq(1 - p_values, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

#' QQ plot for GWAS p-values
#'
#' @param object A \code{LazyGas} object, or a numeric p-value vector.
#' @param pheno Phenotype name or index when \code{object} is \code{LazyGas}.
#' @param max_points Maximum points drawn (full data still used for
#'   \eqn{\lambda_{GC}}). Large genome-wide scans are thinned for plotting.
#' @export
plotQQ <- function(object, pheno = NULL, max_points = 100000L) {
  if (inherits(object, "LazyGas")) {
    pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
    scan_df <- lazyData(object = object, dataset = "scan", pheno = pheno_name)
    if (is.null(scan_df)) {
      stop("No scan data for phenotype: ", pheno_name, call. = FALSE)
    }
    p_values <- scan_df$P.model
    title <- paste0("QQ plot: ", pheno_name)
  } else {
    p_values <- object
    title <- "QQ plot"
  }
  p_values <- p_values[is.finite(p_values) & p_values > 0 & p_values <= 1]
  if (length(p_values) == 0L) {
    stop("No valid p-values to plot.", call. = FALSE)
  }
  n <- length(p_values)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(p_values))
  lambda <- calcGenomicInflation(p_values)

  max_points <- as.integer(max_points)[1L]
  if (is.finite(max_points) && max_points > 0L && n > max_points) {
    # Keep the most extreme observed points; thin the remainder evenly.
    keep_tail <- min(max(2000L, as.integer(max_points * 0.2)), n)
    tail_idx <- seq.int(n - keep_tail + 1L, n)
    body_n <- max_points - length(tail_idx)
    body_idx <- if (body_n > 0L) {
      unique(as.integer(round(seq(1, n - keep_tail, length.out = body_n))))
    } else {
      integer()
    }
    idx <- sort(unique(c(body_idx, tail_idx)))
    expected <- expected[idx]
    observed <- observed[idx]
  }

  ggplot2::ggplot(
    data = data.frame(expected = expected, observed = observed),
    mapping = ggplot2::aes(x = expected, y = observed)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("lambda GC = %.3f (n = %s)", lambda, format(n, big.mark = ",")),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    ggplot2::theme_bw()
}

#' Summarize GWAS quality metrics
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index; if \code{NULL}, all phenotypes.
#' @param small_n_threshold Sample size below which a warning flag is set.
#' @param store If \code{TRUE}, write results to the companion store.
#' @export
summarizeGWASQC <- function(object,
                            pheno = NULL,
                            small_n_threshold = 30L,
                            store = TRUE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  pheno_info <- getPheno(object = object)
  if (is.null(pheno)) {
    pheno_names <- pheno_info$pheno_names
  } else {
    pheno_names <- .determine_phenotype_name(object = object, pheno = pheno)
  }

  out <- lapply(pheno_names, function(pn) {
    scan_df <- lazyData(object = object, dataset = "scan", pheno = pn)
    if (is.null(scan_df)) {
      return(data.frame(
        pheno = pn,
        n_samples = NA_integer_,
        n_markers = NA_integer_,
        n_valid_p = NA_integer_,
        lambda_gc = NA_real_,
        max_negLog10P = NA_real_,
        small_n_warning = NA,
        stringsAsFactors = FALSE
      ))
    }
    p <- scan_df$P.model
    n_valid <- sum(is.finite(p) & p > 0 & p <= 1, na.rm = TRUE)
    n_sam <- sum(!is.na(pheno_info$pheno[[pn]]))
    data.frame(
      pheno = pn,
      n_samples = n_sam,
      n_markers = nrow(scan_df),
      n_valid_p = n_valid,
      lambda_gc = calcGenomicInflation(p),
      max_negLog10P = max(scan_df$negLog10P, na.rm = TRUE),
      small_n_warning = n_sam < small_n_threshold,
      stringsAsFactors = FALSE
    )
  })
  qc_df <- do.call(rbind, out)
  rownames(qc_df) <- NULL

  if (store) {
    .store_write_section_df(object, "qc", "summary", qc_df, pheno_name = NULL)
  }
  qc_df
}

################################################################################
#' Conditional association within a peak block
#'
#' Re-scans markers in a peak while conditioning on the lead variant (and
#' optionally additional independent signals via stepwise conditioning).
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param peak_id Peak block ID.
#' @param recalc Use recalculated peak blocks.
#' @param k_step Maximum number of conditioning steps.
#' @param store If \code{TRUE}, write results to the companion store.
#' @export
conditionalAssoc <- function(object,
                             pheno,
                             peak_id,
                             recalc = TRUE,
                             k_step = 3L,
                             store = TRUE) {
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  block <- .conditional_peak_block(
    object = object,
    pheno_name = pheno_name,
    peak_id = peak_id,
    recalc = recalc
  )
  if (nrow(block) == 0L) {
    stop("No markers found for peak_id = ", peak_id, call. = FALSE)
  }

  pheno_vec <- getPheno(object)$pheno[[pheno_name]]
  geno_format <- .store_read_scan_scalar(object, "geno_format")
  if (is.null(geno_format)) {
    geno_format <- "dosage"
  }
  formula_txt <- .store_read_scan_scalar(object, "formula")
  if (is.null(formula_txt)) {
    formula_txt <- "add"
  }
  conv_fun_txt <- .store_read_scan_scalar(object, "conv_fun")
  conv_fun <- if (!is.null(conv_fun_txt) && any(nzchar(as.character(conv_fun_txt)))) {
    eval(parse(text = paste(as.character(conv_fun_txt), collapse = "\n")))
  } else {
    makeConvFun(geno_format = geno_format, n_levels = 3L)
  }
  binary <- getPheno(object)$pheno_type$binary[
    match(pheno_name, getPheno(object)$pheno_names)
  ]

  conditioned_on <- character()
  results <- list()
  remaining <- block$variant_ID

  for (step in seq_len(k_step)) {
    if (length(remaining) == 0L) {
      break
    }
    fits <- lapply(remaining, function(vid) {
      .conditional_fit_marker(
        object = object,
        variant_id = vid,
        pheno_vec = pheno_vec,
        geno_format = geno_format,
        conv_fun = conv_fun,
        formula_txt = formula_txt,
        binary = binary,
        conditioned_on = conditioned_on
      )
    })
    fit_df <- do.call(rbind, fits)
    fit_df <- fit_df[order(fit_df$P_conditional), , drop = FALSE]
    fit_df$step <- step
    fit_df$conditioned_on <- paste(conditioned_on, collapse = ",")
    results[[step]] <- fit_df

    lead <- fit_df$variant_ID[1]
    if (!is.finite(fit_df$P_conditional[1]) || fit_df$P_conditional[1] >= 0.05) {
      break
    }
    conditioned_on <- c(conditioned_on, as.character(lead))
    remaining <- setdiff(remaining, lead)
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  if (store) {
    .store_write_section_df(
      object, "conditional", paste0("peak_", peak_id), out, pheno_name
    )
  }
  out
}

#' Calculate credible sets within a peak using Wakefield ABF
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param peak_id Peak block ID.
#' @param coverage Target cumulative PIP for the credible set.
#' @param prior_W,prior_V Wakefield ABF prior parameters.
#' @param use_conditional If \code{TRUE}, run [conditionalAssoc()] first.
#' @param recalc Use recalculated peak blocks.
#' @param store If \code{TRUE}, write results to the companion store.
#' @export
calcCredibleSet <- function(object,
                            pheno,
                            peak_id,
                            coverage = 0.95,
                            prior_W = 0.15,
                            prior_V = 0.04,
                            use_conditional = TRUE,
                            recalc = TRUE,
                            store = TRUE) {
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  if (use_conditional) {
    cond <- conditionalAssoc(
      object = object,
      pheno = pheno_name,
      peak_id = peak_id,
      recalc = recalc,
      store = FALSE
    )
    use_df <- cond[cond$step == min(cond$step), , drop = FALSE]
    beta <- use_df$beta
    se <- use_df$se
    variant_ids <- use_df$variant_ID
    p_col <- use_df$P_conditional
  } else {
    block <- .conditional_peak_block(
      object = object,
      pheno_name = pheno_name,
      peak_id = peak_id,
      recalc = recalc
    )
    scan_df <- lazyData(object = object, dataset = "scan", pheno = pheno_name)
    hit <- scan_df$variant_ID %in% block$variant_ID
    use_df <- scan_df[hit, , drop = FALSE]
    coef_cols <- grep("^Coef\\.", names(use_df), value = TRUE)
    if (length(coef_cols) == 0L) {
      stop("Scan results lack Coef.* columns for ABF.", call. = FALSE)
    }
    beta <- use_df[[coef_cols[1]]]
    se <- .abf_se_from_p(beta, use_df$P.model)
    variant_ids <- use_df$variant_ID
    p_col <- use_df$P.model
  }

  abf <- .wakefield_abf(beta = beta, se = se, prior_V = prior_V, prior_W = prior_W)
  pip <- abf / sum(abf, na.rm = TRUE)
  ord <- order(pip, decreasing = TRUE)
  pip_ord <- pip[ord]
  cum_pip <- cumsum(pip_ord)
  in_set <- cum_pip <= coverage
  if (any(in_set)) {
    last_in <- max(which(in_set))
  } else {
    last_in <- 1L
  }
  in_set_idx <- ord[seq_len(last_in)]

  out <- data.frame(
    variant_ID = variant_ids,
    beta = beta,
    se = se,
    P = p_col,
    ABF = abf,
    PIP = pip,
    in_credible_set = seq_along(variant_ids) %in% in_set_idx,
    credible_set_id = 1L,
    cumulative_PIP = pip,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$PIP, decreasing = TRUE), , drop = FALSE]
  attr(out, "summary") <- list(
    n_variants_in_set = sum(out$in_credible_set),
    set_coverage = sum(out$PIP[out$in_credible_set]),
    lead_variant_ID = out$variant_ID[which.max(out$PIP)]
  )

  if (store) {
    .store_write_section_df(
      object, "credible_set", paste0("peak_", peak_id), out, pheno_name
    )
  }
  out
}

#' Plot credible-set PIP values for a peak
#'
#' @param object A \code{LazyGas} object or a data.frame from [calcCredibleSet()].
#' @param pheno,peak_id Required when \code{object} is a \code{LazyGas} object.
#' @param ... Passed to [calcCredibleSet()] when needed.
#' @export
plotCredibleSet <- function(object, pheno = NULL, peak_id = NULL, ...) {
  if (inherits(object, "LazyGas")) {
    cred <- calcCredibleSet(object = object, pheno = pheno, peak_id = peak_id, ...)
  } else {
    cred <- object
  }
  if (nrow(cred) == 0L) {
    stop("No credible-set data to plot.", call. = FALSE)
  }
  cred$variant_ID <- as.factor(cred$variant_ID)
  ggplot2::ggplot(cred, ggplot2::aes(x = variant_ID, y = PIP, fill = in_credible_set)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#E45756", "FALSE" = "#72B7B2")) +
    ggplot2::labs(
      title = "Credible set (PIP)",
      x = "Variant",
      y = "Posterior inclusion probability",
      fill = "In credible set"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

.conditional_peak_block <- function(object, pheno_name, peak_id, recalc) {
  use_recalc <- recalc && .store_section_exists(object, "recalc")
  pk <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = use_recalc)
  if (is.null(pk)) {
    return(data.frame())
  }
  sub <- pk[pk$peak_ID == peak_id, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(data.frame())
  }
  unique(sub[, c("variant_ID", "peak_ID", "Chr", "Pos"), drop = FALSE])
}

.conditional_fit_marker <- function(object,
                                 variant_id,
                                 pheno_vec,
                                 geno_format,
                                 conv_fun,
                                 formula_txt,
                                 binary,
                                 conditioned_on) {
  g <- getGenoPerMarker(
    object = object,
    geno_format = geno_format,
    marker_index = .variant_id_to_marker_index(object, variant_id)
  )
  if (geno_format == "dosage") {
    g <- as.numeric(g)
  } else if (is.matrix(g)) {
    g <- colSums(g == 1, na.rm = TRUE)
  }
  df <- .makeDF(g = g, phe = pheno_vec, conv_fun = conv_fun, formula = formula_txt)
  if (length(conditioned_on) > 0L) {
    cov_mat <- vapply(conditioned_on, function(vid) {
      gv <- getGenoPerMarker(
        object = object,
        geno_format = geno_format,
        marker_index = .variant_id_to_marker_index(object, vid)
      )
      if (geno_format == "dosage") as.numeric(gv) else colSums(gv == 1, na.rm = TRUE)
    }, numeric(length(pheno_vec)))
    if (is.vector(cov_mat)) {
      cov_mat <- matrix(cov_mat, ncol = 1L)
    }
    colnames(cov_mat) <- paste0("cond_", conditioned_on)
    df$df <- cbind(df$df, cov_mat)
    formula_txt <- paste(c(formula_txt, colnames(cov_mat)), collapse = " + ")
    df$fml <- stats::formula(paste0("phe ~ ", formula_txt))
  }
  family <- if (binary) "binomial" else "gaussian"
  fit <- try(stats::glm(df$fml, data = df$df, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(
      variant_ID = variant_id,
      beta = NA_real_,
      se = NA_real_,
      P_conditional = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  sm <- summary(fit)
  coef_nm <- setdiff(rownames(sm$coefficients), "(Intercept)")
  if (length(coef_nm) == 0L) {
    return(data.frame(
      variant_ID = variant_id,
      beta = NA_real_,
      se = NA_real_,
      P_conditional = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  term <- coef_nm[1]
  data.frame(
    variant_ID = variant_id,
    beta = sm$coefficients[term, 1],
    se = sm$coefficients[term, 2],
    P_conditional = sm$coefficients[term, 4],
    stringsAsFactors = FALSE
  )
}

.wakefield_abf <- function(beta, se, prior_V, prior_W) {
  z <- beta / se
  r <- prior_V / (prior_V + se^2)
  sqrt(r) * exp(z^2 * r / 2) / prior_W
}

.abf_se_from_p <- function(beta, p) {
  z <- abs(stats::qnorm(p / 2))
  z[z == 0] <- NA
  abs(beta / z)
}

.variant_id_to_marker_index <- function(object, variant_id) {
  ids <- getMarID(object = object, valid = TRUE)
  idx <- match(variant_id, ids)
  if (is.na(idx)) {
    stop("variant_ID not found in valid markers: ", variant_id, call. = FALSE)
  }
  idx
}

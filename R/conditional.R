################################################################################
#' Conditional association within a peak block
#'
#' Re-scans markers in a peak while conditioning on the lead variant (and
#' optionally additional independent signals via stepwise conditioning).
#' Use [summarizeConditionalSignals()] to count independent associations and
#' [summarizeFineMapping()] for an interpretation-focused peak summary.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param peak_id Peak block ID.
#' @param recalc Use recalculated peak blocks.
#' @param k_step Maximum number of conditioning steps.
#' @param p_threshold Stop adding conditioning leads when the best conditional
#'   \code{P_conditional} is not below this value.
#' @param store If \code{TRUE}, write results to the companion store.
#' @export
conditionalAssoc <- function(object,
                             pheno,
                             peak_id,
                             recalc = TRUE,
                             k_step = 3L,
                             p_threshold = 0.05,
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
    if (!is.finite(fit_df$P_conditional[1]) || fit_df$P_conditional[1] >= p_threshold) {
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
#' Approximate Bayes factors are used to rank variants and form a credible set.
#' Large credible sets indicate limited resolution (often due to LD, sample size,
#' or weak signals) rather than confident causal-variant identification.
#' See [summarizeFineMapping()] for interpretation-ready summaries.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param peak_id Peak block ID.
#' @param coverage Target cumulative PIP for the credible set.
#' @param prior_W,prior_V Wakefield ABF prior parameters.
#' @param use_conditional If \code{TRUE}, use step-1 results from
#'   [conditionalAssoc()] (unconditioned on other leads).
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
  block <- .conditional_peak_block(
    object = object,
    pheno_name = pheno_name,
    peak_id = peak_id,
    recalc = recalc
  )
  n_locus <- nrow(block)

  marginal_lead <- .conditional_marginal_lead(
    object = object,
    pheno_name = pheno_name,
    block = block
  )

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
  cum_pip_ord <- cumsum(pip_ord)
  if (any(cum_pip_ord <= coverage)) {
    last_in <- max(which(cum_pip_ord <= coverage))
  } else {
    last_in <- 1L
  }
  in_set_idx <- ord[seq_len(last_in)]

  cum_pip <- numeric(length(variant_ids))
  cum_pip[ord] <- cum_pip_ord

  out <- data.frame(
    variant_ID = variant_ids,
    beta = beta,
    se = se,
    P = p_col,
    ABF = abf,
    PIP = pip,
    in_credible_set = seq_along(variant_ids) %in% in_set_idx,
    credible_set_id = 1L,
    cumulative_PIP = cum_pip,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$PIP, decreasing = TRUE), , drop = FALSE]

  top_pip_variant <- out$variant_ID[which.max(out$PIP)]
  n_in_set <- sum(out$in_credible_set)
  set_coverage <- sum(out$PIP[out$in_credible_set])
  interp <- .interpret_credible_set(
    n_in_set = n_in_set,
    n_locus = n_locus,
    top_pip_variant = top_pip_variant,
    marginal_lead = marginal_lead
  )

  attr(out, "summary") <- c(
    list(
      n_variants_in_set = n_in_set,
      n_variants_in_locus = n_locus,
      set_coverage = set_coverage,
      lead_variant_ID = top_pip_variant,
      marginal_lead_variant_ID = marginal_lead,
      top_pip_is_marginal_lead = identical(as.character(top_pip_variant), as.character(marginal_lead)),
      credible_set_fraction = if (n_locus > 0L) n_in_set / n_locus else NA_real_,
      resolution = interp$resolution
    ),
    interp
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
  summ <- attr(cred, "summary")
  subtitle <- if (!is.null(summ$message)) summ$message else NULL
  cred$variant_ID <- as.factor(cred$variant_ID)
  p <- ggplot2::ggplot(cred, ggplot2::aes(x = variant_ID, y = PIP, fill = in_credible_set)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#E45756", "FALSE" = "#72B7B2"),
      labels = c("TRUE" = "In credible set", "FALSE" = "Outside set")
    ) +
    ggplot2::labs(
      title = "Credible set (PIP)",
      subtitle = subtitle,
      x = "Variant",
      y = "Posterior inclusion probability",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  p
}

#' Summarize independent signals from conditional association
#'
#' Counts stepwise leads and, when a \code{LazyGas} object is supplied, classifies
#' them with [assessIndependentSignals()] so that \code{n_independent_signals}
#' reflects LD-supported signals rather than raw stepwise picks.
#'
#' @param object A \code{LazyGas} object, or output from [conditionalAssoc()].
#' @param pheno,peak_id Required when \code{object} is a \code{LazyGas} object.
#' @param p_threshold P-value threshold for calling an independent lead.
#' @param r2_threshold Passed to [assessIndependentSignals()] when \code{object}
#'   is a \code{LazyGas} object.
#' @param recalc Use recalculated peak blocks when running [conditionalAssoc()].
#' @param store If \code{TRUE}, write the summary row to the companion store.
#' @export
summarizeConditionalSignals <- function(object,
                                        pheno = NULL,
                                        peak_id = NULL,
                                        p_threshold = 0.05,
                                        r2_threshold = 0.8,
                                        recalc = TRUE,
                                        store = FALSE) {
  if (inherits(object, "LazyGas")) {
    assessed <- assessIndependentSignals(
      object = object,
      pheno = pheno,
      peak_id = peak_id,
      p_threshold = p_threshold,
      r2_threshold = r2_threshold,
      recalc = recalc,
      store = FALSE
    )
    out <- assessed$summary
    pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  } else {
    cond <- object
    if (is.null(peak_id)) {
      stop("peak_id is required when object is a conditional result data.frame.",
           call. = FALSE)
    }
    if (nrow(cond) == 0L) {
      return(data.frame())
    }
    out <- .conditional_summary_from_assessment(
      peak_id = peak_id,
      n_variants_in_locus = length(unique(cond$variant_ID)),
      assessment = .assess_independent_leads(
        object = NULL,
        cond = cond,
        p_threshold = p_threshold,
        r2_threshold = r2_threshold
      ),
      p_threshold = p_threshold,
      r2_threshold = r2_threshold
    )
    pheno_name <- NA_character_
  }

  if (store && inherits(object, "LazyGas")) {
    .store_write_section_df(
      object,
      "fine_mapping",
      paste0("conditional_peak_", peak_id),
      out,
      pheno_name
    )
  }
  out
}

#' Assess whether stepwise conditional leads are independent or LD artifacts
#'
#' Classifies each stepwise lead from [conditionalAssoc()] using pairwise \(r^2\),
#' near-identical step-1 effect estimates, and joint-model p-values. Only leads
#' labelled \code{supported} count toward independent signals in summary APIs.
#'
#' @param object A \code{LazyGas} object (required for genotype-based checks).
#' @param pheno Phenotype name or index.
#' @param peak_id Peak block ID.
#' @param p_threshold P-value threshold for stepwise leads and joint models.
#' @param r2_threshold Markers with pairwise \(r^2\) at or above this value are
#'   treated as LD artifacts of an already-supported lead.
#' @param recalc Use recalculated peak blocks.
#' @param store If \code{TRUE}, write per-lead / pairwise tables to the store.
#' @return A list with \code{pairs}, \code{leads}, and \code{summary}.
#' @export
assessIndependentSignals <- function(object,
                                     pheno = NULL,
                                     peak_id = NULL,
                                     p_threshold = 0.05,
                                     r2_threshold = 0.8,
                                     recalc = TRUE,
                                     store = FALSE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object for independent-signal assessment.",
         call. = FALSE)
  }
  if (is.null(peak_id)) {
    stop("peak_id is required.", call. = FALSE)
  }
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  cond <- conditionalAssoc(
    object = object,
    pheno = pheno_name,
    peak_id = peak_id,
    recalc = recalc,
    p_threshold = p_threshold,
    store = FALSE
  )
  if (nrow(cond) == 0L) {
    empty <- .empty_independent_assessment(peak_id = peak_id,
                                           p_threshold = p_threshold,
                                           r2_threshold = r2_threshold)
    return(empty)
  }

  attr(cond, "pheno_name") <- pheno_name
  assessment <- .assess_independent_leads(
    object = object,
    cond = cond,
    p_threshold = p_threshold,
    r2_threshold = r2_threshold
  )
  summary_df <- .conditional_summary_from_assessment(
    peak_id = peak_id,
    n_variants_in_locus = length(unique(cond$variant_ID)),
    assessment = assessment,
    p_threshold = p_threshold,
    r2_threshold = r2_threshold
  )

  if (store) {
    .store_write_section_df(
      object,
      "fine_mapping",
      paste0("assess_leads_peak_", peak_id),
      assessment$leads,
      pheno_name
    )
    if (nrow(assessment$pairs) > 0L) {
      .store_write_section_df(
        object,
        "fine_mapping",
        paste0("assess_pairs_peak_", peak_id),
        assessment$pairs,
        pheno_name
      )
    }
    .store_write_section_df(
      object,
      "fine_mapping",
      paste0("conditional_peak_", peak_id),
      summary_df,
      pheno_name
    )
  }

  list(pairs = assessment$pairs, leads = assessment$leads, summary = summary_df)
}

#' Run conditional association and credible-set analysis for peaks
#'
#' Stepwise leads are classified with [assessIndependentSignals()] so that
#' \code{n_independent_signals} counts only LD-supported associations.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index. If \code{NULL}, all assigned phenotypes.
#' @param peak_ids Peak IDs to analyse. Default: all peaks for each phenotype.
#' @param recalc Use recalculated peak blocks.
#' @param k_step,p_threshold Passed to [conditionalAssoc()].
#' @param r2_threshold Passed to [assessIndependentSignals()].
#' @param coverage,prior_W,prior_V Passed to [calcCredibleSet()].
#' @param store If \code{TRUE}, write per-peak tables and a summary table.
#' @export
runFineMapping <- function(object,
                           pheno = NULL,
                           peak_ids = NULL,
                           recalc = TRUE,
                           k_step = 3L,
                           p_threshold = 0.05,
                           r2_threshold = 0.8,
                           coverage = 0.95,
                           prior_W = 0.15,
                           prior_V = 0.04,
                           store = TRUE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  pheno_names <- if (is.null(pheno)) {
    getPheno(object)$pheno_names
  } else {
    .determine_phenotype_name(object = object, pheno = pheno)
  }

  summaries <- lapply(pheno_names, function(pn) {
    ids <- peak_ids
    if (is.null(ids)) {
      ids <- .peak_ids_for_pheno(object = object, pheno_name = pn, recalc = recalc)
    }
    if (length(ids) == 0L) {
      return(NULL)
    }
    rows <- lapply(ids, function(pid) {
      cond <- conditionalAssoc(
        object = object,
        pheno = pn,
        peak_id = pid,
        recalc = recalc,
        k_step = k_step,
        p_threshold = p_threshold,
        store = store
      )
      cred <- calcCredibleSet(
        object = object,
        pheno = pn,
        peak_id = pid,
        coverage = coverage,
        prior_W = prior_W,
        prior_V = prior_V,
        use_conditional = TRUE,
        recalc = recalc,
        store = store
      )
      .fine_mapping_summary_row(
        object = object,
        pheno_name = pn,
        peak_id = pid,
        cond = cond,
        cred = cred,
        p_threshold = p_threshold,
        r2_threshold = r2_threshold
      )
    })
    do.call(rbind, rows)
  })

  summary_df <- do.call(rbind, summaries)
  if (is.null(summary_df)) {
    summary_df <- data.frame()
  } else {
    rownames(summary_df) <- NULL
  }

  if (store && nrow(summary_df) > 0L) {
    for (pn in unique(summary_df$pheno)) {
      sub <- summary_df[summary_df$pheno == pn, , drop = FALSE]
      .store_write_section_df(object, "fine_mapping", "summary", sub, pn)
    }
  }

  invisible(list(summary = summary_df, object = object))
}

#' Interpretation-focused fine-mapping summary per peak
#'
#' Combines conditional-signal counts and credible-set resolution metrics.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param peak_ids Peak IDs. Default: all peaks.
#' @param recalc Use recalculated peak blocks.
#' @param run_if_missing If \code{TRUE}, call [runFineMapping()] when stored
#'   summaries are absent.
#' @param store If \code{TRUE}, write results to the companion store.
#' @export
summarizeFineMapping <- function(object,
                                 pheno,
                                 peak_ids = NULL,
                                 recalc = TRUE,
                                 run_if_missing = TRUE,
                                 store = TRUE) {
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  if (is.null(peak_ids)) {
    peak_ids <- .peak_ids_for_pheno(object = object, pheno_name = pheno_name, recalc = recalc)
  }

  stored <- tryCatch(
    lazyData(object = object, dataset = "fine_mapping", pheno = pheno_name, kind = "summary"),
    error = function(e) NULL
  )
  if (!is.null(stored) && nrow(stored) > 0L && !run_if_missing) {
    if (!is.null(peak_ids)) {
      stored <- stored[stored$peak_ID %in% peak_ids, , drop = FALSE]
    }
    return(stored)
  }

  needed <- peak_ids
  if (!is.null(stored) && nrow(stored) > 0L) {
    needed <- setdiff(peak_ids, stored$peak_ID)
    if (length(needed) == 0L) {
      return(stored[stored$peak_ID %in% peak_ids, , drop = FALSE])
    }
  }

  if (length(needed) == 0L) {
    return(data.frame())
  }

  fm <- runFineMapping(
    object = object,
    pheno = pheno_name,
    peak_ids = needed,
    recalc = recalc,
    store = store
  )
  out <- fm$summary
  if (!is.null(stored) && nrow(stored) > 0L) {
    out <- rbind(stored, out)
    out <- out[!duplicated(out[, c("pheno", "peak_ID")]), , drop = FALSE]
    rownames(out) <- NULL
  }
  if (!is.null(peak_ids)) {
    out <- out[out$peak_ID %in% peak_ids, , drop = FALSE]
  }
  out
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

.peak_ids_for_pheno <- function(object, pheno_name, recalc) {
  use_recalc <- recalc && .store_section_exists(object, "recalc")
  if (!use_recalc && !.store_section_exists(object, "peakcall")) {
    return(integer())
  }
  pk <- .get_peakcall(object = object, pheno_name = pheno_name, recalc = use_recalc)
  if (is.null(pk) || nrow(pk) == 0L) {
    return(integer())
  }
  sort(unique(pk$peak_ID))
}

.conditional_independent_leads <- function(cond, p_threshold) {
  steps <- sort(unique(cond$step))
  leads <- character()
  for (step in steps) {
    sub <- cond[cond$step == step, , drop = FALSE]
    sub <- sub[order(sub$P_conditional), , drop = FALSE]
    if (nrow(sub) == 0L) {
      next
    }
    if (is.finite(sub$P_conditional[1]) && sub$P_conditional[1] < p_threshold) {
      leads <- c(leads, as.character(sub$variant_ID[1]))
    }
  }
  unique(leads)
}

.conditional_marginal_lead <- function(object, pheno_name, block) {
  if (nrow(block) == 0L) {
    return(NA_character_)
  }
  scan_df <- lazyData(object = object, dataset = "scan", pheno = pheno_name)
  if (is.null(scan_df)) {
    return(NA_character_)
  }
  hit <- scan_df$variant_ID %in% block$variant_ID
  sub <- scan_df[hit, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(NA_character_)
  }
  as.character(sub$variant_ID[which.min(sub$P.model)])
}

.interpret_credible_set <- function(n_in_set, n_locus, top_pip_variant, marginal_lead) {
  resolution <- if (n_in_set <= 3L) {
    "high"
  } else if (n_in_set <= 10L) {
    "moderate"
  } else {
    "low"
  }

  frac <- if (n_locus > 0L) n_in_set / n_locus else NA_real_
  same_as_marginal <- identical(as.character(top_pip_variant), as.character(marginal_lead))

  message <- switch(
    resolution,
    high = paste0(
      "Credible set contains ", n_in_set, " of ", n_locus,
      " locus variants (well resolved within the LD block)."
    ),
    moderate = paste0(
      "Credible set contains ", n_in_set, " of ", n_locus,
      " variants; resolution is moderate — prioritize peak-neighbourhood variants."
    ),
    low = paste0(
      "Credible set contains ", n_in_set, " of ", n_locus,
      " variants (", sprintf("%.0f%%", 100 * frac),
      "); locus not finely resolved — interpret as uncertainty, not causal SNP identity."
    )
  )

  if (!is.na(marginal_lead) && !same_as_marginal) {
    message <- paste0(
      message,
      " Top-PIP variant (", top_pip_variant,
      ") differs from marginal lead (", marginal_lead, "); both are LD neighbours."
    )
  } else if (same_as_marginal) {
    message <- paste0(message, " Top PIP matches the marginal GWAS lead.")
  }

  list(
    resolution = resolution,
    credible_set_fraction = frac,
    message = message
  )
}

.fine_mapping_summary_row <- function(object,
                                      pheno_name,
                                      peak_id,
                                      cond,
                                      cred,
                                      p_threshold,
                                      r2_threshold = 0.8) {
  attr(cond, "pheno_name") <- pheno_name
  assessment <- .assess_independent_leads(
    object = object,
    cond = cond,
    p_threshold = p_threshold,
    r2_threshold = r2_threshold
  )
  cond_sum <- .conditional_summary_from_assessment(
    peak_id = peak_id,
    n_variants_in_locus = length(unique(cond$variant_ID)),
    assessment = assessment,
    p_threshold = p_threshold,
    r2_threshold = r2_threshold
  )
  cred_sum <- attr(cred, "summary")
  data.frame(
    pheno = pheno_name,
    peak_ID = peak_id,
    n_variants_in_locus = cred_sum$n_variants_in_locus,
    n_stepwise_leads = cond_sum$n_stepwise_leads,
    n_independent_signals = cond_sum$n_independent_signals,
    n_ld_artifacts = cond_sum$n_ld_artifacts,
    n_ambiguous_signals = cond_sum$n_ambiguous_signals,
    stepwise_lead_variants = cond_sum$stepwise_lead_variants,
    independent_lead_variants = cond_sum$independent_lead_variants,
    ld_artifact_variants = cond_sum$ld_artifact_variants,
    n_credible_set = cred_sum$n_variants_in_set,
    credible_set_fraction = cred_sum$credible_set_fraction,
    resolution = cred_sum$resolution,
    top_pip_variant_ID = cred_sum$lead_variant_ID,
    marginal_lead_variant_ID = cred_sum$marginal_lead_variant_ID,
    top_pip_is_marginal_lead = cred_sum$top_pip_is_marginal_lead,
    conditional_message = cond_sum$message,
    credible_set_message = cred_sum$message,
    stringsAsFactors = FALSE
  )
}

.dosage_vector <- function(object, variant_id, geno_format) {
  g <- getGenoPerMarker(
    object = object,
    geno_format = geno_format,
    marker_index = .variant_id_to_marker_index(object, variant_id)
  )
  if (geno_format == "dosage") {
    return(as.numeric(g))
  }
  if (is.matrix(g)) {
    return(colSums(g == 1, na.rm = TRUE))
  }
  as.numeric(g)
}

.variant_pair_r2 <- function(object, id_i, id_j, geno_format) {
  if (identical(as.character(id_i), as.character(id_j))) {
    return(1)
  }
  gi <- .dosage_vector(object, id_i, geno_format)
  gj <- .dosage_vector(object, id_j, geno_format)
  if (geno_format %in% c("genotype", "haplotype", "corrected")) {
    ok <- is.finite(gi) & is.finite(gj)
    if (!any(ok)) {
      return(NA_real_)
    }
    return(mean(gi[ok] == gj[ok], na.rm = TRUE))
  }
  if (stats::sd(gi, na.rm = TRUE) == 0 || stats::sd(gj, na.rm = TRUE) == 0) {
    ok <- is.finite(gi) & is.finite(gj)
    if (!any(ok)) {
      return(NA_real_)
    }
    return(as.numeric(identical(gi[ok], gj[ok])))
  }
  stats::cor(gi, gj, use = "pairwise.complete.obs")^2
}

.stats_nearly_identical <- function(beta_i, se_i, beta_j, se_j, tol = 1e-8) {
  if (!all(is.finite(c(beta_i, se_i, beta_j, se_j)))) {
    return(FALSE)
  }
  scale_b <- max(abs(c(beta_i, beta_j, 1e-12)))
  scale_s <- max(abs(c(se_i, se_j, 1e-12)))
  abs(beta_i - beta_j) <= tol * scale_b && abs(se_i - se_j) <= tol * scale_s
}

.conditional_scan_context <- function(object, pheno_name) {
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
  list(
    pheno_vec = pheno_vec,
    geno_format = geno_format,
    formula_txt = formula_txt,
    conv_fun = conv_fun,
    binary = isTRUE(binary)
  )
}

.joint_conditional_p <- function(object,
                                 lead_target,
                                 lead_vs,
                                 ctx) {
  g_target <- .dosage_vector(object, lead_target, ctx$geno_format)
  g_vs <- .dosage_vector(object, lead_vs, ctx$geno_format)
  df <- .makeDF(
    g = g_target,
    phe = ctx$pheno_vec,
    conv_fun = ctx$conv_fun,
    formula = ctx$formula_txt
  )
  target_cols <- setdiff(names(df$df), "phe")
  vs_df <- .makeDF(
    g = g_vs,
    phe = ctx$pheno_vec,
    conv_fun = ctx$conv_fun,
    formula = ctx$formula_txt
  )
  vs_cols <- setdiff(names(vs_df$df), "phe")
  vs_rename <- paste0("vs_", vs_cols)
  names(vs_df$df)[match(vs_cols, names(vs_df$df))] <- vs_rename
  joint <- cbind(df$df, vs_df$df[, vs_rename, drop = FALSE])
  fml <- stats::formula(
    paste0("phe ~ ", paste(c(target_cols, vs_rename), collapse = " + "))
  )
  family <- if (ctx$binary) "binomial" else "gaussian"
  fit <- try(stats::glm(fml, data = joint, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(NA_real_)
  }
  sm <- summary(fit)
  coef_nm <- rownames(sm$coefficients)
  hit <- intersect(target_cols, coef_nm)
  if (length(hit) == 0L) {
    return(NA_real_)
  }
  sm$coefficients[hit[1], 4]
}

.step1_stats_for_variant <- function(cond, variant_id) {
  step1 <- cond[cond$step == min(cond$step), , drop = FALSE]
  hit <- step1[as.character(step1$variant_ID) == as.character(variant_id), , drop = FALSE]
  if (nrow(hit) == 0L) {
    return(list(beta = NA_real_, se = NA_real_, P = NA_real_))
  }
  list(beta = hit$beta[1], se = hit$se[1], P = hit$P_conditional[1])
}

.lead_step_p <- function(cond, variant_id, p_threshold) {
  steps <- sort(unique(cond$step))
  for (step in steps) {
    sub <- cond[cond$step == step, , drop = FALSE]
    sub <- sub[order(sub$P_conditional), , drop = FALSE]
    if (nrow(sub) == 0L) {
      next
    }
    if (identical(as.character(sub$variant_ID[1]), as.character(variant_id)) &&
        is.finite(sub$P_conditional[1]) &&
        sub$P_conditional[1] < p_threshold) {
      return(list(step = step, P = sub$P_conditional[1]))
    }
  }
  list(step = NA_integer_, P = NA_real_)
}

.assess_independent_leads <- function(object,
                                      cond,
                                      p_threshold = 0.05,
                                      r2_threshold = 0.8) {
  stepwise <- .conditional_independent_leads(cond = cond, p_threshold = p_threshold)
  empty_pairs <- data.frame(
    lead = character(),
    vs_lead = character(),
    r2 = numeric(),
    stats_identical = logical(),
    joint_P = numeric(),
    stringsAsFactors = FALSE
  )
  empty_leads <- data.frame(
    variant_ID = character(),
    step = integer(),
    P_step = numeric(),
    verdict = character(),
    max_r2_to_supported = numeric(),
    min_joint_P = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(stepwise) == 0L) {
    return(list(pairs = empty_pairs, leads = empty_leads, stepwise = character()))
  }

  ctx <- NULL
  geno_format <- "dosage"
  if (!is.null(object) && inherits(object, "LazyGas")) {
    pheno_names <- getPheno(object)$pheno_names
    attr_pheno <- attr(cond, "pheno_name")
    if (!is.null(attr_pheno) && attr_pheno %in% pheno_names) {
      pheno_name <- attr_pheno
    } else if (length(pheno_names) == 1L) {
      pheno_name <- pheno_names[[1]]
    } else {
      pheno_name <- pheno_names[[1]]
      for (pn in pheno_names) {
        scan_df <- tryCatch(
          lazyData(object = object, dataset = "scan", pheno = pn),
          error = function(e) NULL
        )
        if (!is.null(scan_df) && any(scan_df$variant_ID %in% cond$variant_ID)) {
          pheno_name <- pn
          break
        }
      }
    }
    ctx <- .conditional_scan_context(object, pheno_name)
    geno_format <- ctx$geno_format
  }

  pair_rows <- list()
  lead_rows <- list()
  supported <- character()

  for (i in seq_along(stepwise)) {
    lead <- as.character(stepwise[[i]])
    step_info <- .lead_step_p(cond, lead, p_threshold)
    if (i == 1L) {
      lead_rows[[i]] <- data.frame(
        variant_ID = lead,
        step = step_info$step,
        P_step = step_info$P,
        verdict = "supported",
        max_r2_to_supported = NA_real_,
        min_joint_P = NA_real_,
        stringsAsFactors = FALSE
      )
      supported <- lead
      next
    }

    max_r2 <- NA_real_
    min_joint <- Inf
    any_identical <- FALSE
    is_artifact <- FALSE

    for (vs in supported) {
      r2 <- NA_real_
      joint_p <- NA_real_
      st_lead <- .step1_stats_for_variant(cond, lead)
      st_vs <- .step1_stats_for_variant(cond, vs)
      identical_stats <- .stats_nearly_identical(
        st_lead$beta, st_lead$se, st_vs$beta, st_vs$se
      )
      if (!is.null(object) && inherits(object, "LazyGas")) {
        r2 <- tryCatch(
          .variant_pair_r2(object, lead, vs, geno_format),
          error = function(e) NA_real_
        )
        joint_p <- tryCatch(
          .joint_conditional_p(object, lead, vs, ctx),
          error = function(e) NA_real_
        )
      }
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        lead = lead,
        vs_lead = vs,
        r2 = r2,
        stats_identical = identical_stats,
        joint_P = joint_p,
        stringsAsFactors = FALSE
      )
      if (is.finite(r2)) {
        max_r2 <- if (!is.finite(max_r2)) r2 else max(max_r2, r2)
      }
      if (is.finite(joint_p)) {
        min_joint <- min(min_joint, joint_p)
      }
      any_identical <- any_identical || identical_stats
      if (identical_stats || (is.finite(r2) && r2 >= r2_threshold)) {
        is_artifact <- TRUE
      }
    }

    if (!is.finite(min_joint)) {
      min_joint <- NA_real_
    }

    verdict <- if (is_artifact || any_identical) {
      "ld_artifact"
    } else if (is.finite(min_joint) && min_joint < p_threshold) {
      "supported"
    } else {
      "ambiguous"
    }

    lead_rows[[i]] <- data.frame(
      variant_ID = lead,
      step = step_info$step,
      P_step = step_info$P,
      verdict = verdict,
      max_r2_to_supported = max_r2,
      min_joint_P = min_joint,
      stringsAsFactors = FALSE
    )
    if (identical(verdict, "supported")) {
      supported <- c(supported, lead)
    }
  }

  pairs <- if (length(pair_rows) == 0L) empty_pairs else do.call(rbind, pair_rows)
  leads_df <- do.call(rbind, lead_rows)
  rownames(pairs) <- NULL
  rownames(leads_df) <- NULL
  list(pairs = pairs, leads = leads_df, stepwise = stepwise)
}

.assess_message <- function(assessment) {
  leads <- assessment$leads
  if (is.null(leads) || nrow(leads) == 0L) {
    return("No stepwise leads below threshold.")
  }
  n_step <- nrow(leads)
  n_sup <- sum(leads$verdict == "supported")
  n_art <- sum(leads$verdict == "ld_artifact")
  n_amb <- sum(leads$verdict == "ambiguous")
  art_ids <- leads$variant_ID[leads$verdict == "ld_artifact"]
  amb_ids <- leads$variant_ID[leads$verdict == "ambiguous"]

  if (n_step == 1L && n_sup == 1L) {
    return("Single independent association (or none below threshold).")
  }

  msg <- paste0(
    n_step, " stepwise lead", if (n_step == 1L) "" else "s", "; ",
    n_sup, " supported independent signal", if (n_sup == 1L) "" else "s",
    " after LD checks"
  )
  if (n_art > 0L) {
    msg <- paste0(
      msg,
      " (", n_art, " LD artifact", if (n_art == 1L) "" else "s",
      ": ", paste(art_ids, collapse = ","), ")"
    )
  }
  msg <- paste0(msg, ".")
  if (n_amb > 0L) {
    msg <- paste0(
      msg,
      " Ambiguous lead", if (n_amb == 1L) "" else "s",
      " (", paste(amb_ids, collapse = ","),
      "): inspect before claiming multi-signal QTL."
    )
  } else if (n_sup > 1L) {
    msg <- paste0(msg, " Inspect each supported lead before fine-mapping.")
  }
  msg
}

.conditional_summary_from_assessment <- function(peak_id,
                                                 n_variants_in_locus,
                                                 assessment,
                                                 p_threshold,
                                                 r2_threshold) {
  leads <- assessment$leads
  stepwise <- assessment$stepwise
  if (is.null(stepwise)) {
    stepwise <- if (is.null(leads) || nrow(leads) == 0L) {
      character()
    } else {
      as.character(leads$variant_ID)
    }
  }
  supported <- if (is.null(leads) || nrow(leads) == 0L) {
    character()
  } else {
    as.character(leads$variant_ID[leads$verdict == "supported"])
  }
  artifacts <- if (is.null(leads) || nrow(leads) == 0L) {
    character()
  } else {
    as.character(leads$variant_ID[leads$verdict == "ld_artifact"])
  }
  ambiguous <- if (is.null(leads) || nrow(leads) == 0L) {
    character()
  } else {
    as.character(leads$variant_ID[leads$verdict == "ambiguous"])
  }

  data.frame(
    peak_ID = peak_id,
    n_variants_in_locus = n_variants_in_locus,
    n_stepwise_leads = length(stepwise),
    n_independent_signals = length(supported),
    n_ld_artifacts = length(artifacts),
    n_ambiguous_signals = length(ambiguous),
    stepwise_lead_variants = paste(stepwise, collapse = ","),
    independent_lead_variants = paste(supported, collapse = ","),
    ld_artifact_variants = paste(artifacts, collapse = ","),
    p_threshold = p_threshold,
    r2_threshold = r2_threshold,
    message = .assess_message(assessment),
    stringsAsFactors = FALSE
  )
}

.empty_independent_assessment <- function(peak_id, p_threshold, r2_threshold) {
  assessment <- list(
    pairs = data.frame(
      lead = character(),
      vs_lead = character(),
      r2 = numeric(),
      stats_identical = logical(),
      joint_P = numeric(),
      stringsAsFactors = FALSE
    ),
    leads = data.frame(
      variant_ID = character(),
      step = integer(),
      P_step = numeric(),
      verdict = character(),
      max_r2_to_supported = numeric(),
      min_joint_P = numeric(),
      stringsAsFactors = FALSE
    ),
    stepwise = character()
  )
  list(
    pairs = assessment$pairs,
    leads = assessment$leads,
    summary = .conditional_summary_from_assessment(
      peak_id = peak_id,
      n_variants_in_locus = 0L,
      assessment = assessment,
      p_threshold = p_threshold,
      r2_threshold = r2_threshold
    )
  )
}

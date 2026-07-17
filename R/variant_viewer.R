################################################################################
#' Gather variant and gene-structure data for the variant viewer
#'
#' @param object A \code{LazyGas} object.
#' @param gene_id Gene identifier (must match \code{GFF} \code{ID} / \code{gene_id}).
#' @param gff A \code{GRanges} object (e.g. from \code{rtracklayer::import.gff()}).
#' @param pheno Phenotype name or index for [lazyData()].
#' @param snpeff Optional \code{snpeff_gds} object. Used when stored snpeff results
#'   are unavailable.
#' @param peak_id If not \code{NULL}, restrict variants to the given peak block.
#' @param recalc If \code{TRUE}, use recalculated peak blocks.
#' @param ann Optional annotation \code{data.frame} (as in [listCandidate()]) for
#'   header text.
#'
#' @return A list with \code{gene_info}, \code{plot_df}, \code{cds_df},
#'   \code{geno_df}, and \code{ggplot} (a \code{ggplot} object).
#' @export
getVariantViewerData <- function(object,
                                   gene_id,
                                   gff,
                                   pheno,
                                   snpeff = NULL,
                                   peak_id = NULL,
                                   recalc = FALSE,
                                   ann = NULL) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  if (!inherits(gff, "GRanges")) {
    stop("'gff' must be a GRanges object.", call. = FALSE)
  }
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)

  gene_info <- .variant_viewer_gene_info(gff = gff, gene_id = gene_id, ann = ann)
  ann_df <- .variant_viewer_snpeff_df(
    object = object,
    gene_id = gene_id,
    pheno_name = pheno_name,
    snpeff = snpeff,
    peak_id = peak_id,
    recalc = recalc
  )
  if (is.null(ann_df) || nrow(ann_df) == 0L) {
    stop(
      "No SnpEff annotations found for gene '", gene_id, "'. ",
      "Run listCandidate() with snpeff = ..., or provide snpeff.",
      call. = FALSE
    )
  }

  scan_df <- lazyData(object = object, dataset = "scan", pheno = pheno_name)
  if (is.null(scan_df)) {
    stop("No scan data found. Run scanAssoc() first.", call. = FALSE)
  }

  plot_df <- .variant_viewer_build_plot_df(
    ann_df = ann_df,
    scan_df = scan_df,
    gene_id = gene_id
  )
  cds_df <- .variant_viewer_cds_df(gff = gff, gene_id = gene_id, plot_df = plot_df)
  geno_df <- .variant_viewer_geno_table(plot_df = plot_df)

  list(
    gene_id = gene_id,
    pheno = pheno_name,
    gene_info = gene_info,
    plot_df = plot_df,
    cds_df = cds_df,
    geno_df = geno_df,
    ggplot = .variant_viewer_ggplot(
      plot_df = plot_df,
      cds_df = cds_df,
      gene_id = gene_id,
      pheno = pheno_name
    )
  )
}

#' Plot gene structure and variant effects
#'
#' @param object A \code{LazyGas} object, or a list returned by
#'   [getVariantViewerData()].
#' @param gene_id Gene identifier (ignored if \code{object} is a viewer data list).
#' @param gff,gff_fn,pheno,snpeff,peak_id,recalc,ann Passed to [getVariantViewerData()]
#'   when \code{object} is a \code{LazyGas} object.
#' @param scale Vertical spacing factor for transcript tracks (default \code{5}).
#'
#' @return A \code{ggplot} object.
#' @export
plotVariantViewer <- function(object,
                              gene_id = NULL,
                              gff = NULL,
                              pheno = NULL,
                              snpeff = NULL,
                              peak_id = NULL,
                              recalc = FALSE,
                              ann = NULL,
                              scale = 5) {
  if (is.list(object) && !is.null(object$ggplot)) {
    return(object$ggplot)
  }
  data <- getVariantViewerData(
    object = object,
    gene_id = gene_id,
    gff = gff,
    pheno = pheno,
    snpeff = snpeff,
    peak_id = peak_id,
    recalc = recalc,
    ann = ann
  )
  data$ggplot
}

#' Build an interactive dashboard with candidate-gene links to variant viewers
#'
#' Extends [makeInteractiveSummary()] by embedding a variant viewer panel per
#' candidate gene. Clicking \code{Gene_ID} in the candidate table scrolls to and
#' reveals the corresponding gene view.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param gff A \code{GRanges} object for gene models.
#' @param out_fn Output HTML file path.
#' @param snpeff Optional \code{snpeff_gds} object for [getVariantViewerData()].
#' @param what Sections to include (see [makeInteractiveSummary()]).
#' @param genes Character vector of gene IDs to show in the viewer panel. Default
#'   uses unique \code{Gene_ID} values from the candidate list.
#' @param ann Optional gene annotation \code{data.frame}.
#' @param recalc Use recalculated peaks when gathering variants.
#' @param peak_id Restrict variants to a single peak (applied to all genes).
#'
#' @export
makeInteractiveDashboard <- function(object,
                                     pheno,
                                     gff,
                                     out_fn,
                                     snpeff = NULL,
                                     what = c("scan_png", "peakcall", "recalc", "groups",
                                              "recalc_haplo", "candidate"),
                                     genes = NULL,
                                     ann = NULL,
                                     recalc = FALSE,
                                     peak_id = NULL) {
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)

  if (is.null(genes)) {
    candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
    if (is.null(candidate) || nrow(candidate) == 0L) {
      stop("No candidate data. Run listCandidate() first, or set genes = ....",
           call. = FALSE)
    }
    if (!"Gene_ID" %in% names(candidate)) {
      stop("Candidate data has no Gene_ID column.", call. = FALSE)
    }
    genes <- unique(as.character(candidate$Gene_ID))
    genes <- genes[!is.na(genes) & nzchar(genes)]
  }

  tag_list <- .interactive_summary_tags(
    object = object,
    pheno = pheno_name,
    what = what,
    candidate_gene_links = TRUE,
    genes_for_links = genes
  )

  viewer_panels <- tagList(
    div(h1("Variant viewer"), style = "text-align:center; margin-top:2em;"),
    div(
      p("Click a Gene_ID in the candidate table above to open the gene view."),
      style = "text-align:center; color:#555; margin-bottom:1em;"
    )
  )

  for (gene in genes) {
    safe_id <- .variant_viewer_safe_id(gene)
    panel <- tryCatch(
      {
        vdata <- getVariantViewerData(
          object = object,
          gene_id = gene,
          gff = gff,
          pheno = pheno_name,
          snpeff = snpeff,
          peak_id = peak_id,
          recalc = recalc,
          ann = ann
        )
        .variant_viewer_html_tags(data = vdata, visible = FALSE, safe_id = safe_id)
      },
      error = function(e) {
        tagList(
          div(
            id = paste0("lazygas-gene-", safe_id),
            class = "lazygas-gene-panel",
            style = "display:none; margin:2em auto; width:95vw;",
            h2(gene),
            p(conditionMessage(e), style = "color:#a00;")
          )
        )
      }
    )
    viewer_panels <- tagList(viewer_panels, panel)
  }

  tag_list <- tagList(
    .variant_viewer_dashboard_js(),
    tag_list,
    viewer_panels
  )

  save_html(tag_list, file = out_fn)
  invisible(out_fn)
}

# ---- internal helpers ---------------------------------------------------------

.variant_viewer_safe_id <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

.variant_viewer_gene_info <- function(gff, gene_id, ann = NULL) {
  gene_rows <- gff
  if ("type" %in% names(S4Vectors::mcols(gff))) {
    type_hit <- gff$type %in% c("gene", "Gene")
    type_hit[is.na(type_hit)] <- FALSE
    gene_rows <- gff[type_hit]
  }
  id_match <- rep(FALSE, length(gene_rows))
  if ("ID" %in% names(S4Vectors::mcols(gene_rows))) {
    id_match <- id_match | as.character(gene_rows$ID) == gene_id
  }
  if ("gene_id" %in% names(S4Vectors::mcols(gene_rows))) {
    gene_id_hit <- as.character(gene_rows$gene_id) == gene_id
    gene_id_hit[is.na(gene_id_hit)] <- FALSE
    id_match <- id_match | gene_id_hit
  }
  id_match[is.na(id_match)] <- FALSE
  gene_rows <- gene_rows[id_match]
  if (length(gene_rows) == 0L) {
    stop("Gene '", gene_id, "' not found in gff.", call. = FALSE)
  }
  gi <- gene_rows[1L]
  info <- list(
    Gene_ID = gene_id,
    gene_chr = as.character(GenomeInfoDb::seqnames(gi)),
    gene_start = BiocGenerics::start(gi),
    gene_end = BiocGenerics::end(gi),
    strand = as.character(GenomicRanges::strand(gi))
  )
  if (!is.null(ann) && "Gene_ID" %in% names(ann)) {
    hit <- ann[ann$Gene_ID == gene_id, , drop = FALSE]
    if (nrow(hit) > 0L) {
      extra <- as.list(hit[1, , drop = FALSE])
      extra$Gene_ID <- NULL
      info <- c(info, extra)
    }
  }
  info
}

.variant_viewer_snpeff_df <- function(object,
                                      gene_id,
                                      pheno_name,
                                      snpeff,
                                      peak_id,
                                      recalc) {
  out <- lazyData(object = object, dataset = "snpeff", pheno = pheno_name)
  if (!is.null(out) && nrow(out) > 0L) {
    out <- out[out$Gene_ID == gene_id, , drop = FALSE]
  } else {
    out <- NULL
  }

  peakcall <- NULL
  if ((is.null(out) || nrow(out) == 0L) && !is.null(snpeff)) {
    peakcall <- lazyData(
      object = object,
      dataset = if (recalc) "recalc" else "peakcall",
      pheno = pheno_name
    )
    if (!is.null(peakcall) && nrow(peakcall) > 0L) {
      if (!is.null(peak_id)) {
        peakcall <- peakcall[peakcall$peak_ID == peak_id, , drop = FALSE]
      }
      peakcall <- peakcall[peakcall$variant_ID %in% peakcall$peak_variant_ID, , drop = FALSE]
      blocks <- unique(peakcall[, c("Chr", "Pos", "negLog10P")])
      if (nrow(blocks) > 0L) {
        pseudo_block <- data.frame(
          Chr = blocks$Chr,
          Pos = blocks$Pos,
          negLog10P = blocks$negLog10P,
          peak_variant_ID = blocks$Pos
        )
        out <- .addSnpEff(snpeff = snpeff, peakblock = pseudo_block)
        if (!is.null(out) && nrow(out) > 0L) {
          out <- out[out$Gene_ID == gene_id, , drop = FALSE]
        }
      }
    }
  }

  if (!is.null(out) && nrow(out) > 0L && !is.null(peak_id)) {
    if (is.null(peakcall)) {
      peakcall <- lazyData(
        object = object,
        dataset = if (recalc) "recalc" else "peakcall",
        pheno = pheno_name
      )
    }
    if (!is.null(peakcall)) {
      peakcall <- peakcall[peakcall$peak_ID == peak_id, , drop = FALSE]
      if (nrow(peakcall) > 0L) {
        vars <- unique(peakcall$variant_ID)
        out <- out[out$Pos %in% peakcall$Pos[peakcall$variant_ID %in% vars], , drop = FALSE]
      }
    }
  }

  out
}

.variant_viewer_build_plot_df <- function(ann_df, scan_df, gene_id) {
  mar <- scan_df[, c("variant_ID", "Chr", "Pos")]
  if ("P.model" %in% names(scan_df)) {
    mar$P_value <- scan_df$P.model
  } else if ("P_value" %in% names(scan_df)) {
    mar$P_value <- scan_df$P_value
  } else {
    mar$P_value <- NA_real_
  }
  if ("negLog10P" %in% names(scan_df)) {
    mar$negLog10p <- scan_df$negLog10P
  } else {
    mar$negLog10p <- -log10(mar$P_value)
  }

  ann_df$variant_key <- paste(ann_df$Chr, ann_df$Pos, sep = ":")
  mar$variant_key <- paste(mar$Chr, mar$Pos, sep = ":")
  merged <- merge(ann_df, mar, by = "variant_key", all.x = TRUE, suffixes = c("", "_scan"))

  if ("negLog10P" %in% names(merged)) {
    from_ann <- as.numeric(merged$negLog10P)
    from_scan <- as.numeric(merged$negLog10p)
    merged$negLog10p <- ifelse(is.na(from_scan), from_ann, from_scan)
  } else {
    merged$negLog10p <- as.numeric(merged$negLog10p)
  }
  if ("P_value" %in% names(merged)) {
    merged$P_value <- as.numeric(merged$P_value)
  }

  ref_alt <- strsplit(as.character(merged$Allele), ",")
  var_len <- vapply(ref_alt, function(x) {
    if (length(x) < 2L) {
      return(1L)
    }
    max(nchar(x[[1]]), nchar(x[[2]]))
  }, integer(1L))

  plot_df <- data.frame(
    Chr = merged$Chr,
    Start = merged$Pos,
    End = merged$Pos + var_len - 1L,
    Allele = merged$Allele,
    ID = merged$variant_ID,
    P_value = merged$P_value,
    negLog10p = merged$negLog10p,
    Gene_ID = gene_id,
    Transcript_ID = merged$Feature_ID,
    Effect = merged$Annotation_Impact,
    Position_in_CDS = merged$Pos.in.CDS,
    Change_in_DNA = merged$HGVS.c,
    Position_in_AA = merged$Pos.in.AA,
    Change_in_AA = merged$HGVS.p,
    stringsAsFactors = FALSE
  )

  plot_df <- plot_df[!is.na(plot_df$Transcript_ID) & nzchar(plot_df$Transcript_ID), , drop = FALSE]
  plot_df <- .variant_viewer_filter_intergenic(plot_df = plot_df, gene_id = gene_id)
  if (nrow(plot_df) == 0L) {
    stop("No in-gene variants remaining after filtering.", call. = FALSE)
  }
  plot_df
}

.variant_viewer_filter_intergenic <- function(plot_df, gene_id) {
  inter_gene <- grepl(gene_id, plot_df$Transcript_ID, fixed = TRUE)
  if (sum(inter_gene) > 0L && sum(!inter_gene) > 0L) {
    check <- duplicated(rbind(
      plot_df[!inter_gene, c("ID", "Effect"), drop = FALSE],
      plot_df[inter_gene, c("ID", "Effect"), drop = FALSE]
    ))
    check <- check[seq(sum(!inter_gene) + 1L, length.out = sum(inter_gene))]
    inter_gene[inter_gene] <- inter_gene[inter_gene] & check
    plot_df <- plot_df[!inter_gene, , drop = FALSE]
  }
  plot_df
}

.variant_viewer_cds_df <- function(gff, gene_id, plot_df) {
  type_hit <- gff$type == "CDS"
  type_hit[is.na(type_hit)] <- FALSE
  cds <- gff[type_hit]
  if ("gene_id" %in% names(S4Vectors::mcols(cds))) {
    gene_id_hit <- as.character(cds$gene_id) == gene_id
    gene_id_hit[is.na(gene_id_hit)] <- FALSE
    cds <- cds[gene_id_hit]
  } else {
    gene_tx <- unique(plot_df$Transcript_ID)
    if ("Parent" %in% names(S4Vectors::mcols(cds))) {
      parent_hit <- as.character(cds$Parent) %in% gene_tx
      parent_hit[is.na(parent_hit)] <- FALSE
      cds <- cds[parent_hit]
    }
  }
  if (length(cds) == 0L) {
    return(data.frame(
      Transcript_ID = character(),
      cds_start = integer(),
      cds_end = integer(),
      ymin = numeric(),
      ymax = numeric()
    ))
  }
  parent_col <- if ("Parent" %in% names(S4Vectors::mcols(cds))) cds$Parent else cds$ID
  cds_df <- data.frame(
    Transcript_ID = unlist(parent_col),
    cds_start = BiocGenerics::start(cds),
    cds_end = BiocGenerics::end(cds),
    stringsAsFactors = FALSE
  )
  cds_df
}

.variant_viewer_layout_y <- function(plot_df, cds_df, scale = 5) {
  tx_levels <- sort(unique(plot_df$Transcript_ID))
  plot_df$Transcript_ID <- factor(plot_df$Transcript_ID, levels = tx_levels)
  plot_df$ymin <- as.numeric(plot_df$Transcript_ID)
  plot_df$ymin <- max(plot_df$ymin) - plot_df$ymin + scale * 1.2
  plot_df$ymax <- plot_df$ymin + 0.8
  plot_df$y_pos <- plot_df$ymin + (plot_df$ymax - plot_df$ymin) / 2
  plot_df$x_pos <- plot_df$Start + (plot_df$End - plot_df$Start) / 2
  plot_df$Effect <- factor(
    plot_df$Effect,
    levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")
  )

  if (nrow(cds_df) > 0L) {
    cds_df$Transcript_ID <- factor(cds_df$Transcript_ID, levels = tx_levels)
    cds_df$ymin <- as.numeric(cds_df$Transcript_ID)
    cds_df$ymin <- max(cds_df$ymin) - cds_df$ymin + scale * 1.2
    cds_df$ymax <- cds_df$ymin + 0.8
  }

  y_lab <- rev(sort(unique(plot_df$y_pos)))
  list(plot_df = plot_df, cds_df = cds_df, y_lab = y_lab, tx_levels = tx_levels)
}

.variant_viewer_is_coding <- function(position_in_cds, position_in_aa, change_in_aa) {
  cds_chr <- trimws(as.character(position_in_cds))
  aa_pos_chr <- trimws(as.character(position_in_aa))
  aa_chg_chr <- trimws(as.character(change_in_aa))
  cds_ok <- !is.na(position_in_cds) && nzchar(cds_chr) && cds_chr != "."
  aa_pos_ok <- !is.na(position_in_aa) && nzchar(aa_pos_chr) && aa_pos_chr != "."
  aa_chg_ok <- !is.na(change_in_aa) && nzchar(aa_chg_chr) && aa_chg_chr != "."
  cds_ok || aa_pos_ok || aa_chg_ok
}

.variant_viewer_format_tooltip_field <- function(label, value) {
  val <- trimws(as.character(value))
  if (is.na(val) || !nzchar(val) || val == ".") {
    val <- "."
  }
  paste0(label, ": ", val)
}

.variant_viewer_variant_tooltip <- function(df) {
  n <- nrow(df)
  if (n == 0L) {
    return(character())
  }
  vapply(seq_len(n), function(i) {
    row <- df[i, , drop = FALSE]
    allele <- as.character(row$Allele)
    if (is.na(allele) || !nzchar(trimws(allele))) {
      allele <- "."
    }
    if (.variant_viewer_is_coding(
      position_in_cds = row$Position_in_CDS,
      position_in_aa = row$Position_in_AA,
      change_in_aa = row$Change_in_AA
    )) {
      paste(
        .variant_viewer_format_tooltip_field("Allele", allele),
        .variant_viewer_format_tooltip_field("CDS pos", row$Position_in_CDS),
        .variant_viewer_format_tooltip_field("AA pos", row$Position_in_AA),
        .variant_viewer_format_tooltip_field("AA change", row$Change_in_AA),
        sep = "<br>"
      )
    } else {
      .variant_viewer_format_tooltip_field("Allele", allele)
    }
  }, character(1L))
}

.variant_viewer_gwas_tooltip <- function(df) {
  n <- nrow(df)
  if (n == 0L) {
    return(character())
  }
  vapply(seq_len(n), function(i) {
    row <- df[i, , drop = FALSE]
    pos <- if ("Start" %in% names(row)) row$Start else row$Pos
    neg <- as.numeric(row$negLog10p)
    if (!is.finite(neg)) {
      neg <- -log10(as.numeric(row$P_value))
    }
    paste(
      .variant_viewer_format_tooltip_field("Pos", pos),
      .variant_viewer_format_tooltip_field(
        "negLog10P",
        if (is.finite(neg)) signif(neg, digits = 6) else "."
      ),
      sep = "<br>"
    )
  }, character(1L))
}

.variant_viewer_ggplot <- function(plot_df,
                                   cds_df,
                                   scale = 5,
                                   gene_id = NULL,
                                   pheno = NULL) {
  laid <- .variant_viewer_layout_y(plot_df = plot_df, cds_df = cds_df, scale = scale)
  plot_df <- laid$plot_df
  cds_df <- laid$cds_df
  y_lab <- laid$y_lab

  plot_man <- unique(subset(
    plot_df,
    select = c(
      ID, Chr, Start, End, Allele, P_value, negLog10p,
      y_pos, x_pos, Position_in_CDS, Change_in_DNA, Position_in_AA, Change_in_AA
    )
  ))
  max_neg <- max(plot_man$negLog10p, na.rm = TRUE)
  if (!is.finite(max_neg) || max_neg <= 0) {
    max_neg <- 1
  }
  plot_man$scaled_score <- plot_man$negLog10p / max_neg * scale
  plot_df$tooltip <- .variant_viewer_variant_tooltip(plot_df)
  plot_man$tooltip <- .variant_viewer_gwas_tooltip(plot_man)

  p <- ggplot2::ggplot()
  if (nrow(cds_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data = cds_df,
      ggplot2::aes(
        ymin = .data$ymin, ymax = .data$ymax,
        xmin = .data$cds_start, xmax = .data$cds_end
      ),
      fill = "gray99"
    )
  }
  p <- p +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(
        y = .data$y_pos,
        x = .data$x_pos,
        color = .data$Effect,
        text = .data$tooltip
      ),
      size = 0.8
    ) +
    ggplot2::geom_point(
      data = plot_man,
      ggplot2::aes(
        y = .data$scaled_score,
        x = .data$x_pos,
        text = .data$tooltip
      ),
      size = 0.8
    ) +
    ggplot2::scale_color_manual(
      breaks = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
      values = c("magenta", "green", "blue", "gray30"),
      labels = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
      na.value = "gray60"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_lab,
      labels = laid$tx_levels
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("Physical position (bp)")

  title <- .variant_viewer_plot_title(gene_id = gene_id, pheno = pheno)
  if (!is.null(title)) {
    p <- p +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  p
}

.variant_viewer_plot_title <- function(gene_id = NULL, pheno = NULL) {
  gene_id <- if (is.null(gene_id)) NULL else as.character(gene_id)[1L]
  pheno <- if (is.null(pheno)) NULL else as.character(pheno)[1L]
  if (!is.null(gene_id) && !is.na(gene_id) && nzchar(gene_id) &&
      !is.null(pheno) && !is.na(pheno) && nzchar(pheno)) {
    return(paste0(gene_id, " — ", pheno))
  }
  if (!is.null(pheno) && !is.na(pheno) && nzchar(pheno)) {
    return(pheno)
  }
  if (!is.null(gene_id) && !is.na(gene_id) && nzchar(gene_id)) {
    return(gene_id)
  }
  NULL
}

.variant_viewer_geno_table <- function(plot_df) {
  geno_df <- unique(subset(
    plot_df,
    select = c(Chr, Start, Allele, ID, P_value, negLog10p, Effect,
               Position_in_CDS, Change_in_DNA, Position_in_AA, Change_in_AA)
  ))
  geno_df$P_value <- signif(as.numeric(geno_df$P_value), digits = 3)
  geno_df$negLog10p <- signif(as.numeric(geno_df$negLog10p), digits = 3)
  names(geno_df)[names(geno_df) == "negLog10p"] <- "negLog10P"
  geno_df
}

.variant_viewer_plot_height <- function(n_transcripts) {
  h <- n_transcripts * 50 + 100
  max(h, 300)
}

.variant_viewer_html_tags <- function(data,
                                      visible = TRUE,
                                      safe_id = NULL,
                                      title = NULL) {
  if (is.null(safe_id)) {
    safe_id <- .variant_viewer_safe_id(data$gene_id)
  }
  gi <- data$gene_info
  n_tx <- length(unique(data$plot_df$Transcript_ID))
  plot_h <- .variant_viewer_plot_height(n_tx)

  header_tags <- list(
    h2(data$gene_id),
    p(paste0(
      gi$gene_chr, ": ", gi$gene_start, "..", gi$gene_end,
      " (", gi$strand, ")"
    ))
  )
  extra_cols <- setdiff(
    names(gi),
    c("Gene_ID", "gene_chr", "gene_start", "gene_end", "strand")
  )
  for (col in extra_cols) {
    val <- gi[[col]]
    if (!is.null(val) && !is.na(val) && nzchar(as.character(val))) {
      header_tags <- c(header_tags, list(p(paste0(col, ": ", val))))
    }
  }

  p <- data$ggplot
  plot_widget <- plotly::ggplotly(
    p,
    height = plot_h,
    width = NULL,
    tooltip = "text"
  ) |>
    plotly::layout(autosize = TRUE) |>
    plotly::config(responsive = TRUE)

  tbl <- reactable::reactable(
    data = data$geno_df,
    sortable = TRUE,
    resizable = TRUE,
    filterable = TRUE,
    searchable = TRUE,
    showPageSizeOptions = TRUE,
    wrap = TRUE,
    striped = TRUE
  )

  display <- if (visible) "block" else "none"

  tagList(
    div(
      id = paste0("lazygas-gene-", safe_id),
      class = "lazygas-gene-panel",
      style = paste0("display:", display, "; margin:2em auto; width:95vw;"),
      div(header_tags, style = "text-align:center;"),
      div(
        plot_widget,
        class = "lazygas-gene-plot",
        style = "margin:auto;width:100%;max-width:95vw;"
      ),
      div(h3("Variant list"), style = "text-align:center; margin-top:1.5em;"),
      div(tbl, style = "margin:auto;width:95vw;")
    )
  )
}

.variant_viewer_dashboard_js <- function() {
  htmltools::HTML(
    "<style>\n",
    ".lazygas-gene-panel .plotly.html-widget,\n",
    ".lazygas-gene-panel .js-plotly-plot {\n",
    "  width: 100% !important;\n",
    "  max-width: 100%;\n",
    "}\n",
    "</style>\n",
    "<script>\n",
    "function lazyGasShowGene(geneId) {\n",
    "  var panels = document.querySelectorAll('.lazygas-gene-panel');\n",
    "  for (var i = 0; i < panels.length; i++) panels[i].style.display = 'none';\n",
    "  var el = document.getElementById('lazygas-gene-' + geneId);\n",
    "  if (el) {\n",
    "    el.style.display = 'block';\n",
    "    el.scrollIntoView({behavior: 'smooth', block: 'start'});\n",
    "    window.requestAnimationFrame(function() {\n",
    "      if (!window.Plotly) return;\n",
    "      var plots = el.querySelectorAll('.js-plotly-plot');\n",
    "      for (var j = 0; j < plots.length; j++) {\n",
    "        Plotly.Plots.resize(plots[j]);\n",
    "      }\n",
    "    });\n",
    "  }\n",
    "}\n",
    "</script>"
  )
}

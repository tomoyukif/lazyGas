#' @param start The start position on the chromosome. Can be NULL.
#' @param end The end position on the chromosome. Can be NULL.
#' @param recalc Logical. If TRUE, use recalculated scan data.
#' @param ... Additional arguments.
#'
#' @return Generates plots for the specified peaks.
#'
#' @export
#'
setGeneric("plotPeaks", function(object,
                                 pheno = NULL,
                                 chr = NULL,
                                 start = NULL,
                                 end = NULL,
                                 recalc = TRUE,
                                 ...)
  standardGeneric("plotPeaks"))

## Define the plotPeaks method for the LazyGas class
#'
#' @method plotPeaks LazyGas
#'
setMethod("plotPeaks",
          "LazyGas",
          function(object,
                   pheno = NULL,
                   chr = NULL,
                   start = NULL,
                   end = NULL,
                   recalc = FALSE){
            .check_peakcall_data(object = object)
            path <- .determine_path(object = object, recalc = recalc)
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)
            chr_pos_limits <- .get_chromosome_data(object = object)

            p <- .generate_peak_plots(object = object,
                                      pheno_name = pheno_name,
                                      path = path,
                                      chr_pos_limits = chr_pos_limits,
                                      chr = chr,
                                      start = start,
                                      end = end)
            return(p)
          }
)

# Check if peakcall data exists
.check_peakcall_data <- function(object) {
  if (!.store_section_exists(object, "peakcall")) {
    stop("No peakcall data in the input LazyGas object.\n",
         "Run callPeakBlock() to call peaks.")
  }
}

# Determine the path based on recalc parameter
.determine_path <- function(object, recalc) {
  section <- if (recalc) "recalc" else "peakcall"
  if (!.store_section_exists(object, section)) {
    if (recalc) {
      stop(
        "recalc = TRUE was specified but, \n",
        "no recalculated scan data in the input LazyGas object.\n",
        "Run recalcAssoc() to recalculate associations."
      )
    }
    stop("No peakcall data in the input LazyGas object.", call. = FALSE)
  }
  paste0("lazygas/", section)
}

# Get chromosome and position data (order = GDS / scan marker order)
.get_chromosome_data <- function(object) {
  chr_data <- getChromosome(object = object)
  pos_data <- getPosition(object = object)
  chr_lev <- unique(as.character(chr_data))
  chr_fac <- factor(chr_data, levels = chr_lev)
  chr_min <- tapply(pos_data, chr_fac, min)
  chr_max <- tapply(pos_data, chr_fac, max)
  # tapply drops unused levels; restore genome order explicitly
  chr_min <- chr_min[chr_lev]
  chr_max <- chr_max[chr_lev]
  list(chr_min = chr_min, chr_max = chr_max, chr_lev = chr_lev)
}

# Generate peak plots
.generate_peak_plots <- function(object,
                                 pheno_name,
                                 path,
                                 chr_pos_limits,
                                 chr,
                                 start,
                                 end) {
  peakcall <- .get_peakcall(object = object,
                            pheno_name = pheno_name,
                            recalc = grepl(pattern = "recalc", x = path))

  if (is.null(peakcall)) {
    peakcall <- data.frame(peak_variant_ID = numeric(), FDR = numeric(),
                           variant_ID = numeric(), Chr = numeric(),
                           Pos = numeric(), negLog10P = numeric(),
                           peak_ID = numeric())
    start <- NULL
    end <- NULL
    chr <- NULL
  }

  # Generate the peak plot
  p <- .peak_draw(peakcall = peakcall,
                  chr = chr,
                  start = start,
                  end = end,
                  chr_min = chr_pos_limits$chr_min,
                  chr_max = chr_pos_limits$chr_max,
                  chr_lev = chr_pos_limits$chr_lev)
  return(p)
}

# Draw Peak Plot
#' @import ggplot2
#'
.peak_draw <- function(peakcall, chr, start, end, chr_min, chr_max, chr_lev, peak = NULL){
  # Mark peaks in peakcall
  peakcall$is_peak <- peakcall$peak_variant_ID == peakcall$variant_ID
  peakcall$Chr <- factor(peakcall$Chr, levels = chr_lev)

  # Create dummy data for chromosome boundaries
  dummy <- rbind(data.frame(Chr = names(chr_min),
                            Pos = chr_min,
                            negLog10P = 0),
                 data.frame(Chr = names(chr_max),
                            Pos = chr_max,
                            negLog10P = 0))
  dummy$Chr <- factor(dummy$Chr, levels = chr_lev)

  # Subset peakcall and dummy data based on the specified chromosome
  if(!is.null(chr)){
    peakcall <- subset(peakcall, subset = Chr == chr)
    dummy <- subset(dummy, subset = Chr == chr)
  }

  # Subset peakcall and dummy data based on the specified start position
  if(!is.null(start)){
    peakcall <- subset(peakcall, subset = Pos >= start)
    dummy <- subset(dummy, subset = Pos >= start)
  }

  # Subset peakcall and dummy data based on the specified end position
  if(!is.null(end)){
    peakcall <- subset(peakcall, subset = Pos <= end)
    dummy <- subset(dummy, subset = Pos <= end)
  }

  # Set peak ID for dummy data
  dummy$peak_ID <- if (nrow(peakcall)) peakcall$peak_ID[1] else NA
  dummy$is_peak <- FALSE

  # Subset significant peaks
  signif_peak <- subset(peakcall, subset = peak_variant_ID == variant_ID)

  # Align columns before continuous-genome coord transform
  keep_cols <- c("Chr", "Pos", "negLog10P", "peak_ID", "is_peak")
  peak_min <- as.data.frame(peakcall[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  dummy_min <- as.data.frame(dummy[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  signif_min <- if (nrow(signif_peak)) {
    as.data.frame(signif_peak[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  } else {
    peak_min[0, , drop = FALSE]
  }
  plot_df <- rbind(peak_min, dummy_min, signif_min)
  # Keep genome (scan) chromosome order even though peak rows are peak-ordered
  plot_df$Chr <- factor(as.character(plot_df$Chr), levels = as.character(chr_lev))
  coords <- .genome_plot_coords(
    plot_df,
    chr_col = "Chr",
    pos_col = "Pos",
    chr_levels = as.character(chr_lev)
  )
  plot_df <- coords$df
  n_peak <- nrow(peak_min)
  n_dummy <- nrow(dummy_min)
  peakcall <- plot_df[seq_len(n_peak), , drop = FALSE]
  dummy <- plot_df[n_peak + seq_len(n_dummy), , drop = FALSE]
  signif_peak <- plot_df[n_peak + n_dummy + seq_len(nrow(signif_min)), , drop = FALSE]

  # Hover shows chromosome physical position (bp), not concatenated genome_x
  if (nrow(peakcall)) {
    peakcall$hover <- .peak_hover_label(peakcall)
  }
  if (nrow(signif_peak)) {
    signif_peak$hover <- .peak_hover_label(signif_peak)
  }

  # Create ggplot object
  p <- ggplot() +
    geom_line(data = peakcall,
              mapping = aes(
                x = .data$genome_x,
                y = .data$negLog10P,
                color = .data$peak_ID,
                text = .data$hover,
                group = .data$peak_ID
              ),
              linewidth = 1) +
    geom_point(data = signif_peak,
               mapping = aes(
                 x = .data$genome_x,
                 y = .data$negLog10P,
                 group = .data$peak_ID,
                 text = .data$hover
               ),
               color = "magenta",
               size = 2,
               shape = 20) +
    geom_point(data = dummy,
               mapping = aes(x = .data$genome_x, y = .data$negLog10P),
               size = 0,
               color = NA) +
    ylab("-log10(P)") +
    xlab("Chromosome") +
    scale_x_continuous(breaks = unname(coords$breaks), labels = coords$labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.line.x.bottom = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  return(p)
}

#' Hover label for peak plots (physical bp, not genome_x)
#' @noRd
.peak_hover_label <- function(df) {
  pos <- as.numeric(df$Pos)
  pos_lab <- ifelse(
    is.na(pos),
    "NA",
    format(pos, scientific = FALSE, trim = TRUE, big.mark = ",")
  )
  pval <- as.numeric(df$negLog10P)
  pval_lab <- ifelse(is.na(pval), "NA", as.character(signif(pval, 4)))
  paste0(
    "Chr: ", as.character(df$Chr),
    "<br>Pos: ", pos_lab, " bp",
    "<br>-log10(P): ", pval_lab,
    "<br>peak_ID: ", as.character(df$peak_ID)
  )
}

#' ggplotly wrapper for peak plots (hide concatenated genome_x)
#' @noRd
.ggplotly_peak_plot <- function(p, ...) {
  plotly::ggplotly(p, tooltip = "text", ...)
}

################################################################################
#' Generate Haplotype Plots
#'
#' This function generates haplotype plots for the specified phenotype.
#' Plot titles include peak location, \code{-log10P}, and the additive scan
#' effect direction when \code{Coef.add} is available in stored scan results
#' (in original phenotype units for continuous traits).
#'
#' @param object A LazyGas object.
#' @param pheno The phenotype to plot.
#' @param recalc Logical. If TRUE, use recalculated scan data.
#' @param ... Additional arguments.
#'
#' @return A list of ggplot objects representing the haplotype plots.
#' @export
setGeneric("haploPlot", function(object,
                                 pheno,
                                 recalc,
                                 ...)
  standardGeneric("haploPlot"))

#' @rdname haploPlot
#' @method haploPlot LazyGas
#'
setMethod("haploPlot",
          "LazyGas",
          function(object, pheno, recalc){
            # Check if peakcall data exists
            .check_peakcall_data(object = object)

            # Determine the path based on recalc parameter
            path <- .determine_path(object = object, recalc = recalc)

            # Determine phenotype index based on input
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)

            # Get the peak call data for the specified phenotype
            peakcall <- .get_peakcall(object = object,
                                      pheno_name = pheno_name,
                                      recalc = recalc)

            # Generate the haplotype plot
            out <- .draw_haplo_plot(object = object,
                                    peakcall = peakcall,
                                    pheno_name = pheno_name)

            return(out)
          }
)

## Draw Haplotype Plot
#' @importMethodsFrom GBScleanR getMarID getPheno
#' @import ggplot2
#'
.draw_haplo_plot <- function(object, peakcall, pheno_name){
  # Check if there are any peaks
  if(is.null(peakcall)){
    message("No peak")
    return(NULL)
  }

  var_id <- getMarID(object = object)

  # Get the genotype format from the object
  geno_format <- .get_geno_format(object = object)

  peak_info <- .get_peak_info(peakcall = peakcall)

  scan_df <- tryCatch(
    .get_scan(object = object, pheno_name = pheno_name),
    error = function(e) NULL
  )

  out <- NULL
  # Loop through each unique peak ID
  for(i_peak in seq_along(peak_info$peak_height)){
    marker_idx <- match(peak_info$peak_id[i_peak], var_id)
    if (is.na(marker_idx)) {
      warning("Peak variant ", peak_info$peak_id[i_peak],
              " was not found in marker IDs.", call. = FALSE)
      next
    }

    hap <- getGenoPerMarker(
      object = object,
      geno_format = geno_format,
      marker_index = marker_idx
    )

    if (geno_format == "haplotype") {
      if (is.matrix(hap) || is.array(hap)) {
        hap <- apply(X = hap, MARGIN = 2, FUN = function(x) {
          paste(sort(x), collapse = "|")
        })
      }
    } else if (is.matrix(hap) || is.array(hap)) {
      hap <- as.vector(hap)
    }

    hap[is.na(hap)] <- "NA"
    hap[hap == ""] <- "NA"
    hap <- factor(hap, levels = sort(unique(hap)))

    # Get phenotype data
    pheno <- getPheno(object = object)
    df <- data.frame(y = pheno$pheno[, pheno_name], x = hap)

    title <- .haplo_plot_title(
      peak_chr = peak_info$peak_chr[i_peak],
      peak_pos = peak_info$peak_pos[i_peak],
      peak_height = peak_info$peak_height[i_peak],
      peak_variant_id = peak_info$peak_id[i_peak],
      scan_df = scan_df
    )

    # Generate the plot
    p <- ggplot(df) +
      geom_boxplot(aes(x = x, y = y)) +
      labs(title = title) +
      xlab("Haplotypes") +
      ylab(pheno_name)

    out <- c(out, list(p))
  }
  return(out)
}

.peak_scan_additive_coef <- function(scan_df, variant_id) {
  if (is.null(scan_df) || nrow(scan_df) == 0L) {
    return(NA_real_)
  }
  hit <- match(variant_id, scan_df$variant_ID)
  if (is.na(hit)) {
    return(NA_real_)
  }
  if ("Coef.add" %in% names(scan_df)) {
    return(as.numeric(scan_df$Coef.add[hit]))
  }
  coef_cols <- grep("^Coef\\.", names(scan_df), value = TRUE)
  if (length(coef_cols) == 0L) {
    return(NA_real_)
  }
  as.numeric(scan_df[[coef_cols[1]]][hit])
}

.haplo_effect_title_suffix <- function(coef) {
  if (is.na(coef)) {
    return("")
  }
  if (coef > 0) {
    return(sprintf(
      ", alt allele increases phenotype (Coef.add = %s)",
      signif(coef, digits = 3)
    ))
  }
  if (coef < 0) {
    return(sprintf(
      ", alt allele decreases phenotype (Coef.add = %s)",
      signif(coef, digits = 3)
    ))
  }
  ", no additive effect (Coef.add = 0)"
}

.haplo_plot_title <- function(peak_chr,
                              peak_pos,
                              peak_height,
                              peak_variant_id,
                              scan_df) {
  base <- paste(
    "Peak @",
    paste(peak_chr, peak_pos, sep = "_"),
    ", -log10P =",
    signif(peak_height, digits = 3)
  )
  coef <- .peak_scan_additive_coef(scan_df = scan_df,
                                   variant_id = peak_variant_id)
  paste0(base, .haplo_effect_title_suffix(coef))
}

.get_peak_info <- function(peakcall = peakcall){
  peak_height <- tapply(X = peakcall$peak_negLog10P,
                        INDEX = peakcall$peak_ID,
                        FUN = function(x){x[1]})
  peak_id <- tapply(X = peakcall$peak_variant_ID,
                    INDEX = peakcall$peak_ID,
                    FUN = function(x){x[1]})
  peak_chr <- tapply(X = peakcall$peak_Chr,
                     INDEX = peakcall$peak_ID,
                     FUN = function(x){x[1]})
  peak_pos <- tapply(X = peakcall$peak_Pos,
                     INDEX = peakcall$peak_ID,
                     FUN = function(x){x[1]})
  out <- data.frame(peak_height = peak_height,
                    peak_id = peak_id,
                    peak_chr = peak_chr,
                    peak_pos = peak_pos)
  return(out)
}

################################################################################
#' Retrieve LazyGas Data
#'
#' This function retrieves data from the LazyGas object based on the specified dataset and phenotype.
#'
#' @param object A LazyGas object.

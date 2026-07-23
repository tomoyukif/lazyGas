#' Call peak blocks from association scan results
#'
#' @param object A \code{LazyGas} object with scan results.
#' @param signif Significance threshold for peak seed variants (FDR or logical vector).
#' @param threshold LD threshold for extending peak blocks.
#' @param limit_peakcall Maximum number of peaks to call per phenotype.
#' @param n_threads Number of parallel threads (\code{NULL} uses a default).
#' @param ... Not used.
#'
#' @export
setGeneric("callPeakBlock", function(object,
                                     signif = 0.05,
                                     threshold = 0.8,
                                     limit_peakcall = Inf,
                                     n_threads = NULL,
                                     ...)
  standardGeneric("callPeakBlock"))

## Define the callPeakBlock method for the LazyGas class
setMethod("callPeakBlock",
          "LazyGas",
          function(object,
                   signif,
                   threshold,
                   limit_peakcall,
                   n_threads){
            # Check if scan data exists in the GDS object
            .check_scan_data_exists(object = object)

            # Create the necessary folders in the GDS object
            .create_peakcall_folders(object = object, signif = signif, threshold = threshold)

            # Perform peak calling for each phenotype
            .perform_peak_calling(object = object,
                                  signif = signif,
                                  threshold = threshold,
                                  n_threads = n_threads,
                                  limit_peakcall = limit_peakcall)
          }
)

## Function to initialize peakcall storage
.create_peakcall_folders <- function(object, signif, threshold) {
  if (.store_is_gds(object)) {
    peakcall_gdsn <- .create_gdsn(root_node = object$root,
                                  target_node = "lazygas",
                                  new_node = "peakcall",
                                  is_folder = TRUE)
    peaks_gdsn <- .create_gdsn(root_node = object$root,
                               target_node = "lazygas/peakcall",
                               new_node = "peaks",
                               is_folder = TRUE)
    blocks_gdsn <- .create_gdsn(root_node = object$root,
                                target_node = "lazygas/peakcall",
                                new_node = "blocks",
                                is_folder = TRUE)
    gdsfmt::put.attr.gdsn(node = peaks_gdsn,
                          name = "col_names",
                          val = c("peak_ID", "peak_variant_ID"))
    gdsfmt::put.attr.gdsn(node = peaks_gdsn, name = "signif", val = signif)
    gdsfmt::put.attr.gdsn(node = peaks_gdsn, name = "threshold", val = threshold)
    gdsfmt::put.attr.gdsn(
      node = blocks_gdsn,
      name = "col_names",
      val = c("peak_ID", "variant_ID", "dist2peak", "LD2peak", "chr_start", "chr_end")
    )
  } else {
  dir.create(file.path(.store_path(object), "peakcall"),
             recursive = TRUE, showWarnings = FALSE)
    .store_set_peakcall_params(object, signif, threshold)
  }
}

## Function to perform peak calling for each phenotype
.perform_peak_calling <- function(object,
                                  signif,
                                  threshold,
                                  n_threads,
                                  limit_peakcall) {
  for(i in seq_along(object@lazydata$pheno_names)){
    message("Peak calling for the following phenotype: ", object@lazydata$pheno_names[i])
    .peakcaller(object = object,
                signif = signif,
                threshold = threshold,
                pheno_name = object@lazydata$pheno_names[i],
                n_threads = n_threads,
                limit_peakcall = limit_peakcall)
  }
}



## Main peakcaller function
.peakcaller <- function(object,
                        signif,
                        threshold,
                        pheno_name,
                        n_threads,
                        limit_peakcall){
  if (!.store_is_gds(object)) {
    .store_peak_buffer_init(object, "peakcall", pheno_name)
  }

  pvalues <- .get_and_filter_pvalues(object = object,
                                     pheno_name = pheno_name,
                                     signif = signif)

  if(nrow(pvalues) == 0){
    .handle_no_significant_peaks(object = object,
                                 pheno_name = pheno_name)
    return()
  }

  variables <- .peakcall_variables(object = object)

  .call_peaks(object = object,
              pvalues = pvalues,
              variables = variables,
              pheno_name = pheno_name,
              threshold = threshold,
              n_threads = n_threads,
              limit_peakcall = limit_peakcall)
}

## Function to get and filter significant p-values
.get_and_filter_pvalues <- function(object, pheno_name, signif) {
  pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))

  signif <- pvalues$FDR <= signif
  signif[is.na(signif)] <- FALSE
  pvalues <- subset(pvalues, subset = signif)

  return(pvalues)
}

## Function to handle the case where no significant peaks are found
.handle_no_significant_peaks <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    gdsfmt::add.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/peaks"),
      name = pheno_name, storage = "uint32",
      compress = "ZIP_RA", replace = TRUE
    )
    gdsfmt::add.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = "lazygas/peakcall/blocks"),
      name = pheno_name, storage = "double",
      compress = "ZIP_RA", replace = TRUE
    )
  } else {
    .store_write_empty_peak_tables(object, "peakcall", pheno_name)
  }
}

## Function to initialize peak-calling marker table
.peakcall_variables <- function(object) {
  data.frame(
    snp_id = getMarID(object = object),
    chr = getChromosome(object = object),
    pos = getPosition(object = object),
    stringsAsFactors = FALSE
  )
}

## Function to process peaks in a loop
.call_peaks <- function(object,
                        pvalues,
                        variables,
                        pheno_name,
                        threshold,
                        n_threads,
                        limit_peakcall) {
  peak_id <- 0

  # Get the genotype format from the object
  geno_format <- .get_geno_format(object = object)
  geno_cache <- new.env(parent = emptyenv())
  selection_cache <- new.env(parent = emptyenv())
  id_is_rowindex <- .peakcall_snp_id_is_rowindex(variables$snp_id)
  rows_by_chr <- .peakcall_rows_by_chr(variables$chr)

  # Loop until there are no more p-values to process
  while (nrow(pvalues) > 0) {
    peak_id <- peak_id + 1
    if(peak_id > limit_peakcall){
      message("The number of peaks hit the limit.")
      break
    }
    message("Calling peak: ", peak_id)

    # Identify the peak information
    peak_info <- .identify_peak(
      pvalues = pvalues,
      variables = variables,
      rows_by_chr = rows_by_chr,
      id_is_rowindex = id_is_rowindex
    )

    # Set the selection criteria based on the genotype format
    sel_key <- as.character(peak_info$chr)
    if (!exists(sel_key, envir = selection_cache, inherits = FALSE)) {
      assign(
        sel_key,
        .set_selection_criteria_for_peakcall(
          object = object,
          geno_format = geno_format,
          chr = peak_info$chr
        ),
        envir = selection_cache
      )
    }
    selection <- get(sel_key, envir = selection_cache, inherits = FALSE)

    cache_key <- paste(geno_format, peak_info$chr, sep = ":")
    if (!exists(cache_key, envir = geno_cache, inherits = FALSE)) {
      message("Loading genotypes for chr ", peak_info$chr, " ...")
      t_load <- proc.time()[[3L]]
      geno_chr <- .retrieve_geno(object = object,
                                 selection = selection,
                                 geno_format = geno_format)
      if (!is.matrix(geno_chr)) {
        geno_chr <- as.matrix(geno_chr)
      }
      storage.mode(geno_chr) <- "double"
      assign(cache_key, geno_chr, envir = geno_cache)
      message(
        sprintf(
          "Loaded chr %s genotypes: %d markers in %.1fs",
          peak_info$chr, ncol(geno_chr), proc.time()[[3L]] - t_load
        )
      )
    }
    geno <- get(cache_key, envir = geno_cache, inherits = FALSE)

    # Calculate linkage disequilibrium (LD) for the identified peak
    peak_ld <- .calculate_ld(geno = geno,
                             peak_variant_idx = peak_info$peak_index_in_chr,
                             is_categorical = selection$is_categorical,
                             n_threads = n_threads)

    # Identify the peak block based on the LD values
    peak_block <- .identify_peak_block(peak_ld = peak_ld,
                                       variables = variables,
                                       peak_info = peak_info,
                                       threshold = threshold,
                                       geno_format = geno_format)

    # Save the peak data into the GDS object
    .save_peak_data(object = object,
                    pheno_name = pheno_name,
                    peak_info = peak_info,
                    peak_block = peak_block,
                    peak_id = peak_id,
                    peak_ld = peak_ld)

    # Remove markers spanning the peak "mountain" on the chromosome. variant_ID
    # values are sequential integers per lazyGas / GDS import; inclusive range
    # drops all markers between the LD block ends, not only strict LD members.
    rm_lo <- min(peak_block$ids)
    rm_hi <- max(peak_block$ids)
    pvalues <- pvalues[
      pvalues$variant_ID < rm_lo | pvalues$variant_ID > rm_hi,
      ,
      drop = FALSE
    ]
  }

  # Finalize GDS nodes by setting them to read mode
  .finalize_gdsn_peakcall(object = object, pheno_name = pheno_name)
}

## Function to identify the peak
.identify_peak <- function(pvalues,
                           variables,
                           peak_index = NULL,
                           rows_by_chr = NULL,
                           id_is_rowindex = NULL) {
  if (is.null(peak_index)) {
    peak_index <- which.max(pvalues$negLog10P)
  }
  peak_variant_ID <- pvalues$variant_ID[peak_index]

  if (is.null(id_is_rowindex)) {
    id_is_rowindex <- .peakcall_snp_id_is_rowindex(variables$snp_id)
  }
  peak_row <- if (isTRUE(id_is_rowindex)) {
    as.integer(peak_variant_ID)
  } else {
    match(peak_variant_ID, variables$snp_id)
  }
  if (is.na(peak_row)) {
    stop("Peak variant_ID not found in marker table: ", peak_variant_ID, call. = FALSE)
  }
  peak_chr <- variables$chr[[peak_row]]

  if (is.null(rows_by_chr)) {
    chr_with_peak <- which(variables$chr == peak_chr)
  } else {
    key <- as.character(peak_chr)
    chr_with_peak <- rows_by_chr[[key]]
    if (is.null(chr_with_peak)) {
      chr_with_peak <- which(variables$chr == peak_chr)
    }
  }
  id_in_chr <- variables$snp_id[chr_with_peak]
  peak_index_in_chr <- match(peak_variant_ID, id_in_chr)

  list(
    index = peak_index,
    variantID = peak_variant_ID,
    chr = peak_chr,
    chr_with_peak = chr_with_peak,
    id_in_chr = id_in_chr,
    peak_index_in_chr = peak_index_in_chr
  )
}

.peakcall_snp_id_is_rowindex <- function(snp_id) {
  n <- length(snp_id)
  if (!n) {
    return(TRUE)
  }
  is.integer(snp_id) &&
    snp_id[[1L]] == 1L &&
    snp_id[[n]] == n &&
    !is.unsorted(snp_id, strictly = TRUE)
}

.peakcall_rows_by_chr <- function(chr) {
  split(seq_along(chr), chr, drop = FALSE)
}

## Get Genotype Format
.get_geno_format <- function(object) {
  geno_format <- .store_read_scan_scalar(object, "geno_format")
  return(geno_format)
}

# Function to set the selection criteria based on the genotype format
.set_selection_criteria_for_peakcall <- function(object, geno_format, chr) {
  valid_chr <- getChromosome(object = object, valid = FALSE) %in% chr

  if (geno_format == "dosage"){
    # Set selection criteria for genotype or dosage data
    selection <- list(validSam(object = object),
                      validMar(object = object) & valid_chr)
    node <- "annotation/format/EDS/data"
    is_categorical <- FALSE

  } else {
    if(geno_format == "haplotype"){
      # Get the dimensions of the genotype data
      node <- "annotation/format/HAP/data"
      is_categorical <- TRUE

    } else if(geno_format == "corrected"){
      node <- "annotation/format/CGT/data"
      is_categorical <- FALSE

    } else {
      node <- "genotype/data"
      is_categorical <- FALSE
    }

    # Set selection criteria for genotype data
    obj_desp <- objdesp.gdsn(node = index.gdsn(node = object, path = node))
    selection <- list(rep(TRUE, obj_desp$dim[1]),
                      validSam(object = object),
                      validMar(object = object) & valid_chr)
  }

  return(list(node = node, selection = selection, is_categorical = is_categorical))
}

#' @importFrom gdsfmt exist.gdsn index.gdsn objdesp.gdsn
.retrieve_geno <- function(object, selection, geno_format, collapse = TRUE){
  out <- .get_data(object = object,
                   node = selection$node,
                   sel = selection$selection)
  if(geno_format == "haplotype"){
    if(collapse){
      out <- apply(X = out, MARGIN = 3, FUN = c)
    }
    out[out == 0] <- NA

  } else if(geno_format %in% c("corrected", "genotype")){
    out[out == 3] <- NA
    out <- .collapse_geno_alleles(out)

  } else if(geno_format == "dosage"){
    out[out == 63] <- NA
  }
  return(out)
}

.collapse_geno_alleles <- function(out) {
  d <- dim(out)
  if (is.null(d)) {
    return(out)
  }
  if (length(d) == 2L) {
    return(matrix(colSums(out, na.rm = TRUE), nrow = 1L))
  }
  if (length(d) == 3L) {
    # alleles x samples x markers → samples x markers (same as apply(..., sum, na.rm=TRUE))
    a <- out
    a[is.na(a)] <- 0
    if (d[1] == 2L) {
      return(a[1, , ] + a[2, , ])
    }
    res <- a[1, , ]
    if (d[1] > 1L) {
      for (i in seq_len(d[1])[-1]) {
        res <- res + a[i, , ]
      }
    }
    return(res)
  }
  out
}

.peakcall_resolve_n_threads <- function(n_threads) {
  if (is.null(n_threads) || !is.finite(n_threads) || n_threads < 1) {
    cores <- parallel::detectCores()
    if (is.na(cores) || cores < 2) {
      return(1L)
    }
    return(max(1L, as.integer(round(cores / 2))))
  }
  max(1L, as.integer(n_threads))
}

#' Identity / match rate LD for categorical (e.g. haplotype) genotypes.
#'
#' @keywords internal
.ld_identity_categorical <- function(geno, peak_variant_idx, n_threads = NULL) {
  n <- ncol(geno)
  geno_col <- as.integer(geno[, peak_variant_idx])
  n_threads <- .peakcall_resolve_n_threads(n_threads)
  # Chunk columns instead of one mclapply task per marker
  mc_cores <- min(n_threads, 16L, max(1L, as.integer(ceiling(n / 20000L))))
  if (mc_cores > 1L && n >= 2000L) {
    chunk <- ceiling(n / mc_cores)
    parts <- split(seq_len(n), ceiling(seq_len(n) / chunk))
    vals <- parallel::mclapply(
      parts,
      function(idx) {
        vapply(idx, function(j) {
          geno_col_j <- as.integer(geno[, j])
          geno_match <- geno_col == geno_col_j
          sum(geno_match, na.rm = TRUE) / sum(!is.na(geno_match))
        }, numeric(1))
      },
      mc.cores = mc_cores,
      mc.preschedule = TRUE
    )
    return(unlist(vals, use.names = FALSE))
  }
  vapply(seq_len(n), function(j) {
    geno_col_j <- as.integer(geno[, j])
    geno_match <- geno_col == geno_col_j
    sum(geno_match, na.rm = TRUE) / sum(!is.na(geno_match))
  }, numeric(1))
}

#' @import parallel
.calculate_ld <- function(geno, peak_variant_idx, is_categorical, n_threads = NULL) {
  if (!is.matrix(geno)) {
    geno <- as.matrix(geno)
  }
  if (is_categorical) {
    .ld_identity_categorical(geno, peak_variant_idx, n_threads = n_threads)
  } else {
    # Single-vector vs matrix Pearson r^2. pairwise.complete.obs is already
    # fast enough in C for ~1e5–1e6 markers; avoid R-level clean/dirty splits
    # that copy huge submatrices and regress wall time.
    as.numeric(
      stats::cor(
        geno[, peak_variant_idx],
        geno,
        use = "pairwise.complete.obs"
      )^2
    )
  }
}

# Function to identify peak block
.identify_peak_block <- function(peak_ld, variables, peak_info, threshold, geno_format) {
  if(geno_format == "haplotype"){
    ld_block <- peak_ld >= max(peak_ld, na.rm = TRUE) * threshold

  } else {
    ld_block <- peak_ld >= threshold
  }
  one_marker_up <- min(which(ld_block))
  if(one_marker_up > 1){
    ld_block[one_marker_up - 1] <- TRUE
    chr_start <- 0

  } else {
    chr_start <- 1
  }
  one_marker_down <- max(which(ld_block))
  if(!is.na(ld_block[one_marker_down + 1])){
    ld_block[one_marker_down + 1] <- TRUE
    chr_end <- 0

  } else {
    chr_end <- 1
  }
  ld_block <- which(ld_block)
  peak_block <- peak_info$id_in_chr[ld_block]
  ld_to_peak <- peak_ld[ld_block]
  # positions on this chromosome (same order as id_in_chr / peak_ld)
  pos <- variables$pos[peak_info$chr_with_peak]
  dist2peak <- pos[ld_block] - pos[peak_info$peak_index_in_chr]

  return(list(ids = peak_block, ld = ld_to_peak, dist = dist2peak,
              chr_start = chr_start, chr_end = chr_end))
}

# Function to save peak data (buffered for Parquet, GDS append for legacy)
.save_peak_data <- function(object, pheno_name, peak_info, peak_block, peak_id, peak_ld, node = "peakcall") {
  if (.store_is_gds(object)) {
    if (peak_id == 1) {
      gdsfmt::add.gdsn(
        node = gdsfmt::index.gdsn(node = object$root,
                                  path = paste0("lazygas/", node, "/peaks")),
        name = pheno_name, storage = "uint32", compress = "ZIP_RA", replace = TRUE,
        val = rbind(peak_id, peak_info$variantID)
      )
      gdsfmt::add.gdsn(
        node = gdsfmt::index.gdsn(node = object$root,
                                  path = paste0("lazygas/", node, "/blocks")),
        name = pheno_name, storage = "double", compress = "ZIP_RA", replace = TRUE,
        val = rbind(peak_id, peak_block$ids, peak_block$dist, peak_block$ld,
                    peak_block$chr_start, peak_block$chr_end)
      )
    } else {
      peaks_gdsn <- gdsfmt::index.gdsn(
        node = object$root,
        path = paste0("lazygas/", node, "/peaks/", pheno_name)
      )
      gdsfmt::append.gdsn(node = peaks_gdsn, val = rbind(peak_id, peak_info$variantID))
      blocks_gdsn <- gdsfmt::index.gdsn(
        node = object$root,
        path = paste0("lazygas/", node, "/blocks/", pheno_name)
      )
      gdsfmt::append.gdsn(
        node = blocks_gdsn,
        val = rbind(peak_id, peak_block$ids, peak_block$dist, peak_block$ld,
                    peak_block$chr_start, peak_block$chr_end)
      )
    }
  } else {
    .store_peak_buffer_append(
      object, section = node, pheno_name = pheno_name,
      peak_id = peak_id, peak_info = peak_info, peak_block = peak_block
    )
  }
}

## Sub-function to finalize peakcall storage
.finalize_gdsn_peakcall <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/peakcall/peaks/", pheno_name))
    )
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/peakcall/blocks/", pheno_name))
    )
  } else {
    .store_peak_buffer_flush(object, "peakcall", pheno_name)
  }
}


################################################################################
#' Plots the peaks for the specified phenotype and chromosome region.
#'
#' @param object A LazyGas object.
#' @param pheno The phenotype(s) to plot. Can be NULL, numeric, logical, or character.
#' @param chr The chromosome to plot. Can be NULL.

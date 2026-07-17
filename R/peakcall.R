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

  # Loop until there are no more p-values to process
  while (nrow(pvalues) > 0) {
    peak_id <- peak_id + 1
    if(peak_id > limit_peakcall){
      message("The number of peaks hit the limit.")
      break
    }
    message("Calling peak: ", peak_id)

    # Identify the peak information
    peak_info <- .identify_peak(pvalues = pvalues,
                                variables = variables)

    # Set the selection criteria based on the genotype format
    selection <- .set_selection_criteria_for_peakcall(object = object,
                                                      geno_format = geno_format,
                                                      chr = peak_info$chr)

    cache_key <- paste(geno_format, peak_info$chr, sep = ":")
    if (!exists(cache_key, envir = geno_cache, inherits = FALSE)) {
      assign(
        cache_key,
        .retrieve_geno(object = object,
                       selection = selection,
                       geno_format = geno_format),
        envir = geno_cache
      )
    }
    geno <- get(cache_key, envir = geno_cache, inherits = FALSE)

    # Calculate linkage disequilibrium (LD) for the identified peak
    peak_ld <- .calculate_ld(geno = as.matrix(geno),
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
    # values are sequential integers per lazyGas / GDS import; seq(min, max)
    # drops all markers between the LD block ends, not only strict LD members.
    rm_ids <- seq(min(peak_block$ids), max(peak_block$ids))
    pvalues <- subset(pvalues, subset = !variant_ID %in% rm_ids)
  }

  # Finalize GDS nodes by setting them to read mode
  .finalize_gdsn_peakcall(object = object, pheno_name = pheno_name)
}

## Function to identify the peak
.identify_peak <- function(pvalues, variables, peak_index = NULL) {
  if(is.null(peak_index)){
    peak_index <- which.max(pvalues$negLog10P)  # Find the index of the maximum p-value
  }
  peak_variant_ID <- pvalues$variant_ID[peak_index]  # Get the variant ID of the peak
  peak_chr <- variables$chr[variables$snp_id == peak_variant_ID]  # Get the chromosome of the peak

  chr_with_peak <- which(variables$chr == peak_chr)
  id_in_chr <- variables$snp_id[chr_with_peak]
  peak_index_in_chr <- which(id_in_chr == peak_variant_ID)

  list(
    index = peak_index,
    variantID = peak_variant_ID,
    chr = peak_chr,
    chr_with_peak = chr_with_peak,
    id_in_chr = id_in_chr,
    peak_index_in_chr = peak_index_in_chr
  )
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
  if (length(dim(out)) == 2L) {
    return(matrix(colSums(out, na.rm = TRUE), nrow = 1L))
  }
  if (length(dim(out)) == 3L) {
    return(apply(out, c(2L, 3L), sum, na.rm = TRUE))
  }
  out
}

#' @import parallel
.calculate_ld <- function(geno, peak_variant_idx, is_categorical, n_threads = NULL) {
  n <- ncol(geno)
  peak_ld <- numeric(n)

  if (is_categorical) {
    geno_col <- as.integer(geno[, peak_variant_idx])
    levels <- sort(unique(geno_col))
    num_levels <- length(levels)

    # Parallel processing setup
    if(is.null(n_threads)){
      cores <- detectCores()
      if(cores > 1){
        n_threads <- round(cores / 2)
      }
    }

    geno_identity <- mclapply(X = 1:n, mc.cores = n_threads, mc.preschedule = TRUE,
                              FUN = function(j) {
                                geno_col_j <- as.integer(geno[, j])
                                geno_match <- geno_col == geno_col_j
                                geno_identity <- sum(geno_match, na.rm = TRUE) / sum(!is.na(geno_match))
                                return(geno_identity)
                              })
    peak_ld <- unlist(geno_identity)

  } else {
    # Parallel processing setup
    corr_values <- cor(geno[, peak_variant_idx], geno, use = "pairwise.complete.obs") ^ 2
    peak_ld <- corr_values
  }

  return(peak_ld)
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
  pos <- variables$pos[variables$snp_id %in% peak_info$id_in_chr]
  dist2peak <- na.omit(pos[ld_block] - pos[peak_info$peak_index_in_chr])

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

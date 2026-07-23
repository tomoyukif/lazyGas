#' @param n_threads description
#' @param refine_position
#'
#' @return None. The function modifies the LazyGas object in place.
#' @export
#'
setGeneric("recalcAssoc", function(object,
                                   n_threads = NULL,
                                   refine_position = FALSE,
                                   grouping_threshold = 0.05,
                                   ...)
  standardGeneric("recalcAssoc"))

#'
#' @rdname recalcAssoc
#' @method recalcAssoc LazyGas
#' @export
#'
setMethod("recalcAssoc",
          "LazyGas",
          function(object, n_threads, refine_position, grouping_threshold){
            if (!.store_section_exists(object, "peakcall")) {
              stop("No peakcall data in the input LazyGas object.\n",
                   "Run callPeakBlock() to call peaks.")
            }

            if (is.null(n_threads)) {
              cores <- parallel::detectCores()
              if (cores > 1) {
                n_threads <- round(cores / 2)
              }
            }

            kruskal <- .store_read_scan_scalar(object, "kruskal")
            if (is.null(kruskal)) {
              kruskal <- character()
            }

            ## Function to create the necessary folders in the GDS object
            .create_recalc_folders(object = object)

            pheno <- getPheno(object = object)

            for(i in seq_along(pheno$pheno_names)){
              pheno_name <- pheno$pheno_names[i]
              dokruskal <- pheno_name %in% kruskal

              if (dokruskal) {
                message("Skip recalculation for non-parametric data: ", pheno_name)
                .store_write_recalc_empty(object, pheno_name)

              } else {
                .peakrecalculator(object = object,
                                  pheno = pheno,
                                  pheno_name = pheno_name,
                                  binary = pheno$pheno_type$binary[i],
                                  n_threads = n_threads,
                                  refine_position = refine_position,
                                  grouping_threshold = grouping_threshold)
              }

              .finalize_gdsn_recalc(object = object, pheno_name = pheno_name)
            }
          }
)

.create_recalc_folders <- function(object) {
  if (.store_is_gds(object)) {
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas",
                 new_node = "recalc",
                 is_folder = TRUE)
    peaks_gdsn <- .create_gdsn(root_node = object$root,
                               target_node = "lazygas/recalc",
                               new_node = "peaks",
                               is_folder = TRUE)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc",
                 new_node = "blocks",
                 is_folder = TRUE)
    groups_gdsn <- .create_gdsn(root_node = object$root,
                                target_node = "lazygas/recalc",
                                new_node = "groups",
                                is_folder = TRUE)
    att <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(object$root, path = "lazygas/peakcall/peaks")
    )
    gdsfmt::put.attr.gdsn(node = peaks_gdsn, name = "signif", val = att$signif)
    gdsfmt::put.attr.gdsn(node = peaks_gdsn, name = "threshold", val = att$threshold)
  } else {
    dir.create(file.path(.store_path(object), "recalc"),
               recursive = TRUE, showWarnings = FALSE)
  }
}

## Sub-function to finalize recalc storage
.finalize_gdsn_recalc <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/recalc/peaks/", pheno_name))
    )
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/recalc/blocks/", pheno_name))
    )
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/recalc/groups/", pheno_name))
    )
  }
}

# Recalculate peaks for a given phenotype.
#' @importFrom dplyr setequal
.peakrecalculator <- function(object, pheno, pheno_name, binary, n_threads, refine_position, grouping_threshold){
  message("Processing: ", pheno_name)

  peakcall <- .get_peakcall(object = object, pheno_name = pheno_name)
  att <- attributes(peakcall)
  # peak_grp_summary <- NULL

  if (is.null(peakcall)) {
    .store_write_recalc_empty(object, pheno_name)

  } else {
    peak_obj <- .make_peakobj(object = object,
                              threshold = att$threshold,
                              signif = att$signif,
                              pheno = pheno,
                              pheno_name = pheno_name,
                              peakcall = peakcall,
                              binary = binary)
    message("Grouping associations...")

    # Cap workers for tiny peak-vs-peak GLM jobs (forking 60+ procs is slower)
    n_peaks <- length(peak_obj$peak_variant_id)
    grp_threads <- min(.peakcall_resolve_n_threads(n_threads), max(1L, n_peaks - 1L), 8L)

    cov_scan <- lapply(X = peak_obj$peak_variant_id,
                       FUN = .composite,
                       peak_obj = peak_obj,
                       n_threads = grp_threads,
                       mode = "excl")

    peak_grp <- .group_peaks(cov_scan = cov_scan,
                             peak_variant_id = peak_obj$peak_variant_id,
                             grouping_threshold = grouping_threshold)

    new_peaks <- .get_newpeaks(object = object,
                               peak_grp = peak_grp$group,
                               peak_obj = peak_obj,
                               n_threads = n_threads)

    peak_obj <- .remake_peakobj(object = object,
                                peak_obj = peak_obj,
                                peak_variant_id = new_peaks$newpeaks)

    message("Recalling peak blocks...")
    variables <- .peakcall_variables(object = object)
    id_is_rowindex <- .peakcall_snp_id_is_rowindex(variables$snp_id)
    rows_by_chr <- .peakcall_rows_by_chr(variables$chr)
    geno_cache <- new.env(parent = emptyenv())
    selection_cache <- new.env(parent = emptyenv())
    peak_block <- lapply(
      X = peak_obj$peak_variant_id,
      FUN = .recallPeak,
      peak_obj = peak_obj,
      object = object,
      peakcall = peakcall,
      n_threads = n_threads,
      geno_cache = geno_cache,
      selection_cache = selection_cache,
      variables = variables,
      rows_by_chr = rows_by_chr,
      id_is_rowindex = id_is_rowindex
    )
    peak_obj$peak_block <- lapply(seq_along(peak_block),
                                  function(i) {peak_block[[i]]$variant_ID} )

    if(refine_position){
      message("Refine peakpositions...")
      new_peak_blocks <- .refine_peakpositions(peak_obj = peak_obj)
      out1 <- unique(subset(new_peak_blocks,
                            select = c(peak_ID, peak_variant_ID)))
      out2 <- unique(subset(new_peak_blocks,
                            select = c(peak_ID, variant_ID,
                                       dist2peak, LD2peak, chr_start, chr_end,
                                       P.model:negLog10P)))
      out3 <- subset(peak_obj$peakcall,
                     subset = variant_ID %in% unlist(new_peaks$member),
                     select = c(peak_ID, variant_ID,
                                P.model:negLog10P))
      peak_id_map <- unique(subset(peak_obj$peakcall,
                                   subset = variant_ID == peak_variant_ID,
                                   select = c(peak_ID, variant_ID)))
      grouped_with <- lapply(seq_along(new_peaks$member), function(i){
        return(data.frame(
          new_peak = peak_id_map$peak_ID[match(new_peaks$newpeaks[i],
                                               peak_id_map$variant_ID)],
          member = new_peaks$member[[i]]
        ))
      })
      grouped_with <- do.call("rbind", grouped_with)
      hit <- match(out3$variant_ID, grouped_with$member)
      out3$grouped_with <- grouped_with$new_peak[hit]
      hit <- match(out3$variant_ID, peak_grp$pvalue$variant_ID)
      out3$full_vs_reduce_pvalue <- peak_grp$pvalue$p_values[hit]

    } else {
      out1 <- unique(subset(peak_obj$peakcall,
                            subset = peak_variant_ID %in% new_peaks$newpeaks,
                            select = c(peak_ID, peak_variant_ID)))
      out2 <- unique(subset(peak_obj$peakcall,
                            subset = peak_variant_ID %in% new_peaks$newpeaks,
                            select = c(peak_ID, variant_ID,
                                       dist2peak, LD2peak, chr_start, chr_end,
                                       FDR:negLog10P)))
      out3 <- subset(peak_obj$peakcall,
                     subset = variant_ID %in% unlist(new_peaks$member) & variant_ID == peak_variant_ID,
                     select = c(peak_ID, variant_ID,
                                FDR:negLog10P))
      peak_id_map <- unique(subset(peak_obj$peakcall,
                                   subset = variant_ID == peak_variant_ID,
                                   select = c(peak_ID, variant_ID)))
      grouped_with <- lapply(seq_along(new_peaks$member), function(i){
        return(data.frame(
          new_peak = peak_id_map$peak_ID[match(new_peaks$newpeaks[i],
                                               peak_id_map$variant_ID)],
          member = new_peaks$member[[i]]
        ))
      })
      grouped_with <- do.call("rbind", grouped_with)
      hit <- match(out3$variant_ID, grouped_with$member)
      out3$grouped_with <- grouped_with$new_peak[hit]
      hit <- match(out3$variant_ID, peak_grp$pvalue$variant_ID)
      out3$full_vs_reduce_pvalue <- peak_grp$pvalue$p_values[hit]
    }

    .store_save_recalc_results(object, pheno_name, out1, out2, out3)
  }
}

# Create a peak object for recalculating associations.
.make_peakobj <- function(object, threshold, signif, pheno, pheno_name, peakcall, binary){
  peak_variant_id <- unique(peakcall$peak_variant_ID)

  geno_format <- .get_geno_format(object = object)

  selection <- .set_selection_criteria_for_recalc(object = object,
                                                  geno_format = geno_format,
                                                  peak_variant_id = peak_variant_id)

  geno <- .retrieve_geno(object = object,
                         selection = selection,
                         geno_format = geno_format,
                         collapse = FALSE)

  if(geno_format == "haplotype"){
    if(length(dim(geno)) == 3){
      geno[, , order(peak_variant_id)] <- geno
    }

  } else if(length(dim(geno)) == 2) {
    geno[, order(peak_variant_id)] <- geno

  } else {
    geno <- matrix(data = geno, ncol = 1)
  }

  if(geno_format == "haplotype"){
    na_val <- 0

  } else if(geno_format == "dosage"){
    na_val <- 63

  } else {
    na_val <- 3
  }

  null_formula <- .store_read_scan_scalar(object, "null_formula")
  if (identical(null_formula, "NULL") || is.null(null_formula)) {
    null_formula <- NULL
  }

  fixed_effect <- .store_read_fixed_effect(object)

  conv_fun_txt <- .store_read_scan_scalar(object, "conv_fun")
  formula_txt <- .store_read_scan_scalar(object, "formula")

  out <- list(
    geno = geno,
    geno_format = geno_format,
    pheno = .standardize(val = pheno$pheno[, pheno_name], binary = binary),
    threshold = threshold,
    signif = signif,
    conv_fun = eval(parse(text = conv_fun_txt)),
    formula = formula_txt,
    null_formula = null_formula,
    fixed_effect = fixed_effect,
    na_val = na_val,
    peakcall = peakcall,
    peak_variant_id = peak_variant_id,
    peak_block = tapply(peakcall$variant_ID, peakcall$peak_variant_ID, c),
    snp_id = selection$snp_id,
    n_sample = nsam(object = object),
    binary = binary
  )
  return(out)
}

# Function to set the selection criteria based on the genotype format
.set_selection_criteria_for_recalc <- function(object,
                                               geno_format,
                                               peak_variant_id) {
  snp_id <- getMarID(object = object, valid = FALSE)

  if (geno_format == "dosage"){
    selection <- list(validSam(object = object),
                      snp_id %in% peak_variant_id)
    node <- "annotation/format/EDS/data"
    is_categorical <- FALSE

  } else {
    if(geno_format == "haplotype"){
      node <- "annotation/format/HAP/data"
      is_categorical <- TRUE

    } else if(geno_format == "corrected"){
      node <- "annotation/format/CGT/data"
      is_categorical <- FALSE

    } else {
      node <- "genotype/data"
      is_categorical <- FALSE
    }

    # Set selection criteria for haplotype data
    obj_desp <- objdesp.gdsn(node = index.gdsn(node = object, path = node))
    selection <- list(rep(TRUE, obj_desp$dim[1]),
                      validSam(object = object),
                      snp_id %in% peak_variant_id)
  }

  return(list(node = node,
              selection = selection,
              snp_id = snp_id,
              is_categorical = is_categorical))
}

.composite <- function(query_peak,
                       subject_peak = NULL,
                       peak_obj,
                       output_all = FALSE,
                       n_threads,
                       mode){
  if(mode == "composite"){
    formula_terms <- unlist(strsplit(peak_obj$formula, split = "\\+|\\*|\\:"))
    formula_terms <- gsub("\\s", "", formula_terms)
    null_df <- NULL
    null_fml <- NULL
    for(j in seq_along(query_peak)){
      index <- peak_obj$peak_variant_id %in% query_peak[j]
      tmp_df <- .makeDF(g = peak_obj$geno[, index],
                        phe = peak_obj$pheno,
                        conv_fun = peak_obj$conv_fun,
                        formula = peak_obj$formula)
      if(is.null(null_df)){
        names(tmp_df$df)[-1] <- paste(names(tmp_df$df)[-1], j, sep = "_")
        null_df <- tmp_df

      } else {
        names(tmp_df$df)[-1] <- paste(names(tmp_df$df)[-1], j, sep = "_")
        null_df$df <- cbind(null_df$df, subset(tmp_df$df, select = -phe))
      }
      add_fml <- peak_obj$formula
      for(k in seq_along(formula_terms)){
        add_fml <- gsub(formula_terms[k], paste(formula_terms[k], j, sep = "_"), add_fml)
      }
      if(is.null(null_fml)){
        null_fml <- paste0("phe ~ ", add_fml)
      } else {
        null_fml <- paste(null_fml, add_fml, sep = "+")
      }
    }
    null_df$fml <- formula(null_fml)

    if(is.null(subject_peak)){
      index <- peak_obj$peak_variant_id %in% query_peak
      subject_peak <- peak_obj$peak_variant_id[!index]

    } else {
      index <- !peak_obj$peak_variant_id %in% subject_peak
    }

  }  else if(mode == "excl"){
    index <- peak_obj$peak_variant_id %in% query_peak
    subject_peak <- peak_obj$peak_variant_id[!index]

    if(length(subject_peak) == 0){
      return(data.frame(variant_ID = NA, FDR = NA))
    }

    if(peak_obj$geno_format == "haplotype"){
      g <- peak_obj$geno[, , index]
    } else {
      g <- peak_obj$geno[, index]
    }

    fml <- peak_obj$formula
    null_df <- .makeDF(g = g,
                       phe = peak_obj$pheno,
                       conv_fun = peak_obj$conv_fun,
                       formula = fml)
  }


  if(peak_obj$geno_format == "haplotype"){
    target_geno <- peak_obj$geno[, , !index, drop = FALSE]
    if(length(dim(target_geno)) == 2){
      target_geno <- array(data = target_geno, dim = c(dim(target_geno), 1))
    }
    n_target <- dim(target_geno)[3]
    target_list <- lapply(seq_len(n_target), function(j) target_geno[, , j])

  } else {
    target_geno <- peak_obj$geno[, !index, drop = FALSE]
    if(is.null(dim(target_geno))){
      target_geno <- matrix(data = target_geno, ncol = 1)
    }
    n_target <- ncol(target_geno)
    target_list <- lapply(seq_len(n_target), function(j) target_geno[, j])
  }

  n_threads <- .peakcall_resolve_n_threads(n_threads)
  # Small peak-count jobs: serial is faster than mass fork
  n_threads <- min(n_threads, max(1L, n_target), 8L)

  na_val <- peak_obj$na_val

  glm_one <- function(g) {
    g[g == na_val] <- NA
    if (all(is.na(g))) {
      return(NA)
    }
    if (length(unique(stats::na.omit(as.vector(g)))) == 1) {
      return(NA)
    }
    df <- .makeDF(g = g,
                  phe = peak_obj$pheno,
                  conv_fun = peak_obj$conv_fun,
                  formula = peak_obj$formula)
    if (length(query_peak) != 0) {
      tmp <- subset(null_df$df, select = -phe)
      names(tmp) <- paste0("qtl_", names(tmp))
      df$df <- cbind(df$df, tmp)
      df$fml <- paste(c(df$fml, names(tmp)), collapse = " + ")
    }

    if (peak_obj$binary) {
      family <- "binomial"
    } else {
      family <- "gaussian"
    }

    if (!all(is.na(peak_obj$fixed_effect))) {
      df$df <- cbind(df$df, peak_obj$fixed_effect)
      null_df$df <- cbind(null_df$df, peak_obj$fixed_effect)
    }

    .doGLM(df = df, family = family, null_df = null_df)
  }

  if (n_threads <= 1L || n_target <= 1L) {
    p_values <- lapply(target_list, glm_one)
  } else {
    p_values <- parallel::mclapply(
      X = target_list,
      mc.cores = n_threads,
      mc.preschedule = TRUE,
      FUN = glm_one
    )
  }

  if(is.list(p_values)){
    p_values <- data.frame(do.call("rbind", p_values))

  } else {
    p_values <- data.frame(t(p_values))
  }

  if(output_all){
    p_values$FDR <- p.adjust(p_values$P.model, "fdr")
    return(data.frame(variant_ID = subject_peak, p_values))

  } else {
    return(data.frame(variant_ID = subject_peak, p_values = p_values$P.model))
  }
}

.group_peaks <- function(cov_scan, peak_variant_id, grouping_threshold){
  peak_list <- peak_variant_id
  out2 <- out1 <- NULL
  for(j in seq_along(cov_scan)){
    if(!peak_variant_id[j] %in% peak_list){
      next
    }
    grp <- c(peak_variant_id[j],
             cov_scan[[j]]$variant_ID[cov_scan[[j]]$p_values > grouping_threshold])
    grp <- grp[grp %in% peak_list]
    out1 <- c(out1, list(grp))
    out2 <- rbind(out2, data.frame(cov_scan[[j]][cov_scan[[j]]$variant_ID %in% grp, ]))
    peak_list <- peak_list[!peak_list %in% grp]
  }
  return(list(group = out1, pvalue = out2))
}

.get_newpeaks <- function(object, peak_grp, peak_obj, n_threads){
  out1 <- NULL
  out2 <- NULL
  for(j in seq_along(peak_grp)){
    grp <- peak_grp[[j]]
    query_peak <- grp[1]
    out1 <- c(out1, query_peak)
    out2 <- c(out2, list(grp))
  }
  return(list(newpeaks = out1, member = out2))
}

.remake_peakobj <- function(object, peak_obj, peak_variant_id){
  peak_variant_id <- unique(peak_variant_id)
  old_ids <- peak_obj$peak_variant_id

  if (!is.null(peak_obj$geno) && all(peak_variant_id %in% old_ids)) {
    keep <- match(peak_variant_id, old_ids)
    if (peak_obj$geno_format == "haplotype" && length(dim(peak_obj$geno)) == 3L) {
      geno <- peak_obj$geno[, , keep, drop = FALSE]
    } else {
      geno <- peak_obj$geno[, keep, drop = FALSE]
    }
  } else {
    selection <- .set_selection_criteria_for_recalc(object = object,
                                                    geno_format = peak_obj$geno_format,
                                                    peak_variant_id = peak_variant_id)

    geno <- .retrieve_geno(object = object,
                           selection = selection,
                           geno_format = peak_obj$geno_format,
                           collapse = FALSE)

    if(is.vector(geno)){
      geno <- matrix(geno, ncol = 1)
    }
  }

  peak_obj$peak_variant_id <- peak_variant_id

  if(peak_obj$geno_format == "haplotype"){
    if(length(dim(geno)) == 3){
      geno[, , order(peak_variant_id)] <- geno
    }

  } else {
    geno[, order(peak_variant_id)] <- geno
  }

  peak_obj$geno <- geno
  return(peak_obj)
}

.recallPeak <- function(peak_id,
                        peak_obj,
                        object,
                        peakcall,
                        n_threads,
                        geno_cache = NULL,
                        selection_cache = NULL,
                        variables = NULL,
                        rows_by_chr = NULL,
                        id_is_rowindex = NULL) {
  if (is.null(variables)) {
    variables <- .peakcall_variables(object = object)
  }
  if (is.null(id_is_rowindex)) {
    id_is_rowindex <- .peakcall_snp_id_is_rowindex(variables$snp_id)
  }
  if (is.null(rows_by_chr)) {
    rows_by_chr <- .peakcall_rows_by_chr(variables$chr)
  }

  peak_row <- if (isTRUE(id_is_rowindex)) {
    as.integer(peak_id)
  } else {
    match(peak_id, variables$snp_id)
  }
  if (is.na(peak_row)) {
    stop("Peak variant_ID not found in marker table: ", peak_id, call. = FALSE)
  }
  peak_chr <- variables$chr[[peak_row]]
  key <- as.character(peak_chr)
  chr_with_peak <- rows_by_chr[[key]]
  if (is.null(chr_with_peak)) {
    chr_with_peak <- which(variables$chr == peak_chr)
  }
  id_in_chr <- variables$snp_id[chr_with_peak]
  peak_index_in_chr <- match(peak_id, id_in_chr)
  peak_info <- list(
    index = peak_row,
    variantID = peak_id,
    chr = peak_chr,
    chr_with_peak = chr_with_peak,
    id_in_chr = id_in_chr,
    peak_index_in_chr = peak_index_in_chr
  )

  if (is.null(selection_cache)) {
    selection_cache <- new.env(parent = emptyenv())
  }
  if (!exists(key, envir = selection_cache, inherits = FALSE)) {
    assign(
      key,
      .set_selection_criteria_for_peakcall(
        object = object,
        geno_format = peak_obj$geno_format,
        chr = peak_info$chr
      ),
      envir = selection_cache
    )
  }
  selection <- get(key, envir = selection_cache, inherits = FALSE)

  cache_key <- paste(peak_obj$geno_format, peak_info$chr, sep = ":")
  if (!is.null(geno_cache) && exists(cache_key, envir = geno_cache, inherits = FALSE)) {
    geno <- get(cache_key, envir = geno_cache, inherits = FALSE)
  } else {
    geno <- .retrieve_geno(object = object,
                           selection = selection,
                           geno_format = peak_obj$geno_format)
    if (!is.matrix(geno)) {
      geno <- as.matrix(geno)
    }
    storage.mode(geno) <- "double"
    if (!is.null(geno_cache)) {
      assign(cache_key, geno, envir = geno_cache)
    }
  }

  peak_ld <- .calculate_ld(
    geno = geno,
    peak_variant_idx = peak_info$peak_index_in_chr,
    is_categorical = selection$is_categorical,
    n_threads = n_threads
  )

  peak_block <- .identify_peak_block(
    peak_ld = peak_ld,
    variables = variables,
    peak_info = peak_info,
    threshold = peak_obj$threshold,
    geno_format = peak_obj$geno_format
  )
  data.frame(
    dist2peak = peak_block$dist,
    LD2peak = peak_block$ld,
    variant_ID = peak_block$ids,
    chr_start = peak_block$chr_start,
    chr_end = peak_block$chr_end
  )
}

.make_newblocks <- function(i, object, peak_obj, n_threads){
  query_peak <- peak_obj$peak_variant_id[-i]
  hit_id <- peak_obj$peakcall$variant_ID %in% peak_obj$peak_variant_id[i]
  i_peak <- peak_obj$peakcall[hit_id, ]
  variant_id <- getMarID(object = object, valid = TRUE, chr = i_peak$peak_Chr)
  peak_obj <- .remake_peakobj(object = object,
                              peak_obj = peak_obj,
                              peak_variant_id = sort(c(variant_id, query_peak)))
  peakcall <- .composite(query_peak = query_peak,
                         peak_obj = peak_obj,
                         output_all = TRUE,
                         n_threads = n_threads,
                         mode = "composite")
  peakcall$P.model <- peakcall$P.model
  peakcall$P.model[peakcall$P.model %in% c(0, 1)] <- NA
  peak_id <- peakcall$variant_ID[which.min(peakcall$P.model)]
  peakcall$peak_ID <- i
  recall_peak <- .recallPeak(peak_id = peak_id,
                             peak_obj = peak_obj,
                             object = object,
                             peakcall = peakcall,
                             n_threads = n_threads)

  peak_obj <- .remake_peakobj(object = object,
                              peak_obj = peak_obj,
                              peak_variant_id = c(query_peak,
                                                  recall_peak$variant_ID))
  peakcall <- .composite(query_peak = query_peak,
                         subject_peak = recall_peak$variant_ID,
                         peak_obj = peak_obj,
                         output_all = TRUE,
                         n_threads = n_threads,
                         mode = "composite")
  peakcall$negLog10P <- -log10(peakcall$P.model)
  out <- data.frame(peak_ID = i,
                    peak_variant_ID = peak_id,
                    peakNegLog10P = peakcall$negLog10P[peakcall$variant_ID == peak_id],
                    recall_peak, subset(peakcall, select = -variant_ID))
  return(out)
}

.refine_peakpositions <- function(peak_obj){
  new_peak_blocks <- lapply(X = seq_along(peak_obj$peak_variant_id),
                            FUN = .make_newblocks,
                            object = object,
                            peak_obj = peak_obj,
                            n_threads = n_threads)

  peak_order <- sapply(new_peak_blocks, function(x) {max(x$peakNegLog10P)})
  new_peak_blocks <- new_peak_blocks[order(peak_order, decreasing = TRUE)]
  header <- sapply(new_peak_blocks, names)
  common <- apply(header, 1, function(x){length(unique(x)) == 1})
  common_header <- header[common]
  uncommon_header <- sort(unique(header[!common]))
  uncommon_header <- c(uncommon_header[grepl("^P.", uncommon_header)],
                       uncommon_header[grepl("^Coef.", uncommon_header)])
  new_peak_blocks <- lapply(seq_along(new_peak_blocks), function(k){
    new_peak_blocks[[k]]$peak_ID <- k
    common_df <- subset(new_peak_blocks[[k]], select = common_header)
    not_exist <- uncommon_header[!uncommon_header %in% names(new_peak_blocks[[k]])]
    not_exist_df <- matrix(NA, nrow = nrow(new_peak_blocks[[k]]), ncol = length(not_exist))
    not_exist_df <- data.frame(not_exist_df)
    names(not_exist_df) <- not_exist
    new_peak_blocks[[k]] <- cbind(new_peak_blocks[[k]], not_exist_df)
    uncommon_df <- subset(new_peak_blocks[[k]], select = uncommon_header)
    new_peak_blocks[[k]] <- cbind(common_df, uncommon_df)
    return(new_peak_blocks[[k]])
  })
  new_peak_blocks <- do.call("rbind", new_peak_blocks)
  return(new_peak_blocks)
}

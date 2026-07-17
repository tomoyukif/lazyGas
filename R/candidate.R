################################################################################

#' @export
setGeneric("listCandidate", function(object,
                                     gff,
                                     snpeff = NULL,
                                     ann = NULL,
                                     recalc = FALSE,
                                     ...)
  standardGeneric("listCandidate"))

setMethod("listCandidate",
          "LazyGas",
          function(object,
                   gff,
                   snpeff = NULL,
                   ann = NULL,
                   recalc = FALSE){
            if (recalc) {
              if (!.store_section_exists(object, "recalc")) {
                stop("No recalc data in the input LazyGas object.\n",
                     "Run recalcAssoc() to recalculate associations.")
              }
            } else {
              if (!.store_section_exists(object, "peakcall")) {
                stop("No peakcall data in the input LazyGas object.\n",
                     "Run callPeakBlock() to call peak blocks.")
              }
            }

            chr <- unique(getChromosome(object))
            .validateGFF(chr = chr, gff = gff)
            if(!is.null(snpeff)){
              .validateSnpEff(chr = chr, snpeff = snpeff)
            }
            .validateANN(gff = gff, ann = ann)

            ## Function to create the necessary folders in the GDS object
            .create_candidate_folders(object = object)

            pheno <- getPheno(object = object)

            for(i in seq_along(pheno$pheno_names)){
              pheno_name <- pheno$pheno_names[i]

              .candidatelistor(object = object,
                               ann = ann,
                               gff = gff,
                               snpeff = snpeff,
                               pheno_name = pheno_name,
                               recalc = recalc)

              .finalize_gdsn_candidate(object = object, pheno_name = pheno_name)
            }
          })

#' @importFrom GenomeInfoDb seqlevels
.validateGFF <- function(chr, gff){
  if(!inherits(gff, "GRanges")){
    stop("The input gff must be a GRanges class object imported by rtracklayer::import.gff3().",
         call. = FALSE)
  }

  hit <- chr %in% seqlevels(gff)
  if(all(!hit)){
    stop("Any of chromosome IDs do not match the chromosomes IDs in the given GFF data.",
         call. = FALSE)
  }
  if(!all(hit)){
    warning("The following chromosome IDs do not appeare in the given GFF data: ",
            chr[!hit])
  }
}

.validateSnpEff <- function(chr, snpeff){
  if(!inherits(snpeff, "snpeff_gds")){
    stop("The input snpeff must be a snpeff_gds class object loaded by lazyGas::open_snpeff().",
         call. = FALSE)
  }

  snpeff_chr <- unique(read.gdsn(index.gdsn(snpeff$root, "chromosome")))
  hit <- chr %in% snpeff_chr
  if(all(!hit)){
    stop("Any of chromosome IDs do not match the chromosomes IDs in the given snpeff GDS data.",
         call. = FALSE)
  }
  if(!all(hit)){
    warning("The following chromosome IDs do not appeare in the given snpeff GDS data: ",
            chr[!hit])
  }
}

.validateANN <- function(gff, ann){
  if(!inherits(ann, "data.frame")){
    stop("The input ann must be a data.frame class object.",
         call. = FALSE)
  }

  if(is.null(ann$Gene_ID)){
    stop("The input ann should have a column named 'Gene_ID' that should match the IDs in the input GFF data.",
         call. = FALSE)
  }

  hit <- ann$Gene_ID %in% gff$ID
  if(all(!hit)){
    stop("Any of Gene_IDs in ann do not match the IDs in the given GFF data.",
         call. = FALSE)
  }
  if(!all(hit)){
    warning("The following Gene_IDs in ann do not appeare in the given GFF data: ",
            unique(ann$Gene_ID[!hit]))
  }
}

.create_candidate_folders <- function(object) {
  if (.store_is_gds(object)) {
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas",
                 new_node = "candidate",
                 is_folder = TRUE)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas",
                 new_node = "snpeff",
                 is_folder = TRUE)
  } else {
    dir.create(file.path(.store_path(object), "candidate"),
               recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(.store_path(object), "snpeff"),
               recursive = TRUE, showWarnings = FALSE)
  }
}

## Sub-function to finalize candidate storage (GDS legacy)
.finalize_gdsn_candidate <- function(object, pheno_name) {
  if (.store_is_gds(object)) {
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/candidate/", pheno_name))
    )
    gdsfmt::readmode.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0("lazygas/snpeff/", pheno_name))
    )
  }
}

#' @importFrom dplyr left_join
.candidatelistor <- function(object,
                             ann,
                             gff,
                             snpeff,
                             pheno_name,
                             recalc){
  message("Processing: ", pheno_name)
  peakcall <- .get_peakcall(object = object,
                            pheno_name = pheno_name,
                            recalc = recalc)

  if (is.null(peakcall)) {
    .store_write_candidate(object, pheno_name, candidate = NULL, snpeff = NULL)

  } else {
    candidate_list <- list()
    snpeff_list <- list()
    snpeff_index <- if (!is.null(snpeff)) .snpeff_index(snpeff) else NULL
    for(i_peak in unique(peakcall$peak_ID)){
      peakblock <- peakcall[peakcall$peak_ID == i_peak, ]
      tmp <- .getCandidate(peakblock = peakblock,
                           gff = gff,
                           snpeff = snpeff,
                           snpeff_index = snpeff_index)
      candidate_list[[length(candidate_list) + 1L]] <- tmp$candidate_list
      if (!is.null(tmp$snpeff_out)) {
        snpeff_list[[length(snpeff_list) + 1L]] <- tmp$snpeff_out
      }
    }
    candidate_list <- if (length(candidate_list) > 0L) {
      do.call(rbind, candidate_list)
    } else {
      NULL
    }
    snpeff_list <- if (length(snpeff_list) > 0L) {
      do.call(rbind, snpeff_list)
    } else {
      NULL
    }

    if(!is.null(candidate_list)){
      if(!is.null(ann)){
        candidate_list <- left_join(candidate_list, ann, by = "Gene_ID")

      }
      candidate_list[is.na(candidate_list)] <- ""
    }

    snpeff_out <- NULL
    if (!is.null(snpeff_list)) {
      snpeff_list[is.na(snpeff_list)] <- ""
      snpeff_out <- snpeff_list
    }
    .store_write_candidate(object, pheno_name, candidate_list, snpeff_out)
  }
}

#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
.getCandidate <- function(peakblock, gff, snpeff, snpeff_index = NULL){
  peak_variant_id <- peakblock$peak_variant_ID[1]
  is_peak <- peakblock$variant_ID == peak_variant_id
  if(peakblock$chr_start[1] == 1){
    peak_start <- 1

  } else {
    peak_start <- peakblock$Pos[1]
  }
  if(peakblock$chr_end[1] == 1){
    peak_end <- 2^30

  } else {
    peak_end <- tail(peakblock$Pos, 1)
  }
  peak_gff <- GRanges(seqnames = peakblock$Chr[is_peak],
                      ranges = IRanges(start = peak_start,
                                       end = peak_end))
  type_hit <- gff$type %in% "gene"
  type_hit[is.na(type_hit)] <- FALSE
  gene_gff <- gff[type_hit]
  hit <- gene_gff[queryHits(findOverlaps(gene_gff, peak_gff))]
  hit <- hit[order(start(hit))]
  hit$dist2peak <- start(hit) - peakblock$Pos[is_peak]
  nearest_p <- vapply(start(hit), function(x) {
    peakblock$negLog10P[which.min(abs(peakblock$Pos - x))]
  }, numeric(1))
  if(length(hit) == 0){
    out <- data.frame(peak_ID = peakblock$peak_ID[1],
                      Gene_ID = NA,
                      Gene_chr = NA,
                      Gene_start = NA,
                      dist2peak = NA,
                      negLog10P = NA)

  } else {
    out <- data.frame(peak_ID = peakblock$peak_ID[1],
                      Gene_ID = hit$ID,
                      Gene_chr = as.character(seqnames(hit)),
                      Gene_start = start(hit),
                      dist2peak = hit$dist2peak,
                      negLog10P = nearest_p)
  }

  if(!is.null(snpeff)){
    snpeff_out <- .addSnpEff(snpeff = snpeff,
                             peakblock = peakblock,
                             snpeff_index = snpeff_index)
    collapse_snpeff_out <- .collapseSnpEff(snpeff_out = snpeff_out)
    out <- left_join(out, collapse_snpeff_out, by ="Gene_ID")
    out <- list(candidate_list = out, snpeff_out = snpeff_out)

  } else {
    out <- list(candidate_list = out, snpeff_out = NULL)
  }
  return(out)
}

#' @importFrom vcfR getFIX getINFO
.snpeff_index <- function(snpeff) {
  list(
    chr = read.gdsn(index.gdsn(snpeff$root, "chromosome")),
    pos = read.gdsn(index.gdsn(snpeff$root, "position")),
    at_ann = read.gdsn(index.gdsn(snpeff$root, "annotation/info/@ANN")),
    ann_node = index.gdsn(snpeff$root, "annotation/info/ANN")
  )
}

.addSnpEff <- function(snpeff, peakblock, snpeff_index = NULL){
  if (is.null(snpeff_index)) {
    snpeff_index <- .snpeff_index(snpeff)
  }
  snpeff_chr <- snpeff_index$chr
  snpeff_pos <- snpeff_index$pos
  at_ann <- snpeff_index$at_ann
  target_chr <- snpeff_chr %in% peakblock$Chr
  target_pos <-  snpeff_pos >= min(peakblock$Pos) & snpeff_pos <= max(peakblock$Pos)
  peak_ann <- target_chr & target_pos

  if(length(peak_ann) == 0){
    out <- data.frame(t(rep(NA, 19)))
    names(out) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                    "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                    "Rank", "HGVS.c", "HGVS.p", "Pos.in.tx", "Pos.in.CDS",
                    "Pos.in.AA", "Distance", "INFO", "Chr", "Pos", "negLog10P")

  } else {
    cum_at_ann <- cumsum(at_ann)
    peak_ann_cum_at_ann <- cum_at_ann[peak_ann]
    peak_ann_cum_at_ann_start <- peak_ann_cum_at_ann - at_ann[peak_ann] + 1
    target_index <- unlist(
      mapply(
        function(start_i, end_i) seq.int(start_i, end_i),
        peak_ann_cum_at_ann_start,
        peak_ann_cum_at_ann,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      ),
      use.names = FALSE
    )
    target_ann <- rep(FALSE, max(cum_at_ann))
    target_ann[target_index] <- TRUE
    out <- readex.gdsn(node = snpeff_index$ann_node, sel = list(target_ann))
    out[grepl("\\|$", out)] <- paste(out[grepl("\\|$", out)], "")
    out <- strsplit(out, "\\|")
    len <- sapply(out, length)
    out <- do.call("rbind", out)
    out <- as.data.frame(out)
    colnames(out) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                       "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                       "Rank", "HGVS.c", "HGVS.p", "Pos.in.tx", "Pos.in.CDS",
                       "Pos.in.AA", "Distance", "INFO")
    idx <- which(peak_ann)
    out$Chr <- rep(snpeff_chr[idx], at_ann[idx])
    out$Pos <- rep(snpeff_pos[idx], at_ann[idx])
    hit <- match(out$Pos, peakblock$Pos)
    out$negLog10P <- peakblock$negLog10P[hit]
    out <- out[!out[, 2] %in% c("custom", "intergenic_region"), ]
    if(nrow(out) == 0){
      out <- NULL
    }

    intergenic <- grep("&", out$Gene_ID)
    if(length(intergenic) > 0){
      intergenic_split <- lapply(intergenic, function(i){
        x <- unlist(strsplit(out$Gene_ID[i], "&"))
        i_out <- data.frame(split_id = x, out[i, ], row.names = NULL)
        return(i_out)
      })
      intergenic_split <- do.call("rbind", intergenic_split)
      intergenic_split$Gene_ID <- intergenic_split$split_id
      intergenic_split$split_id <- NULL
      out <- rbind(out[-intergenic, ], intergenic_split)
    }
  }
  return(out)
}

.collapseSnpEff <- function(snpeff_out){
  if(all(is.na(snpeff_out[1, ]))){
    out <- data.frame(t(rep(NA, 9)))
    names(out) <- c("Gene_ID", "HIGH", "MODERATE", "LOW", "MODIFIER",
                    "HIGH_at_var", "MODERATE_at_var", "LOW_at_var", "MODIFIER_at_var")

  } else {
    impacts <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
    out <- tapply(seq_along(snpeff_out$Annotation_Impact), snpeff_out$Gene_ID, function(i){
      data.frame(t(as.vector(table(factor(snpeff_out$Annotation_Impact[i], impacts)))),
                 t(as.vector(table(factor(snpeff_out$Annotation_Impact[i][!is.na(snpeff_out$negLog10P[i])], impacts)))))
    })
    Gene_ID = names(out)
    names(out) <- NULL
    out <- data.frame(Gene_ID, do.call("rbind", out))
    colnames(out) <- c("Gene_ID", "HIGH", "MODERATE", "LOW", "MODIFIER",
                       "HIGH_at_var", "MODERATE_at_var", "LOW_at_var", "MODIFIER_at_var")
  }
  return(out)
}

################################################################################
#' Generate interactive table for a candidate gene list
#'
#' @param object QTLscan object
#' @param pehno Phenotype names to be drawn
#' @param out_fn Prefix of output file
#'
#' @import reactable
#' @import htmltools
#' @export
#'
makeCanditeList <- function(object, pheno, out_fn, peak_id = NULL){
  candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno)
  table <- reactable(data = candidate, sortable = TRUE,
                     resizable = TRUE, filterable = TRUE,
                     searchable = TRUE, showPageSizeOptions = TRUE,
                     wrap = TRUE)
  save_html(tagList(table), file = out_fn)
}

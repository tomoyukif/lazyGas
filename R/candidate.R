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
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas",
                 new_node = "simple_candidate",
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
      dplyr::bind_rows(candidate_list)
    } else {
      NULL
    }
    snpeff_list <- if (length(snpeff_list) > 0L) {
      dplyr::bind_rows(snpeff_list)
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

    # E0: simple candidate list + Gene↔transcript/protein map (Parquet sidecar)
    simple <- .simple_candidate_from_wide(candidate_list, gff = gff)
    protein_map <- .gff_gene_protein_map(gff)
    .store_write_simple_candidate(object, pheno_name, simple)
    .store_write_gene_protein_map(object, protein_map)
  }
}

#' Minimal Gene_ID + position table for E (E0 §2.1)
#' @keywords internal
.simple_candidate_from_wide <- function(candidate_list, gff = NULL) {
  empty <- data.frame(
    peak_ID = character(),
    Gene_ID = character(),
    Gene_chr = character(),
    Gene_start = integer(),
    Gene_end = integer(),
    dist2peak = numeric(),
    negLog10P = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(candidate_list) || !nrow(candidate_list)) {
    return(empty)
  }
  need <- c("peak_ID", "Gene_ID", "Gene_chr", "Gene_start", "dist2peak", "negLog10P")
  for (cn in need) {
    if (!cn %in% names(candidate_list)) {
      candidate_list[[cn]] <- NA
    }
  }
  out <- candidate_list[, need, drop = FALSE]
  out$Gene_ID <- as.character(out$Gene_ID)
  out$Gene_chr <- as.character(out$Gene_chr)
  out$Gene_start <- as.integer(out$Gene_start)
  out$Gene_end <- as.integer(NA)
  if (!is.null(gff) && inherits(gff, "GRanges") && length(gff)) {
    md <- S4Vectors::mcols(gff)
    typ <- as.character(md$type %||% "")
    gene_i <- grepl("^gene$", typ, ignore.case = TRUE)
    if (any(gene_i)) {
      gg <- gff[gene_i]
      gid <- as.character(S4Vectors::mcols(gg)$ID %||% S4Vectors::mcols(gg)$gene_id)
      gend <- as.integer(BiocGenerics::end(gg))
      out$Gene_end <- gend[match(out$Gene_ID, gid)]
    }
  }
  # drop empty Gene_ID rows
  out <- out[!is.na(out$Gene_ID) & nzchar(out$Gene_ID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
.getCandidate <- function(peakblock, gff, snpeff, snpeff_index = NULL){
  peak_variant_id <- peakblock$peak_variant_ID[1]
  is_peak <- peakblock$variant_ID == peak_variant_id
  # One genomic window per chromosome present in the (possibly merged) block.
  chr_levels <- unique(as.character(peakblock$Chr))
  chr_levels <- chr_levels[!is.na(chr_levels) & nzchar(chr_levels)]
  if (length(chr_levels) == 0L) {
    chr_levels <- as.character(peakblock$Chr[is_peak][1])
  }
  ir_list <- vector("list", length(chr_levels))
  for (i in seq_along(chr_levels)) {
    ch <- chr_levels[[i]]
    sub <- peakblock[as.character(peakblock$Chr) == ch, , drop = FALSE]
    pos_num <- as.numeric(sub$Pos)
    if (nrow(sub) && !is.na(sub$chr_start[1]) && sub$chr_start[1] == 1) {
      peak_start <- 1
    } else {
      peak_start <- suppressWarnings(min(pos_num, na.rm = TRUE))
    }
    if (nrow(sub) && !is.na(sub$chr_end[1]) && sub$chr_end[1] == 1) {
      peak_end <- 2^30
    } else {
      peak_end <- suppressWarnings(max(pos_num, na.rm = TRUE))
    }
    if (!is.finite(peak_start) || !is.finite(peak_end) || peak_end < peak_start) {
      lead_pos <- as.numeric(peakblock$Pos[is_peak][1])
      peak_start <- if (is.finite(lead_pos)) lead_pos else 1
      peak_end <- peak_start
    }
    ir_list[[i]] <- IRanges::IRanges(
      start = as.integer(peak_start),
      end = as.integer(peak_end)
    )
  }
  peak_gff <- GenomicRanges::GRanges(
    seqnames = chr_levels,
    ranges = do.call(c, ir_list)
  )
  type_hit <- gff$type %in% "gene"
  type_hit[is.na(type_hit)] <- FALSE
  gene_gff <- gff[type_hit]
  hit <- gene_gff[queryHits(findOverlaps(gene_gff, peak_gff))]
  hit <- hit[order(start(hit))]
  hit$dist2peak <- start(hit) - as.numeric(peakblock$Pos[is_peak][1])
  nearest_p <- vapply(seq_along(hit), function(i) {
    x <- BiocGenerics::start(hit)[i]
    ch <- as.character(GenomeInfoDb::seqnames(hit)[i])
    same <- as.character(peakblock$Chr) == ch
    if (!any(same)) {
      return(NA_real_)
    }
    vals <- as.numeric(peakblock$negLog10P[same])
    pos_s <- as.numeric(peakblock$Pos[same])
    vals[which.min(abs(pos_s - x))]
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
  empty <- data.frame(
    Allele = NA_character_,
    Annotation = NA_character_,
    Annotation_Impact = NA_character_,
    Gene_Name = NA_character_,
    Gene_ID = NA_character_,
    Feature_Type = NA_character_,
    Feature_ID = NA_character_,
    Transcript_BioType = NA_character_,
    Rank = NA_character_,
    `HGVS.c` = NA_character_,
    `HGVS.p` = NA_character_,
    Pos.in.tx = NA_character_,
    Pos.in.CDS = NA_character_,
    Pos.in.AA = NA_character_,
    Distance = NA_character_,
    INFO = NA_character_,
    Chr = NA_character_,
    Pos = NA_real_,
    negLog10P = NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (is.null(snpeff_index)) {
    snpeff_index <- .snpeff_index(snpeff)
  }
  snpeff_chr <- snpeff_index$chr
  snpeff_pos <- snpeff_index$pos
  at_ann <- snpeff_index$at_ann
  target_chr <- snpeff_chr %in% unique(as.character(peakblock$Chr))
  if (!any(target_chr)) {
    return(empty)
  }
  pos_lo <- suppressWarnings(min(as.numeric(peakblock$Pos), na.rm = TRUE))
  pos_hi <- suppressWarnings(max(as.numeric(peakblock$Pos), na.rm = TRUE))
  if (!is.finite(pos_lo) || !is.finite(pos_hi)) {
    return(empty)
  }
  target_pos <- snpeff_pos >= pos_lo & snpeff_pos <= pos_hi
  peak_ann <- target_chr & target_pos
  if (!any(peak_ann)) {
    return(empty)
  }

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
  if (!length(target_index)) {
    return(empty)
  }
  target_ann <- rep(FALSE, max(cum_at_ann))
  target_ann[target_index] <- TRUE
  out <- readex.gdsn(node = snpeff_index$ann_node, sel = list(target_ann))
  if (!length(out)) {
    return(empty)
  }
  out[grepl("\\|$", out)] <- paste(out[grepl("\\|$", out)], "")
  out <- strsplit(out, "\\|")
  # Pad ANN fields to a common width before binding rows.
  n_fields <- 16L
  out <- lapply(out, function(x) {
    if (length(x) < n_fields) {
      c(x, rep("", n_fields - length(x)))
    } else if (length(x) > n_fields) {
      x[seq_len(n_fields)]
    } else {
      x
    }
  })
  out <- do.call("rbind", out)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
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
  if (nrow(out) == 0L) {
    return(empty)
  }

  intergenic <- grep("&", out$Gene_ID)
  if (length(intergenic) > 0L) {
    intergenic_split <- lapply(intergenic, function(i) {
      x <- unlist(strsplit(out$Gene_ID[i], "&"))
      base <- out[i, , drop = FALSE]
      base$Gene_ID <- NULL
      data.frame(Gene_ID = x, base, row.names = NULL, stringsAsFactors = FALSE)
    })
    out <- dplyr::bind_rows(out[-intergenic, , drop = FALSE], dplyr::bind_rows(intergenic_split))
  }
  out
}

#' Keep worst SnpEff ANN per gene × variant site
#'
#' Site key is Gene_ID + Chr + Pos + Allele (Allele omitted when absent).
#' Rows without a usable Pos are not merged across rows. Impact order:
#' HIGH > MODERATE > LOW > MODIFIER; ties prefer a non-empty HGVS.p.
#'
#' @keywords internal
.snpeff_worst_ann_per_site <- function(snpeff) {
  if (is.null(snpeff) || !is.data.frame(snpeff) || !nrow(snpeff)) {
    return(snpeff)
  }
  n <- nrow(snpeff)
  gene <- if ("Gene_ID" %in% names(snpeff)) {
    as.character(snpeff$Gene_ID)
  } else {
    rep("", n)
  }
  gene[is.na(gene)] <- ""

  impact_col <- if ("Annotation_Impact" %in% names(snpeff)) {
    "Annotation_Impact"
  } else if ("Impact" %in% names(snpeff)) {
    "Impact"
  } else {
    NULL
  }
  rank <- rep(0L, n)
  if (!is.null(impact_col)) {
    lab <- toupper(as.character(snpeff[[impact_col]]))
    rank <- match(lab, c("MODIFIER", "LOW", "MODERATE", "HIGH"), nomatch = 0L)
  }

  if (!("Pos" %in% names(snpeff))) {
    site <- paste0("row", seq_len(n))
  } else {
    pos <- suppressWarnings(as.numeric(snpeff$Pos))
    chr <- if ("Chr" %in% names(snpeff)) as.character(snpeff$Chr) else rep("", n)
    chr[is.na(chr)] <- ""
    allele <- if ("Allele" %in% names(snpeff)) {
      as.character(snpeff$Allele)
    } else {
      rep("", n)
    }
    allele[is.na(allele)] <- ""
    ok <- is.finite(pos)
    site <- ifelse(
      ok,
      paste(chr, pos, allele, sep = "\t"),
      paste0("row", seq_len(n))
    )
  }
  key <- paste(gene, site, sep = "\t")

  hgvs <- rep("", n)
  if ("HGVS.p" %in% names(snpeff)) {
    hgvs <- as.character(snpeff[["HGVS.p"]])
  } else if ("HGVS_P" %in% names(snpeff)) {
    hgvs <- as.character(snpeff$HGVS_P)
  }
  hgvs[is.na(hgvs)] <- ""
  hgvs_score <- as.integer(nzchar(trimws(hgvs)))

  ord <- order(key, -rank, -hgvs_score, seq_len(n))
  keep_idx <- sort(ord[!duplicated(key[ord])])
  snpeff[keep_idx, , drop = FALSE]
}

.collapseSnpEff <- function(snpeff_out){
  empty_collapse <- data.frame(
    Gene_ID = NA_character_,
    HIGH = NA_real_,
    MODERATE = NA_real_,
    LOW = NA_real_,
    MODIFIER = NA_real_,
    HIGH_at_var = NA_real_,
    MODERATE_at_var = NA_real_,
    LOW_at_var = NA_real_,
    MODIFIER_at_var = NA_real_,
    stringsAsFactors = FALSE
  )
  if (is.null(snpeff_out) || !is.data.frame(snpeff_out) || !nrow(snpeff_out)) {
    return(empty_collapse)
  }
  if (all(is.na(snpeff_out[1, ]))) {
    return(empty_collapse)
  }
  impacts <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  gene_ids <- unique(as.character(snpeff_out$Gene_ID))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  if (!length(gene_ids)) {
    return(empty_collapse)
  }
  rows <- lapply(gene_ids, function(gid) {
    i <- which(as.character(snpeff_out$Gene_ID) == gid)
    np <- suppressWarnings(as.numeric(as.character(snpeff_out$negLog10P[i])))
    linked <- is.finite(np)
    sub <- snpeff_out[i[linked], , drop = FALSE]
    sub <- .snpeff_worst_ann_per_site(sub)
    if (!nrow(sub)) {
      counts <- rep(0L, 4L)
    } else {
      imp <- factor(sub$Annotation_Impact, levels = impacts)
      counts <- as.integer(table(imp))
    }
    data.frame(
      Gene_ID = gid,
      HIGH = counts[[1L]],
      MODERATE = counts[[2L]],
      LOW = counts[[3L]],
      MODIFIER = counts[[4L]],
      HIGH_at_var = counts[[1L]],
      MODERATE_at_var = counts[[2L]],
      LOW_at_var = counts[[3L]],
      MODIFIER_at_var = counts[[4L]],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Peak-block markers for one peak (recalc, else peakcall)
#' @keywords internal
.peak_block_markers <- function(object, pheno_name, peak_id) {
  empty <- data.frame(
    Chr = character(),
    Pos = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(object) || is.null(pheno_name) || is.null(peak_id)) {
    return(empty)
  }
  pk <- tryCatch(
    lazyData(object = object, dataset = "recalc", pheno = pheno_name),
    error = function(e) NULL
  )
  if (is.null(pk) || !nrow(pk)) {
    pk <- tryCatch(
      lazyData(object = object, dataset = "peakcall", pheno = pheno_name),
      error = function(e) NULL
    )
  }
  if (is.null(pk) || !nrow(pk) || !"peak_ID" %in% names(pk)) {
    return(empty)
  }
  b <- pk[as.character(pk$peak_ID) == as.character(peak_id)[1L], , drop = FALSE]
  if (!nrow(b) || !"Pos" %in% names(b)) {
    return(empty)
  }
  chr <- if ("Chr" %in% names(b)) as.character(b$Chr) else NA_character_
  data.frame(
    Chr = chr,
    Pos = as.numeric(b$Pos),
    stringsAsFactors = FALSE
  )
}

#' SnpEff rows for one peak: store load + peak-block-marker link only
#'
#' Keeps annotations whose Chr/Pos match peak-block markers and that carry a
#' finite \code{negLog10P} (the listCandidate link to GWAS peak markers).
#' Falls back to unique Chr/Pos matches if all \code{negLog10P} are blank
#' (legacy store rows from other peaks).
#'
#' @keywords internal
.snpeff_for_peak <- function(snpeff = NULL,
                             object = NULL,
                             pheno_name = NULL,
                             peak_id = NULL,
                             peak_markers = NULL) {
  empty <- if (!is.null(snpeff) && is.data.frame(snpeff)) {
    snpeff[0, , drop = FALSE]
  } else {
    data.frame()
  }
  if (is.null(snpeff) || !is.data.frame(snpeff) || !nrow(snpeff)) {
    if (!is.null(object) && !is.null(pheno_name)) {
      snpeff <- tryCatch(
        lazyData(object = object, dataset = "snpeff", pheno = pheno_name),
        error = function(e) NULL
      )
    }
  }
  if (is.null(snpeff) || !is.data.frame(snpeff) || !nrow(snpeff)) {
    return(empty)
  }
  if (is.null(peak_markers)) {
    peak_markers <- .peak_block_markers(object, pheno_name, peak_id)
  }
  if (is.null(peak_markers) || !nrow(peak_markers) || !"Pos" %in% names(peak_markers)) {
    return(snpeff[0, , drop = FALSE])
  }
  pos <- unique(as.numeric(peak_markers$Pos))
  pos <- pos[is.finite(pos)]
  if (!length(pos) || !"Pos" %in% names(snpeff)) {
    return(snpeff[0, , drop = FALSE])
  }
  sp <- as.numeric(snpeff$Pos)
  keep <- is.finite(sp) & sp %in% pos
  if ("Chr" %in% names(snpeff) && "Chr" %in% names(peak_markers)) {
    chrs <- unique(as.character(peak_markers$Chr))
    chrs <- chrs[!is.na(chrs) & nzchar(chrs)]
    if (length(chrs)) {
      keep <- keep & as.character(snpeff$Chr) %in% chrs
    }
  }
  out <- snpeff[keep, , drop = FALSE]
  if (!nrow(out)) {
    return(out)
  }
  if ("negLog10P" %in% names(out)) {
    np <- suppressWarnings(as.numeric(as.character(out$negLog10P)))
    linked <- is.finite(np)
    if (any(linked)) {
      out <- out[linked, , drop = FALSE]
    } else {
      # Drop blank-negLog10P duplicates from other peaks' rbind
      key_cols <- intersect(
        c("Gene_ID", "Chr", "Pos", "Allele", "Annotation", "Feature_ID", "HGVS.p", "HGVS_P"),
        names(out)
      )
      if (length(key_cols)) {
        out <- out[!duplicated(out[, key_cols, drop = FALSE]), , drop = FALSE]
      }
    }
  }
  rownames(out) <- NULL
  out
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

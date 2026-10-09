################################################################################
# GFF gene windows and gene↔protein maps (shared by candidate listing / finemap)

.gff_gene_windows <- function(gff,
                              upstream_bp = 3000L,
                              downstream_bp = 500L) {
  empty <- data.frame(
    Gene_ID = character(),
    Chr = character(),
    window_start = integer(),
    window_end = integer(),
    mid = numeric(),
    strand = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(gff) || !inherits(gff, "GRanges") || !length(gff)) {
    return(empty)
  }
  md <- S4Vectors::mcols(gff)
  typ <- as.character(md$type %||% md$Type %||% "")
  cds_i <- grepl("CDS", typ, ignore.case = TRUE)
  if (!any(cds_i)) {
    # fall back to gene features
    gene_i <- grepl("^gene$", typ, ignore.case = TRUE)
    if (!any(gene_i)) {
      return(empty)
    }
    gg <- gff[gene_i]
    gid <- as.character(
      S4Vectors::mcols(gg)$ID %||% S4Vectors::mcols(gg)$gene_id %||% ""
    )
    st <- as.character(BiocGenerics::strand(gg))
    start <- as.integer(BiocGenerics::start(gg))
    end <- as.integer(BiocGenerics::end(gg))
    chr <- as.character(GenomeInfoDb::seqnames(gg))
    w_start <- ifelse(st == "-", start - downstream_bp, start - upstream_bp)
    w_end <- ifelse(st == "-", end + upstream_bp, end + downstream_bp)
    return(data.frame(
      Gene_ID = gid,
      Chr = chr,
      window_start = pmax(1L, as.integer(w_start)),
      window_end = as.integer(w_end),
      mid = (as.numeric(w_start) + as.numeric(w_end)) / 2,
      strand = st,
      stringsAsFactors = FALSE
    ))
  }
  cds <- gff[cds_i]
  cds_parent <- as.character(S4Vectors::mcols(cds)$Parent %||% "")
  # Parent may be transcript; map transcript → gene via Parent of mRNA
  mrna_i <- grepl("mRNA|transcript", typ, ignore.case = TRUE)
  tx_to_gene <- character()
  if (any(mrna_i)) {
    mr <- gff[mrna_i]
    tx_id <- as.character(S4Vectors::mcols(mr)$ID %||% "")
    tx_parent <- trimws(sub(",.*$", "", as.character(S4Vectors::mcols(mr)$Parent %||% "")))
    keep <- nzchar(tx_id)
    tx_to_gene <- setNames(tx_parent[keep], tx_id[keep])
  }
  # Vectorized first-Parent → gene (avoid per-CDS vapply/strsplit)
  first_parent <- trimws(sub(",.*$", "", cds_parent))
  mapped <- unname(tx_to_gene[first_parent])
  gene_for_cds <- ifelse(
    !is.na(mapped) & nzchar(mapped),
    mapped,
    first_parent
  )
  gene_for_cds[!nzchar(gene_for_cds)] <- NA_character_

  start <- as.integer(BiocGenerics::start(cds))
  end <- as.integer(BiocGenerics::end(cds))
  chr <- as.character(GenomeInfoDb::seqnames(cds))
  strand <- as.character(BiocGenerics::strand(cds))
  ok <- !is.na(gene_for_cds) & nzchar(gene_for_cds)
  if (!any(ok)) {
    return(empty)
  }
  gene_for_cds <- gene_for_cds[ok]
  start <- start[ok]
  end <- end[ok]
  chr <- chr[ok]
  strand <- strand[ok]

  # Aggregate per gene without split()/lapply()
  cds_start <- tapply(start, gene_for_cds, min, na.rm = TRUE)
  cds_end <- tapply(end, gene_for_cds, max, na.rm = TRUE)
  # strand / Chr: first occurrence in original order
  gid_levels <- names(cds_start)
  idx_first <- match(gid_levels, gene_for_cds)
  st <- strand[idx_first]
  chr_out <- chr[idx_first]
  cds_start <- as.integer(cds_start[gid_levels])
  cds_end <- as.integer(cds_end[gid_levels])

  w_start <- ifelse(st == "-", cds_start - downstream_bp, cds_start - upstream_bp)
  w_end <- ifelse(st == "-", cds_end + upstream_bp, cds_end + downstream_bp)
  data.frame(
    Gene_ID = gid_levels,
    Chr = chr_out,
    window_start = pmax(1L, as.integer(w_start)),
    window_end = as.integer(w_end),
    mid = (as.numeric(w_start) + as.numeric(w_end)) / 2,
    strand = st,
    stringsAsFactors = FALSE
  )
}

#' Gene_ID ↔ transcript / protein IDs from GFF attributes
#' @keywords internal
.gff_gene_protein_map <- function(gff) {
  empty <- data.frame(
    Gene_ID = character(),
    transcript_id = character(),
    protein_id = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(gff) || !inherits(gff, "GRanges") || !length(gff)) {
    return(empty)
  }
  md <- S4Vectors::mcols(gff)
  typ <- as.character(md$type %||% "")
  id <- as.character(md$ID %||% "")
  parent <- as.character(md$Parent %||% "")
  protein_id <- as.character(md$protein_id %||% md$proteinId %||% "")
  first_parent <- trimws(sub(",.*$", "", parent))

  rows <- list()
  # transcripts (vectorized)
  mrna_i <- grepl("mRNA|transcript", typ, ignore.case = TRUE)
  if (any(mrna_i)) {
    gid <- first_parent[mrna_i]
    tid <- id[mrna_i]
    pid <- protein_id[mrna_i]
    pid <- ifelse(nzchar(pid), pid, tid)
    keep <- nzchar(gid) & nzchar(tid)
    if (any(keep)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Gene_ID = gid[keep],
        transcript_id = tid[keep],
        protein_id = pid[keep],
        stringsAsFactors = FALSE
      )
    }
  }
  # CDS with protein_id (vectorized)
  cds_i <- grepl("CDS", typ, ignore.case = TRUE) & nzchar(protein_id)
  if (any(cds_i)) {
    tx_to_gene <- character()
    if (any(mrna_i)) {
      tx_to_gene <- setNames(first_parent[mrna_i], id[mrna_i])
    }
    p <- first_parent[cds_i]
    mapped <- unname(tx_to_gene[p])
    gid <- ifelse(!is.na(mapped) & nzchar(mapped), mapped, p)
    keep <- nzchar(gid) & nzchar(protein_id[cds_i])
    if (any(keep)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Gene_ID = gid[keep],
        transcript_id = p[keep],
        protein_id = protein_id[cds_i][keep],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

################################################################################
#' Cluster peaks across traits by distance and/or peak-marker genotype correlation
#'
#' @param object A \code{LazyGas} object.
#' @param recalc Use recalculated peaks when available.
#' @param dist_threshold Maximum physical distance (bp) on the same chromosome.
#' @param r2_threshold Minimum r-squared between peak representative markers.
#' @param store If \code{TRUE}, write cluster assignments to the companion store.
#' @export
clusterCrossTraitPeaks <- function(object,
                                   recalc = TRUE,
                                   dist_threshold = 500000L,
                                   r2_threshold = 0.6,
                                   store = TRUE) {
  peaks <- .multitrait_collect_peaks(object = object, recalc = recalc)
  if (nrow(peaks) == 0L) {
    return(data.frame())
  }
  clusters <- .multitrait_cluster_peaks(
    object = object,
    peaks = peaks,
    dist_threshold = dist_threshold,
    r2_threshold = r2_threshold
  )
  if (store) {
    .store_write_section_df(object, "multitrait", "clusters", clusters, pheno_name = NULL)
  }
  clusters
}

#' Summarize cross-trait peak clusters
#'
#' @param object A \code{LazyGas} object, or a cluster table from
#'   [clusterCrossTraitPeaks()].
#' @param recalc Use recalculated peaks when collecting peaks from \code{object}.
#' @export
summarizeCrossTraitPeaks <- function(object, recalc = TRUE) {
  if (inherits(object, "LazyGas")) {
    clusters <- clusterCrossTraitPeaks(object = object, recalc = recalc, store = FALSE)
  } else {
    clusters <- object
  }
  if (nrow(clusters) == 0L) {
    return(data.frame())
  }
  split_clusters <- split(clusters, clusters$cluster_id)
  out <- lapply(names(split_clusters), function(cid) {
    cl <- split_clusters[[cid]]
    traits <- unique(cl$trait)
    data.frame(
      cluster_id = as.integer(cid),
      traits = paste(sort(traits), collapse = ","),
      n_traits = length(traits),
      chr = cl$chr[1],
      cluster_start = min(cl$pos, na.rm = TRUE),
      cluster_end = max(cl$pos, na.rm = TRUE),
      lead_variant_ID = cl$peak_variant_ID[which.max(cl$negLog10P)],
      max_negLog10P = max(cl$negLog10P, na.rm = TRUE),
      merge_reason = paste(unique(cl$merge_reason), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, out)
  rownames(summary_df) <- NULL
  if (inherits(object, "LazyGas")) {
    .store_write_section_df(object, "multitrait", "summary", summary_df, pheno_name = NULL)
  }
  summary_df
}

#' Plot multi-trait peak overview
#'
#' @param object A \code{LazyGas} object.
#' @param recalc Use recalculated peaks when clustering.
#' @param dist_threshold,r2_threshold Passed to [clusterCrossTraitPeaks()].
#'
#' @return A list with \code{heatmap} (genome-wide trait-by-chromosome tile plot),
#'   \code{clusters}, \code{summary}, and \code{shared_genes}. The heatmap shows
#'   all chromosomes in the marker set; colored tiles mark chromosomes with a
#'   peak for that trait (fill = \code{cluster_id}). \code{shared_genes} lists
#'   candidate genes that appear in two or more traits within the same cluster
#'   window; it is empty when [listCandidate()] has not been run or no gene is
#'   shared across traits in that region.
#' @export
plotMultiTraitOverview <- function(object,
                                   recalc = TRUE,
                                   dist_threshold = 500000L,
                                   r2_threshold = 0.6) {
  clusters <- clusterCrossTraitPeaks(
    object = object,
    recalc = recalc,
    dist_threshold = dist_threshold,
    r2_threshold = r2_threshold,
    store = TRUE
  )
  summary_df <- summarizeCrossTraitPeaks(clusters)
  if (nrow(clusters) == 0L) {
    stop("No peaks found across traits.", call. = FALSE)
  }

  heat_df <- .multitrait_heatmap_df(object = object, clusters = clusters)
  cluster_levels <- sort(unique(clusters$cluster_id))
  cluster_cols <- c(
    "#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B",
    "#EECA3B", "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC"
  )
  cluster_cols <- stats::setNames(
    cluster_cols[seq_along(cluster_levels)],
    as.character(cluster_levels)
  )

  p_heat <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = .data$chr, y = .data$trait, fill = .data$cluster_id)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(
      values = cluster_cols,
      na.value = "gray97",
      name = "Cluster",
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Cross-trait peak presence",
      subtitle = paste0(
        length(cluster_levels), " cluster(s); colored tiles = peak on chromosome"
      ),
      x = "Chromosome",
      y = "Trait"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  shared_genes <- .multitrait_shared_genes(object = object, clusters = clusters, recalc = recalc)
  if (nrow(shared_genes) > 0L) {
  .store_write_section_df(object, "multitrait", "shared_genes", shared_genes, pheno_name = NULL)
  }

  list(
    heatmap = p_heat,
    clusters = clusters,
    summary = summary_df,
    shared_genes = shared_genes
  )
}

.multitrait_empty_shared_genes <- function() {
  data.frame(
    cluster_id = integer(),
    Gene_ID = character(),
    trait = character(),
    n_traits = integer(),
    stringsAsFactors = FALSE
  )
}

.multitrait_heatmap_df <- function(object, clusters) {
  chr_lev <- as.character(.get_chromosome_data(object = object)$chr_lev)
  traits <- sort(unique(clusters$trait))
  presence <- stats::aggregate(
    cluster_id ~ trait + chr,
    data = clusters,
    FUN = function(x) {
      ux <- unique(x)
      if (length(ux) == 1L) {
        ux[[1L]]
      } else {
        paste(ux, collapse = ",")
      }
    },
    simplify = TRUE
  )
  presence$chr <- as.character(presence$chr)
  grid <- expand.grid(
    trait = traits,
    chr = chr_lev,
    stringsAsFactors = FALSE
  )
  merged <- merge(grid, presence, by = c("trait", "chr"), all.x = TRUE)
  merged$chr <- factor(merged$chr, levels = chr_lev)
  merged$trait <- factor(merged$trait, levels = traits)
  merged$cluster_id <- factor(as.character(merged$cluster_id), levels = as.character(sort(unique(clusters$cluster_id))))
  merged
}

.multitrait_collect_peaks <- function(object, recalc = TRUE) {
  pheno_names <- getPheno(object)$pheno_names
  rows <- lapply(pheno_names, function(pn) {
    use_recalc <- recalc && .store_section_exists(object, "recalc")
    pk <- .get_peakcall(object = object, pheno_name = pn, recalc = use_recalc)
    if (is.null(pk) || nrow(pk) == 0L) {
      return(NULL)
    }
    lead <- pk[pk$peak_variant_ID == pk$variant_ID, , drop = FALSE]
    if (nrow(lead) == 0L) {
      lead <- unique(pk[, c("peak_ID", "peak_variant_ID"), drop = FALSE])
      lead <- lead[order(lead$peak_ID), , drop = FALSE]
      pos_col <- if ("peak_Pos" %in% names(pk)) "peak_Pos" else "Pos"
      chr_col <- if ("peak_Chr" %in% names(pk)) "peak_Chr" else "Chr"
      p_col <- if ("peak_negLog10P" %in% names(pk)) "peak_negLog10P" else "negLog10P"
      hit <- match(lead$peak_variant_ID, pk$variant_ID)
      data.frame(
        trait = pn,
        peak_ID = lead$peak_ID,
        peak_variant_ID = lead$peak_variant_ID,
        chr = pk[[chr_col]][hit],
        pos = pk[[pos_col]][hit],
        negLog10P = pk[[p_col]][hit],
        stringsAsFactors = FALSE
      )
    } else {
      pos_col <- if ("peak_Pos" %in% names(lead)) "peak_Pos" else "Pos"
      chr_col <- if ("peak_Chr" %in% names(lead)) "peak_Chr" else "Chr"
      p_col <- if ("peak_negLog10P" %in% names(lead)) "peak_negLog10P" else "negLog10P"
      data.frame(
        trait = pn,
        peak_ID = lead$peak_ID,
        peak_variant_ID = lead$peak_variant_ID,
        chr = lead[[chr_col]],
        pos = lead[[pos_col]],
        negLog10P = lead[[p_col]],
        stringsAsFactors = FALSE
      )
    }
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame())
  }
  rownames(out) <- NULL
  out
}

.multitrait_cluster_peaks <- function(object, peaks, dist_threshold, r2_threshold) {
  n <- nrow(peaks)
  parent <- seq_len(n)

  find_root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }

  union_nodes <- function(i, j) {
    ri <- find_root(i)
    rj <- find_root(j)
    if (ri != rj) {
      parent[rj] <<- ri
    }
  }

  reasons <- matrix(NA_character_, n, n)
  diag(reasons) <- "self"

  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        if (as.character(peaks$chr[i]) != as.character(peaks$chr[j])) {
          next
        }
        dist_ok <- abs(peaks$pos[i] - peaks$pos[j]) <= dist_threshold
        r2 <- .multitrait_peak_r2(
          object = object,
          variant_i = peaks$peak_variant_ID[i],
          variant_j = peaks$peak_variant_ID[j],
          chr = peaks$chr[i]
        )
        r2_ok <- is.finite(r2) && r2 >= r2_threshold
        if (dist_ok || r2_ok) {
          union_nodes(i, j)
          reasons[i, j] <- reasons[j, i] <- if (dist_ok && r2_ok) {
            "both"
          } else if (dist_ok) {
            "distance"
          } else {
            "correlation"
          }
        }
      }
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1L))
  cluster_id <- match(roots, unique(roots))
  merge_reason <- vapply(seq_len(n), function(i) {
    rs <- reasons[i, roots == roots[i]]
    rs <- rs[!is.na(rs) & rs != "self"]
    if (length(rs) == 0L) "singleton" else paste(unique(rs), collapse = ";")
  }, character(1L))

  cbind(
    data.frame(cluster_id = cluster_id, merge_reason = merge_reason),
    peaks
  )[, c(
    "cluster_id", "merge_reason", "trait", "peak_ID", "peak_variant_ID",
    "chr", "pos", "negLog10P"
  ), drop = FALSE]
}

.multitrait_peak_r2 <- function(object, variant_i, variant_j, chr) {
  if (identical(as.character(variant_i), as.character(variant_j))) {
    return(1)
  }
  geno_format <- .store_read_scan_scalar(object, "geno_format")
  if (is.null(geno_format)) {
    geno_format <- "dosage"
  }
  selection <- .set_selection_criteria_for_peakcall(
    object = object,
    geno_format = geno_format,
    chr = chr
  )
  geno <- as.matrix(.retrieve_geno(
    object = object,
    selection = selection,
    geno_format = geno_format
  ))
  variables <- .multitrait_marker_table(object = object)
  chr_ids <- variables$snp_id[variables$chr == chr]
  idx_i <- match(variant_i, chr_ids)
  idx_j <- match(variant_j, chr_ids)
  if (is.na(idx_i) || is.na(idx_j)) {
    return(NA_real_)
  }
  if (selection$is_categorical) {
    gi <- as.integer(geno[, idx_i])
    gj <- as.integer(geno[, idx_j])
    mean(gi == gj, na.rm = TRUE)
  } else {
    cor(geno[, idx_i], geno[, idx_j], use = "pairwise.complete.obs")^2
  }
}

.multitrait_shared_genes <- function(object, clusters, recalc) {
  pheno_names <- unique(clusters$trait)
  cand_list <- lapply(pheno_names, function(pn) {
    cand <- lazyData(object = object, dataset = "candidate", pheno = pn)
    if (is.null(cand) || nrow(cand) == 0L || !"Gene_ID" %in% names(cand)) {
      return(NULL)
    }
    cand$trait <- pn
    cand
  })
  cand_all <- do.call(rbind, cand_list)
  if (is.null(cand_all) || nrow(cand_all) == 0L) {
    return(.multitrait_empty_shared_genes())
  }

  multi_clusters <- unique(clusters$cluster_id[duplicated(clusters$cluster_id) |
                                                 duplicated(clusters$cluster_id, fromLast = TRUE)])
  if (length(multi_clusters) == 0L) {
    return(.multitrait_empty_shared_genes())
  }

  rows <- lapply(multi_clusters, function(cid) {
    cl <- clusters[clusters$cluster_id == cid, , drop = FALSE]
    chr <- cl$chr[1]
    pad <- 500000L
    start <- min(cl$pos, na.rm = TRUE) - pad
    end <- max(cl$pos, na.rm = TRUE) + pad
    traits <- unique(cl$trait)
    sub <- cand_all[cand_all$trait %in% traits, , drop = FALSE]
    if (!"Gene_chr" %in% names(sub)) {
      return(NULL)
    }
    hit <- as.character(sub$Gene_chr) == as.character(chr) &
      sub$Gene_start <= end &
      sub$Gene_end >= start
    sub <- sub[hit, , drop = FALSE]
    if (nrow(sub) == 0L) {
      return(NULL)
    }
    gene_tab <- stats::aggregate(
      trait ~ Gene_ID,
      data = sub,
      FUN = function(x) paste(sort(unique(x)), collapse = ",")
    )
    gene_tab$cluster_id <- cid
    gene_tab$n_traits <- vapply(
      strsplit(gene_tab$trait, ","),
      length,
      integer(1L)
    )
    gene_tab[gene_tab$n_traits >= 2, c("cluster_id", "Gene_ID", "trait", "n_traits")]
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(.multitrait_empty_shared_genes())
  }
  rownames(out) <- NULL
  out
}

.multitrait_marker_table <- function(object) {
  data.frame(
    snp_id = getMarID(object = object),
    chr = getChromosome(object = object),
    pos = getPosition(object = object),
    stringsAsFactors = FALSE
  )
}

################################################################################
#' Compare cross-trait peaks by credible-set overlap
#'
#' For multi-trait peak clusters from [clusterCrossTraitPeaks()], compute
#' Wakefield ABF credible sets per peak (via [calcCredibleSet()]) and label
#' whether pairs of peaks from different traits are compatible with a shared
#' signal. This is an **optional** follow-up to the main pipeline — it is not a
#' default [runLazyGas()] step.
#'
#' Labels:
#' \itemize{
#'   \item \code{shared-compatible}: overlap \eqn{\ge} 1, lead of either peak in
#'     the other CS, both resolutions not \code{low}
#'   \item \code{possible-shared}: overlap but weaker evidence
#'   \item \code{independent-nearby}: same cluster, no CS overlap, both resolved
#'   \item \code{inconclusive}: both \code{low}, or overlap driven by large CS
#'     fractions, or unresolved non-overlap
#' }
#'
#' @param object A \code{LazyGas} object.
#' @param recalc Use recalculated peaks / blocks.
#' @param clusters Optional cluster table from [clusterCrossTraitPeaks()]. If
#'   \code{NULL}, clusters are computed with \code{dist_threshold} /
#'   \code{r2_threshold}.
#' @param dist_threshold,r2_threshold Passed to [clusterCrossTraitPeaks()] when
#'   \code{clusters} is \code{NULL}.
#' @param coverage,prior_W,prior_V,use_conditional Passed to [calcCredibleSet()].
#' @param large_cs_fraction Both sides at or above this locus fraction with
#'   overlap → \code{inconclusive} when either resolution is \code{low}
#'   (default 0.8).
#' @param store If \code{TRUE}, write \code{multitrait/cs_overlap} to the
#'   companion store.
#'
#' @return A list with \code{pairs} (data.frame) and \code{clusters} (the
#'   cluster table used). Retrieve stored pairs with
#'   \code{lazyData(object, dataset = "multitrait", kind = "cs_overlap")}.
#' @seealso [clusterCrossTraitPeaks()], [calcCredibleSet()], [plotMultiTraitOverview()]
#' @export
compareCrossTraitCredibleSets <- function(object,
                                          recalc = TRUE,
                                          clusters = NULL,
                                          dist_threshold = 500000L,
                                          r2_threshold = 0.6,
                                          coverage = 0.95,
                                          prior_W = 0.15,
                                          prior_V = 0.04,
                                          use_conditional = TRUE,
                                          large_cs_fraction = 0.8,
                                          store = TRUE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  if (is.null(clusters)) {
    clusters <- clusterCrossTraitPeaks(
      object = object,
      recalc = recalc,
      dist_threshold = dist_threshold,
      r2_threshold = r2_threshold,
      store = store
    )
  }
  if (nrow(clusters) == 0L) {
    pairs <- .multitrait_empty_cs_overlap()
    if (store) {
      .store_write_section_df(object, "multitrait", "cs_overlap", pairs, pheno_name = NULL)
    }
    return(list(pairs = pairs, clusters = clusters))
  }

  multi_ids <- unique(clusters$cluster_id[
    duplicated(clusters$cluster_id) | duplicated(clusters$cluster_id, fromLast = TRUE)
  ])
  if (length(multi_ids) == 0L) {
    pairs <- .multitrait_empty_cs_overlap()
    if (store) {
      .store_write_section_df(object, "multitrait", "cs_overlap", pairs, pheno_name = NULL)
    }
    return(list(pairs = pairs, clusters = clusters))
  }

  cs_cache <- new.env(parent = emptyenv())
  get_cs <- function(trait, peak_id) {
    key <- paste(trait, peak_id, sep = "\r")
    if (exists(key, envir = cs_cache, inherits = FALSE)) {
      return(get(key, envir = cs_cache, inherits = FALSE))
    }
    cred <- tryCatch(
      calcCredibleSet(
        object = object,
        pheno = trait,
        peak_id = peak_id,
        coverage = coverage,
        prior_W = prior_W,
        prior_V = prior_V,
        use_conditional = use_conditional,
        recalc = recalc,
        store = store
      ),
      error = function(e) {
        warning(
          "calcCredibleSet failed for trait=", trait, " peak_id=", peak_id,
          ": ", conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    info <- .multitrait_cs_info(cred)
    assign(key, info, envir = cs_cache)
    info
  }

  pair_rows <- list()
  for (cid in multi_ids) {
    cl <- clusters[clusters$cluster_id == cid, , drop = FALSE]
    n <- nrow(cl)
    if (n < 2L) {
      next
    }
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        if (identical(as.character(cl$trait[i]), as.character(cl$trait[j]))) {
          next
        }
        # Stable trait order for reproducible pair orientation
        if (as.character(cl$trait[i]) <= as.character(cl$trait[j])) {
          a <- i
          b <- j
        } else {
          a <- j
          b <- i
        }
        cs_a <- get_cs(cl$trait[a], cl$peak_ID[a])
        cs_b <- get_cs(cl$trait[b], cl$peak_ID[b])
        pair_rows[[length(pair_rows) + 1L]] <- .multitrait_cs_overlap_row(
          cluster_id = cid,
          trait_a = cl$trait[a],
          peak_ID_a = cl$peak_ID[a],
          lead_a = cl$peak_variant_ID[a],
          cs_a = cs_a,
          trait_b = cl$trait[b],
          peak_ID_b = cl$peak_ID[b],
          lead_b = cl$peak_variant_ID[b],
          cs_b = cs_b,
          large_cs_fraction = large_cs_fraction
        )
      }
    }
  }

  pairs <- if (length(pair_rows) == 0L) {
    .multitrait_empty_cs_overlap()
  } else {
    out <- do.call(rbind, pair_rows)
    rownames(out) <- NULL
    out
  }

  if (store) {
    .store_write_section_df(object, "multitrait", "cs_overlap", pairs, pheno_name = NULL)
  }
  list(pairs = pairs, clusters = clusters)
}

.multitrait_empty_cs_overlap <- function() {
  data.frame(
    cluster_id = integer(),
    trait_a = character(),
    peak_ID_a = integer(),
    lead_a = character(),
    trait_b = character(),
    peak_ID_b = integer(),
    lead_b = character(),
    n_cs_a = integer(),
    n_cs_b = integer(),
    n_overlap = integer(),
    jaccard = numeric(),
    lead_a_in_cs_b = logical(),
    lead_b_in_cs_a = logical(),
    resolution_a = character(),
    resolution_b = character(),
    cs_fraction_a = numeric(),
    cs_fraction_b = numeric(),
    share_label = character(),
    stringsAsFactors = FALSE
  )
}

.multitrait_cs_info <- function(cred) {
  empty <- list(
    members = character(),
    resolution = NA_character_,
    cs_fraction = NA_real_,
    n_in_set = 0L
  )
  if (is.null(cred) || !is.data.frame(cred) || nrow(cred) == 0L) {
    return(empty)
  }
  if (!"in_credible_set" %in% names(cred)) {
    return(empty)
  }
  members <- as.character(cred$variant_ID[cred$in_credible_set %in% TRUE])
  summ <- attr(cred, "summary")
  resolution <- if (!is.null(summ$resolution)) {
    as.character(summ$resolution)
  } else {
    n <- length(members)
    if (n <= 3L) "high" else if (n <= 10L) "moderate" else "low"
  }
  cs_fraction <- if (!is.null(summ$credible_set_fraction)) {
    as.numeric(summ$credible_set_fraction)
  } else {
    NA_real_
  }
  list(
    members = unique(members[!is.na(members)]),
    resolution = resolution,
    cs_fraction = cs_fraction,
    n_in_set = length(unique(members[!is.na(members)]))
  )
}

.multitrait_cs_jaccard <- function(set_a, set_b) {
  set_a <- unique(as.character(set_a))
  set_b <- unique(as.character(set_b))
  set_a <- set_a[!is.na(set_a) & nzchar(set_a)]
  set_b <- set_b[!is.na(set_b) & nzchar(set_b)]
  if (length(set_a) == 0L && length(set_b) == 0L) {
    return(NA_real_)
  }
  n_int <- length(intersect(set_a, set_b))
  n_uni <- length(union(set_a, set_b))
  if (n_uni == 0L) {
    return(NA_real_)
  }
  n_int / n_uni
}

.multitrait_label_cs_overlap <- function(n_overlap,
                                         lead_a_in_cs_b,
                                         lead_b_in_cs_a,
                                         resolution_a,
                                         resolution_b,
                                         cs_fraction_a,
                                         cs_fraction_b,
                                         large_cs_fraction = 0.8) {
  res_a <- as.character(resolution_a)
  res_b <- as.character(resolution_b)
  both_low <- identical(res_a, "low") && identical(res_b, "low")
  either_low <- identical(res_a, "low") || identical(res_b, "low")
  neither_low <- !either_low && !is.na(res_a) && !is.na(res_b) &&
    nzchar(res_a) && nzchar(res_b)
  large_both <- is.finite(cs_fraction_a) && is.finite(cs_fraction_b) &&
    cs_fraction_a >= large_cs_fraction && cs_fraction_b >= large_cs_fraction

  if (is.na(n_overlap)) {
    return("inconclusive")
  }
  if (both_low) {
    return("inconclusive")
  }
  if (n_overlap >= 1L && large_both && either_low) {
    return("inconclusive")
  }
  if (n_overlap == 0L) {
    if (neither_low) {
      return("independent-nearby")
    }
    return("inconclusive")
  }
  lead_ok <- isTRUE(lead_a_in_cs_b) || isTRUE(lead_b_in_cs_a)
  if (lead_ok && neither_low) {
    return("shared-compatible")
  }
  "possible-shared"
}

.multitrait_cs_overlap_row <- function(cluster_id,
                                       trait_a,
                                       peak_ID_a,
                                       lead_a,
                                       cs_a,
                                       trait_b,
                                       peak_ID_b,
                                       lead_b,
                                       cs_b,
                                       large_cs_fraction = 0.8) {
  members_a <- cs_a$members
  members_b <- cs_b$members
  n_overlap <- length(intersect(members_a, members_b))
  jaccard <- .multitrait_cs_jaccard(members_a, members_b)
  lead_a_chr <- as.character(lead_a)
  lead_b_chr <- as.character(lead_b)
  lead_a_in_b <- length(members_b) > 0L && lead_a_chr %in% members_b
  lead_b_in_a <- length(members_a) > 0L && lead_b_chr %in% members_a
  share_label <- .multitrait_label_cs_overlap(
    n_overlap = n_overlap,
    lead_a_in_cs_b = lead_a_in_b,
    lead_b_in_cs_a = lead_b_in_a,
    resolution_a = cs_a$resolution,
    resolution_b = cs_b$resolution,
    cs_fraction_a = cs_a$cs_fraction,
    cs_fraction_b = cs_b$cs_fraction,
    large_cs_fraction = large_cs_fraction
  )
  data.frame(
    cluster_id = as.integer(cluster_id),
    trait_a = as.character(trait_a),
    peak_ID_a = as.integer(peak_ID_a),
    lead_a = lead_a_chr,
    trait_b = as.character(trait_b),
    peak_ID_b = as.integer(peak_ID_b),
    lead_b = lead_b_chr,
    n_cs_a = as.integer(cs_a$n_in_set),
    n_cs_b = as.integer(cs_b$n_in_set),
    n_overlap = as.integer(n_overlap),
    jaccard = as.numeric(jaccard),
    lead_a_in_cs_b = as.logical(lead_a_in_b),
    lead_b_in_cs_a = as.logical(lead_b_in_a),
    resolution_a = as.character(cs_a$resolution),
    resolution_b = as.character(cs_b$resolution),
    cs_fraction_a = as.numeric(cs_a$cs_fraction),
    cs_fraction_b = as.numeric(cs_b$cs_fraction),
    share_label = as.character(share_label),
    stringsAsFactors = FALSE
  )
}

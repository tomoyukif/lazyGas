# Phase 2 E0: WD config, ann column map, GFF windows / ID map, CS gene selection

#' Write explore-pipeline config under a working directory
#'
#' @param work_dir Working directory path.
#' @param mapping_mode \code{"gwas"} or \code{"qtl"}.
#' @param top_n Max genes for E narration (\code{Inf} = all).
#' @param ... Extra YAML fields (thresholds etc.).
#' @return Invisibly, path to \code{config.yaml}.
#' @export
writeLazyGasExploreConfig <- function(work_dir,
                                      mapping_mode = c("gwas", "qtl"),
                                      top_n = Inf,
                                      ...) {
  mapping_mode <- match.arg(mapping_mode)
  work_dir <- path.expand(as.character(work_dir)[1L])
  if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  }
  cfg <- c(
    list(
      mapping_mode = mapping_mode,
      top_n = if (is.infinite(top_n)) "Inf" else as.numeric(top_n),
      tau_high = 0.8,
      tau_moderate = 0.4,
      relative_high = 0.8,
      relative_moderate = 0.4,
      row_sd_z_warn_lo = 0.8,
      row_sd_z_warn_hi = 1.2,
      # abundance (default) | row_zscore — Z forbids Yanai tau / relative height
      matrix_scale = "abundance",
      # NULL = auto (axes with any non-NA coarse label in sample_map)
      query_axes = NULL
    ),
    list(...)
  )
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for explore config.", call. = FALSE)
  }
  path <- file.path(work_dir, "config.yaml")
  yaml::write_yaml(cfg, path)
  invisible(path)
}

#' Read explore-pipeline config
#'
#' @param work_dir Working directory path.
#' @return Named list. Errors if \code{mapping_mode} missing/invalid.
#' @export
readLazyGasExploreConfig <- function(work_dir) {
  work_dir <- path.expand(as.character(work_dir)[1L])
  path <- file.path(work_dir, "config.yaml")
  if (!file.exists(path)) {
    stop("Missing config.yaml in work_dir: ", work_dir, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for explore config.", call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  mode <- cfg$mapping_mode
  if (is.null(mode) || !as.character(mode)[1L] %in% c("gwas", "qtl")) {
    stop(
      "config mapping_mode must be 'gwas' or 'qtl' (no silent default).",
      call. = FALSE
    )
  }
  cfg$mapping_mode <- as.character(mode)[1L]
  tn <- cfg$top_n
  if (is.null(tn) || identical(tn, "Inf") || (is.character(tn) && tn == "Inf")) {
    cfg$top_n <- Inf
  } else {
    cfg$top_n <- as.numeric(tn)[1L]
  }
  cfg
}

#' Classify annotation columns for E1/E5–E7
#'
#' @param ann Annotation \code{data.frame}.
#' @param use_llm If \code{FALSE}, all columns become \code{free_description}.
#'   If \code{TRUE}, Gemma classifies from column names only (no cell examples).
#' @param llm LLM settings from \code{.ai_config_llm_settings()} (model / base_url /
#'   timeout). Required when \code{use_llm = TRUE}.
#' @param overrides Named character vector: column → class (wins over LLM).
#' @return \code{data.frame} with \code{column}, \code{class}, \code{note}.
#' @export
classifyAnnColumns <- function(ann,
                               use_llm = FALSE,
                               llm = NULL,
                               overrides = NULL) {
  if (is.null(ann) || !is.data.frame(ann) || !ncol(ann)) {
    return(data.frame(
      column = character(),
      class = character(),
      note = character(),
      stringsAsFactors = FALSE
    ))
  }
  cols <- names(ann)
  classes <- rep("free_description", length(cols))
  notes <- rep("use_llm=FALSE default", length(cols))
  allowed <- c(
    "free_description", "go", "kegg", "pfam_interpro", "numeric", "other"
  )
  if (isTRUE(use_llm)) {
    if (is.null(llm) || is.null(llm$model)) {
      llm <- .ai_config_llm_settings()
    }
    classes <- .e0_llm_classify_columns(cols, llm = llm, allowed = allowed)
    notes <- rep("llm column-name classify", length(cols))
  }
  if (!is.null(overrides) && length(overrides)) {
    ov <- as.character(overrides)
    names(ov) <- names(overrides)
    for (i in seq_along(cols)) {
      if (cols[i] %in% names(ov)) {
        classes[i] <- ov[[cols[i]]]
        notes[i] <- "user override"
      }
    }
  }
  # mixed GO+description in one column → go wins (plan E0 §4.4): handled by
  # classifier / override; no post-hoc cell inspection.
  bad <- !classes %in% allowed
  if (any(bad)) {
    stop(
      "Invalid ann column class: ",
      paste(unique(classes[bad]), collapse = ", "),
      call. = FALSE
    )
  }
  data.frame(
    column = cols,
    class = classes,
    note = notes,
    stringsAsFactors = FALSE
  )
}

#' Normalize a single ann-column class token toward the allowed set
#' @keywords internal
.e0_normalize_ann_class <- function(cl, allowed) {
  cl <- as.character(cl)[1L]
  if (!nzchar(cl)) {
    return(cl)
  }
  if (cl %in% allowed) {
    return(cl)
  }
  low <- tolower(trimws(cl))
  hit <- allowed[tolower(allowed) == low]
  if (length(hit)) {
    return(hit[1L])
  }
  # Common drift: "class-kegg", "class_kegg", "class:kegg"
  stripped <- sub("^(class[-_.:\\s]+)", "", low, perl = TRUE)
  hit2 <- allowed[tolower(allowed) == stripped]
  if (length(hit2)) {
    return(hit2[1L])
  }
  cl
}

#' Parse classify JSON → named character vector (column → class)
#' @keywords internal
.e0_parse_classify_json <- function(raw) {
  parsed <- tryCatch(
    jsonlite::fromJSON(raw, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(parsed) || !is.list(parsed)) {
    return(NULL)
  }
  rows <- parsed$columns %||% parsed$Classes %||% parsed
  if (!is.list(rows)) {
    return(NULL)
  }
  class_by <- character()
  for (row in rows) {
    if (!is.list(row)) {
      next
    }
    cn <- as.character(row$column %||% row$name %||% "")[1L]
    cl <- as.character(row$class %||% "")[1L]
    if (!nzchar(cn) || !nzchar(cl)) {
      next
    }
    class_by[[cn]] <- cl
  }
  class_by
}

#' Validate column→class map; return problems + optional normalized out
#' @keywords internal
.e0_validate_classify_map <- function(cols, class_by, allowed, normalize = TRUE) {
  cols <- as.character(cols)
  out <- character(length(cols))
  problems <- character()
  nms <- names(class_by)
  for (i in seq_along(cols)) {
    cn <- cols[i]
    # Character [[ lookup errors on missing names; use %in% + [ ].
    cl <- if (!is.null(nms) && cn %in% nms) {
      as.character(class_by[[cn]])[1L]
    } else {
      NA_character_
    }
    if (is.na(cl) || !nzchar(cl)) {
      problems <- c(problems, paste0("missing class for '", cn, "'"))
      out[i] <- NA_character_
      next
    }
    if (isTRUE(normalize)) {
      cl <- .e0_normalize_ann_class(cl, allowed)
    }
    if (!cl %in% allowed) {
      problems <- c(
        problems,
        paste0("invalid class '", cl, "' for '", cn, "'")
      )
      out[i] <- cl
      next
    }
    out[i] <- cl
  }
  list(ok = !length(problems), out = out, problems = problems)
}

#' One LLM scrutiny pass: re-check prior JSON against allowed labels
#' @keywords internal
.e0_llm_classify_scrutiny <- function(cols, prior_raw, problems, llm, allowed) {
  sys <- paste(
    "You audit and correct annotation-column class assignments.",
    "Allowed class labels exactly (use these strings only):",
    paste(allowed, collapse = ", "),
    ".",
    "Do not invent labels such as class-kegg; use kegg.",
    "Fix every listed problem. Keep correct assignments unchanged.",
    "Return JSON only: {\"columns\":[{\"column\":\"...\",\"class\":\"...\"}, ...]}",
    "Include every input column name exactly once."
  )
  user <- paste0(
    "Column names (JSON array):\n",
    jsonlite::toJSON(as.character(cols), auto_unbox = TRUE),
    "\n\nAllowed classes (JSON array):\n",
    jsonlite::toJSON(as.character(allowed), auto_unbox = TRUE),
    "\n\nProblems found:\n",
    paste(paste0("- ", problems), collapse = "\n"),
    "\n\nPrevious model JSON to correct:\n",
    as.character(prior_raw)[1L]
  )
  tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = sys),
        list(role = "user", content = user)
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = llm$timeout %||% 120
    ),
    error = function(e) e
  )
}

#' Gemma column-name-only classify (one scrutiny pass if labels invalid)
#' @keywords internal
.e0_llm_classify_columns <- function(cols, llm, allowed) {
  cols <- as.character(cols)
  sys <- paste(
    "You classify genome annotation table column names for a plant genetics report.",
    "Allowed classes exactly:",
    paste(allowed, collapse = ", "),
    ".",
    "free_description = prose gene function descriptions;",
    "go = Gene Ontology terms/IDs; kegg = KEGG pathway/KO/EC;",
    "pfam_interpro = Pfam/InterPro/SMART/CDD domains;",
    "numeric = numeric scores/counts; other = IDs/symbols/misc.",
    "If a column mixes GO IDs and prose, choose go.",
    "The JSON field is named class; its value must be one allowed label",
    "(e.g. kegg), never a compound like class-kegg.",
    "Return JSON only: {\"columns\":[{\"column\":\"...\",\"class\":\"...\"}, ...]}",
    "Include every input column name exactly once. No other keys."
  )
  user <- paste0(
    "Column names (JSON array):\n",
    jsonlite::toJSON(cols, auto_unbox = TRUE)
  )
  raw <- tryCatch(
    llmChat(
      messages = list(
        list(role = "system", content = sys),
        list(role = "user", content = user)
      ),
      model = llm$model,
      base_url = llm$base_url,
      json_mode = TRUE,
      timeout = llm$timeout %||% 120
    ),
    error = function(e) e
  )
  if (inherits(raw, "error")) {
    stop(
      "Ann column classify LLM failed: ", conditionMessage(raw),
      call. = FALSE
    )
  }
  class_by <- .e0_parse_classify_json(raw)
  if (is.null(class_by)) {
    stop("Ann column classify returned invalid JSON.", call. = FALSE)
  }
  checked <- .e0_validate_classify_map(cols, class_by, allowed, normalize = TRUE)
  if (isTRUE(checked$ok)) {
    return(checked$out)
  }
  # One scrutiny pass: feed prior answer + allowed labels back to the LLM
  message(
    "Ann column classify needs scrutiny: ",
    paste(checked$problems, collapse = "; ")
  )
  raw2 <- .e0_llm_classify_scrutiny(
    cols = cols,
    prior_raw = raw,
    problems = checked$problems,
    llm = llm,
    allowed = allowed
  )
  if (inherits(raw2, "error")) {
    stop(
      "Ann column classify scrutiny LLM failed: ", conditionMessage(raw2),
      call. = FALSE
    )
  }
  class_by2 <- .e0_parse_classify_json(raw2)
  if (is.null(class_by2)) {
    stop("Ann column classify scrutiny returned invalid JSON.", call. = FALSE)
  }
  checked2 <- .e0_validate_classify_map(cols, class_by2, allowed, normalize = TRUE)
  if (isTRUE(checked2$ok)) {
    return(checked2$out)
  }
  # Soft fill: keep valid labels; missing/invalid → free_description.
  message(
    "Ann column classify scrutiny incomplete; filling gaps with free_description: ",
    paste(checked2$problems, collapse = "; ")
  )
  out <- checked2$out
  bad <- is.na(out) | !nzchar(out) | !out %in% allowed
  out[bad] <- "free_description"
  out
}

#' @keywords internal
.write_ann_column_map <- function(work_dir, map_df) {
  work_dir <- path.expand(as.character(work_dir)[1L])
  if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(work_dir, "ann_column_map.yaml")
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required.", call. = FALSE)
  }
  yaml::write_yaml(map_df, path)
  path
}

#' @keywords internal
.read_ann_column_map <- function(work_dir) {
  path <- file.path(path.expand(work_dir), "ann_column_map.yaml")
  if (!file.exists(path)) {
    return(NULL)
  }
  raw <- yaml::read_yaml(path)
  if (is.null(raw)) {
    return(NULL)
  }
  if (is.data.frame(raw)) {
    return(raw)
  }
  # yaml may return named list of columns
  if (is.list(raw) && all(c("column", "class") %in% names(raw))) {
    return(data.frame(
      column = as.character(raw$column),
      class = as.character(raw$class),
      note = as.character(raw$note %||% ""),
      stringsAsFactors = FALSE
    ))
  }
  # or list of row-lists
  if (is.list(raw) && length(raw) && is.list(raw[[1L]])) {
    return(data.frame(
      column = vapply(raw, function(r) as.character(r$column %||% ""), character(1L)),
      class = vapply(raw, function(r) as.character(r$class %||% ""), character(1L)),
      note = vapply(raw, function(r) as.character(r$note %||% ""), character(1L)),
      stringsAsFactors = FALSE
    ))
  }
  NULL
}

#' Longest-CDS GFF window per gene (−3 kb / +0.5 kb, strand-aware)
#'
#' @param gff \code{GRanges} with gene / CDS features.
#' @param upstream_bp Upstream flank (default 3000).
#' @param downstream_bp Downstream flank (default 500).
#' @return \code{data.frame}: Gene_ID, Chr, window_start, window_end, mid, strand.
#' @keywords internal
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

#' Session cache for GFF → windows / protein_map (avoids per-peak rebuild)
#' @keywords internal
.lazygas_gff_e_cache <- new.env(parent = emptyenv())

#' @keywords internal
.gff_e_cache_key <- function(gff) {
  if (is.null(gff) || !inherits(gff, "GRanges") || !length(gff)) {
    return("empty")
  }
  n <- length(gff)
  paste(
    n,
    as.character(GenomeInfoDb::seqnames(gff))[1L],
    as.character(GenomeInfoDb::seqnames(gff))[n],
    BiocGenerics::start(gff)[1L],
    BiocGenerics::end(gff)[n],
    sep = "|"
  )
}

#' Build (or reuse) gene windows + protein map for E
#' @keywords internal
.gff_e_tables <- function(gff,
                          upstream_bp = 3000L,
                          downstream_bp = 500L,
                          protein_map = NULL) {
  key <- .gff_e_cache_key(gff)
  cached <- .lazygas_gff_e_cache[[key]]
  if (!is.null(cached)) {
    if (!is.null(protein_map) && is.data.frame(protein_map) && nrow(protein_map)) {
      cached$protein_map <- protein_map
    }
    return(cached)
  }
  message("E: building GFF gene windows / protein map (one-time)...")
  t0 <- proc.time()[["elapsed"]]
  windows <- .gff_gene_windows(
    gff,
    upstream_bp = upstream_bp,
    downstream_bp = downstream_bp
  )
  if (is.null(protein_map) || !is.data.frame(protein_map) || !nrow(protein_map)) {
    protein_map <- .gff_gene_protein_map(gff)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  message(
    "E: GFF tables ready — ",
    nrow(windows), " gene windows, ",
    nrow(protein_map), " protein rows (",
    round(elapsed, 1), "s)"
  )
  out <- list(windows = windows, protein_map = protein_map)
  .lazygas_gff_e_cache[[key]] <- out
  out
}

#' Attach Chr/Pos to a credible-set table via peak-block variant map
#'
#' Compact CS stores may omit coordinates. Peak-call tables from
#' \code{.get_peakcall()} carry \code{variant_ID} + \code{Chr} + \code{Pos}.
#'
#' @keywords internal
.e_enrich_cs_coords <- function(cred,
                                object = NULL,
                                pheno_name = NULL,
                                peak_id = NULL) {
  if (is.null(cred) || !nrow(cred)) {
    return(cred)
  }
  has_pos <- "Pos" %in% names(cred) &&
    any(is.finite(as.numeric(cred$Pos)), na.rm = TRUE)
  has_chr <- any(c("Chr", "chr") %in% names(cred)) &&
    any(nzchar(as.character(cred$Chr %||% cred$chr)), na.rm = TRUE)
  if (has_pos && has_chr) {
    return(cred)
  }
  if (is.null(object) || is.null(pheno_name) || is.null(peak_id)) {
    return(cred)
  }
  if (!"variant_ID" %in% names(cred)) {
    return(cred)
  }
  pk <- tryCatch(
    .get_peakcall(object = object, pheno_name = pheno_name, recalc = TRUE),
    error = function(e) NULL
  )
  if (is.null(pk) || !nrow(pk)) {
    pk <- tryCatch(
      .get_peakcall(object = object, pheno_name = pheno_name, recalc = FALSE),
      error = function(e) NULL
    )
  }
  if (is.null(pk) || !nrow(pk)) {
    return(cred)
  }
  pk <- pk[as.character(pk$peak_ID) == as.character(peak_id), , drop = FALSE]
  pos_cols <- intersect(c("variant_ID", "Chr", "Pos"), names(pk))
  if (!all(c("variant_ID", "Chr", "Pos") %in% pos_cols)) {
    return(cred)
  }
  pos_map <- unique(pk[, pos_cols, drop = FALSE])
  drop_old <- intersect(c("Chr", "Pos", "chr", "pos"), names(cred))
  cred2 <- if (length(drop_old)) {
    cred[, setdiff(names(cred), drop_old), drop = FALSE]
  } else {
    cred
  }
  out <- merge(cred2, pos_map, by = "variant_ID", all.x = TRUE)
  if ("PIP" %in% names(out)) {
    out <- out[order(-as.numeric(out$PIP), out$variant_ID), , drop = FALSE]
  }
  rownames(out) <- NULL
  summ <- attr(cred, "summary")
  if (!is.null(summ)) {
    attr(out, "summary") <- summ
  }
  out
}

#' @keywords internal
.e_cs_chr_minmax <- function(cred) {
  if (is.null(cred) || !nrow(cred)) {
    return(NULL)
  }
  in_cs <- if ("in_credible_set" %in% names(cred)) {
    as.logical(cred$in_credible_set)
  } else {
    rep(TRUE, nrow(cred))
  }
  in_cs[is.na(in_cs)] <- FALSE
  sub <- cred[in_cs, , drop = FALSE]
  if (!nrow(sub)) {
    return(NULL)
  }
  if (!"Pos" %in% names(sub) || all(!is.finite(as.numeric(sub$Pos)))) {
    return(NULL)
  }
  chr <- as.character(sub$Chr %||% sub$chr)
  chr <- chr[nzchar(chr) & !is.na(chr)]
  if (!length(chr)) {
    return(NULL)
  }
  list(
    Chr = chr[1L],
    start = min(as.numeric(sub$Pos), na.rm = TRUE),
    end = max(as.numeric(sub$Pos), na.rm = TRUE),
    n = nrow(sub),
    sub = sub
  )
}

#' Select E genes for one peak (CS + GFF rules)
#' @keywords internal
.report_e_select_genes <- function(mapping_mode,
                                   simple_candidates,
                                   credible_set,
                                   gff_windows,
                                   snpeff = NULL,
                                   top_n = Inf) {
  empty <- data.frame(
    Gene_ID = character(),
    Name = character(),
    Chr = character(),
    dist2peak = numeric(),
    max_PIP = numeric(),
    nearest_cs_PIP = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(simple_candidates) || !nrow(simple_candidates)) {
    return(empty)
  }
  cand <- simple_candidates
  cand$Gene_ID <- as.character(cand$Gene_ID)
  win <- gff_windows
  if (is.null(win) || !nrow(win)) {
    return(empty)
  }
  csinfo <- .e_cs_chr_minmax(credible_set)
  if (is.null(csinfo)) {
    return(empty)
  }
  mode <- as.character(mapping_mode)[1L]
  ids <- character()
  if (identical(mode, "qtl")) {
    w <- win[win$Chr == csinfo$Chr, , drop = FALSE]
    hit <- w$window_end >= csinfo$start & w$window_start <= csinfo$end
    ids <- unique(w$Gene_ID[hit])
  } else {
    # gwas
    if (!is.null(snpeff) && nrow(snpeff) && "Gene_ID" %in% names(snpeff) &&
        "Pos" %in% names(snpeff)) {
      sub <- csinfo$sub
      snp <- snpeff
      snp$Gene_ID <- as.character(snp$Gene_ID)
      snp$Pos <- as.numeric(snp$Pos)
      snp$Chr <- as.character(snp$Chr %||% snp$chr)
      key_cs <- paste(as.character(sub$Chr %||% sub$chr), as.numeric(sub$Pos), sep = ":")
      key_sn <- paste(snp$Chr, snp$Pos, sep = ":")
      ids <- unique(snp$Gene_ID[key_sn %in% key_cs])
    } else {
      w <- win[win$Chr == csinfo$Chr, , drop = FALSE]
      pos <- as.numeric(csinfo$sub$Pos)
      hit <- vapply(seq_len(nrow(w)), function(i) {
        any(pos >= w$window_start[i] & pos <= w$window_end[i], na.rm = TRUE)
      }, logical(1L))
      ids <- unique(w$Gene_ID[hit])
    }
  }
  ids <- intersect(ids, unique(cand$Gene_ID))
  if (!length(ids)) {
    return(empty)
  }
  out <- cand[cand$Gene_ID %in% ids, , drop = FALSE]
  # attach PIP columns
  out$max_PIP <- NA_real_
  out$nearest_cs_PIP <- NA_real_
  pip_map <- NULL
  if ("PIP" %in% names(csinfo$sub) && "Pos" %in% names(csinfo$sub)) {
    pip_map <- csinfo$sub
  }
  if (!is.null(snpeff) && nrow(snpeff) && !is.null(pip_map) &&
      identical(mode, "gwas")) {
    # max PIP among CS variants assigned to gene
    snp <- snpeff
    snp$Gene_ID <- as.character(snp$Gene_ID)
    snp$Pos <- as.numeric(snp$Pos)
    snp$Chr <- as.character(snp$Chr %||% snp$chr)
    for (gid in unique(out$Gene_ID)) {
      gpos <- snp$Pos[snp$Gene_ID == gid]
      if (!length(gpos)) next
      hit <- as.numeric(pip_map$Pos) %in% gpos
      if (any(hit)) {
        out$max_PIP[out$Gene_ID == gid] <- max(as.numeric(pip_map$PIP[hit]), na.rm = TRUE)
      }
    }
  }
  if (identical(mode, "qtl") && !is.null(pip_map)) {
    for (i in seq_len(nrow(out))) {
      gid <- out$Gene_ID[i]
      wrow <- win[win$Gene_ID == gid, , drop = FALSE]
      if (!nrow(wrow)) next
      mid <- wrow$mid[1L]
      d <- abs(as.numeric(pip_map$Pos) - mid)
      j <- which.min(d)
      out$nearest_cs_PIP[i] <- as.numeric(pip_map$PIP[j])
    }
  }
  if (identical(mode, "gwas") && !is.null(pip_map)) {
    for (i in seq_len(nrow(out))) {
      if (is.finite(out$max_PIP[i])) next
      gid <- out$Gene_ID[i]
      wrow <- win[win$Gene_ID == gid & win$Chr == csinfo$Chr, , drop = FALSE]
      if (!nrow(wrow)) next
      pos <- as.numeric(pip_map$Pos)
      hit <- pos >= wrow$window_start[1L] & pos <= wrow$window_end[1L]
      if (any(hit)) {
        out$max_PIP[i] <- max(as.numeric(pip_map$PIP[hit]), na.rm = TRUE)
      }
    }
  }
  # sort + top_n
  if (identical(mode, "qtl")) {
    ord <- order(-out$nearest_cs_PIP, out$dist2peak, out$Gene_ID, na.last = TRUE)
  } else {
    ord <- order(-out$max_PIP, out$dist2peak, out$Gene_ID, na.last = TRUE)
  }
  out <- out[ord, , drop = FALSE]
  if (is.finite(top_n) && top_n >= 1 && nrow(out) > top_n) {
    out <- out[seq_len(as.integer(top_n)), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' @keywords internal
.write_run_metadata <- function(work_dir, meta) {
  path <- file.path(path.expand(work_dir), "run_metadata.json")
  jsonlite::write_json(meta, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  path
}

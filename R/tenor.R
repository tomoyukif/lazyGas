################################################################################
#' RAP-DB / TENOR expression download and evaluation (Phase 1)
#'
#' Production pipelines are expected to use a pre-built gene x condition
#' expression spreadsheet. For Phase 1 testing, these helpers download TENOR
#' RNA-seq TPM values from RAP-DB for genes in candidate regions, cache them
#' under \code{lazygas/ai_report/tenor_download/}, and score expression against
#' a phenotype query.
#'
#' @name tenor-rapdb
#' @keywords internal
NULL

.TENOR_RAPDB_BASE <- "https://rapdb.dna.naro.go.jp"

#' Built-in TENOR dataset catalog used by RAP-DB Expression (TENOR) UI
#'
#' @return A data.frame with columns \code{dataset} and \code{publication}.
#' @export
tenorDatasetCatalog <- function() {
  data.frame(
    dataset = c(
      "Strand-specific RNA-Seq (8 tissues, such as flower, leaf, root, seed, etc.)",
      "Total RNA-Seq (8 tissues, such as flower, leaf, root, seed, etc.)",
      "Strand-specific RNA-Seq (shoot, seed, anther, pistil)",
      "Strand-specific mRNA-Seq (leaf)",
      "Strand-specific mRNA-Seq (seedling under salinity, cold and drought stress conditions)",
      "Strand-specific mRNA-Seq (shoot under dehydration and rehydration stress conditions)",
      "Strand-specific mRNA-Seq (shoot and root under nitrogen starvation)",
      "Strand-specific mRNA-Seq (lamina joint, OsDCL3a RNAi mutant)",
      "Strand-specific total RNA-Seq (root under phosphate starvation and recovery conditions, dcl3a RNAi mutant)",
      "Strand-specific RNA-Seq (root and shoot under iron-deficient conditions, osbhlh156-1 mutant)",
      "Strand-specific RNA-Seq (root after inoculation of arbuscular mycorrhizal fungus, phr2 mutant)",
      "Total RNA-Seq (root under mineral deficiency conditions)",
      "mRNA-Seq (7 tissues, such as leaf, root, panicle, callus, etc.)",
      "mRNA-Seq (9 tissues, such as leaf, flower, anther, pstill, etc)",
      "mRNA-Seq (leaf)",
      "mRNA-Seq (leaf and callus)",
      "mRNA-Seq (anther)",
      "mRNA-Seq (egg, sperm, zygote)",
      "mRNA-Seq (aleurone)",
      "mRNA-Seq (root meristematic, elongation, differentiation zone)",
      "mRNA-Seq (culm tissue containing shoot apical meristem)",
      "mRNA-Seq (root and root tips in early development)",
      "mRNA-Seq under salinity stress conditions (TENOR)",
      "mRNA-Seq under phosphate stress conditions (TENOR)",
      "mRNA-Seq under phosphate stress conditions (TENOR)",
      "mRNA-Seq under cadmium stress conditions (TENOR)",
      "mRNA-Seq under cadmium stress conditions (TENOR)",
      "mRNA-Seq under cold, flood and osmotic stress conditions (TENOR)",
      "mRNA-Seq under drought stress conditions (TENOR)",
      "mRNA-Seq under JA treatment (TENOR)",
      "mRNA-Seq under ABA treatment (TENOR)",
      "mRNA-Seq under no treatment (control data in TENOR)",
      "mRNA-Seq (root and shoot under phosphate starvation and recovery conditions)",
      "mRNA-Seq (root under phosphate starvation and recovery conditions)",
      "mRNA-Seq (shoot and root under arsenic stress)",
      "mRNA-Seq (coleoptile and the first leaf under shade)",
      "mRNA-Seq (leaf under drought stress condition in long day and short day)",
      "mRNA-Seq (shoot, MET1-2 mutant)",
      "mRNA-Seq (shoot, osino80 RNAi mutant)",
      "mRNA-Seq (whole tissue of above-ground seedling, oscmt3a mutant)",
      "mRNA-Seq (leaf after inoculation of Magnaporthe oryzae strains)",
      "mRNA-Seq (leaf after inoculation of Xanthomonas oryzae pv. oryzicola (Xoc) strains)",
      "mRNA-Seq (leaf after inoculation of Xanthomonas oryzae pv. oryzicola (Xoo) strains)",
      "mRNA-Seq (root after inoculation of arbuscular mycorrhizal fungus)",
      "mRNA-Seq (root tip and root after inoculation of root knot and root rot nematodes)",
      "mRNA-Seq (root vascular cell after inoculation of root knot nematode)",
      "mRNA-Seq (root under AM fungus-secreted substances treatment, oscerk1 mutant)"
    ),
    publication = c(
      "Wang et al. 2015 Plant J.",
      "Wang et al. 2015 Plant J.",
      "Zhang YC et al. 2014 Genome Biol.",
      "Lu et al. 2019 Nat Plants",
      "Lu et al. 2012 BMC Genomics",
      "Park et al. 2023 Int J Mol Sci.",
      "Shin et al. 2018 BMC Genomics",
      "Wei et al. 2014 Proc Natl Acad Sci U S A.",
      "Secco et al. 2015 Elife.",
      "Wang et al. 2019 New Phytol",
      "Das et al. 2022 Nat Commun.",
      "Dong et al. 2018 Plant Cell",
      "Sakai H et al. 2011 Genome Biol Evol.",
      "Davidson RM et al. 2012 Plant J.",
      "Yang et al. 2019 BMC Plant Biol",
      "Wu et al. 2011 Plant Cell",
      "Komiya et al. 2014 Plant J.",
      "Rahman et al. 2019 Plant Cell Physiol.",
      "Watanabe et al. 2014 Genomics.",
      "Huang et al. 2015 Plant Cell.",
      "Song et al. 2015 Proc Natl Acad Sci U S A.",
      "Kyndt et al. 2012 J Exp Bot.",
      "Mizuno H et al. 2010 BMC Genomics",
      "Oono et al. 2011 Rice",
      "Oono et al. 2013 Plant Mol Biol.",
      "Oono et al. 2014 PLoS One",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Kawahara et al. 2016 Plant Cell Physiol.",
      "Secco et al. 2013 Plant Cell",
      "Secco et al. 2015 Elife.",
      "Yu et al. 2012 New Phytol.",
      "Liu et al. 2016 J Integr Plant Biol.",
      "Galbiati et al. 2016 Plant Cell Environ.",
      "Hu et al. 2014 Proc Natl Acad Sci U S A.",
      "Li et al. 2018 J Integr Plant Biol.",
      "Cheng et al. 2015 Plant J.",
      "Kawahara et al. 2012 PLoS ONE",
      "Wilkins et al. 2015 Front Plant Sci.",
      "Hummel et al. 2017 Mol Plant Pathol.",
      "Wang et al. 2020 Proc Natl Acad Sci U S A.",
      "Kyndt et al. 2012 New Phytol",
      "Ji et al. 2013 J Exp Bot.",
      "Miyata et al. 2014 Plant Cell Physiol."
    ),
    stringsAsFactors = FALSE
  )
}

#' Resolve cache directory for TENOR downloads
#'
#' @param object Optional \code{LazyGas} object. When given, files are stored
#'   under \code{<dirname(store)>/lazygas/ai_report/tenor_download/}.
#' @param out_dir Explicit output directory (overrides \code{object}).
#' @return Character path.
#' @export
tenorCacheDir <- function(object = NULL, out_dir = NULL) {
  if (!is.null(out_dir) && nzchar(out_dir)) {
    return(normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  }
  if (!is.null(object) && inherits(object, "LazyGas")) {
    root <- dirname(.store_path(object))
    return(file.path(root, "lazygas", "ai_report", "tenor_download"))
  }
  normalizePath(
    file.path("lazygas", "ai_report", "tenor_download"),
    winslash = "/",
    mustWork = FALSE
  )
}

#' Download TENOR expression profiles from RAP-DB for gene IDs
#'
#' For each gene ID, resolves RAP-DB locus/transcripts via
#' \code{/tools/Feature}, then queries \code{/api/rnaseq} for each TENOR
#' dataset. Results are written to one CSV per gene under
#' \code{lazygas/ai_report/tenor_download/} with columns
#' \code{data_description} and \code{tpm}. Non-TPM fields are joined with
#' \code{":"} into \code{data_description}. Existing files are skipped unless
#' \code{overwrite = TRUE}.
#'
#' Bulk matrix download is not provided by RAP-DB; this per-gene API is the
#' supported Phase 1 path for candidate-region testing.
#'
#' @param gene_ids Character vector of RAP gene / transcript IDs
#'   (e.g. \code{Os01g0911700} or \code{Os01t0911700-01}).
#' @param object Optional \code{LazyGas} object for default cache location.
#' @param out_dir Cache directory (see [tenorCacheDir()]).
#' @param datasets Optional data.frame with \code{dataset} and
#'   \code{publication} columns. Default: [tenorDatasetCatalog()].
#' @param primary_transcript_only If \code{TRUE}, download only the first
#'   transcript per locus (typically \code{*-01}).
#' @param overwrite Re-download even if CSV already exists.
#' @param sleep_sec Pause between API calls (be polite to RAP-DB).
#' @param timeout HTTP timeout in seconds.
#' @param quiet Suppress progress messages.
#'
#' @return A data.frame summarizing download status per gene
#'   (\code{gene_id}, \code{locus_id}, \code{transcript_id}, \code{path},
#'   \code{status}, \code{n_rows}).
#' @export
#'
#' @examples
#' \dontrun{
#' downloadTenorExpression("Os01g0911700")
#' }
downloadTenorExpression <- function(gene_ids,
                                    object = NULL,
                                    out_dir = NULL,
                                    datasets = NULL,
                                    primary_transcript_only = TRUE,
                                    overwrite = FALSE,
                                    sleep_sec = 0.15,
                                    timeout = 60,
                                    quiet = FALSE) {
  gene_ids <- unique(as.character(gene_ids))
  gene_ids <- gene_ids[nzchar(gene_ids) & !is.na(gene_ids)]
  if (!length(gene_ids)) {
    stop("'gene_ids' must contain at least one ID.", call. = FALSE)
  }
  if (is.null(datasets)) {
    datasets <- tenorDatasetCatalog()
  }
  if (!all(c("dataset", "publication") %in% names(datasets))) {
    stop("'datasets' must have columns dataset and publication.", call. = FALSE)
  }

  cache_dir <- tenorCacheDir(object = object, out_dir = out_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- vector("list", length(gene_ids))
  for (i in seq_along(gene_ids)) {
    gid <- gene_ids[[i]]
    path <- file.path(cache_dir, paste0(.store_safe_name(gid), ".csv"))
    if (file.exists(path) && !isTRUE(overwrite)) {
      if (!quiet) {
        message("Skip (cached): ", gid)
      }
      n_rows <- tryCatch(nrow(utils::read.csv(path, stringsAsFactors = FALSE)),
                         error = function(e) NA_integer_)
      rows[[i]] <- data.frame(
        gene_id = gid,
        locus_id = NA_character_,
        transcript_id = NA_character_,
        path = path,
        status = "cached",
        n_rows = n_rows,
        stringsAsFactors = FALSE
      )
      next
    }

    feat <- tryCatch(
      .rapdb_feature(gid, timeout = timeout),
      error = function(e) {
        if (!quiet) {
          warning("Feature lookup failed for ", gid, ": ",
                  conditionMessage(e), call. = FALSE)
        }
        NULL
      }
    )
    transcripts <- .rapdb_transcript_ids(feat, gene_id = gid)
    if (!length(transcripts)) {
      rows[[i]] <- data.frame(
        gene_id = gid,
        locus_id = NA_character_,
        transcript_id = NA_character_,
        path = path,
        status = "no_transcript",
        n_rows = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    if (isTRUE(primary_transcript_only)) {
      transcripts <- transcripts[[1L]]
    }

    locus_id <- if (!is.null(feat$locus_name)) {
      as.character(feat$locus_name)
    } else {
      .rapdb_guess_locus(gid)
    }

    chunks <- list()
    for (tid in transcripts) {
      for (j in seq_len(nrow(datasets))) {
        ds <- datasets$dataset[[j]]
        pub <- datasets$publication[[j]]
        recs <- tryCatch(
          .rapdb_rnaseq(ds, pub, tid, timeout = timeout),
          error = function(e) {
            if (!quiet) {
              message("API miss ", tid, " / dataset ", j, ": ",
                      conditionMessage(e))
            }
            NULL
          }
        )
        if (!is.null(recs) && length(recs)) {
          df <- .rapdb_rnaseq_to_df(recs)
          df$gene_id <- gid
          df$locus_id <- locus_id
          df$transcript_id <- tid
          df$dataset <- ds
          df$publication <- pub
          chunks[[length(chunks) + 1L]] <- df
        }
        if (sleep_sec > 0) {
          Sys.sleep(sleep_sec)
        }
      }
    }

    if (!length(chunks)) {
      rows[[i]] <- data.frame(
        gene_id = gid,
        locus_id = locus_id,
        transcript_id = paste(transcripts, collapse = ","),
        path = path,
        status = "empty",
        n_rows = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    out <- do.call(rbind, chunks)
    out <- .tenor_export_table(out)
    utils::write.csv(out, path, row.names = FALSE)
    if (!quiet) {
      message("Wrote ", nrow(out), " rows -> ", path)
    }
    rows[[i]] <- data.frame(
      gene_id = gid,
      locus_id = locus_id,
      transcript_id = paste(transcripts, collapse = ","),
      path = path,
      status = "downloaded",
      n_rows = nrow(out),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Download TENOR for genes in a LazyGas candidate region
#'
#' Reads candidate genes for \code{pheno} and calls [downloadTenorExpression()].
#' Intended to run after peak calling / \code{recalcAssoc} / candidate listing.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno Phenotype name or index.
#' @param ... Passed to [downloadTenorExpression()].
#' @return Download status data.frame.
#' @export
downloadTenorForCandidates <- function(object, pheno, ...) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }
  pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
  candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno_name)
  if (is.null(candidate) || nrow(candidate) == 0L) {
    stop("No candidate data for phenotype '", pheno_name, "'.", call. = FALSE)
  }
  if (!"Gene_ID" %in% names(candidate)) {
    stop("Candidate table lacks Gene_ID.", call. = FALSE)
  }
  downloadTenorExpression(
    gene_ids = unique(as.character(candidate$Gene_ID)),
    object = object,
    ...
  )
}

#' Evaluate TENOR expression for one gene against a phenotype query
#'
#' Reads a cached CSV from [downloadTenorExpression()] (or downloads if
#' missing) and returns an evidence-like record with a 0--1 score.
#'
#' Scoring (deterministic Phase 1 default):
#' \enumerate{
#'   \item Match rows whose \code{data_description} overlaps query keywords
#'     (case-insensitive).
#'   \item If matches exist, score =
#'     \code{min(1, mean(TPM_match) / max(mean(TPM_all), eps)) * presence},
#'     where presence is 1 if any matched TPM > 0.
#'   \item If no keyword match, fall back to a weak score from genome-wide
#'     max TPM (\code{min(1, log1p(max_tpm) / log1p(50))}).
#' }
#' LLM-based tissue/phenotype relevance can be layered later via
#' \code{use_llm = TRUE} (currently ignored with a message).
#'
#' @param gene_id Gene identifier.
#' @param query Phenotype query text, character keywords, or
#'   \code{PhenotypeQuery}.
#' @param object,out_dir Cache location (see [tenorCacheDir()]).
#' @param download_if_missing Call [downloadTenorExpression()] when CSV absent.
#' @param use_llm Reserved for Phase 1 LLM tissue relevance (not yet applied).
#' @param ... Passed to [downloadTenorExpression()] when downloading.
#'
#' @return A list with \code{source}, \code{score}, \code{snippets},
#'   \code{details} (compatible with [collectGeneEvidence()] records).
#' @export
evaluateTenorExpression <- function(gene_id,
                                    query = NULL,
                                    object = NULL,
                                    out_dir = NULL,
                                    download_if_missing = TRUE,
                                    use_llm = FALSE,
                                    ...) {
  gene_id <- as.character(gene_id)[[1L]]
  cache_dir <- tenorCacheDir(object = object, out_dir = out_dir)
  path <- file.path(cache_dir, paste0(.store_safe_name(gene_id), ".csv"))

  if (!file.exists(path)) {
    if (!isTRUE(download_if_missing)) {
      return(list(
        source = "expression",
        score = 0,
        snippets = character(),
        details = list(status = "missing_cache", path = path)
      ))
    }
    downloadTenorExpression(
      gene_ids = gene_id,
      object = object,
      out_dir = out_dir,
      ...
    )
  }
  if (!file.exists(path)) {
    return(list(
      source = "expression",
      score = 0,
      snippets = character(),
      details = list(status = "download_failed", path = path)
    ))
  }

  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(dat) || !"tpm" %in% names(dat)) {
    return(list(
      source = "expression",
      score = 0,
      snippets = character(),
      details = list(status = "empty", path = path)
    ))
  }
  # Accept legacy wide CSVs by collapsing non-TPM columns if needed
  if (!"data_description" %in% names(dat)) {
    dat <- .tenor_export_table(dat)
  }
  dat$tpm <- as.numeric(dat$tpm)
  dat <- dat[is.finite(dat$tpm), , drop = FALSE]
  dat$data_description <- as.character(dat$data_description)

  keywords <- .tenor_query_keywords(query)
  row_text <- dat$data_description

  matched <- rep(FALSE, nrow(dat))
  if (length(keywords)) {
    row_text_l <- tolower(row_text)
    matched <- Reduce(`|`, lapply(keywords, function(kw) {
      grepl(tolower(kw), row_text_l, fixed = TRUE)
    }))
  }

  max_tpm <- max(dat$tpm, na.rm = TRUE)
  mean_all <- mean(dat$tpm, na.rm = TRUE)
  eps <- 1e-6

  if (any(matched)) {
    mean_m <- mean(dat$tpm[matched], na.rm = TRUE)
    presence <- as.numeric(any(dat$tpm[matched] > 0))
    enrich <- mean_m / max(mean_all, eps)
    score <- presence * min(1, enrich)
    top <- dat[matched, , drop = FALSE]
    top <- top[order(top$tpm, decreasing = TRUE), , drop = FALSE]
    snippets <- utils::head(
      sprintf("TPM=%.3f | %s", top$tpm, top$data_description),
      5L
    )
    method <- "keyword_enrichment"
  } else {
    score <- min(1, log1p(max_tpm) / log1p(50))
    snippets <- sprintf("max_tpm=%.3f (no phenotype keyword match)", max_tpm)
    method <- "max_tpm_fallback"
  }

  if (isTRUE(use_llm)) {
    message("evaluateTenorExpression(use_llm=TRUE): LLM tissue scoring ",
            "not yet applied; using deterministic score.")
  }

  list(
    source = "expression",
    score = as.numeric(score),
    snippets = as.character(snippets),
    details = list(
      path = path,
      method = method,
      n_rows = nrow(dat),
      n_matched = sum(matched),
      max_tpm = max_tpm,
      mean_tpm = mean_all,
      keywords = keywords
    )
  )
}

# ---- internals -------------------------------------------------------------

.tenor_query_keywords <- function(query) {
  if (is.null(query)) {
    return(character())
  }
  if (inherits(query, "PhenotypeQuery")) {
    kw <- unique(c(
      query$trait_keywords,
      query$tissues,
      query$developmental_stage,
      query$conditions
    ))
    return(as.character(kw[nzchar(kw) & !is.na(kw)]))
  }
  if (is.character(query)) {
    toks <- unlist(strsplit(paste(query, collapse = " "), "[[:space:],;]+"))
    return(unique(toks[nzchar(toks)]))
  }
  character()
}

.rapdb_http_get_json <- function(url, timeout = 60) {
  # RAP-DB often returns JSON with Content-Type: text/html
  parse_body <- function(txt) {
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
  }
  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    req <- httr2::req_timeout(req, timeout)
    req <- httr2::req_headers(
      req,
      Accept = "application/json",
      `User-Agent` = "lazyGas/TENOR"
    )
    resp <- httr2::req_perform(req)
    txt <- httr2::resp_body_string(resp)
    return(parse_body(txt))
  }
  if (requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    curl::handle_setheaders(
      h,
      Accept = "application/json",
      `User-Agent` = "lazyGas/TENOR"
    )
    curl::handle_setopt(h, timeout = timeout)
    raw <- curl::curl_fetch_memory(url, handle = h)$content
    return(parse_body(rawToChar(raw)))
  }
  stop("HTTP requests require 'httr2' or 'curl'.", call. = FALSE)
}

.rapdb_feature <- function(name, timeout = 60) {
  url <- paste0(
    .TENOR_RAPDB_BASE,
    "/tools/Feature?name=",
    utils::URLencode(as.character(name), reserved = TRUE)
  )
  .rapdb_http_get_json(url, timeout = timeout)
}

.rapdb_rnaseq <- function(dataset, publication, transcript, timeout = 60) {
  q <- paste0(
    "dataset=", utils::URLencode(dataset, reserved = TRUE),
    "&publication=", utils::URLencode(publication, reserved = TRUE),
    "&transcript=", utils::URLencode(transcript, reserved = TRUE)
  )
  url <- paste0(.TENOR_RAPDB_BASE, "/api/rnaseq?", q)
  .rapdb_http_get_json(url, timeout = timeout)
}

.rapdb_transcript_ids <- function(feat, gene_id) {
  if (is.null(feat) || !length(feat)) {
    # If gene_id already looks like a transcript, try it directly
    if (grepl("t[0-9]", gene_id, ignore.case = TRUE)) {
      return(gene_id)
    }
    return(character())
  }
  attrs <- feat$transcripts_attributes
  if (is.null(attrs) || !length(attrs)) {
    if (grepl("t[0-9]", gene_id, ignore.case = TRUE)) {
      return(gene_id)
    }
    return(character())
  }
  ids <- vapply(attrs, function(a) {
    id <- a$ID %||% a$Name %||% NA_character_
    as.character(id)
  }, character(1))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  # Prefer *-01 first
  ord <- order(!grepl("-01$", ids), ids)
  ids[ord]
}

.rapdb_guess_locus <- function(gene_id) {
  if (grepl("g[0-9]", gene_id, ignore.case = TRUE) &&
      !grepl("t[0-9]", gene_id, ignore.case = TRUE)) {
    return(gene_id)
  }
  sub("t", "g", sub("-\\d+$", "", gene_id), ignore.case = TRUE)
}

.rapdb_rnaseq_to_df <- function(recs) {
  if (!length(recs)) {
    return(data.frame())
  }
  data.frame(
    name = vapply(recs, function(x) as.character(x$name %||% NA_character_), character(1)),
    tissue_and_condition = vapply(
      recs,
      function(x) as.character(x$tissue_and_condition %||% x$name %||% NA_character_),
      character(1)
    ),
    plant_ontology = vapply(
      recs,
      function(x) as.character(x$plant_ontology %||% NA_character_),
      character(1)
    ),
    plant_experimental_conditions_ontology = vapply(
      recs,
      function(x) {
        as.character(x$plant_experimental_conditions_ontology %||% NA_character_)
      },
      character(1)
    ),
    experiment = vapply(
      recs,
      function(x) as.character(x$experiment %||% NA_character_),
      character(1)
    ),
    tpm = vapply(recs, function(x) as.numeric(x$tpm %||% NA_real_), numeric(1)),
    alignment_rate = vapply(
      recs,
      function(x) as.numeric(x$alignment_rate %||% NA_real_),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

#' Collapse all non-TPM columns into colon-separated data_description
#'
#' @param df Wide TENOR table including a \code{tpm} column.
#' @return data.frame with \code{data_description} and \code{tpm}.
#' @noRd
.tenor_export_table <- function(df) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(data.frame(
      data_description = character(),
      tpm = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  if (!"tpm" %in% names(df)) {
    stop("TENOR table must contain a 'tpm' column.", call. = FALSE)
  }
  desc_cols <- setdiff(names(df), "tpm")
  desc <- if (!length(desc_cols)) {
    rep("", nrow(df))
  } else {
    apply(df[desc_cols], 1L, function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(trimws(x))]
      # Avoid embedding colons that break round-trips awkwardly: keep as-is
      # (ontology strings already contain colons like PO:0009010).
      paste(x, collapse = ":")
    })
  }
  data.frame(
    data_description = as.character(desc),
    tpm = as.numeric(df$tpm),
    stringsAsFactors = FALSE
  )
}

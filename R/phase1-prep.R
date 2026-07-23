################################################################################
#' Prepare Phase 1 peaked LazyGas object (JRC + WRC)
#'
#' Runs the Phase 1 preprocessing pipeline against paths from
#' \code{inst/config/ai-lazygas-data.yaml}:
#' \enumerate{
#'   \item Open existing genotype GDS (reuse; no VCF reconversion required)
#'   \item Reset sample filter; apply Phase 1 marker INFO / QC filters
#'   \item Assign JRC + WRC phenotypes (sample IDs \code{JRC*}/\code{WRC*})
#'   \item \code{scanAssoc(..., method = "mlm")} (LMM; \code{"lmm"} accepted)
#'   \item \code{callPeakBlock} → \code{recalcAssoc} → \code{listCandidate}
#'     with nb_combined GFF/ann and Phase 1 SnpEff GDS
#' }
#'
#' @param work_dir Working directory (default from config).
#' @param phenos Phenotype names to assign/scan. Default:
#'   \code{Amylose content}, \code{Leaf width}, \code{Heading date}.
#' @param method Scan method: \code{"mlm"} (default) or alias \code{"lmm"}.
#' @param overwrite Replace existing companion results.
#' @param n_threads Threads for peak/recalc steps.
#' @param limit_peakcall Max peaks per phenotype for \code{callPeakBlock}
#'   (default 50 for Phase 1 final retest).
#' @param skip_if_ready If \code{TRUE} and candidates already exist, skip the
#'   heavy scan/peak/recalc/candidate steps and return the opened object.
#' @param apply_qc_filter If \code{TRUE} (default), apply Phase 1 GBScleanR
#'   INFO / marker filters after resetting sample filters.
#' @param skip_scan If \code{TRUE}, skip \code{scanAssoc} and reuse existing
#'   scan results in the companion store.
#' @param skip_peakcall If \code{TRUE}, skip \code{callPeakBlock} and reuse
#'   existing peakcall results (implies scan already present).
#' @param omit_outlier If \code{TRUE} (Phase 1 default), pass
#'   \code{omit_outlier = TRUE} to \code{scanAssoc} so boxplot outliers are
#'   set to NA before the LMM scan for each continuous phenotype.
#'
#' @return A \code{LazyGas} object (opened).
#' @export
prepPhase1JrcWrc <- function(work_dir = NULL,
                             phenos = NULL,
                             method = NULL,
                             overwrite = FALSE,
                             n_threads = NULL,
                             limit_peakcall = 50L,
                             skip_if_ready = FALSE,
                             apply_qc_filter = TRUE,
                             skip_scan = FALSE,
                             skip_peakcall = FALSE,
                             omit_outlier = TRUE) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read phenotype xlsx files.",
         call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required to import GFF.", call. = FALSE)
  }

  work_dir <- if (!is.null(work_dir) && nzchar(work_dir)) {
    work_dir
  } else {
    phase1DataPath(
      "work_dir",
      default = file.path(getwd(), "demo_output", "phase1_jrc_wrc")
    )
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- normalizePath(work_dir, winslash = "/")

  method <- .phase1_resolve_scan_method(method)
  gds_src <- phase1DataPath(
    "genotype_gds",
    env_var = "LAZYGAS_PHASE1_GENOTYPE_GDS",
    default = "/home/ftom/01_wd/galis/nhan_gwas/input/jrc_wrc_wgs_on_nb_genome.gds",
    must_exist = TRUE
  )
  gff_path <- phase1DataPath("gff", env_var = "LAZYGAS_PHASE1_GFF", must_exist = TRUE)
  ann_path <- phase1DataPath("ann", env_var = "LAZYGAS_PHASE1_ANN", must_exist = TRUE)
  jrc_xlsx <- phase1DataPath("phenotype_jrc", must_exist = TRUE)
  wrc_xlsx <- phase1DataPath("phenotype_wrc", must_exist = TRUE)
  snpeff_path <- phase1SnpEffGds(must_exist = TRUE)

  gds_link <- file.path(work_dir, "jrc_wrc_wgs_on_nb_genome.gds")
  if (!file.exists(gds_link)) {
    message("Linking genotype GDS into work_dir...")
    ok <- file.symlink(gds_src, gds_link)
    if (!isTRUE(ok) && !file.exists(gds_link)) {
      stop("Failed to symlink GDS into ", gds_link, call. = FALSE)
    }
  }

  companion <- file.path(work_dir, "jrc_wrc_wgs_on_nb_genome.lazygas")
  message("Opening LazyGas: ", gds_link)
  # Stored GDS filters may drop JRC samples; Phase 1 needs the full JRC+WRC panel.
  lg <- buildLazyGas(
    gds_fn = gds_link,
    load_filter = FALSE,
    companion_path = companion,
    overwrite = overwrite,
    lazygas_store = "parquet"
  )
  lg <- GBScleanR::resetSamFilter(lg)
  lg <- GBScleanR::resetMarFilter(lg)
  message(
    "Samples after filter reset: ",
    length(GBScleanR::getSamID(lg)),
    " / markers: ", sum(GBScleanR::validMar(lg))
  )
  if (isTRUE(apply_qc_filter)) {
    lg <- .phase1_apply_qc_filter(lg)
  }
  lg <- tryCatch(restorePhenoFromStore(lg, warn_stub = FALSE), error = function(e) lg)

  pheno_all <- .phase1_merge_jrc_wrc_pheno(jrc_xlsx = jrc_xlsx, wrc_xlsx = wrc_xlsx)
  if (is.null(phenos)) {
    phenos <- c("Amylose content", "Leaf width", "Heading date")
  }
  phenos <- as.character(phenos)
  miss <- setdiff(phenos, setdiff(names(pheno_all), "id"))
  if (length(miss)) {
    stop("Requested phenos not in JRC/WRC tables: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  pheno <- pheno_all[, c("id", phenos), drop = FALSE]
  message("Assigning phenotypes: ", paste(phenos, collapse = ", "))
  lg <- assignPheno(lg, pheno = pheno)

  if (isTRUE(skip_if_ready) && !isTRUE(overwrite)) {
    ready <- all(vapply(phenos, function(ph) {
      !is.null(lazyData(lg, dataset = "candidate", pheno = ph))
    }, logical(1)))
    if (ready) {
      message("Candidates already present; skipping scan/peak/recalc/candidate.")
      return(invisible(lg))
    }
  }

  message("Importing GFF (may take a while)...")
  gff <- rtracklayer::import.gff(gff_path)
  ann <- utils::read.delim(ann_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Gene_ID" %in% names(ann)) {
    if ("gene_id" %in% names(ann)) {
      ann$Gene_ID <- ann$gene_id
    } else if ("ID" %in% names(ann)) {
      ann$Gene_ID <- ann$ID
    } else {
      stop("Annotation table needs Gene_ID / gene_id / ID column.", call. = FALSE)
    }
  }
  snpeff <- open_snpeff(snpeff_path)

  message(
    "scanAssoc method=", method,
    " omit_outlier=", isTRUE(omit_outlier),
    " for: ", paste(phenos, collapse = ", ")
  )
  # scanAssoc / callPeakBlock / recalcAssoc mutate the store in place and may
  # return NULL; do not replace `lg` with their return values.
  if (!isTRUE(skip_scan)) {
    scanAssoc(object = lg, method = method, omit_outlier = isTRUE(omit_outlier))
  } else {
    message("skip_scan=TRUE; reusing existing scan results.")
  }

  if (!isTRUE(skip_peakcall)) {
    message("callPeakBlock (limit_peakcall=", limit_peakcall, ")...")
    callPeakBlock(
      object = lg,
      n_threads = n_threads,
      limit_peakcall = as.integer(limit_peakcall)
    )
  } else {
    message("skip_peakcall=TRUE; reusing existing peakcall results.")
  }

  message("recalcAssoc...")
  recalcAssoc(object = lg, n_threads = n_threads)

  message("listCandidate...")
  listCandidate(
    object = lg,
    gff = gff,
    ann = ann,
    snpeff = snpeff,
    recalc = TRUE
  )

  message("Phase 1 peaked object ready under ", work_dir)
  invisible(lg)
}

#' Phase 1 genotype QC filters (GBScleanR).
#'
#' Applies INFO filters, depth via \code{countRead}/\code{setMarFilter}, then
#' missing / MAF using VCF INFO \code{AN}/\code{AF} (equivalent thresholds to
#' \code{setMarFilter(missing = 0.2, maf = 0.05)}). \code{countGenotype()} is
#' avoided because it segfaults on this ~17M-variant GDS.
#'
#' @param object A \code{LazyGas} / \code{GbsrGenotypeData} object.
#' @return The filtered object.
#' @keywords internal
.phase1_apply_qc_filter <- function(object) {
  message("Applying Phase 1 setInfoFilter...")
  object <- GBScleanR::setInfoFilter(
    object,
    mq = 40,
    fs = 60,
    qd = 2,
    sor = 3,
    mqranksum = c(-12.5, Inf),
    readposranksum = c(-8, Inf)
  )
  message("Markers after setInfoFilter: ", sum(GBScleanR::validMar(object)))

  need_read <- inherits(
    try(GBScleanR::getCountRead(object, target = "marker", valid = TRUE),
        silent = TRUE),
    "try-error"
  )
  if (need_read) {
    message("countRead() for dp filter (may take a while)...")
    object <- GBScleanR::countRead(object)
  }
  message("Applying setMarFilter(dp = c(0, 5000))...")
  object <- GBScleanR::setMarFilter(object, dp = c(0, 5000))
  message("Markers after dp filter: ", sum(GBScleanR::validMar(object)))

  message("Applying missing/maf via INFO AN/AF (missing<=0.2, maf>=0.05)...")
  object <- .phase1_filter_missing_maf_from_info(
    object,
    missing = 0.2,
    maf = 0.05
  )
  message(
    "Samples after QC: ",
    length(GBScleanR::getSamID(object)),
    " / markers: ",
    sum(GBScleanR::validMar(object))
  )
  object
}

#' Filter markers by missing rate and MAF using INFO AN / AF.
#'
#' @param object GbsrGenotypeData / LazyGas
#' @param missing Max missing genotype rate (default 0.2)
#' @param maf Min minor allele frequency (default 0.05)
#' @keywords internal
.phase1_filter_missing_maf_from_info <- function(object,
                                                 missing = 0.2,
                                                 maf = 0.05) {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("Package 'SeqArray' is required for INFO-based missing/maf filter.",
         call. = FALSE)
  }
  n_sam <- length(GBScleanR::getSamID(object))
  cur <- GBScleanR::validMar(object)
  an <- as.numeric(SeqArray::seqGetData(object, "annotation/info/AN"))
  af <- as.numeric(SeqArray::seqGetData(object, "annotation/info/AF"))
  # Diploid: missing ≈ 1 - AN / (2 * n_samples); maf = min(AF, 1-AF)
  pass_vec <- function(an_v, af_v) {
    miss_rate <- 1 - an_v / (2 * n_sam)
    maf_val <- pmin(af_v, 1 - af_v)
    is.finite(miss_rate) & is.finite(maf_val) &
      miss_rate <= missing & maf_val >= maf
  }
  if (length(an) == length(cur) && length(af) == length(cur)) {
    keep <- pass_vec(an, af)
    GBScleanR::validMar(object) <- cur & keep
  } else if (length(an) == sum(cur) && length(af) == sum(cur)) {
    keep_valid <- pass_vec(an, af)
    cur[cur] <- keep_valid
    GBScleanR::validMar(object) <- cur
  } else {
    stop(
      "INFO AN/AF length unexpected (AN=", length(an),
      ", AF=", length(af), ", nmar=", length(cur),
      ", valid=", sum(cur), ").",
      call. = FALSE
    )
  }
  object
}

.phase1_resolve_scan_method <- function(method) {
  if (is.null(method) || !nzchar(as.character(method)[[1L]])) {
    cfg <- .ai_config_phase1_data()
    method <- cfg$scan_method %||% "mlm"
  }
  method <- tolower(as.character(method)[[1L]])
  if (identical(method, "lmm")) {
    method <- "mlm"
  }
  if (!method %in% c("glm", "mlm")) {
    stop("scan method must be 'mlm' (LMM) or 'glm' (got: ", method, ").",
         call. = FALSE)
  }
  method
}

.phase1_merge_jrc_wrc_pheno <- function(jrc_xlsx, wrc_xlsx) {
  jrc <- as.data.frame(readxl::read_excel(jrc_xlsx), stringsAsFactors = FALSE)
  wrc <- as.data.frame(readxl::read_excel(wrc_xlsx), stringsAsFactors = FALSE)

  names(wrc) <- gsub("\\s+", " ", gsub("\n", " ", names(wrc)))
  names(jrc) <- gsub("\\s+", " ", gsub("\n", " ", names(jrc)))

  jrc_id_col <- intersect(names(jrc), c("JRC_number", "JRC number"))[1]
  wrc_id_col <- intersect(names(wrc), c("WRC_number", "WRC number"))[1]
  if (is.na(jrc_id_col) || is.na(wrc_id_col)) {
    stop("Could not find JRC_number / WRC_number columns.", call. = FALSE)
  }

  jrc_traits <- c("Amylose content", "Leaf width", "Heading date")
  wrc_traits <- c(
    "Heading date", "High latitude flowering", "Grain color",
    "Amylose content", "Grain length"
  )
  all_traits <- unique(c(jrc_traits, wrc_traits))

  sanitize_id <- function(x) {
    x <- trimws(as.character(x))
    # Drop non-ASCII trailing junk sometimes present in spreadsheet cells
    x <- gsub("[^A-Za-z0-9_-]+$", "", x)
    x
  }
  jrc_df <- data.frame(
    id = sanitize_id(jrc[[jrc_id_col]]),
    stringsAsFactors = FALSE
  )
  wrc_df <- data.frame(
    id = sanitize_id(wrc[[wrc_id_col]]),
    stringsAsFactors = FALSE
  )
  for (tr in all_traits) {
    jrc_df[[tr]] <- if (tr %in% names(jrc)) {
      suppressWarnings(as.numeric(jrc[[tr]]))
    } else {
      NA_real_
    }
    wrc_df[[tr]] <- if (tr %in% names(wrc)) {
      suppressWarnings(as.numeric(wrc[[tr]]))
    } else {
      NA_real_
    }
  }

  out <- rbind(jrc_df, wrc_df)
  out <- out[!is.na(out$id) & nzchar(out$id), , drop = FALSE]
  out
}

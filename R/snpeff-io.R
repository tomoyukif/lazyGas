################################################################################
#'
#' @importFrom SeqArray seqVCF2GDS
#' @export
snpeff2gds <- function(vcf_fn, out_fn, verbose = TRUE, overwrite = FALSE){
  if (file.exists(out_fn)) {
    if (!overwrite) {
      stop(out_fn, " already exists. Set overwrite = TRUE to replace.", call. = FALSE)
    }
    file.remove(out_fn)
  }

  # Create GDS file formatted in the SeqArray style.
  fmt.import <- ""
  info.import <- "ANN"
  out_fn <- seqVCF2GDS(vcf_fn, out_fn, ignore.chr.prefix = "",
                       info.import = info.import,
                       fmt.import = fmt.import,
                       verbose = verbose)
  return(out_fn)
}

#'
#' @importFrom gdsfmt exist.gdsn delete.gdsn objdesp.gdsn setdim.gdsn add.gdsn
#' @export
open_snpeff <- function(gds_fn) {
  # Check input data type
  if(file.exists(gds_fn)){
    gds <- openfn.gds(filename = gds_fn, readonly = FALSE, allow.duplicate = TRUE)

  } else {
    stop("gds_fn must be a valid file path", call. = FALSE)
  }
  class(gds) <- c(class(gds), "snpeff_gds")
  return(gds)
}

#' Phase 1 evaluation SnpEff GDS path (JRC/WRC, pre-annotated)
#'
#' Returns the existing annotated GDS used for SnpEff impact scoring in Phase 1.
#' Does **not** re-run SnpEff on the raw VCF.
#'
#' Resolution order:
#' 1. \code{LAZYGAS_PHASE1_SNPEFF_GDS}
#' 2. \code{snpeff_gds} in \code{inst/config/ai-lazygas-data.yaml}
#' 3. Built-in default under \code{/home/ftom/01_wd/galis/...}
#'
#' @param must_exist If \code{TRUE} (default), stop when the file is missing.
#' @return Character scalar path.
#' @export
#' @seealso [open_snpeff()], [snpeff2gds()]
phase1SnpEffGds <- function(must_exist = TRUE) {
  phase1DataPath(
    key = "snpeff_gds",
    env_var = "LAZYGAS_PHASE1_SNPEFF_GDS",
    default = "/home/ftom/01_wd/galis/nhan_gwas/input/jrc_wrc_wgs_on_nb_genome.snpeff.gds",
    must_exist = must_exist
  )
}

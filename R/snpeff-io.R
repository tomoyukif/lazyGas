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

################################################################################
# Define the LazyGas class object
#' Class `LazyGas`
#'
#' The `LazyGas` class is the main class of `lazyGas` and
#' user work with this class object.
#'
#' @details

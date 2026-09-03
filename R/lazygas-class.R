#' The `LazyGas` class is an extention of `GbsrGenotypeData` in the
#' [GBScleanR] package.
#'
#' @importClassesFrom GBScleanR GbsrGenotypeData
#' @aliases  LazyGas-class LazyGas
#' @slot lazydata A [list] object (phenotype data).
#' @slot store Companion storage configuration (`mode`, `path`, `gds_fn`).
#'
#' @examples
#' # `buildLazyGas()` initialize the `LazyGas` object.
#'
#' # Load a GDS file and instantiate a `LazyGas` object.
#' gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")
#' lgas <- buildLazyGas(gds_fn)
#'
#' # Close connection to the GDS file.
#' closeGDS(lgas)
#'
#' @exportClass LazyGas
#' @importFrom methods setClass slot setAs setMethod S3Part
#' @import GBScleanR
#'
setClass(
  Class = "LazyGas",
  contains = "GbsrGenotypeData",
  slots = c(lazydata = "list", store = "list"),
  prototype = list(lazydata = list(), store = list())
)

# SeqArray/GBScleanR hybrid: LazyGas extends SeqVarGDSClass + gds.class (oldClass).
# Under devtools::load_all(), setClass regenerates LazyGas→SeqVarGDSClass coerce as a
# hollow SeqVarGDSClass (copies only .S3Class). SeqArray validity then calls
# index.gdsn(object) → object$root and fails with "$ operator not defined for this
# S4 class". Installed packages get S3Part-based coerce instead. Force S3Part here
# so load_all matches install behavior for new()/validObject().
setAs("LazyGas", "SeqVarGDSClass", function(from) {
  methods::S3Part(from, strictS3 = TRUE)
})
setAs("LazyGas", "gds.class", function(from) {
  methods::S3Part(from, strictS3 = TRUE)
})

################################################################################
# Inherited methods
#' @importMethodsFrom GBScleanR closeGDS
setMethod("closeGDS",
          "LazyGas",
          function(object, verbose){
            object <- .store_close(object)
            closefn.gds(gdsfile = object)
            if(verbose){
              message('The connection to the GDS file was closed.')
            }
            invisible(object)
          }
)

################################################################################
# Functions to be exported

#' Build a LazyGas object
#'
#' @param gds_fn A path to a GDS file
#' @param load_filter A logical value to indicate whether apply filrering stored in the input GDS file to samples and markers.
#' @param overwrite If `TRUE`, replace existing lazyGas results in the companion store (or GDS `lazygas/` folder when `lazygas_store = "gds"`).
#' @param lazygas_store Where to store association results: `"parquet"` (default companion folder `{stem}.lazygas/`), `"sqlite"`, `"gds"` (legacy GDS subtree), or `"auto"`.
#' @param companion_path Optional path to an existing companion store when it
#'   does not sit next to \code{gds_fn} (e.g. Shiny file uploads). See
#'   [resolveLazyGasPaths()].
#' @param create_gds A named list to create a new GDS file with specified genotype and marker information. See the Details section.
#'
#' @details
#' As the default, the `buildLazyGas()` function loads the genotype and marker
#' information from a GDS file specified to the `gds_fn` argument.
#' If a named list storing genotype and marker information was specified to the
#' `create_gds` argument, a GDS file will be created with supplied genotype and
#' marker information. the `create_gds` list should have the following elements.
#' - genotype:
#' - sample.id:
#' - snp.id:
#' - snp.rs.id:
#' - snp.chromosome:
#' - snp.position:
#' - snp.allele:
#' - haplotype:
#' - dosage:
#' Either genotype, haplotype, or dosage should
#' be supplied in the `create_gds` list. If `sample.id`, `snp.id`, and
#' `snp.rs.id` were left as `NULL`, serial numbers will be assigned as sample
#' and SNP IDs. Random alleles will be assigned if `snp.allele = NULL`.
#' `snp.chromosome` and `snp.marker` must be supplied.
#'
#' @importFrom GBScleanR loadGDS countGenotype
#' @importClassesFrom GBScleanR GbsrGenotypeData
#' @importFrom gdsfmt closefn.gds exist.gdsn objdesp.gdsn
#'
#' @export
#'
buildLazyGas <- function(gds_fn = "",
                         load_filter = TRUE,
                         overwrite = FALSE,
                         create_gds = NULL,
                         lazygas_store = c("parquet", "sqlite", "gds", "auto"),
                         companion_path = NULL){
  if(!all(sapply(X = create_gds, FUN = is.null))){
    raw_gds <- .createGDS(gds_fn = gds_fn, create_gds = create_gds)
    closefn.gds(gdsfile = raw_gds)
    out <- loadGDS(x = gds_fn, load_filter = load_filter, verbose = FALSE)

  } else {
    # SeqArray connection required by GBScleanR accessors (e.g. getSamID)
    out <- loadGDS(x = gds_fn, load_filter = load_filter, verbose = FALSE)
  }

  lazygas_store <- match.arg(lazygas_store)

  if (lazygas_store == "gds") {
    if (exist.gdsn(node = out$root, path = "lazygas")) {
      if (overwrite) {
        .create_gdsn(root_node = out$root, target_node = "",
                     new_node = "lazygas", is_folder = TRUE)
      } else {
        message(
          "Some LazyGas data has already been recorded in the input GDS file.",
          "\nCreate a LazyGas object with the recorded LazyGas data."
        )
      }
    } else {
      .create_gdsn(root_node = out$root, target_node = "",
                   new_node = "lazygas", is_folder = TRUE)
    }
  }

  out <- new(Class = "LazyGas", out, lazydata = list(), store = list())
  out <- .store_init(
    object = out,
    gds_fn = gds_fn,
    store = lazygas_store,
    overwrite = overwrite,
    companion_path = companion_path
  )
  return(out)
}

#' @importFrom SNPRelate snpgdsCreateGeno
#' @importFrom SeqArray seqSNP2GDS
#' @importFrom gdsfmt openfn.gds closefn.gds addfolder.gdsn add.gdsn index.gdsn createfn.gds put.attr.gdsn
.createGDS <- function(gds_fn, create_gds){
  check <- list(genotype = TRUE,
                haplotype = TRUE,
                dosage = TRUE)
  check$genotype <- is.null(create_gds$genotype)
  check$haplotype <- is.null(create_gds$haplotype)
  check$dosage <- is.null(create_gds$dosage)

  if(check$genotype & check$haplotype & check$dosage){
    stop("Either genotype, haplotype, or dosage should be supplied at least.",
         call. = FALSE)

  } else {
    if(!check$genotype){
      n_sample <- nrow(create_gds$genotype)
      n_snp <- ncol(create_gds$genotype)

    } else if(!check$dosage){
      n_sample <- nrow(create_gds$dosage)
      n_snp <- ncol(create_gds$dosage)

    } else {
      if(length(dim(create_gds$haplotype)) == 2){
        n_sample <- nrow(create_gds$haplotype)
        n_snp <- ncol(create_gds$haplotype)

      } else {
        n_sample <- dim(create_gds$haplotype)[2]
        n_snp <- dim(create_gds$haplotype)[3]
      }
    }
  }

  if(is.null(create_gds$sample.id)){
    message("No sample ID was supplied. ",
            "Assign serial numbers to samples as sample IDs")
    create_gds$sample.id <- seq_len(n_sample)
  }

  if(is.null(create_gds$snp.id)){
    message("No SNP ID was supplied. ",
            "Assign serial numbers to SNPs as SNP IDs")
    create_gds$snp.id <- seq_len(n_snp)
  }

  if(is.null(create_gds$snp.rs.id)){
    message("No SNP RS ID was supplied. ",
            "Assign serial numbers to SNPs as SNP RS IDs")
    create_gds$snp.rs.id <- seq_len(n_snp)
  }

  if(is.null(create_gds$snp.chromosome)){
    stop("Chromosomes in which SNPs locate should be supplied as snp.chromosome.",
         call. = FALSE)
  }

  if(is.null(create_gds$snp.position)){
    stop("SNP positions should be supplied as snp.position.",
         call. = FALSE)
  }

  if(is.null(create_gds$snp.allele)){
    message("No SNP allele was supplied. ",
            "Assign random alleles to SNPs.")
    nuc <- c("A", "T", "G", "C")
    allele_comb <- expand.grid(nuc, nuc)
    allele_comb <- subset(allele_comb, subset = Var1 != Var2)
    allele_comb <- apply(X = allele_comb, MARGIN = 1, FUN = paste,
                         collapse = "/")
    create_gds$snp.allele <- sample(x = allele_comb,
                                    size = n_snp,
                                    replace = TRUE)
  } else {
    create_gds$snp.allele <- sub(",", "/",  create_gds$snp.allele)
  }

  gds <- createfn.gds(filename = gds_fn)

  put.attr.gdsn(node = index.gdsn(gds, ""),
                name = "FileFormat", val = "SEQ_ARRAY")
  put.attr.gdsn(node = index.gdsn(gds, ""),
                name = "FileVersion", val = "v1.0")

  addfolder.gdsn(node = index.gdsn(gds, ""), name = "description")

  .create_gdsn(root_node = gds$root,
               target_node = "",
               new_node = "sample.id",
               val = create_gds$sample.id,
               storage = "string")

  .create_gdsn(root_node = gds$root,
               target_node = "",
               new_node = "variant.id",
               val = create_gds$snp.id,
               storage = "int32")

  .create_gdsn(root_node = gds$root,
               target_node = "",
               new_node = "position",
               val = create_gds$snp.position,
               storage = "int32")

  .create_gdsn(root_node = gds$root,
               target_node = "",
               new_node = "chromosome",
               val = create_gds$snp.chromosome,
               storage = "string")

  .create_gdsn(root_node = gds$root,
               target_node = "",
               new_node = "allele",
               val = create_gds$snp.allele,
               storage = "string")

  if(!check$genotype){
    addfolder.gdsn(node = index.gdsn(gds, ""), name = "genotype")
    genotype <- array(0L, dim = c(2L, dim(create_gds$genotype)))
    gc(); gc()
    genotype[1L, , ] <- as.integer(create_gds$genotype != 0L)
    gc(); gc()
    genotype[2L, , ] <- as.integer(create_gds$genotype == 2L)
    gc(); gc()
    .create_gdsn(root_node = gds,
                 target_node = "genotype",
                 new_node = "data",
                 val = genotype,
                 storage = "bit2")
  }

  # GBScleanR::loadGDS() runs seqOptimize when genotype/~data is absent;
  # SeqArray requires a phase/data node (sample x SNP, unphased = 0).
  if (!exist.gdsn(node = index.gdsn(gds, ""), path = "phase/data")) {
    if (!exist.gdsn(node = index.gdsn(gds, ""), path = "phase")) {
      addfolder.gdsn(node = index.gdsn(gds, ""), name = "phase")
    }
    phase <- matrix(0L, nrow = n_sample, ncol = n_snp)
    .create_gdsn(root_node = gds,
                 target_node = "phase",
                 new_node = "data",
                 val = phase,
                 storage = "bit1")
  }

  if(!check$haplotype){
    create_gds$haplotype[is.na(create_gds$haplotype)] <- 63L
    gc(); gc()
    addfolder.gdsn(node = index.gdsn(gds, ""), name = "annotation")
    addfolder.gdsn(node = index.gdsn(gds, "annotation"), name = "format")
    addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "HAP")
    .create_gdsn(root_node = gds,
                 target_node = "annotation/format/HAP",
                 new_node = "data",
                 val = create_gds$haplotype,
                 storage = "bit6")
  }

  if(!check$dosage){
    create_gds$dosage[is.na(create_gds$dosage)] <- 63L
    gc(); gc()
    if (!exist.gdsn(node = index.gdsn(gds, ""), path = "annotation")) {
      addfolder.gdsn(node = index.gdsn(gds, ""), name = "annotation")
    }
    if (!exist.gdsn(node = index.gdsn(gds, ""), path = "annotation/format")) {
      addfolder.gdsn(node = index.gdsn(gds, "annotation"), name = "format")
    }
    addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "EDS")
    .create_gdsn(root_node = gds,
                 target_node = "annotation/format/EDS",
                 new_node = "data",
                 val = create_gds$dosage,
                 storage = "bit6")
  }

  return(gds)
}

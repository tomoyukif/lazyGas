################################################################################
# Utility functions to handle a GDS file
## Create a GDS node
#' @importFrom gdsfmt exist.gdsn addfolder.gdsn index.gdsn add.gdsn
.create_gdsn <- function(root_node,
                         target_node,
                         new_node,
                         val = NULL,
                         storage = "string32",
                         valdim = dim(val),
                         replace = TRUE,
                         attr = NULL,
                         is_folder = FALSE){

  # Check if the target node exists in the root node
  check <- exist.gdsn(node = root_node, path = target_node)
  if(!check){ stop(target_node, " does not exist!") }

  # Get the index of the target node
  target_node_index <- index.gdsn(node = root_node, path = target_node)

  # Create a folder node if specified
  if(is_folder){
    new_node_index <- addfolder.gdsn(node = target_node_index,
                                     name = new_node,
                                     replace = TRUE)
  } else {
    # Create a new node with the specified value, storage type, and dimensions
    new_node_index <- add.gdsn(node = target_node_index, name = new_node,
                               val = val, storage = storage, valdim = valdim,
                               replace = replace)
    # Assign attributes to the new node if provided
    .putAttrGdsn(node = new_node_index, attr = attr)
    # Compress the new node
    .compressNodes(node = new_node_index)
  }
  invisible(new_node_index)
}

# Assign attributes to a GDS node
#' @importFrom gdsfmt put.attr.gdsn
.putAttrGdsn <- function(node, attr){
  # Check if attributes are provided
  if(!is.null(attr)){
    # Loop through each attribute and assign it to the node
    for(i in seq_along(attr)){
      put.attr.gdsn(node = node,
                    name = names(attr)[i], val = attr[[i]])
    }
  }
}

# Compress GDS nodes
#' @importFrom gdsfmt index.gdsn compression.gdsn readmode.gdsn
.compressNodes <- function(node){
  # Loop through each node
  for(i in seq_along(node)){
    # Compress the node using ZIP_RA compression
    compression.gdsn(node = node, compress = "ZIP_RA")
    # Set the node to read mode
    readmode.gdsn(node = node)
  }
}

## Retrieve data from a GDS node
#' @importFrom gdsfmt index.gdsn read.gdsn readex.gdsn
.get_data <- function(object, node, start = NULL, count, sel = NULL){
  # Check if selection criteria (sel) is provided
  if(is.null(sel)){
    # If selection criteria is not provided
    if(is.null(start)){
      # If start is not provided, read the entire node
      out <- read.gdsn(node = index.gdsn(node = object$root, path = node))
    } else {
      # If start is provided, read the node from the specified start position and count
      out <- read.gdsn(node = index.gdsn(node = object$root, path = node),
                       start = start,
                       count = count)
    }
  } else {
    # If selection criteria is provided, read the node with the selection criteria
    out <- readex.gdsn(node = index.gdsn(node = object$root, path = node), sel = sel)
  }
  # Return the retrieved data
  return(out)
}

# Retrieve a scan result for a phenotype
#' @importFrom gdsfmt get.attr.gdsn objdesp.gdsn read.gdsn readex.gdsn index.gdsn
.get_scan <- function(object, pheno_name, sel = NULL){
  # Construct the scan node path for the specified phenotype
  scan_node <- paste0("lazygas/scan/", pheno_name)

  # Get the attribute names from the scan node
  col_names <- get.attr.gdsn(node = index.gdsn(node = object$root, path = scan_node))
  col_names <- col_names$colnames

  # Check if selection criteria (sel) is provided
  if(is.null(sel)){
    # If selection criteria is not provided, read the entire scan node
    out <- data.frame(getMarID(object = object),
                      getChromosome(object = object),
                      getPosition(object = object),
                      read.gdsn(node = index.gdsn(node = object$root, path = scan_node)))
    # Set the column names for the output data frame
    colnames(out) <- c("variant_ID", "Chr", "Pos", col_names)
  } else {
    # If selection criteria is provided, read the scan node with the selection criteria
    dim <- objdesp.gdsn(node = index.gdsn(node = object$root, path = scan_node))
    sel <- list(rep(TRUE, dim$dim[1]), col_names %in% sel)
    out <- data.frame(getMarID(object = object),
                      getChromosome(object = object),
                      getPosition(object = object),
                      readex.gdsn(node = index.gdsn(node = object$root, path = scan_node), sel = sel))
    # Set the column names for the output data frame based on the selection criteria
    colnames(out) <- c("variant_ID", "Chr", "Pos", col_names[sel[[2]]])
  }

  # Return the retrieved scan data
  return(out)
}


# Retrieve a peakcall result for a phenotype
#' @importFrom dplyr right_join left_join
#'
.get_peakcall <- function(object, pheno_name, recalc = FALSE){
  # Determine the paths for peaks and blocks based on whether recalculation is required
  if(recalc){
    peaks_gdsn <- "lazygas/recalc/peaks"
    blocks_gdsn <- "lazygas/recalc/blocks"
    peaks_col_names <- get.attr.gdsn(node = index.gdsn(node = object$root, path =  paste0(peaks_gdsn, "/", pheno_name)))
    blocks_col_names <- get.attr.gdsn(node = index.gdsn(node = object$root, path =  paste0(blocks_gdsn, "/", pheno_name)))

  } else {
    peaks_gdsn <- "lazygas/peakcall/peaks"
    blocks_gdsn <- "lazygas/peakcall/blocks"
    peaks_col_names <- get.attr.gdsn(node = index.gdsn(node = object$root, path = peaks_gdsn))
    blocks_col_names <- get.attr.gdsn(node = index.gdsn(node = object$root, path = blocks_gdsn))
  }

  # Retrieve the peaks data for the specified phenotype
  peaks <- .get_data(object = object, node = paste0(peaks_gdsn, "/", pheno_name))
  if(length(peaks) == 0){
    return(NULL)
  }

  # Ensure the peaks data is in a data frame format
  if(!is.matrix(peaks)){
    peaks <- matrix(peaks, nrow = 1)
    peaks <- data.frame(peaks)
  } else {
    if(recalc){
      peaks <- data.frame(peaks)
    } else {
      peaks <- data.frame(t(peaks))
    }
  }
  names(peaks) <- peaks_col_names$col_names

  # Retrieve the blocks data for the specified phenotype
  blocks <- .get_data(object = object, node = paste0(blocks_gdsn, "/", pheno_name))
  if(!is.matrix(blocks)){
    blocks <- matrix(blocks, nrow = 1)
    blocks <- data.frame(blocks)

  } else {
    if(recalc){
      blocks <- data.frame(blocks)

    } else {
      blocks <- data.frame(t(blocks))
    }
  }
  names(blocks) <-blocks_col_names$col_names

  # Add variant_ID to peaks data
  peaks$variant_ID <- peaks$peak_variant_ID

  # Retrieve the scan results for the specified phenotype
  pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))

  if(recalc){
    # For recalculated data, join peaks and pvalues data
    pvalues <- subset(pvalues, select = variant_ID:Pos)
    peaks_join <- left_join(x = peaks, y = pvalues, "variant_ID")
    names(peaks_join)[-(1:3)] <- paste0("peak_", names(peaks_join)[-(1:3)])
    peaks_join <- subset(peaks_join, select = -variant_ID)

    # Join blocks and pvalues data
    blocks_join <- left_join(x = blocks, y = pvalues, by = "variant_ID")
    hit <- match(peaks_join$peak_variant_ID, blocks_join$variant_ID)
    peaks_join$peak_negLog10P <- blocks_join$negLog10P[hit]

    # Right join peaks and blocks data
    out <- right_join(x = peaks_join, y = blocks_join, by = "peak_ID")
  } else {
    # For non-recalculated data, join peaks and pvalues data
    peaks_join <- left_join(x = peaks, y = pvalues, "variant_ID")
    names(peaks_join)[-(1:3)] <- paste0("peak_", names(peaks_join)[-(1:3)])
    peaks_join <- subset(peaks_join, select = -c(variant_ID, peak_FDR))

    # Join blocks and pvalues data
    blocks_join <- left_join(x = blocks, y = pvalues, by = "variant_ID")

    # Right join peaks and blocks data
    out <- right_join(x = peaks_join, y = blocks_join, by = "peak_ID")
  }

  # Return the final joined data
  return(out)
}


################################################################################
# Define the LazyGas class object
#' Class `LazyGas`
#'
#' The `LazyGas` class is the main class of [lazyGas] and
#' user work with this class object.
#'
#' @details
#' The `LazyGas` class is an extention of `GbsrGenotypeData` in the
#' [GBScleanR] package.
#'
#' @importClassesFrom GBScleanR GbsrGenotypeData
#' @aliases  LazyGas-class LazyGas
#' @slot lazydata A [list] object.
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
#' @importFrom methods setClass slot
#' @import GBScleanR
#'
setClass(Class = "LazyGas",
         contains = "GbsrGenotypeData",
         slots = c(lazydata = "list"))

################################################################################
# Inherited methods
#' @importMethodsFrom GBScleanR closeGDS
setMethod("closeGDS",
          "LazyGas",
          function(object, verbose){
            closefn.gds(gdsfile = object)
            if(verbose){
              message('The connection to the GDS file was closed.')
            }
          }
)

################################################################################
# Functions to be exported

#' Build a LazyGas object
#'
#' @param gds_fn A path to a GDS file
#' @param load_filter A logical value to indicate whether apply filrering stored in the input GDS file to samples and markers.
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
#' @importFrom gdsfmt  exist.gdsn
#'
#' @export
#'
buildLazyGas <- function(gds_fn = "",
                         load_filter = TRUE,
                         overwrite = FALSE,
                         create_gds = NULL){
  if(!all(sapply(X = create_gds, FUN = is.null))){
    .createGbsrGDS(gds_fn = gds_fn, create_gds = create_gds)
  }

  # Load the GDS file and apply the load filter if specified
  out <- loadGDS(x = gds_fn, load_filter = load_filter)

  # Count the genotypes in the GDS object
  out <- countGenotype(object = out)

  # Check if the "lazygas" node already exists in the GDS file
  if(exist.gdsn(node = out$root, path = "lazygas")){
    if(overwrite){
      # If overwrite is TRUE, create (or replace) the "lazygas" folder node
      .create_gdsn(root_node = out$root, target_node = "",
                   new_node = "lazygas", is_folder = TRUE)
    } else {
      # If overwrite is FALSE, notify the user that LazyGas data already exists
      message("Some LazyGas data has already been recorded in the input GDS file.",
              "\nCreate a LazyGas object with the recorded LazyGas data.")
    }
  } else {
    # If the "lazygas" node does not exist, create it
    .create_gdsn(root_node = out$root, target_node = "",
                 new_node = "lazygas", is_folder = TRUE)
  }

  # Create a new LazyGas object and return it
  out <- new(Class = "LazyGas", out)
  return(out)
}

#' @importFrom SNPRelate snpgdsCreateGeno
#' @importFrom SeqArray seqSNP2GDS
#' @importFrom gdsfmt openfn.gds closefn.gds addfolder.gdsn add.gdsn index.gdsn
.createGbsrGDS <- function(gds_fn, create_gds){
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
        n_sample <- nrow(create_gds$dosage)
        n_snp <- ncol(create_gds$dosage)

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

  create_gds$genotype[is.na(create_gds$genotype)] <- 3
  snpgdsCreateGeno(gds.fn = paste0(gds_fn, ".temp"),
                   genmat = create_gds$genotype,
                   sample.id = create_gds$sample.id,
                   snp.id = create_gds$snp.id,
                   snp.rs.id = create_gds$snp.rs.id,
                   snp.chromosome = create_gds$snp.chromosome,
                   snp.position = create_gds$snp.position,
                   snp.allele = create_gds$snp.allele,
                   snpfirstdim = FALSE,
                   compress.annotation = "LZMA_RA",
                   compress.geno = "LZMA_RA")

  seqSNP2GDS(gds.fn = paste0(gds_fn, ".temp"), out.fn = gds_fn)
  unlink(x = paste0(gds_fn, ".temp"), force = TRUE)

  gds <- openfn.gds(gds_fn, readonly = FALSE)

  if(!check$haplotype){
    create_gds$haplotype[is.na(create_gds$haplotype)] <- 63
    if(length(dim(create_gds$haplotype)) == 2){
      rep_row <- rep(seq_len(nrow(create_gds$haplotype)), each = 2)
      create_gds$haplotype <- create_gds$haplotype[rep_row, ]
      create_gds$haplotype <- array(data = create_gds$haplotype,
                                    dim = c(2,
                                            nrow(create_gds$haplotype) / 2,
                                            ncol(create_gds$haplotype)))

    }
    addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "HAP")
    add.gdsn(index.gdsn(gds, "annotation/format/HAP"), name = "data",
             val = create_gds$haplotype,
             storage = "bit6",
             compress = "LZMA_RA")
  }

  if(!check$dosage){
    create_gds$dosage[is.na(create_gds$dosage)] <- 63
    addfolder.gdsn(node = index.gdsn(gds, "annotation/format"), name = "EDS")
    add.gdsn(index.gdsn(gds, "annotation/format/EDS"), name = "data",
             val = create_gds$dosage,
             storage = "bit6",
             compress = "LZMA_RA")
  }

  closefn.gds(gds)
}

#' Assign phenotype data to a LazyGas object
#'
#' @param object A LazyGas object
#' @param pheno A data.frame of phenotype data that must contain a sample id
#' column named `id` or `ID`
#'
#' @importMethodsFrom GBScleanR getSamID
#'
#' @export
#'
setGeneric("assignPheno", function(object, pheno, rename = NULL, ...)
  standardGeneric("assignPheno"))

setMethod("assignPheno",
          "LazyGas",
          function(object, pheno, rename){
            # Get index for either 'ID' or 'id' column
            sampleID_pheno_index <- grep("^ID$|^id$", colnames(pheno))

            # Extract phenotype names
            pheno_names <- colnames(pheno)[-sampleID_pheno_index]

            # Check if the sample ID column exists and is valid
            if(length(sampleID_pheno_index) != 1){
              stop("The input `pheno` data.frame has an invalid structure.",
                   "\n It must contain a sample id column named `id` or `ID`.")
            }

            # Extract sample IDs from the phenotype data
            sampleID_pheno <- pheno[, sampleID_pheno_index]

            # Extract sample IDs from the GDS object
            sampleID_geno <- getSamID(object = object)

            # Identify samples with missing phenotype or genotype data
            nopheno <- sampleID_geno[!sampleID_geno %in% sampleID_pheno]
            nogeno <- sampleID_pheno[!sampleID_pheno %in% sampleID_geno]

            # Report samples with missing phenotype or genotype data
            if(length(nopheno) > 0){
              message("The following samples have no phenotype info: \n",
                      paste(nopheno, collapse = " "))
            }
            if(length(nogeno) > 0){
              message("The following samples have no genotype info: \n",
                      paste(nogeno, collapse = " "))
            }

            # Reorder phenotype data to match genotype data
            pheno <- subset(pheno, subset = !sampleID_pheno %in% nogeno)
            sampleID_pheno <- sampleID_pheno[!sampleID_pheno %in% nogeno]
            id_match <- match(sampleID_geno, sampleID_pheno)
            pheno <- pheno[id_match, -sampleID_pheno_index]

            # Convert to data frame if pheno is a vector
            if(is.vector(pheno)){
              pheno <- data.frame(pheno)
              colnames(pheno) <- pheno_names
            }

            # Rename phenotypes
            if(!is.null(rename)){
              colnames(pheno) <- rename
            }

            # Identify phenotypes with less than 3 observations
            invalid_pheno <- sapply(lapply(pheno, is.na), sum)
            invalid_pheno <- invalid_pheno > nrow(pheno) - 3

            if(sum(invalid_pheno) > 0){
              message("The following phenotype(s) have less than 3 observations: \n",
                      paste(pheno_names[invalid_pheno], collapse = " "))
              message("These phenotype(s) were omitted.")
              pheno <- subset(pheno, select = !invalid_pheno)
            }

            # Check phenotype type
            pheno_type <- .checkPhenoType(pheno = pheno)

            # Add phenotype data to the GDS object
            object@lazydata$pheno <- pheno
            object@lazydata$pheno_type <- pheno_type
            object@lazydata$pheno_names <- colnames(pheno)

            return(object)
          }
)

# Function to check the phenotype type
.checkPhenoType <- function(pheno) {
  # Determine the unique levels for each phenotype
  pheno_levels <- lapply(pheno, unique)

  # Remove any NA values from the list of phenotype levels
  pheno_levels <- lapply(pheno_levels, na.omit)

  # Count the number of unique levels for each phenotype
  n_levels <- sapply(pheno_levels, length)

  # Determine if the phenotype is binary (i.e., has exactly 2 unique levels)
  binary <- n_levels == 2

  # Perform the Shapiro-Wilk test for normality on each phenotype
  pheno_normality <- lapply(pheno, shapiro.test)

  # Extract the p-value from each Shapiro-Wilk test result
  pheno_normality <- sapply(pheno_normality, function(x) x[[2]])

  # Set the normality test result to NA for binary phenotypes
  pheno_normality[binary] <- NA

  # Create a data frame to store the binary status and normality p-values
  out <- data.frame(binary = binary, normality = pheno_normality)

  # Return the resulting data frame
  return(out)
}

#' Get Phenotype Data from LazyGas Object
#'
#' @param object An object of class \code{LazyGas}.
#' @return A list containing phenotype data from the \code{LazyGas} object.
#' @export
setGeneric("getPheno", function(object, ...)
  standardGeneric("getPheno"))

#' @rdname getPheno
setMethod("getPheno",
          "LazyGas",
          function(object){
            return(object@lazydata[grepl("pheno", names(object@lazydata))])
          }
)

#' Assign phenotype data to the samples in the LazyGas object
#'
#' @param object A LazyGas object
#' @param pheno A data.frame of phenotype data that must contain a sample id
#' column named `id` or `ID`
#'
#' @importFrom cowplot plot_grid
#'
#' @export
#'
setGeneric("plotPheno", function(object,
                                 pheno = 1,
                                 xlab = "phenotype",
                                 axis_title_size = 14,
                                 axis_text_size = 12,
                                 fill = "skyblue",
                                 color = "darkblue",
                                 boxplot = TRUE,
                                 ...)
  standardGeneric("plotPheno"))

setMethod("plotPheno",
          "LazyGas",
          function(object,
                   pheno,
                   xlab,
                   axis_title_size,
                   axis_text_size,
                   fill,
                   color,
                   boxplot){
            # Check if pheno is not NULL and process accordingly
            if(!is.null(pheno)){
              if(is.numeric(pheno)){
                # If pheno is numeric, use it directly as an index
                phe_index <- pheno
              } else if(is.logical(pheno)){
                # If pheno is logical, find the indices where it is TRUE
                phe_index <- which(pheno)
              } else if(is.character(pheno)){
                # If pheno is a character, find the matching phenotype names in the object
                phe_index <- which(object@lazydata$pheno_names %in% pheno)
              }

              # If multiple phenotypes are selected, use the first one and notify the user
              if(length(phe_index) > 1){
                message("Multiple phenotypes were selected.",
                        "\nUse the first phenotype.")
              }
            } else {
              # If pheno is NULL, use the first phenotype and notify the user
              message("Use the first phenotype.")
              phe_index <- 1
            }

            # Create a data frame with the selected phenotype values
            df <- data.frame(value = object@lazydata$pheno[, phe_index])

            # Create a histogram plot of the phenotype values
            p1 <- ggplot(df, aes(x = value)) +
              geom_histogram(fill = fill, color = color) +
              ylab('Count') +
              theme(axis.title.y = element_text(size = axis_title_size),
                    axis.text.y = element_text(size = axis_text_size),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank())

            # Create a boxplot of the phenotype values
            p2 <- ggplot(df, aes(x = value)) +
              geom_boxplot(fill = fill, color = color) +
              xlab(xlab) +
              theme(axis.title.x = element_text(size = axis_title_size),
                    axis.text.x = element_text(size = axis_text_size),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    panel.grid = element_blank())

            # Combine the histogram and boxplot into a single plot
            if(boxplot){
              p <- plot_grid(p1, p2, ncol = 1,
                             rel_heights = c(3, 1),
                             align = 'v', axis = 'lr')
            } else {
              p <- p1
            }

            # Return the combined plot
            return(p)
          }
)

#'
#' @importFrom methods show
#' @importFrom GBScleanR nsam nmar
#'
setMethod("show",
          "LazyGas",
          function(object){
            # Display the number of samples
            message("Number of samples")
            print(nsam(object = object))

            # Display the number of markers
            message("Number of markers")
            print(nmar(object = object))

            # Display phenotype information
            message("Phenotype")
            if(is.null(object@lazydata$pheno_names)){
              # If phenotype names are not assigned, print a message
              print("Not assigned")
            } else {
              # Print the assigned phenotype names
              print(object@lazydata$pheno_names)

              # Display binary phenotype information
              message("Binary phenotype")
              if(any(object@lazydata$pheno_type$binary)){
                # If there are binary phenotypes, print their names
                print(object@lazydata$pheno_names[object@lazydata$pheno_type$binary])
              } else {
                # If no binary phenotypes are assigned, print a message
                print("Not assigned")
              }

              # Display non-normally distributing phenotype information
              message("Non-normally distributing phenotype")
              check <- object@lazydata$pheno_type$normality <= 0.05
              check[is.na(check)] <- FALSE
              if(any(check)){
                # If there are non-normally distributing phenotypes, print their names
                print(object@lazydata$pheno_names[check])
              } else {
                # If no non-normally distributing phenotypes are assigned, print a message
                print("Not assigned")
              }
            }
          }
)

################################################################################
#' Confirm data structure of genotype data per marker
#' @export
#'
setGeneric("getGenoPerMarker", function(object, geno_format, marker_id, ...)
  standardGeneric("getGenoPerMarker"))

setMethod("getGenoPerMarker",
          "LazyGas",
          function(object,
                   geno_format = c("genotype", "corrected", "dosage", "haplotype"),
                   marker_index = 1){

            geno_format <- match.arg(arg = geno_format, c("genotype", "corrected", "dosage", "haplotype"))

            path <- switch(geno_format,
                           "genotype" = "genotype/data",
                           "corrected" = "annotation/format/CGT/data",
                           "dosage" = "annotation/format/EDS/data",
                           "haplotype" = "annotation/format/HAP/data")

            if(geno_format == "dosage"){
              selection <- list(validSam(object = object),
                                validMar(object = object))
              selection[[2]][selection[[2]]][-marker_index] <- FALSE

            } else {
              # Get the index of the target node in the GDS object
              target_node_index <- index.gdsn(node = object, path = path)
              obj_desp <- objdesp.gdsn(node = target_node_index)
              selection <- list(rep(TRUE, obj_desp$dim[1]),
                                validSam(object = object),
                                validMar(object = object))
              selection[[3]][selection[[3]]][-marker_index] <- FALSE
            }

            out <- .get_data(object = object, node = path, sel = selection)
            return(out)
          }
)

################################################################################
#' Make a conversion function
#' @export
#'
makeConvFun <-  function(geno_format = c("genotype", "corrected", "dosage", "haplotype"),
                         n_levels = 3,
                         as_categorical = FALSE){

  geno_format <- match.arg(arg = geno_format, c("genotype", "corrected", "dosage", "haplotype"))

  if(geno_format == "dosage"){
    out <- .makeConvFunDosage(n_levels = n_levels,
                              as_categorical = as_categorical)

  } else if(geno_format %in% c("genotype", "corrected")){
    out <- .makeConvFunGenotype(as_categorical = as_categorical)

  } else {
    out <- .makeConvFunHaplotype(n_levels = n_levels)
  }

  return(out)
}

.makeConvFunDosage <- function(n_levels, as_categorical){
  if(as_categorical){
    out <- .categoricalConvFun(n_levels = n_levels, prefix = "plex")

  } else {
    out <- .numericalConvFun(n_levels = n_levels)
  }
  return(out)
}

.makeConvFunGenotype <- function(as_categorical){
  if(as_categorical){
    out <- .categoricalConvFun(n_levels = 3, prefix = "plex")

  } else {
    out <- .numericalConvFun(n_levels = 3)
  }
  return(out)
}


.makeConvFunHaplotype <- function(n_levels){
  out <- .categoricalConvFun(n_levels = n_levels, prefix = "hap")
  return(out)
}

.categoricalConvFun <- function(n_levels, prefix){
  out <- "function(g){"
  for(i in seq_len(n_levels)){
    if(i == 2){
      out <- paste0(out,
                    paste0(prefix, i, " <- colSums(g == ", i, ")"))

    } else if(i > 2){
      out <- paste(out,
                   paste0(prefix, i, " <- colSums(g == ", i, ")"),
                   sep = ";")
    }
  }

  out <- paste0(out,
                paste0("; out <- data.frame(", prefix, "2 = ", prefix, "2"))

  for(i in seq_len(n_levels)){
    if(i > 2){
      out <- paste0(out,
                    paste0(", ", prefix, i, " = ", prefix, i))

    }
  }

  out <- paste0(out,
                paste0("); return(out)}"))

  out <- eval(parse(text = out))

  fml <- paste0(prefix, "2")

  for(i in seq_len(n_levels)){
    if(i > 2){
      fml <- paste(fml, paste0(prefix, i), sep = " + ")

    }
  }

  message("Use the following formula: \n", fml)
  return(out)
}

.numericalConvFun <- function(n_levels){
  if(n_levels < 3){
    stop("n_levels should be greater than 2", call. = FALSE)
  }

  if(n_levels == 3){
    out <- function(g){
      add <- g
      dom <- as.numeric(g == 1)
      out <- data.frame(add = add, dom = dom)
      return(out)
    }
    message("Use the following formula: \n",
            "add + dom")
  }

  if(n_levels > 3){
    out <- function(g){
      add <- g
      dom <- as.numeric(g != c(0, n_levels - 1))
      if(n_levels %% 2 == 0){
        hetadd <- g - {(n_levels - 1) / 2}
        hetadd[hetadd > 0] <- ceiling(hetadd[hetadd > 0])
        hetadd[hetadd < 0] <- floor(hetadd[hetadd < 0])
        hetadd[g == c(0, n_levels - 1)] <- 0

      } else {
        hetadd <- g - {(n_levels - 1) / 2}
        hetadd[g == c(0, n_levels - 1)] <- 0
      }
      out <- data.frame(add = add, dom = dom, hetadd = hetadd)
      return(out)
    }
    message("Use the following formula: \n",
            "add + dom + hetadd")
  }

  return(out)
}

################################################################################
#' Scan QTL
#'
#' @param object A LazyGas object
#' @param out_fn Prefix of output file name
#' @param formula The formula of the regression model
#' @param conv_fun The function to convert genotype data to a model matrix for the regression
#'
#' @importFrom gdsfmt apply.gdsn objdesp.gdsn
#'
#' @export
#'
setGeneric("scanAssoc", function(object,
                                 formula,
                                 conv_fun = NULL,
                                 geno_format = c("genotype", "dosage", "haplotype"),
                                 kruskal = NULL,
                                 ...)
  standardGeneric("scanAssoc"))

## Define the scanAssoc method for the LazyGas class
setMethod("scanAssoc",
          "LazyGas",
          function(object,
                   formula,
                   conv_fun = NULL,
                   geno_format = c("genotype", "dosage", "haplotype"),
                   kruskal = NULL){

            # Notify the user of the phenotypes being analyzed
            message("Going to analyze the following phenotypes: \n",
                    paste(object@lazydata$pheno_names, collapse = ", "))

            # Add a folder to store scan results in the GDS object
            .create_scan_folder(object)

            # Loop through each phenotype name and perform analysis
            for(i in seq_along(object@lazydata$pheno_names)){
              i_pheno_names <- object@lazydata$pheno_names[i]
              dokruskal <- i_pheno_names %in% kruskal  # Determine if Kruskal-Wallis test is needed

              # Notify the user of the phenotype being processed
              message("Processing: \n",
                      paste(i_pheno_names, collapse = ", "))

              # Retrieve the phenotype data from the GDS object
              binary <- getPheno(object = object)$pheno_type$binary[i]
              i_pheno <- getPheno(object = object)$pheno[, i_pheno_names]
              i_pheno <- .standardize(val = i_pheno, binary = binary)

              # Perform regression analysis and store results
              .perform_regression(object = object,
                                  i_pheno = i_pheno,
                                  geno_format = geno_format,
                                  conv_fun = conv_fun,
                                  formula = formula,
                                  dokruskal = dokruskal,
                                  i_pheno_names = i_pheno_names,
                                  binary = binary)
            }

            # Add additional information to the root node in the GDS object
            .store_additional_info(object = object,
                                   kruskal = kruskal,
                                   formula = formula,
                                   conv_fun = conv_fun,
                                   geno_format = geno_format)
          }
)

## Sub-function to create the scan folder in the GDS object
.create_scan_folder <- function(object) {
  .create_gdsn(root_node = object$root, target_node = "lazygas",
               new_node = "scan", is_folder = TRUE)
}

.standardize <- function(val, binary){
  if(binary){
    return(val)

  } else {
    return((val - mean(val, na.rm = TRUE)) / sd(val, na.rm = TRUE))
  }
}

## Sub-function to perform regression analysis and store results
#' @importFrom gdsfmt apply.gdsn
.perform_regression <- function(object,
                                i_pheno,
                                geno_format,
                                conv_fun,
                                formula,
                                dokruskal,
                                i_pheno_names,
                                binary) {
  geno_format <- match.arg(arg = geno_format, c("genotype", "corrected", "dosage", "haplotype"))

  margin <- switch(geno_format,
                   "genotype" = 3,
                   "corrected" = 3,
                   "dosage" = 2,
                   "haplotype" = 3)

  path <- switch(geno_format,
                 "genotype" = "genotype/data",
                 "corrected" = "annotation/format/CGT/data",
                 "dosage" = "annotation/format/EDS/data",
                 "haplotype" = "annotation/format/HAP/data")

  # Get the index of the target node in the GDS object
  target_node_index <- index.gdsn(node = object, path = path)

  # Set the selection criteria based on the genotype format
  selection <- .set_selection_criteria_for_scan(object = object,
                                                geno_format = geno_format,
                                                target_node_index = target_node_index)

  if(geno_format == "haplotype"){
    na_val <- 0

  } else if(geno_format == "corrected"){
    na_val <- 3


  } else if(geno_format == "dosage"){
    na_val <- 63

  } else {
    na_val <- 2
  }

  # Perform regression analysis on genotype data
  p_values <- apply.gdsn(node = target_node_index,
                         margin = margin,
                         as.is = "list",
                         selection = selection,
                         FUN = .regression,
                         pheno = i_pheno,
                         binary = binary,
                         na_val = na_val,
                         geno_format = geno_format,
                         conv_fun = conv_fun,
                         formula = formula,
                         dokruskal = dokruskal)

  # Fill NA values
  p_values <- .fill.na(p_values = p_values)

  # Combine results into a single data frame
  p_values <- do.call("rbind", p_values)

  # Calculate FDR and -log10(p-values)
  p_values <- cbind(p_values,
                    FDR = p.adjust(p = p_values[, "P.model"], method = "fdr"),
                    negLog10P = -log10(p_values[, "P.model"]))

  # Add results to the scan node in the GDS object
  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = i_pheno_names,
               val = p_values,
               storage = "float",
               valdim = dim(p_values),
               attr = list(colnames = colnames(p_values)))
}

## Sub-function to set selection criteria based on the genotype format
.set_selection_criteria_for_scan <- function(object, geno_format, target_node_index) {
  if(geno_format == "dosage"){
    selection <- list(validSam(object = object),
                      validMar(object = object))

  } else {
    obj_desp <- objdesp.gdsn(node = target_node_index)
    selection <- list(rep(TRUE, obj_desp$dim[1]),
                      validSam(object = object),
                      validMar(object = object))
  }
  return(selection)
}

## Sub-function to store additional information in the GDS object
.store_additional_info <- function(object,
                                   kruskal,
                                   formula,
                                   conv_fun,
                                   geno_format) {
  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "kruskal",
               val = kruskal,
               storage = "string32",
               valdim = dim(kruskal))

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "formula",
               val = formula,
               storage = "string32",
               valdim = dim(formula))

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "conv_fun",
               val = deparse(conv_fun),
               storage = "string32",
               valdim = dim(deparse(conv_fun)))

  .create_gdsn(root_node = object$root,
               target_node = "lazygas/scan",
               new_node = "geno_format",
               val = geno_format,
               storage = "string32",
               valdim = dim(geno_format))
}

## Perform regression analysis on genotype data
.regression <- function(g,
                        pheno,
                        binary,
                        na_val,
                        geno_format,
                        conv_fun,
                        formula,
                        dokruskal){
    g[g == na_val] <- NA
    if(all(is.na(g))){
      return(NA)
    }

  # Check if the genotype data has only one unique value after removing NAs
  if(length(unique(na.omit(as.vector(g)))) == 1){
    return(NA)  # Return NA if there is no variability in the genotype data
  }

  if(dokruskal){
    # If Kruskal-Wallis test is required
    df <- data.frame(phe = pheno, group = g)  # Create a data frame with phenotype and genotype data
    out <- .doKruskal(df)  # Perform Kruskal-Wallis test

  } else {
    # If GLM is required
    df <- .makeDF(g = g,
                  phe = pheno,
                  conv_fun = conv_fun,
                  formula = formula)  # Create a data frame for GLM

    # Set the family for GLM based on whether the phenotype is binary or continuous
    if(binary){
      family <- "binomial"

    } else {
      family <- "gaussian"
    }
    out <- .doGLM(df = df, family = family)  # Perform GLM
  }

  return(out)  # Return the result of the regression analysis
}

## Create a data frame for regression analysis
.makeDF <- function(g, phe, conv_fun, formula){
  if(is.null(conv_fun)){
    # If no custom function is provided, convert genotype data to numeric and create a simple data frame
    g <- as.numeric(g)
    df <- data.frame(add = g)  # Create a data frame with the genotype data
    fml <- formula("phe ~ add")  # Define the formula for regression

  } else {
    # If a custom function is provided, use it to create the data frame
    df <- conv_fun(g)  # Apply the custom function to the genotype data
    if(formula == ""){
      stop("Provide formula if you specified conv_fun",
           call. = FALSE)  # Ensure that a formula is provided if a custom function is used
    }
    fml <- formula(paste0("phe ~ ", formula))  # Define the formula for regression based on the provided string
  }

  return(list(df = data.frame(phe = phe, df), fml = fml))  # Return the data frame and formula as a list
}

## Perform Generalized Linear Model (GLM) analysis
#' @importFrom parameters p_value
.doGLM <- function(df, family, null_df = NULL, model = FALSE){
  # Try to fit the GLM to the data
  res <- try(glm(formula = df$fml, data = df$df, family = family))

  # If there is an error in fitting the model, return NA
  if(inherits(res, "try-error")){ return(NA) }

  # If the model parameter is TRUE, return the fitted model object
  if(model){
    return(res)
  }

  # Calculate p-value based on the model's deviance
  if(is.null(null_df)){
    # If no null model is provided, calculate p-value using the null and residual deviance of the model
    p <- pchisq(res$null.deviance - res$deviance,
                res$df.null - res$df.residual,
                lower.tail = FALSE)

  } else {
    # If a null model is provided, fit the null model and calculate p-value
    null_res <- try(glm(formula = null_df$fml,
                        data = null_df$df,
                        family = family))
    p <- pchisq(null_res$deviance - res$deviance,
                null_res$df.residual - res$df.residual,
                lower.tail = FALSE)
  }

  # Calculate the proportion of variance explained by the model
  pervar <- 1 - res$deviance / res$null.deviance
  s <- summary(res)
  coef <- s$coefficients

  # If there is only one row in the coefficients, return NA
  if(nrow(coef) == 1){ return(NA) }

  # Get the terms of the model
  att <- attributes(res$terms)

  # Create the output with the p-value, coefficients, proportion of variance, and term labels
  return(.makeOut(p, coef, pervar, att$term.labels))
}

## Perform Kruskal-Wallis test
.doKruskal <- function(df){
  # Try to perform the Kruskal-Wallis test using the formula and data provided
  res <- try(kruskal.test(formula = "phe ~ group", data = df$df))

  # If there is an error in performing the test, return NA
  if(inherits(res, "try-error")){ return(NA) }

  # Return the p-value of the test
  return(c(P.model = res$p.value))
}

## Create output for regression results
.makeOut <- function(p, coef, pervar, terms){
  # Define the labels for the output
  labs <- c("P.model", paste("P", terms, sep = "."),
            paste("Coef", terms, sep = "."),
            "PVE")

  # Initialize the output vector with NA values
  out <- rep(NA, length(labs))

  # Assign the p-value to the first position in the output
  out[1] <- p

  # Get the row names of the coefficients
  coef_row <- rownames(coef)

  # Loop through each term to populate p-values and coefficients in the output
  for(j in seq_along(terms)){
    if(terms[j] %in% coef_row){
      out[1 + j] <- coef[terms[j], 4]  # p-value for the term
      out[1 + j + length(terms)] <- coef[terms[j], 1]  # coefficient for the term
    }
  }

  # Assign the proportion of variance explained to the output
  out[2 + length(terms) * 2] <- pervar

  # Assign names to the output
  names(out) <- labs

  return(out)  # Return the output vector
}

## Fill NA values in a list of p-values
.fill.na <- function(p_values){
  # Determine the length of each element in the list
  v_len <- sapply(X = p_values, FUN = function(x){
    if(all(is.na(x))){
      return(NA)  # Return NA if all elements in the list are NA

    } else {
      return(length(x))  # Return the length of the element if it is not all NA
    }
  })

  # Replace NA elements in the list with a vector of NA values of appropriate length
  p_values[is.na(v_len)] <- list(rep(NA, v_len[!is.na(v_len)][1]))

  return(p_values)  # Return the modified list
}

################################################################################
#' Draw a Manhattan plot
#'
#' @param object QTLscan object
#' @param pehno Phenotype names to be drawn
#' @param chr Chromosome ID to be drawn
#' @param start Start position of the range to be drawn
#' @param end End position of the range to be drawn
#' @param signif Expression to define the significant markers
#' @param out_fn Prefix of output file
#'
#' @importFrom gdsfmt ls.gdsn
#' @import ggplot2
#' @export
#'
setGeneric("plotManhattan", function(object,
                                     pheno = NULL,
                                     chr = NULL,
                                     start = NULL,
                                     end = NULL,
                                     signif = "FDR<0.05",
                                     out_fn = "",
                                     ...)
  standardGeneric("plotManhattan"))

setMethod("plotManhattan",
          "LazyGas",
          function(object,
                   pheno = NULL,
                   chr = NULL,
                   start = NULL,
                   end = NULL,
                   signif = "FDR<0.05"){
            # Check if scan data exists in the GDS object
            .check_scan_data_exists(object = object)

            # Determine the phenotype index based on user input
            pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)

            # List the scan nodes in the GDS object
            scan_node <- ls.gdsn(node = index.gdsn(node = object, path = "lazygas/scan"))

            # Retrieve the scan data for the specified phenotype
            pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))

            # Generate the Manhattan plot
            p <- .mhplot_draw(x = pvalues, chr = chr, start = start, end = end, signif = signif)

            return(p)  # Return the plot
          }
)

## Function to check if scan data exists in the GDS object
.check_scan_data_exists <- function(object) {
  if(!exist.gdsn(node = object$root, path = "lazygas/scan")){
    stop("No scan data in the input LazyQTL object.\n",
         "Run scanAssoc() to scan associations.")
  }
}

## Function to determine the phenotype index based on user input
.determine_phenotype_name <- function(object, pheno) {
  if(!is.null(pheno)){
    if(is.numeric(pheno)){
      phe_index <- pheno

    } else if(is.logical(pheno)){
      phe_index <- which(pheno)

    } else if(is.character(pheno)){
      phe_index <- which(getPheno(object = object)$pheno_names %in% pheno)
    }

    # Notify the user if multiple phenotypes are selected and use the first one
    if(length(phe_index) != 1){
      message("Multiple (or no) phenotypes were selected.",
              "\nUse the first phenotype.")
      phe_index <- 1
    }
  } else {
    # If no phenotype is specified, use the first one and notify the user
    message("Use the first phenotype.")
    phe_index <- 1
  }
  return(getPheno(object = object)$pheno_names[phe_index])
}

## Draw Manhattan plot
#' @import ggplot2
## Draw Manhattan plot
.mhplot_draw <- function(x,
                         chr = NULL,
                         start = NULL,
                         end = NULL,
                         signif = NULL){
  chr_lev <- unique(x$Chr)
  x$Chr <- factor(x$Chr, levels = chr_lev)

  # Filter data by chromosome, start, and end positions
  x <- .filter_data(x, chr, start, end)

  # Evaluate significant points
  signif <- .evaluate_signif(x, signif)
  x_signif <- subset(x, subset = signif)
  x <- subset(x, subset = !signif)

  # Create ggplot object
  p <- .create_ggplot(x, x_signif)

  # Return the plot
  return(p)
}

## Function to filter data by chromosome, start, and end positions
.filter_data <- function(x, chr, start, end) {
  if(!is.null(chr)){
    x <- subset(x, subset = Chr == chr)
  }
  if(!is.null(start)){
    x <- subset(x, subset = Pos >= start)
  }
  if(!is.null(end)){
    x <- subset(x, subset = Pos <= end)
  }
  return(x)
}

## Function to evaluate significant points
.evaluate_signif <- function(x, signif) {
  if(is.null(signif)){
    signif <- rep(TRUE, nrow(x))
  }
  if(!is.logical(signif)){
    if(is.character(signif)){
      signif <- eval(parse(text = paste0("x$", signif)))
    } else {
      stop("signif should be a string or a vector of logical values.", call. = FALSE)
    }
  }
  signif[is.na(signif)] <- FALSE
  return(signif)
}

## Function to create ggplot object
.create_ggplot <- function(x, x_signif) {
  p <- ggplot() +
    geom_point(data = x,
               mapping = aes(x = Pos, y = negLog10P),
               color = "darkgray",
               size = 1,
               shape = 20)

  if(nrow(x_signif) != 0){
    p <- p + geom_point(data = x_signif,
                        mapping = aes(x = Pos, y = negLog10P),
                        color = "magenta", size = 1, shape = 18)
  }

  p <- p + facet_wrap(~ Chr,
                      nrow = 1,
                      scales = "free_x",
                      strip.position = "bottom") +
    ylab("-log10(P)") +
    xlab("Chromosome") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          strip.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.line.x.bottom = element_line(colour = "black"),
          panel.spacing.x = unit(0.2, "lines"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_rect(fill = "white", colour = "white"))

  return(p)
}

################################################################################
#' Call peack blocks
#'
#' @param object QTLscan object
#' @param signif expression to define significant markers
#' @param rsquare threshold on squared R values to define peak blocks
#'
#' @export
#'
setGeneric("callPeakBlock", function(object,
                                     signif = 0.05,
                                     threshold = 0.8,
                                     n_threads = NULL,
                                     ...)
  standardGeneric("callPeakBlock"))

## Define the callPeakBlock method for the LazyGas class
setMethod("callPeakBlock",
          "LazyGas",
          function(object,
                   signif,
                   threshold,
                   n_threads){
            # Check if scan data exists in the GDS object
            .check_scan_data_exists(object = object)

            # Create the necessary folders in the GDS object
            .create_peakcall_folders(object = object)

            # Perform peak calling for each phenotype
            .perform_peak_calling(object = object,
                                  signif = signif,
                                  threshold = threshold,
                                  n_threads = n_threads)
          }
)

## Function to create the necessary folders in the GDS object
.create_peakcall_folders <- function(object) {
  peakcall_gdsn <- .create_gdsn(root_node = object$root,
                                target_node = "lazygas",
                                new_node = "peakcall",
                                is_folder = TRUE)
  peaks_gdsn <- .create_gdsn(root_node = object$root,
                             target_node = "lazygas/peakcall",
                             new_node = "peaks",
                             is_folder = TRUE)
  blocks_gdsn <- .create_gdsn(root_node = object$root,
                              target_node = "lazygas/peakcall",
                              new_node = "blocks",
                              is_folder = TRUE)

  put.attr.gdsn(node = peaks_gdsn,
                name = "col_names",
                val = c("peak_ID", "peak_variant_ID"))
  put.attr.gdsn(node = blocks_gdsn,
                name = "col_names",
                val = c("peak_ID", "variant_ID", "dist2peak", "LD2peak"))
}

## Function to perform peak calling for each phenotype
.perform_peak_calling <- function(object, signif, threshold, n_threads) {
  for(i in seq_along(object@lazydata$pheno_names)){
    message("Peak calling for the following phenotype: ", object@lazydata$pheno_names[i])
    .peakcaller(object = object,
                signif = signif,
                threshold = threshold,
                pheno_name = object@lazydata$pheno_names[i],
                n_threads = n_threads)
  }
}



## Main peakcaller function
.peakcaller <- function(object, signif, threshold, pheno_name, n_threads){
  pvalues <- .get_and_filter_pvalues(object = object,
                                     pheno_name = pheno_name,
                                     signif = signif)

  if(nrow(pvalues) == 0){
    .handle_no_significant_peaks(object = object,
                                 pheno_name = pheno_name)
    return()
  }

  variables <- .initialize_variables(object = object)

  .call_peaks(object = object,
              pvalues = pvalues,
              variables = variables,
              pheno_name = pheno_name,
              threshold = threshold,
              n_threads = n_threads)
}

## Function to get and filter significant p-values
.get_and_filter_pvalues <- function(object, pheno_name, signif) {
  pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))

  signif <- pvalues$FDR <= signif
  signif[is.na(signif)] <- FALSE
  pvalues <- subset(pvalues, subset = signif)

  return(pvalues)
}

## Function to handle the case where no significant peaks are found
.handle_no_significant_peaks <- function(object, pheno_name) {
  peaks_gdsn <- add.gdsn(node = index.gdsn(node = object$root,
                                           path = "lazygas/peakcall/peaks"),
                         name = pheno_name, storage = "uint32",
                         compress = "ZIP_RA", replace = TRUE)
  blocks_gdsn <- add.gdsn(node = index.gdsn(node = object$root,
                                            path = "lazygas/peakcall/blocks"),
                          name = pheno_name, storage = "double",
                          compress = "ZIP_RA", replace = TRUE)
}

## Function to initialize variables
.initialize_variables <- function(object) {
  list(
    n_sample = nsam(object = object),
    snp_id = getMarID(object = object),
    chr = getChromosome(object = object),
    pos = getPosition(object = object)
  )
}

## Function to process peaks in a loop
.call_peaks <- function(object, pvalues, variables, pheno_name, threshold, n_threads) {
  peak_id <- 0

  # Get the genotype format from the object
  geno_format <- .get_geno_format(object = object)

  # Loop until there are no more p-values to process
  while (nrow(pvalues) > 0) {
    peak_id <- peak_id + 1
    message("Calling peak: ", peak_id)

    # Identify the peak information
    peak_info <- .identify_peak(pvalues = pvalues,
                                variables = variables)

    # Set the selection criteria based on the genotype format
    selection <- .set_selection_criteria_for_peakcall(object = object,
                                                      geno_format = geno_format,
                                                      chr = peak_info$chr)

    # Retrieve genotype data based on the selection criteria
    geno <- .retrieve_geno(object = object,
                           selection = selection,
                           geno_format = geno_format)

    # Calculate linkage disequilibrium (LD) for the identified peak
    peak_ld <- .calculate_ld(geno = as.matrix(geno),
                             peak_variant_idx = peak_info$peak_index_in_chr,
                             is_categorical = selection$is_categorical,
                             n_threads = n_threads)

    # Identify the peak block based on the LD values
    peak_block <- .identify_peak_block(peak_ld = peak_ld,
                                       variables = variables,
                                       peak_info = peak_info,
                                       threshold = threshold,
                                       geno_format = geno_format)

    # Save the peak data into the GDS object
    .save_peak_data(object = object,
                    pheno_name = pheno_name,
                    peak_info = peak_info,
                    peak_block = peak_block,
                    peak_id = peak_id,
                    peak_ld = peak_ld)

    # Update pvalues by removing the variants that are in the peak block
    rm_ids <- seq(min(peak_block$ids), max(peak_block$ids))
    pvalues <- subset(pvalues, subset = !variant_ID %in% rm_ids)
  }

  # Finalize GDS nodes by setting them to read mode
  .finalize_gdsn_peakcall(object = object, pheno_name = pheno_name)
}

## Function to identify the peak
.identify_peak <- function(pvalues, variables, peak_index = NULL) {
  if(is.null(peak_index)){
    peak_index <- which.max(pvalues$negLog10P)  # Find the index of the maximum p-value
  }
  peak_variant_ID <- pvalues$variant_ID[peak_index]  # Get the variant ID of the peak
  peak_chr <- variables$chr[variables$snp_id == peak_variant_ID]  # Get the chromosome of the peak

  chr_with_peak <- which(variables$chr == peak_chr)
  id_in_chr <- variables$snp_id[chr_with_peak]
  peak_index_in_chr <- which(id_in_chr == peak_variant_ID)

  list(
    index = peak_index,
    variantID = peak_variant_ID,
    chr = peak_chr,
    chr_with_peak = chr_with_peak,
    id_in_chr = id_in_chr,
    peak_index_in_chr = peak_index_in_chr
  )
}

## Get Genotype Format
.get_geno_format <- function(object) {
  geno_format <- .get_data(object = object, node = "lazygas/scan/geno_format")
  return(geno_format)
}

# Function to set the selection criteria based on the genotype format
.set_selection_criteria_for_peakcall <- function(object, geno_format, chr) {
  if (geno_format == "dosage"){
    # Identify valid chromosomes
    valid_chr <- getChromosome(object = object, valid = FALSE) %in% chr

    # Set selection criteria for genotype or dosage data
    selection <- list(validSam(object = object),
                      validMar(object = object) & valid_chr)
    node <- "annotation/format/EDS/data"
    is_categorical <- TRUE

  } else {
    # Identify valid chromosomes
    valid_chr <- getChromosome(object = object, valid = FALSE) %in% chr

    # Get the dimensions of the haplotype data
    obj_desp <- objdesp.gdsn(node = index.gdsn(node = object,
                                               path = "annotation/format/HAP/data"))

    # Set selection criteria for haplotype data
    selection <- list(rep(TRUE, obj_desp$dim[1]),
                      validSam(object = object),
                      validMar(object = object) & valid_chr)

    if(geno_format == "haplotype"){
      node <- "annotation/format/HAP/data"
      is_categorical <- TRUE

    } else {
      node <- ifelse(geno_format == "corrected",
                     "annotation/format/CGT/data",
                     "genotype/data")
      is_categorical <- FALSE
    }
  }
  return(list(node = node, selection = selection, is_categorical = is_categorical))
}

#' @importFrom gdsfmt exist.gdsn index.gdsn objdesp.gdsn
.retrieve_geno <- function(object, selection, geno_format, collapse = TRUE){
  out <- .get_data(object = object,
                   node = selection$node,
                   sel = selection$selection)
  if(geno_format == "haplotype"){
    if(collapse){
      out <- apply(X = out, MARGIN = 3, FUN = c)
    }
    out[out == 0] <- NA

  } else if(geno_format %in% c("corrected", "genotype")){
    out[out == 2] <- NA

  } else if(geno_format == "dosage"){
    out[out == 63] <- NA
  }
  return(out)
}

#' @import parallel
.calculate_ld <- function(geno, peak_variant_idx, is_categorical, n_threads = NULL) {
  n <- ncol(geno)
  peak_ld <- numeric(n)

  if (is_categorical) {
    geno_col <- as.integer(geno[, peak_variant_idx])
    levels <- sort(unique(geno_col))
    num_levels <- length(levels)

    # Parallel processing setup
    if(is.null(n_threads)){
      cores <- detectCores()
      if(cores > 1){
        n_threads <- round(cores / 2)
      }
    }

    geno_identity <- mclapply(X = 1:n, mc.cores = n_threads, mc.preschedule = TRUE,
                              FUN = function(j) {
                                geno_col_j <- as.integer(geno[, j])
                                geno_match <- geno_col == geno_col_j
                                geno_identity <- sum(geno_match, na.rm = TRUE) / sum(!is.na(geno_match))
                                return(geno_identity)
                                # cont_table <- matrix(0, num_levels, num_levels)
                                #
                                # for (k in seq_len(nrow(geno))) {
                                #   row_idx <- match(geno_col[k], levels)
                                #   col_idx <- match(geno_col_j[k], levels)
                                #   cont_table[row_idx, col_idx] <- cont_table[row_idx, col_idx] + 1
                                # }
                                #
                                # chi_sq <- 0
                                # row_sums <- rowSums(cont_table)
                                # col_sums <- colSums(cont_table)
                                # total_sum <- sum(cont_table)
                                #
                                # for (r in seq_len(num_levels)) {
                                #   for (c in seq_len(num_levels)) {
                                #     expected <- (row_sums[r] * col_sums[c]) / total_sum
                                #     chi_sq <- chi_sq + ((cont_table[r, c] - expected) ^ 2) / expected
                                #   }
                                # }
                                # return(chi_sq / nrow(geno))
                              })
    peak_ld <- unlist(geno_identity)

  } else {
    # Parallel processing setup
    corr_values <- cor(geno[, peak_variant_idx], geno, use = "pairwise.complete.obs") ^ 2
    peak_ld <- corr_values
  }

  return(peak_ld)
}

# Function to identify peak block
.identify_peak_block <- function(peak_ld, variables, peak_info, threshold, geno_format) {
  if(geno_format == "haplotype"){
    ld_block <- peak_ld >= max(peak_ld, na.rm = TRUE) * threshold

  } else {
    ld_block <- peak_ld >= threshold
  }
  peak_block <- na.omit(peak_info$id_in_chr[ld_block])
  ld_to_peak <- na.omit(peak_ld[ld_block])
  pos <- variables$pos[variables$snp_id %in% peak_info$id_in_chr]
  dist2peak <- na.omit(pos[ld_block] - pos[peak_info$peak_index_in_chr])

  return(list(ids = peak_block, ld = ld_to_peak, dist = dist2peak))
}

# Function to save the peak data into the GDS object
#' @importFrom gdsfmt append.gdsn
.save_peak_data <- function(object, pheno_name, peak_info, peak_block, peak_id, peak_ld, node = "peakcall") {
  # If this is the first peak, create new nodes for peaks and blocks
  if (peak_id == 1) {
    # Create a new node for storing peak information
    add.gdsn(node = index.gdsn(node = object$root,
                               path = paste0("lazygas/", node, "/peaks")),
             name = pheno_name, storage = "uint32", compress = "ZIP_RA", replace = TRUE,
             val = rbind(peak_id, peak_info$variantID))

    # Create a new node for storing block information
    add.gdsn(node = index.gdsn(node = object$root,
                               path = paste0("lazygas/", node, "/blocks")),
             name = pheno_name, storage = "double", compress = "ZIP_RA", replace = TRUE,
             val = rbind(peak_id, peak_block$ids, peak_block$dist, peak_block$ld))

  } else {
    # Append the peak information to the existing node
    peaks_gdsn <- index.gdsn(node = object$root,
                             path = paste0("lazygas/", node, "/peaks/", pheno_name))
    append.gdsn(node = peaks_gdsn, val = rbind(peak_id, peak_info$variantID))

    # Append the block information to the existing node
    blocks_gdsn <- index.gdsn(node = object$root,
                              path = paste0("lazygas/", node, "/blocks/", pheno_name))
    append.gdsn(node = blocks_gdsn,
                val = rbind(peak_id, peak_block$ids, peak_block$dist, peak_block$ld))
  }
}

## Sub-function to finalize GDS nodes
.finalize_gdsn_peakcall <- function(object, pheno_name) {
  readmode.gdsn(node = index.gdsn(node = object$root,
                                  path = paste0("lazygas/peakcall/peaks/", pheno_name)))
  readmode.gdsn(node = index.gdsn(node = object$root,
                                  path = paste0("lazygas/peakcall/blocks/", pheno_name)))
}


################################################################################
#' Plots the peaks for the specified phenotype and chromosome region.
#'
#' @param object A LazyGas object.
#' @param pheno The phenotype(s) to plot. Can be NULL, numeric, logical, or character.
#' @param chr The chromosome to plot. Can be NULL.
#' @param start The start position on the chromosome. Can be NULL.
#' @param end The end position on the chromosome. Can be NULL.
#' @param recalc Logical. If TRUE, use recalculated scan data.
#' @param ... Additional arguments.
#'
#' @return Generates plots for the specified peaks.
#'
#' @export
#'
setGeneric("plotPeaks", function(object,
                                 pheno = NULL,
                                 chr = NULL,
                                 start = NULL,
                                 end = NULL,
                                 recalc = TRUE,
                                 ...)
  standardGeneric("plotPeaks"))

## Define the plotPeaks method for the LazyGas class
#'
#' @method plotPeaks LazyGas
#'
setMethod("plotPeaks",
          "LazyGas",
          function(object,
                   pheno = NULL,
                   chr = NULL,
                   start = NULL,
                   end = NULL,
                   recalc = TRUE){
            .check_peakcall_data(object = object)
            path <- .determine_path(object = object, recalc = recalc)
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)
            chr_pos_limits <- .get_chromosome_data(object = object)

            p <- .generate_peak_plots(object = object,
                                      pheno_name = pheno_name,
                                      path = path,
                                      chr_pos_limits = chr_pos_limits,
                                      chr = chr,
                                      start = start,
                                      end = end)
            return(p)
          }
)

# Check if peakcall data exists
.check_peakcall_data <- function(object) {
  if (!exist.gdsn(node = object$root, path = "lazygas/peakcall")) {
    stop("No peakcall data in the input LazyQTL object.\n",
         "Run callPeakBlock() to call peaks.")
  }
}

# Determine the path based on recalc parameter
.determine_path <- function(object, recalc) {
  if (recalc) {
    if (!exist.gdsn(node = object$root, path = "lazygas/recalc")) {
      stop("recalc = TRUE was specified but, \n",
           "no recalculated scan data in the input LazyQTL object.\n",
           "Run recalcAssoc() to recalculate associations.")
    }
    return("lazygas/recalc")

  } else {
    return("lazygas/peakcall")
  }
}

# Get chromosome and position data
.get_chromosome_data <- function(object) {
  chr_data <- getChromosome(object = object)
  pos_data <- getPosition(object = object)
  chr_min <- tapply(pos_data, chr_data, min)
  chr_max <- tapply(pos_data, chr_data, max)
  return(list(chr_min = chr_min, chr_max = chr_max, chr_lev = unique(chr_data)))
}

# Generate peak plots
.generate_peak_plots <- function(object,
                                 pheno_name,
                                 path,
                                 chr_pos_limits,
                                 chr,
                                 start,
                                 end) {
  scan_node <- ls.gdsn(node = index.gdsn(node = object,
                                         path = paste0(path, "/peaks")))
  peakcall <- .get_peakcall(object = object,
                            pheno_name = pheno_name,
                            recalc = grepl(pattern = "recalc", x = path))

  if (is.null(peakcall)) {
    peakcall <- data.frame(peak_variant_ID = numeric(), FDR = numeric(),
                           variant_ID = numeric(), Chr = numeric(),
                           Pos = numeric(), negLog10P = numeric(),
                           peak_ID = numeric())
    start <- NULL
    end <- NULL
    chr <- NULL
  }

  # Generate the peak plot
  p <- .peak_draw(peakcall = peakcall,
                  chr = chr,
                  start = start,
                  end = end,
                  chr_min = chr_pos_limits$chr_min,
                  chr_max = chr_pos_limits$chr_max,
                  chr_lev = chr_pos_limits$chr_lev)
  return(p)
}

# Draw Peak Plot
#' @import ggplot2
#'
.peak_draw <- function(peakcall, chr, start, end, chr_min, chr_max, chr_lev, peak = NULL){
  # Mark peaks in peakcall
  peakcall$is_peak <- peakcall$peak_variant_ID == peakcall$variant_ID
  peakcall$Chr <- factor(peakcall$Chr, levels = chr_lev)

  # Create dummy data for chromosome boundaries
  dummy <- rbind(data.frame(Chr = names(chr_min),
                            Pos = chr_min,
                            negLog10P = 0),
                 data.frame(Chr = names(chr_max),
                            Pos = chr_max,
                            negLog10P = 0))
  dummy$Chr <- factor(dummy$Chr, levels = chr_lev)

  # Subset peakcall and dummy data based on the specified chromosome
  if(!is.null(chr)){
    peakcall <- subset(peakcall, subset = Chr == chr)
    dummy <- subset(dummy, subset = Chr == chr)
  }

  # Subset peakcall and dummy data based on the specified start position
  if(!is.null(start)){
    peakcall <- subset(peakcall, subset = Pos >= start)
    dummy <- subset(dummy, subset = Pos >= start)
  }

  # Subset peakcall and dummy data based on the specified end position
  if(!is.null(end)){
    peakcall <- subset(peakcall, subset = Pos <= end)
    dummy <- subset(dummy, subset = Pos <= end)
  }

  # Set peak ID for dummy data
  dummy$peak_ID <- peakcall$peak_ID[1]

  # Subset significant peaks
  signif_peak <- subset(peakcall, subset = FDR <= 0.05 & peak_variant_ID == variant_ID)

  # Create ggplot object
  p <- ggplot() +
    geom_line(data = peakcall,
              mapping = aes(x = Pos, y = negLog10P, color = peak_ID, group = peak_ID),
              linewidth = 1) +
    geom_point(data = signif_peak,
               mapping = aes(x = Pos, y = negLog10P),
               color = "magenta",
               size = 2,
               shape = 20) +
    geom_point(data = dummy,
               mapping = aes(x = Pos, y = negLog10P),
               size = 0,
               color = NA) +
    facet_wrap(~ Chr, nrow = 1, scales = "free_x", strip.position = "bottom") +
    ylab("-log10(P)") +
    xlab("Chromosome") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          strip.text.x = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.line.x.bottom = element_line(colour = "black"),
          panel.spacing.x = unit(0.2, "lines"),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.placement = "outside",
          strip.background = element_rect(fill = "white", colour = "white"))
  return(p)
}

################################################################################
#' Generate Haplotype Plots
#'
#' This function generates haplotype plots for the specified phenotype.
#'
#' @param object A LazyGas object.
#' @param pheno The phenotype to plot.
#' @param recalc Logical. If TRUE, use recalculated scan data.
#' @param ... Additional arguments.
#'
#' @return A list of ggplot objects representing the haplotype plots.
#' @export
setGeneric("haploPlot", function(object,
                                 pheno,
                                 recalc,
                                 ...)
  standardGeneric("haploPlot"))

#' @rdname haploPlot
#' @method haploPlot LazyGas
#'
setMethod("haploPlot",
          "LazyGas",
          function(object, pheno, recalc){
            # Check if peakcall data exists
            .check_peakcall_data(object = object)

            # Determine the path based on recalc parameter
            path <- .determine_path(object = object, recalc = recalc)

            # Determine phenotype index based on input
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)

            # Get the list of scan nodes
            scan_node <- ls.gdsn(node = index.gdsn(node = object,
                                                   path = paste0(path, "/peaks")))

            # Get the peak call data for the specified phenotype
            peakcall <- .get_peakcall(object = object,
                                      pheno_name = pheno_name,
                                      recalc = grepl(pattern = "recalc", x = path))

            # Generate the haplotype plot
            out <- .draw_haplo_plot(object = object,
                                    peakcall = peakcall,
                                    pheno_name = pheno_name)

            return(out)
          }
)

## Draw Haplotype Plot
#' @importMethodsFrom GBScleanR getMarID getHaplotype getGenotype
#' @import ggplot2
#'
.draw_haplo_plot <- function(object, peakcall, pheno_name){
  # Check if there are any peaks
  if(is.null(peakcall)){
    message("No peak")
    return(NULL)
  }

  haplo <- NULL
  var_id <- getMarID(object = object)

  # Get the genotype format from the object
  geno_format <- .get_geno_format(object = object)

  peak_info <- .get_peak_info(peakcall = peakcall)

  out <- NULL
  # Loop through each unique peak ID
  for(i_peak in seq_along(peak_info$peak_height)){
    valid <- var_id %in% peak_info$peak_id[i_peak]

    # Get the haplotype or genotype data based on the format
    if(geno_format == "haplotype"){
      hap <- getHaplotype(object = object)[, , valid]
      hap <- apply(X = hap, MARGIN = 2, FUN = function(x){
        return(paste(sort(x), collapse= "|"))
      })

    } else {
      node <- switch (geno_format,
                      dosage = "dosage",
                      corrected = "cor",
                      genotype = "raw")
      hap <- getGenotype(object = object, node = node)[, valid]
      hap <- factor(hap, levels = sort(unique(hap)))
    }

    # Get phenotype data
    pheno <- getPheno(object = object)
    df <- data.frame(y = pheno$pheno[, pheno_name], x = hap)

    # Generate the plot
    p <- ggplot(df) +
      geom_boxplot(aes(x = x, y = y)) +
      labs(title = paste("Peak @",
                         paste(peak_info$peak_chr[i_peak],
                               peak_info$peak_pos[i_peak],
                               sep = "_"),
                         ", -log10P = ",
                         signif(x = peak_info$peak_height[i_peak],
                                digits = 3))) +
      xlab("Haplotypes") +
      ylab(pheno_name)

    out <- c(out, list(p))
  }
  return(out)
}

.get_peak_info <- function(peakcall = peakcall){
  peak_height <- tapply(X = peakcall$peak_negLog10P,
                        INDEX = peakcall$peak_ID,
                        FUN = function(x){x[1]})
  peak_id <- tapply(X = peakcall$peak_variant_ID,
                    INDEX = peakcall$peak_ID,
                    FUN = function(x){x[1]})
  peak_chr <- tapply(X = peakcall$peak_Chr,
                     INDEX = peakcall$peak_ID,
                     FUN = function(x){x[1]})
  peak_pos <- tapply(X = peakcall$peak_Pos,
                     INDEX = peakcall$peak_ID,
                     FUN = function(x){x[1]})
  out <- data.frame(peak_height = peak_height,
                    peak_id = peak_id,
                    peak_chr = peak_chr,
                    peak_pos = peak_pos)
  return(out)
}

################################################################################
#' Retrieve LazyGas Data
#'
#' This function retrieves data from the LazyGas object based on the specified dataset and phenotype.
#'
#' @param object A LazyGas object.
#' @param dataset The dataset to retrieve. One of "scan", "peakcall", or "recalc".
#' @param pheno The phenotype to retrieve.
#' @param ... Additional arguments.
#'
#' @return The requested data from the LazyGas object.
#' @export
#'
setGeneric("lazyData", function(object,
                                dataset = c("scan", "peakcall", "recalc", "candidate", "snpeff"),
                                pheno,
                                ...)
  standardGeneric("lazyData"))

#'
#' @rdname lazyData
#' @method lazyData LazyGas
#' @export
#'
setMethod("lazyData",
          "LazyGas",
          function(object, dataset, pheno){
            # Match the dataset argument to one of the allowed choices
            dataset <- match.arg(arg = dataset,
                                 choices = c("scan", "peakcall", "recalc", "candidate", "snpeff"))

            # Determine phenotype index based on input
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)

            # Retrieve the data based on the dataset type
            if(dataset == "scan"){
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = "lazygas"),
                                       path = dataset)
              if(!check_node){
                return(NULL)
              }
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = paste0("lazygas/",
                                                                       dataset)),
                                       path = pheno_name)
              if(!check_node){
                return(NULL)
              }
              out <- .get_scan(object = object,
                               pheno_name = pheno_name)

            } else if(dataset %in% c("candidate", "snpeff")){
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = "lazygas"),
                                       path = dataset)
              if(!check_node){
                return(NULL)
              }
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = paste0("lazygas/",
                                                                       dataset)),
                                       path = pheno_name)
              if(!check_node){
                return(NULL)
              }
              out <- .get_data(object = object,
                               node = paste0("lazygas/",
                                             dataset, "/",
                                             pheno_name))
              if(length(out) == 0){
                return(NULL)
              }
              if(!is.matrix(out)){
                out <- matrix(out, nrow = 1)
              }

              col_names <- get.attr.gdsn(node = index.gdsn(node = object$root,
                                                           path = paste0("lazygas/", dataset)))
              colnames(out) <- col_names$col_names

            } else {
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = "lazygas"),
                                       path = dataset)
              if(!check_node){
                return(NULL)
              }
              check_node <- exist.gdsn(node = index.gdsn(node = object,
                                                         path = paste0("lazygas/",
                                                                       dataset,
                                                                       "/peaks")),
                                       path = pheno_name)
              if(!check_node){
                return(NULL)
              }
              out <- .get_peakcall(object = object,
                                   pheno_name = pheno_name,
                                   recalc = (dataset == "recalc"))
            }

            return(out)
          }
)

################################################################################
#' Recalculate Associations
#'
#' Recalculate associations using the LazyGas object and specified r-square threshold.
#'
#' @param object A LazyGas object.
#'
#' @return None. The function modifies the LazyGas object in place.
#' @export
#'
setGeneric("recalcAssoc", function(object,
                                   signif = 0.05,
                                   threshold = 0.6,
                                   n_threads = NULL,
                                   ...)
  standardGeneric("recalcAssoc"))

#'
#' @rdname recalcAssoc
#' @method recalcAssoc LazyGas
#' @export
#'
setMethod("recalcAssoc",
          "LazyGas",
          function(object, signif, threshold, n_threads){
            if(!exist.gdsn(node = object$root, path = "lazygas/peakcall")){
              stop("No peakcall data in the input LazyGas object.\n",
                   "Run callPeakBlock() to scan associations.")
            }

            # Parallel processing setup
            if(is.null(n_threads)){
              cores <- detectCores()
              if(cores > 1){
                n_threads <- round(cores / 2)
              }
            }

            kruskal <- .get_data(object = object, node = "lazygas/scan/kruskal")

            ## Function to create the necessary folders in the GDS object
            .create_recalc_folders(object = object)

            pheno <- getPheno(object = object)

            for(i in seq_along(pheno$pheno_names)){
              pheno_name <- pheno$pheno_names[i]
              dokruskal <- pheno_name %in% kruskal

              if(dokruskal){
                message("Skip recalculation for non-parametric data: ", pheno_name)
                .create_gdsn(root_node = object$root,
                             target_node = "lazygas/recalc/peaks",
                             new_node = pheno_name,
                             storage = "uint32",
                             replace = TRUE)
                .create_gdsn(root_node = object$root,
                             target_node = "lazygas/recalc/blocks",
                             new_node = pheno_name,
                             storage = "double",
                             replace = TRUE)

              } else {
                .peakrecalculator(object = object,
                                  threshold = threshold,
                                  signif = signif,
                                  pheno = pheno,
                                  pheno_name = pheno_name,
                                  binary = pheno$pheno_type$binary[i],
                                  n_threads = n_threads)
              }

              .finalize_gdsn_recalc(object = object, pheno_name = pheno_name)
            }
          }
)

.create_recalc_folders <- function(object) {
  peakcall_gdsn <- .create_gdsn(root_node = object$root,
                                target_node = "lazygas",
                                new_node = "recalc",
                                is_folder = TRUE)
  peaks_gdsn <- .create_gdsn(root_node = object$root,
                             target_node = "lazygas/recalc",
                             new_node = "peaks",
                             is_folder = TRUE)
  blocks_gdsn <- .create_gdsn(root_node = object$root,
                              target_node = "lazygas/recalc",
                              new_node = "blocks",
                              is_folder = TRUE)
  # group_gdsn <- .create_gdsn(root_node = object$root,
  #                            target_node = "lazygas/recalc",
  #                            new_node = "group",
  #                            is_folder = TRUE)
  # put.attr.gdsn(node = group_gdsn, name = "col_names",
  #               val = c("round", "peak_variant_ID", "member_peak", "new_peak"))
}

## Sub-function to finalize GDS nodes
.finalize_gdsn_recalc <- function(object, pheno_name) {
  readmode.gdsn(node = index.gdsn(node = object$root,
                                  path = paste0("lazygas/recalc/peaks/", pheno_name)))
  readmode.gdsn(node = index.gdsn(node = object$root,
                                  path = paste0("lazygas/recalc/blocks/", pheno_name)))
  # readmode.gdsn(node = index.gdsn(node = object$root,
  #                                 path = paste0("lazygas/recalc/group/", pheno_name)))
}

# Recalculate peaks for a given phenotype.
#' @importFrom dplyr setequal
.peakrecalculator <- function(object, signif, threshold, pheno, pheno_name, binary, n_threads){
  message("Processing: ", pheno_name)

  peakcall <- .get_peakcall(object = object, pheno_name = pheno_name)
  # peak_grp_summary <- NULL

  if(is.null(peakcall)){
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc/peaks",
                 new_node = pheno_name,
                 storage = "uint32",
                 replace = TRUE)
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/recalc/blocks",
                 new_node = pheno_name,
                 storage = "double",
                 replace = TRUE)
    # .create_gdsn(root_node = object$root,
    #              target_node = "lazygas/recalc/group",
    #              new_node = pheno_name,
    #              storage = "uint32",
    #              replace = TRUE)

  } else {
    peak_obj <- .make_peakobj(object = object,
                              threshold = threshold,
                              signif = signif,
                              pheno = pheno,
                              pheno_name = pheno_name,
                              peakcall = peakcall,
                              binary = binary)
    message("Recalculating associations...")
    count <- 0

    while(TRUE){
      count <- count + 1
      message("Round ", count)

      cov_scan <- lapply(X = peak_obj$peak_variant_id,
                         FUN = .composite,
                         peak_obj = peak_obj,
                         n_threads = n_threads,
                         mode = "excl")

      peak_grp <- .group_peaks(cov_scan = cov_scan,
                               peak_variant_id = peak_obj$peak_variant_id)

      new_peaks <- .get_newpeaks(object = object, peak_grp = peak_grp, peak_obj = peak_obj, n_threads = n_threads)
      new_peaks <- unique(unlist(new_peaks))

      if(setequal(new_peaks, peak_obj$peak_variant_id)){
        # peak_grp_summary <-  do.call("rbind",
        #                              lapply(seq_along(peak_grp), function(i){
        #                                data.frame(round = count,
        #                                           peak_variant_ID = peak_grp[[i]][1],
        #                                           member_peak = peak_grp[[i]],
        #                                           new_peak = new_peaks[[i]])
        #                              }))
        break

      } else {
        peak_obj <- .remake_peakobj(object = object,
                                    peak_obj = peak_obj,
                                    peak_variant_id = new_peaks)
      }

      peak_block <- lapply(X = peak_obj$peak_variant_id,
                           FUN = .recallPeak,
                           peak_obj = peak_obj,
                           object = object,
                           peakcall = peakcall,
                           values = FALSE,
                           n_threads = n_threads)
      names(peak_block) <- peak_obj$peak_variant_id
      peak_obj$peak_block <- peak_block
      if(length(new_peaks) == 1){
        break
      }
    }

    new_peak_blocks <- lapply(X = seq_along(peak_obj$peak_variant_id),
                              FUN = .make_newblocks,
                              object = object,
                              peak_obj = peak_obj,
                              n_threads = n_threads)

    peak_order <- sapply(new_peak_blocks, function(x) {max(x$peakNegLog10P)})
    new_peak_blocks <- new_peak_blocks[order(peak_order, decreasing = TRUE)]
    header <- sapply(new_peak_blocks, names)
    common <- apply(header, 1, function(x){length(unique(x)) == 1})
    common_header <- header[common]
    uncommon_header <- sort(unique(header[!common]))
    uncommon_header <- c(uncommon_header[grepl("^P.", uncommon_header)],
                         uncommon_header[grepl("^Coef.", uncommon_header)])
    new_peak_blocks <- lapply(seq_along(new_peak_blocks), function(k){
      new_peak_blocks[[k]]$peak_ID <- k
      common_df <- subset(new_peak_blocks[[k]], select = common_header)
      not_exist <- uncommon_header[!uncommon_header %in% names(new_peak_blocks[[k]])]
      not_exist_df <- matrix(NA, nrow = nrow(new_peak_blocks[[k]]), ncol = length(not_exist))
      not_exist_df <- data.frame(not_exist_df)
      names(not_exist_df) <- not_exist
      new_peak_blocks[[k]] <- cbind(new_peak_blocks[[k]], not_exist_df)
      uncommon_df <- subset(new_peak_blocks[[k]], select = uncommon_header)
      new_peak_blocks[[k]] <- cbind(common_df, uncommon_df)
      return(new_peak_blocks[[k]])
    })
    new_peak_blocks <- do.call("rbind", new_peak_blocks)
    out1 <- unique(subset(new_peak_blocks,
                          select = c(peak_ID, peak_variant_ID)))
    out2 <- unique(subset(new_peak_blocks,
                          select = c(peak_ID, variant_ID,
                                     Dist2peak, LD2peak,
                                     P.model:negLog10P)))

    out1_gdsn <- .create_gdsn(root_node = object$root,
                              target_node = "lazygas/recalc/peaks",
                              new_node = pheno_name,
                              val = as.matrix(out1),
                              storage = "uint32",
                              replace = TRUE)
    out2_gdsn <- .create_gdsn(root_node = object$root,
                              target_node = "lazygas/recalc/blocks",
                              new_node = pheno_name,
                              val = as.matrix(out2),
                              storage = "double",
                              replace = TRUE)
    # .create_gdsn(root_node = object$root,
    #              target_node = "lazygas/recalc/group",
    #              new_node = pheno_name,
    #              val = as.matrix(peak_grp_summary),
    #              storage = "uint32",
    #              replace = TRUE)

    put.attr.gdsn(node = out1_gdsn,
                  name = "col_names",
                  val = names(out1))
    put.attr.gdsn(node = out2_gdsn,
                  name = "col_names",
                  val = names(out2))
  }
}

# Create a peak object for recalculating associations.
.make_peakobj <- function(object, threshold, signif, pheno, pheno_name, peakcall, binary){
  peak_variant_id <- unique(peakcall$peak_variant_ID)

  geno_format <- .get_geno_format(object = object)

  selection <- .set_selection_criteria_for_recalc(object = object,
                                                  geno_format = geno_format,
                                                  peak_variant_id = peak_variant_id)

  geno <- .retrieve_geno(object = object,
                         selection = selection,
                         geno_format = geno_format,
                         collapse = FALSE)

  if(geno_format == "haplotype"){
    if(length(dim(geno)) == 3){
      geno[, , order(peak_variant_id)] <- geno
    }

  } else if(length(dim(geno)) == 2) {
    geno[, order(peak_variant_id)] <- geno

  } else {
    geno <- matrix(data = geno, ncol = 1)
  }

  if(geno_format == "haplotype"){
    na_val <- 0

  } else if(geno_format == "corrected"){
    na_val <- 3


  } else if(geno_format == "dosage"){
    na_val <- 63

  } else {
    na_val <- 2
  }

  out <- list(
    geno = geno,
    geno_format = geno_format,
    pheno = .standardize(val = pheno$pheno[, pheno_name], binary = binary),
    threshold = threshold,
    signif = signif,
    conv_fun = eval(parse(text = .get_data(object,
                                           node = "lazygas/scan/conv_fun"))),
    formula = .get_data(object, node = "lazygas/scan/formula"),
    na_val = na_val,
    peakcall = peakcall,
    peak_variant_id = peak_variant_id,
    peak_block = tapply(peakcall$variant_ID, peakcall$peak_variant_ID, c),
    snp_id = selection$snp_id,
    n_sample = nsam(object = object),
    binary = binary
  )
  return(out)
}

# Function to set the selection criteria based on the genotype format
.set_selection_criteria_for_recalc <- function(object,
                                               geno_format,
                                               peak_variant_id) {
  snp_id <- getMarID(object = object, valid = FALSE)

  if (geno_format == "dosage"){
    selection <- list(validSam(object = object),
                      snp_id %in% peak_variant_id)
    node <- "annotation/format/EDS/data"
    is_categorical <- FALSE

  } else {
    # Get the dimensions of the haplotype data
    obj_desp <- objdesp.gdsn(node = index.gdsn(node = object,
                                               path = "annotation/format/HAP/data"))

    # Set selection criteria for haplotype data
    selection <- list(rep(TRUE, obj_desp$dim[1]),
                      validSam(object = object),
                      snp_id %in% peak_variant_id)

    if(geno_format == "haplotype"){
      node <- "annotation/format/HAP/data"
      is_categorical <- TRUE

    } else {
      node <- ifelse(geno_format == "corrected",
                     "annotation/format/CGT/data",
                     "genotype/data")
      is_categorical <- FALSE
    }
  }

  return(list(node = node,
              selection = selection,
              snp_id = snp_id,
              is_categorical = is_categorical))
}

.composite <- function(query_peak,
                       subject_peak = NULL,
                       peak_obj,
                       output_all = FALSE,
                       n_threads,
                       mode){
  if(mode == "all"){
    null_df <- NULL
    subject_peak <- peak_obj$peak_variant_id
    target_geno <- peak_obj$geno

    if(peak_obj$geno_format == "haplotype"){
      target_geno <- apply(X = target_geno, MARGIN = 3, FUN = list)
      target_geno <- lapply(X = target_geno, FUN = "[[", 1)

    } else {
      target_geno <- apply(X = target_geno, MARGIN = 2, FUN = list)
      target_geno <- lapply(X = target_geno, FUN = "[[", 1)
    }


  } else if(mode == "composite"){
    null_df <- NULL
    for(j in query_peak){
      index <- peak_obj$peak_variant_id %in% j
      tmp_df <- .makeDF(g = peak_obj$geno[, index],
                        phe = peak_obj$pheno,
                        conv_fun = peak_obj$conv_fun,
                        formula = peak_obj$formula)
      if(is.null(null_df)){
        names(tmp_df$df)[-1] <- paste(names(tmp_df$df)[-1], j, sep = "_")
        null_df <- tmp_df

      } else {
        names(tmp_df$df)[-1] <- paste(names(tmp_df$df)[-1], j, sep = "_")
        null_df$df <- cbind(null_df$df, subset(tmp_df$df, select = -phe))
      }
      null_df$fml <- paste0("phe ~ ",
                            paste(names(null_df$df)[-1], collapse = " + "))
    }

    if(is.null(subject_peak)){
      index <- peak_obj$peak_variant_id %in% query_peak
      subject_peak <- peak_obj$peak_variant_id[!index]

    } else {
      index <- !peak_obj$peak_variant_id %in% subject_peak
    }

  }  else if(mode == "excl")  {
    index <- peak_obj$peak_variant_id %in% query_peak
    subject_peak <- peak_obj$peak_variant_id[!index]

    if(length(subject_peak) == 0){
      return(data.frame(variant_ID = NA, FDR = NA))
    }

    if(peak_obj$geno_format == "haplotype"){
      g <- peak_obj$geno[, , index]
    } else {
      g <- peak_obj$geno[, index]
    }
    null_df <- .makeDF(g = g,
                       phe = peak_obj$pheno,
                       conv_fun = peak_obj$conv_fun,
                       formula = peak_obj$formula)
  }

  if(peak_obj$geno_format == "haplotype"){
    target_geno <- peak_obj$geno[, , !index]
    if(is.vector(target_geno)){
      target_geno <- matrix(data = target_geno, ncol = 1)
    }
    target_geno <- apply(X = target_geno, MARGIN = 3, FUN = list)
    target_geno <- lapply(X = target_geno, FUN = "[[", 1)

  } else {
    target_geno <- peak_obj$geno[, !index]
    if(is.vector(target_geno)){
      target_geno <- matrix(data = target_geno, ncol = 1)
    }
    target_geno <- apply(X = target_geno, MARGIN = 2, FUN = list)
    target_geno <- lapply(X = target_geno, FUN = "[[", 1)
  }

  # Parallel processing setup
  if(is.null(n_threads)){
    cores <- detectCores()
    if(cores > 1){
      n_threads <- round(cores / 2)
    }
  }

  na_val <- peak_obj$na_val

  p_values <- mclapply(X = target_geno, mc.cores = n_threads, mc.preschedule = TRUE,
                       function(g){
                         g[g == na_val] <- NA
                         if(all(is.na(g))){
                           return(NA)
                         }
                         if(length(unique(na.omit(as.vector(g)))) == 1){
                           return(NA)  # Return NA if there is no variability in the genotype data
                         }
                         df <- .makeDF(g = g,
                                       phe = peak_obj$pheno,
                                       conv_fun = peak_obj$conv_fun,
                                       formula = peak_obj$formula)
                         if(length(query_peak) != 0){
                           tmp <- subset(null_df$df, select = -phe)
                           names(tmp) <- paste0("qtl_", names(tmp))
                           df$df <- cbind(df$df, tmp)
                           df$fml <- paste(c(df$fml, names(tmp)),
                                           collapse = " + ")
                         }

                         if(peak_obj$binary){
                           family <- "binomial"

                         } else {
                           family <- "gaussian"
                         }

                         out <- .doGLM(df = df,
                                       family = family,
                                       null_df = null_df)
                         return(out)
                       })

  if(is.list(p_values)){
    p_values <- data.frame(do.call("rbind", p_values))

  } else {
    p_values <- data.frame(t(p_values))
  }
  p_values$FDR <- p.adjust(p_values$P.model, "fdr")

  if(output_all){
    return(data.frame(variant_ID = subject_peak, p_values))

  } else {
    return(data.frame(variant_ID = subject_peak, FDR = p_values$FDR))
  }
}

.group_peaks <- function(cov_scan, peak_variant_id){
  peak_list <- peak_variant_id
  out <- NULL
  for(j in seq_along(cov_scan)){
    if(!peak_variant_id[j] %in% peak_list){
      next
    }
    grp <- c(peak_variant_id[j],
             cov_scan[[j]]$variant_ID[cov_scan[[j]]$FDR > 0.05])
    grp <- grp[grp %in% peak_list]
    out <- c(out, list(grp))
    peak_list <- peak_list[!peak_list %in% grp]
  }
  return(out)
}

.get_newpeaks <- function(object, peak_grp, peak_obj, n_threads){
  out <- NULL
  for(j in seq_along(peak_grp)){
    grp <- peak_grp[[j]]
    query_peak <- grp[1]
    out <- c(out, list(query_peak))
    subject_peak <- as.numeric(names(peak_obj$peak_block)) %in% grp
    subject_peak <- sort(unique(unlist(peak_obj$peak_block[subject_peak])))

    if(length(subject_peak) == 0){
      next
    }

    peak_obj <- .remake_peakobj(object = object,
                                peak_obj = peak_obj,
                                peak_variant_id = subject_peak)

    peak_scan <- .composite(query_peak = query_peak,
                            peak_obj = peak_obj,
                            n_threads = n_threads,
                            mode = "composite")

    near_peak <- peak_obj$peakcall$peak_variant_ID == query_peak
    max_cor <- max(peak_obj$peakcall$LD2peak[near_peak], na.rm = TRUE)
    near_peak_cor <- peak_obj$peakcall$LD2peak[near_peak] == max_cor
    cor1 <- peak_obj$peakcall$variant_ID[near_peak][near_peak_cor]
    peak_scan$FDR[peak_scan$variant_ID %in% cor1] <- 1

    other_peaks <- peak_scan[peak_scan$FDR < peak_obj$signif, ]
    other_peaks_hit <- peak_obj$peakcall$variant_ID %in% other_peaks$variant_ID
    other_peaks_id <- unique(peak_obj$peakcall$peak_variant_ID[other_peaks_hit])
    if(query_peak %in% other_peaks_id){
      query_peak_member <- peak_obj$peakcall$peak_variant_ID %in% query_peak
      in_query_peak <- other_peaks$variant_ID %in% peak_obj$peakcall$variant_ID[query_peak_member]
      in_query_peak <- other_peaks[in_query_peak, ]
      in_query_peak_id <- in_query_peak$variant_ID[which.min(in_query_peak$FDR)]
      other_peaks_id[other_peaks_id %in% query_peak] <- in_query_peak_id
    }
    if(length(other_peaks_id) != 0){
      out <- c(out, list(other_peaks_id))
    }
  }
  return(out)
}

.remake_peakobj <- function(object, peak_obj, peak_variant_id){
  peak_variant_id <- unique(peak_variant_id)
  peak_obj$peak_variant_id <- unique(peak_variant_id)

  selection <- .set_selection_criteria_for_recalc(object = object,
                                                  geno_format = peak_obj$geno_format,
                                                  peak_variant_id = peak_variant_id)

  geno <- .retrieve_geno(object = object,
                         selection = selection,
                         geno_format = peak_obj$geno_format,
                         collapse = FALSE)

  if(is.vector(geno)){
    geno <- matrix(geno, ncol = 1)
  }

  if(peak_obj$geno_format == "haplotype"){
    if(length(dim(geno)) == 3){
      geno[, , order(peak_variant_id)] <- geno
    }

  } else {
    geno[, order(peak_variant_id)] <- geno
  }

  peak_obj$geno <- geno
  return(peak_obj)
}

.recallPeak <- function(peak_id, peak_obj, object, peakcall, values, n_threads){
  target_variants <- peakcall$peak_ID %in% peakcall$peak_ID[peakcall$variant_ID %in% peak_id]
  out <- target_variants <- sort(peakcall$variant_ID[target_variants])

  if(length(target_variants) == 1){
    if(values){
      out <- data.frame(Dist2peak = 0,
                        LD2peak = 1,
                        variant_ID = out)
    }

  } else{
    variables <- .initialize_variables(object = object)
    index <- which(variables$snp_id == peak_id)
    chr <- variables$chr[index]
    chr_with_peak <- which(variables$chr == chr)
    id_in_chr <- variables$snp_id[chr_with_peak]
    peak_index_in_chr <- which(id_in_chr == peak_id)
    peak_info <-   list(index = index,
                        variantID = peak_id,
                        chr = chr,
                        chr_with_peak = chr_with_peak,
                        id_in_chr = id_in_chr,
                        peak_index_in_chr = peak_index_in_chr)

    # Set the selection criteria based on the genotype format
    selection <- .set_selection_criteria_for_peakcall(object = object,
                                                      geno_format = peak_obj$geno_format,
                                                      chr = peak_info$chr)

    # Retrieve genotype data based on the selection criteria
    geno <- .retrieve_geno(object = object,
                           selection = selection,
                           geno_format = peak_obj$geno_format)

    # Calculate linkage disequilibrium (LD) for the identified peak
    peak_ld <- .calculate_ld(geno = as.matrix(geno),
                             peak_variant_idx = peak_info$peak_index_in_chr,
                             is_categorical = selection$is_categorical,
                             n_threads = n_threads)

    # Identify the peak block based on the LD values
    peak_block <- .identify_peak_block(peak_ld = peak_ld,
                                       variables = variables,
                                       peak_info = peak_info,
                                       threshold = peak_obj$threshold,
                                       geno_format = peak_obj$geno_format)
    if(values){
      out <- data.frame(Dist2peak = peak_block$dist,
                        LD2peak = peak_block$ld,
                        variant_ID = peak_block$ids)

    } else {
      out <- peak_block$ids
    }
  }

  return(out)
}

.make_newblocks <- function(i, object, peak_obj, n_threads){
  query_peak <- peak_obj$peak_variant_id[-i]
  hit_id <- peak_obj$peakcall$variant_ID %in% peak_obj$peak_variant_id[i]
  i_peak <- peak_obj$peakcall[hit_id, ]
  variant_id <- getMarID(object = object, valid = TRUE, chr = i_peak$peak_Chr)
  peak_obj <- .remake_peakobj(object = object,
                              peak_obj = peak_obj,
                              peak_variant_id = sort(c(variant_id, query_peak)))
  peakcall <- .composite(query_peak = query_peak,
                         peak_obj = peak_obj,
                         output_all = TRUE,
                         n_threads = n_threads,
                         mode = "composite")
  peak_id <- peakcall$variant_ID[which.min(peakcall$P.model)]
  peakcall$peak_ID <- i
  recall_peak <- .recallPeak(peak_id = peak_id,
                             peak_obj = peak_obj,
                             object = object,
                             peakcall = peakcall,
                             values = TRUE,
                             n_threads = n_threads)

  peak_obj <- .remake_peakobj(object = object,
                              peak_obj = peak_obj,
                              peak_variant_id = c(query_peak,
                                                  recall_peak$variant_ID))
  peakcall <- .composite(query_peak = query_peak,
                         subject_peak = recall_peak$variant_ID,
                         peak_obj = peak_obj,
                         output_all = TRUE,
                         n_threads = n_threads,
                         mode = "composite")
  peakcall$negLog10P <- -log10(peakcall$P.model)
  out <- data.frame(peak_ID = i,
                    peak_variant_ID = peak_id,
                    peakNegLog10P = peakcall$negLog10P[peakcall$variant_ID == peak_id],
                    recall_peak, subset(peakcall, select = -variant_ID))
  return(out)
}

################################################################################

#' @export
setGeneric("listCandidate", function(object,
                                     gff_fn,
                                     snpeff_fn = NULL,
                                     ann = NULL,
                                     recalc = FALSE,
                                     numerize_chr = FALSE,
                                     ...)
  standardGeneric("listCandidate"))

setMethod("listCandidate",
          "LazyGas",
          function(object,
                   gff_fn,
                   snpeff_fn = NULL,
                   ann = NULL,
                   recalc = FALSE,
                   numerize_chr = FALSE){
            if(recalc){
              if(!exist.gdsn(node = object$root, path = "lazygas/recalc")){
                stop("No peakcall data in the input LazyGas object.\n",
                     "Run recalcAssoc() to recalculate associations.")
              }

            } else {
              if(!exist.gdsn(node = object$root, path = "lazygas/peakcall")){
                stop("No peakcall data in the input LazyGas object.\n",
                     "Run callPeakBlock() to call peak blocks.")
              }
            }

            gff <- .loadGFF(gff_fn = gff_fn, numerize_chr = numerize_chr)
            snpeff <- .loadSNPEff(snpeff_fn = snpeff_fn,
                                  numerize_chr = numerize_chr)

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

              # .finalize_gdsn_candidate(object = object,
              #                          pheno_name = pheno_name)
            }
          })

.create_candidate_folders <- function(object) {
  candidate_gdsn <- .create_gdsn(root_node = object$root,
                                 target_node = "lazygas",
                                 new_node = "candidate",
                                 is_folder = TRUE)

  snpeff_gdsn <- .create_gdsn(root_node = object$root,
                              target_node = "lazygas",
                              new_node = "snpeff",
                              is_folder = TRUE)
}

#' @importFrom rtracklayer import.gff
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
.loadGFF <- function(gff_fn, numerize_chr){
  if(!file.exists(gff_fn)){
    stop("The specified GFF file does not exist!",
         call. = FALSE)
  }
  gff <- import.gff(gff_fn)

  if(numerize_chr){
    seq_lev <- seqlevels(gff)
    num_seq_lev <- suppressWarnings(as.numeric(seq_lev))
    if(!all(!is.na(num_seq_lev))){
      seq_lev <- gsub("[^0-9]", "", seq_lev)
      missing <- seq_lev == ""
    }
    num_seq_lev <- suppressWarnings(as.numeric(seq_lev))
    if(!all(!is.na(num_seq_lev[!missing]))){
      seq_lev <- gsub("[^0-9]", "", seq_lev)
    }
    new_seq_lev <- as.numeric(seq_lev)
    new_seq_lev[missing] <- seqlevels(gff)[missing]
    seqlevels(gff) <- new_seq_lev
  }

  return(gff)
}

#' @importFrom vcfR read.vcfR
.loadSNPEff <- function(snpeff_fn, numerize_chr){
  if(is.null(snpeff_fn)){
    return(NULL)
  } else {
    if(!file.exists(snpeff_fn)){
      stop("The specified SNPEff file does not exist!",
           call. = FALSE)
    }
    snpeff <- read.vcfR(snpeff_fn)
  }

  if(numerize_chr){
    snpeff_fix <- as.data.frame(snpeff@fix)
    seq_lev <- unique(snpeff_fix$CHROM)
    num_seq_lev <- suppressWarnings(as.numeric(seq_lev))
    if(!all(!is.na(num_seq_lev))){
      seq_lev <- gsub("[^0-9]", "", seq_lev)
      missing <- seq_lev == ""
    }
    num_seq_lev <- suppressWarnings(as.numeric(seq_lev))
    if(!all(!is.na(num_seq_lev[!missing]))){
      seq_lev <- gsub("[^0-9]", "", seq_lev)
    }
    new_seq_lev <- as.numeric(seq_lev)
    new_seq_lev[missing] <- seqlevels(gff)[missing]
    seq_lev <- unique(snpeff_fix$CHROM)
    hit <- match(snpeff_fix$CHROM, seq_lev)
    snpeff_fix$CHROM <- new_seq_lev[hit]
    snpeff@fix <- as.matrix(snpeff_fix)
  }
  return(snpeff)
}

.loadANN <- function(ann){
  if(!is.null(ann)){
    if(!file.exists(ann)){
      stop("The specified annotation file does not exist!",
           call. = FALSE)
    }
    ann <- read.csv(ann)
  }
  return(ann)
}

#' ## Sub-function to finalize GDS nodes
#' #' @importFrom gdsfmt index.gdsn readmode.gdsn
#' .finalize_gdsn_candidate <- function(object, pheno_name) {
#'   readmode.gdsn(node = index.gdsn(node = object$root,
#'                                   path = paste0("lazygas/candidate/", pheno_name)))
#'   readmode.gdsn(node = index.gdsn(node = object$root,
#'                                   path = paste0("lazygas/snpeff/", pheno_name)))
#' }

#' @importFrom dplyr left_join
.candidatelistor <- function(object,
                             ann,
                             gff,
                             snpeff,
                             pheno_name,
                             recalc){
  peakcall <- .get_peakcall(object = object,
                            pheno_name = pheno_name,
                            recalc = recalc)

  if(is.null(peakcall)){
    .create_gdsn(root_node = object$root,
                 target_node = "lazygas/candidate",
                 new_node = pheno_name,
                 storage = "string",
                 replace = TRUE)

  } else {
    candidate_list <- NULL
    snpeff_list <- NULL
    for(i_peak in unique(peakcall$peak_ID)){
      peakblock <- peakcall[peakcall$peak_ID == i_peak, ]
      tmp <- .getCandidate(peakblock = peakblock,
                           gff = gff,
                           snpeff = snpeff)
      candidate_list <- rbind(candidate_list, tmp$candidate_list)
      snpeff_list <- rbind(snpeff_list, tmp$snpeff_out)
    }

    if(!is.null(candidate_list)){
      if(!is.null(ann)){
        candidate_list <- left_join(candidate_list, ann, by = "Gene_ID")

      }
      candidate_list[is.na(candidate_list)] <- ""
    }

    if(!is.null(candidate_list)){
      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/candidate",
                   new_node = pheno_name,
                   val = as.matrix(candidate_list),
                   storage = "string",
                   replace = TRUE)

      put.attr.gdsn(node = index.gdsn(node = object$root,
                                      path = "lazygas/candidate"),
                    name = "col_names",
                    val = colnames(candidate_list))

      readmode.gdsn(node = index.gdsn(node = object$root,
                                      path = paste0("lazygas/candidate/", pheno_name)))
    }
    if(!is.null(snpeff_list)){
      snpeff_list[is.na(snpeff_list)] <- ""

      .create_gdsn(root_node = object$root,
                   target_node = "lazygas/snpeff",
                   new_node = pheno_name,
                   val = as.matrix(snpeff_list),
                   storage = "string",
                   replace = TRUE)

      put.attr.gdsn(node = index.gdsn(node = object$root,
                                      path = "lazygas/snpeff"),
                    name = "col_names",
                    val = colnames(snpeff_list))

      readmode.gdsn(node = index.gdsn(node = object$root,
                                      path = paste0("lazygas/snpeff/", pheno_name)))
    }
  }
}

#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
.getCandidate <- function(peakblock, gff, snpeff){
  peak_variant_id <- peakblock$peak_variant_ID[1]
  is_peak <- peakblock$variant_ID == peak_variant_id
  peak_gff <- GRanges(seqnames = peakblock$Chr[is_peak],
                      ranges = IRanges(start = peakblock$Pos[1],
                                       end = tail(peakblock$Pos, 1)))
  gene_gff <- gff[gff$type %in% "gene"]
  hit <- gene_gff[queryHits(findOverlaps(gene_gff, peak_gff))]
  hit <- hit[order(start(hit))]
  hit$dist2peak <- start(hit) - peakblock$Pos[is_peak]
  nearest_p <- sapply(start(hit), function(x){
    return(peakblock$negLog10P[which.min(abs(peakblock$Pos - x))])
  })
  if(length(hit) == 0){
    out <- NULL

  } else {
    out <- data.frame(peak_ID = peakblock$peak_ID[1],
                      Gene_ID = hit$ID,
                      Gene_chr = as.character(seqnames(hit)),
                      Gene_start = start(hit),
                      Dist2peak = hit$dist2peak,
                      negLog10P = nearest_p)
  }

  if(!is.null(snpeff)){
    snpeff_out <- .addSnpEff(snpeff = snpeff, peakblock = peakblock)
    collapse_snpeff_out <- .collapseSnpEff(snpeff_out = snpeff_out)
    out <- left_join(out, collapse_snpeff_out, by ="Gene_ID")
    out <- list(candidate_list = out, snpeff_out = snpeff_out)

  } else {
    out <- list(candidate_list = out, snpeff_out = NULL)
  }
  return(out)
}

#' @importFrom vcfR getFIX getINFO
.addSnpEff <- function(snpeff, peakblock){
  snpeff_fix <- as.data.frame(snpeff@fix)
  target_chr <- snpeff_fix$CHROM %in% peakblock$Chr
  target_pos <-  snpeff_fix$POS >= min(peakblock$Pos) &  snpeff_fix$POS <= max(peakblock$Pos)
  peak_ann <- snpeff_fix[target_chr & target_pos, ]

  if(nrow(peak_ann) == 0){
    out <- data.frame(t(rep(NA, 19)))
    names(out) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                    "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                    "Rank", "HGVS.c", "HGVS.p", "Pos.in.tx", "Pos.in.CDS",
                    "Pos.in.AA", "Distance", "INFO", "Chr", "Pos", "negLog10P")

  } else {
    var_id <- paste0(peak_ann$CHROM, "_", peak_ann$POS, "_")
    out <- unlist(peak_ann$INFO)
    out <- sub(".+ANN=", "", out)
    names(out) <- var_id
    out <- unlist(strsplit(out, ","))
    out <- sub(";.+", "", out)
    out[grepl("\\|$", out)] <- paste(out[grepl("\\|$", out)], "")
    out <- strsplit(out, "\\|")
    len <- sapply(out, length)
    out <- out[len %in% 16]
    if(length(out) == 0){
      return(NULL)
    }
    out <- do.call("rbind", out)
    out <- out[!out[, 2] %in% c("custom", "intergenic_region"), ]
    if(length(out) == 0){
      return(NULL)
    }
    out <- as.data.frame(out)
    colnames(out) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name",
                       "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
                       "Rank", "HGVS.c", "HGVS.p", "Pos.in.tx", "Pos.in.CDS",
                       "Pos.in.AA", "Distance", "INFO")
    chr_pos <- sub("_[0-9]$", "", rownames(out))
    out$Chr <- sub("_.+", "", chr_pos)
    out$Pos <- sub(".+_", "", chr_pos)
    hit <- match(out$Pos, peakblock$Pos)
    out$negLog10P <- peakblock$negLog10P[hit]
    return(out)
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
                 t(as.vector(table(factor(snpeff_out$Annotation_Impact[i][!is.na(snpeff_out$negLog10P)], impacts)))))
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

################################################################################
#' Generate interactive summary
#'
#' @param object QTLscan object
#' @param pehno Phenotype names to be drawn
#' @param out_fn Prefix of output file
#'
#' @import reactable
#' @import htmltools
#' @import plotly
#' @export
#'
makeInteractiveSummary <- function(object, pheno, out_fn){
  plot_pheno <- plotPheno(object = object, pheno = pheno, xlab = pheno, boxplot = FALSE)
  plot_pheno <- ggplotly(plot_pheno)
  plot_man <- plotManhattan(object = object, pheno = pheno) + labs(title = pheno)
  plot_man <- ggplotly(plot_man)
  tag_list <- tagList(div(plot_pheno, style = "margin:auto;width:80vw;"),
                      div(plot_man, style = "margin:auto;width:80vw;"))

  peak1 <- lazyData(object = object, dataset = "peakcall", pheno = pheno)
  if(!is.null(peak1)){
    plot_peak1 <- plotPeaks(object = object, pheno = pheno, recalc = FALSE) + labs(title = pheno)
    tag_list <- tagList(tag_list,
                        div(ggplotly(plot_peak1), style = "margin:auto;width:80vw;"))
  }

  peak2 <- lazyData(object = object, dataset = "recalc", pheno = pheno)
  if(!is.null(peak2)){
    plot_peak2 <- plotPeaks(object = object, pheno = pheno, recalc = TRUE) + labs(title = pheno)
    tag_list <- tagList(tag_list,
                        div(ggplotly(plot_peak2), style = "margin:auto;width:80vw;"))
  }

  plot_hap <- haploPlot(object = object, pheno = pheno, recalc = TRUE)
  if(!is.null(plot_hap)){
    for(i in seq_along(plot_hap)){
      tag_list <- tagList(tag_list,
                          div(ggplotly(plot_hap[[i]]), style = "margin:auto;width:60vw;"))
    }
  }

  candidate <- lazyData(object = object, dataset = "candidate", pheno = pheno)
  if(length(candidate) == 0){
    candidate <- NULL
  }
  if(!is.null(candidate)){
    col_names <- colnames(candidate)
    col_def <- NULL
    for(i in seq_along(col_names)){
      min_width <- nchar(col_names[i]) * 16
      if(min_width < 100){
        min_width <- 100
      }
      col_def <- c(col_def, list(colDef(minWidth = min_width)))
    }
    names(col_def) <- col_names
    table <- reactable(data = candidate, columns = col_def, sortable = TRUE,
                       resizable = TRUE, filterable = TRUE,
                       searchable = TRUE, showPageSizeOptions = TRUE,
                       wrap = TRUE, striped = TRUE, onClick = JS())
    tag_list <- tagList(tag_list,
                        div(table, style = "margin:auto;width:90vw;"))
  }

  save_html(tag_list, file = out_fn)
}

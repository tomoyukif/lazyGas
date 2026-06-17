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

#' Assign phenotype data to a LazyGas object
#'
#' @param object A LazyGas object.
#' @param pheno_name A character vector.
#' @param p_values A numeric matrix or data.frame with the number of rows matching the number of markers.
#' @importMethodsFrom GBScleanR nmar getMarID getChromosome getPosition
#'
#' @export
#'
setGeneric("assignPvalues", function(object, pheno_name, p_values,
                                     geno_format = c("genotype", "dosage", "haplotype"),
                                     conv_fun = NULL,
                                     formula = "phe ~ add",
                                     null_formula = NULL,
                                     kruskal = NULL, ...)
  standardGeneric("assignPvalues"))

setMethod("assignPvalues",
          "LazyGas",
          function(object, pheno_name, p_values,
                   geno_format = c("genotype", "dosage", "haplotype"),
                   conv_fun = NULL, formula = "phe ~ add", null_formula = NULL, kruskal = NULL){
            geno_format <- match.arg(arg = geno_format, choices = c("genotype", "dosage", "haplotype"))

            if(is.null(conv_fun)){
              formula <- formula("phe ~ add")  # Define the formula for regression
              null_formula <- NULL
            }

            hit <- pheno_name %in% object@lazydata$pheno_names
            if(any(!hit)){
              stop('Following phenotype name(s) was not found in the input LazyGas object: ',
                   paste(pheno_name[!hit], collapse = ", "),
                   call. = FALSE)
            }
            if(is.vector(p_values)){
              p_values <- matrix(p_values, ncol = 1)
            }

            nrow_p_values <- nrow(p_values)
            lg_nmar <- nmar(object)
            if(lg_nmar != nrow_p_values){
              stop('The input p_values has ',
                   nrow_p_values, ' rows',
                   ', but the marker number assigned in the lazyGas object is ',
                   lg_nmar)
            }

            # Add a folder to store scan results in the GDS object
            .create_scan_folder(object)

            # Loop through each phenotype name and perform analysis
            for(i in seq_along(pheno_name)){
              i_pheno_names <- pheno_name[i]

              p_values_i <- p_values[, i]
              p_values_i <- cbind(P.model = p_values_i,
                                  FDR = p.adjust(p = p_values_i, method = "fdr"),
                                  negLog10P = -log10(p_values_i))

              .store_write_scan(
                object = object,
                pheno_name = i_pheno_names,
                mat = p_values_i,
                colnames_stats = colnames(p_values_i)
              )
            }

            .store_additional_info(object = object,
                                   kruskal = kruskal,
                                   formula = formula,
                                   null_formula = null_formula,
                                   fixed_effect = NULL,
                                   conv_fun = conv_fun,
                                   geno_format = geno_format)
            return(object)
          }
)

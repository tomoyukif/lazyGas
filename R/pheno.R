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

            .store_save_pheno_snapshot(object)

            return(object)
          }
)

#' Restore phenotype metadata from the companion store
#'
#' Re-opening a GDS file does not reload phenotype data assigned by
#' [assignPheno()] in a previous session. This function restores phenotype
#' names (and values when a snapshot was saved) from the companion Parquet/SQLite
#' store so that [getPheno()], [lazyData()], and explorer tools work again.
#'
#' @param object A \code{LazyGas} object.
#' @param warn_stub If \code{TRUE}, warn when only trait names are restored
#'   without phenotype values.
#'
#' @return The \code{LazyGas} object with phenotype metadata attached.
#' @export
restorePhenoFromStore <- function(object, warn_stub = TRUE) {
  if (!inherits(object, "LazyGas")) {
    stop("'object' must be a LazyGas object.", call. = FALSE)
  }

  existing <- getPheno(object)$pheno_names
  if (!is.null(existing) && length(existing)) {
    return(object)
  }

  snap <- .store_read_pheno_snapshot(object)
  if (!is.null(snap) && nrow(snap) > 0L && "sample_id" %in% names(snap)) {
    sam <- getSamID(object)
    id_match <- match(sam, snap$sample_id)
    if (any(is.na(id_match))) {
      warning(
        "Phenotype snapshot sample IDs do not fully match the GDS; ",
        "missing samples will be NA.",
        call. = FALSE
      )
    }
    pheno_cols <- setdiff(names(snap), "sample_id")
    pheno <- snap[id_match, pheno_cols, drop = FALSE]
    rownames(pheno) <- sam

    meta <- if (.store_is_gds(object)) list() else .store_read_meta(object)
    if (!is.null(meta$pheno_type)) {
      pheno_type <- meta$pheno_type
    } else {
      pheno_type <- .checkPhenoType(pheno = pheno)
    }

    object@lazydata$pheno <- pheno
    object@lazydata$pheno_type <- pheno_type
    object@lazydata$pheno_names <- pheno_cols
    return(object)
  }

  pheno_names <- .store_list_pheno_names(object)
  if (!length(pheno_names)) {
    warning(
      "No phenotype information found in the companion store.",
      call. = FALSE
    )
    return(object)
  }

  if (warn_stub) {
    message(
      "Restored phenotype names from companion store (values not available). ",
      "Re-run assignPheno() to attach phenotype data for new analyses."
    )
  }

  n <- nsam(object)
  pheno <- as.data.frame(
    setNames(
      rep(list(rep(NA_real_, n)), length(pheno_names)),
      pheno_names
    ),
    stringsAsFactors = FALSE
  )
  pheno_type <- data.frame(
    binary = rep(FALSE, length(pheno_names)),
    normality = rep(NA_real_, length(pheno_names)),
    row.names = pheno_names,
    stringsAsFactors = FALSE
  )
  object@lazydata$pheno <- pheno
  object@lazydata$pheno_type <- pheno_type
  object@lazydata$pheno_names <- pheno_names
  object
}

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

#' Register externally computed association results in the companion store
#'
#' Use this when p-values (and optionally coefficients or other per-marker
#' columns) come from another GWAS tool, not from [scanAssoc()]. Each call
#' writes a **new** scan table for \code{pheno_name}; it does not merge into
#' an existing \code{scanAssoc()} result. To keep downstream plots such as
#' [haploPlot()] informative, pass \code{coef} and/or \code{any_data} together
#' with \code{p_values}.
#'
#' @param object A \code{LazyGas} object.
#' @param pheno_name Character vector of phenotype name(s).
#' @param p_values Numeric vector, matrix, or data.frame. Rows must equal
#'   \code{nmar(object)}. A vector is \code{P.model} for one phenotype; matrix
#'   columns map to \code{pheno_name} in order. \code{FDR} and \code{negLog10P}
#'   are derived from \code{P.model}.
#' @param coef Optional effect sizes to store with the scan. A numeric vector
#'   (length \code{nmar}) becomes \code{Coef.add}. A matrix or \code{data.frame}
#'   with \code{nrow = nmar} uses its column names (e.g. \code{Coef.add},
#'   \code{Coef.dom}). When \code{p_values} has multiple columns, \code{coef}
#'   may have the same number of columns (one per phenotype) or a single
#'   coefficient block shared across phenotypes.
#' @param any_data Optional \code{data.frame} with \code{nrow = nmar}. Its
#'   columns are appended to the stored scan table (e.g. external test
#'   statistics). Column names must not duplicate \code{P.model}, \code{FDR},
#'   \code{negLog10P}, or \code{coef} columns.
#' @param geno_format Genotype format metadata for downstream steps.
#' @param conv_fun,formula,null_formula,kruskal Scan metadata recorded for
#'   peak calling and recalculation.
#' @param ... Unused.
#'
#' @return The updated \code{LazyGas} object (invisibly).
#'
#' @seealso [scanAssoc()], [lazyData()]
#'
#' @export
#'
setGeneric("assignPvalues", function(object, pheno_name, p_values,
                                     coef = NULL,
                                     any_data = NULL,
                                     geno_format = c("genotype", "dosage", "haplotype"),
                                     conv_fun = NULL,
                                     formula = "phe ~ add",
                                     null_formula = NULL,
                                     kruskal = NULL, ...)
  standardGeneric("assignPvalues"))

.assignpvalues_as_matrix <- function(x, n_markers, label = "p_values") {
  if (is.data.frame(x)) {
    mat <- as.matrix(x)
  } else if (is.matrix(x)) {
    mat <- x
  } else {
    mat <- matrix(as.numeric(x), ncol = 1L)
  }
  if (nrow(mat) != n_markers) {
    stop(
      "The input ", label, " has ", nrow(mat), " rows, but nmar(object) is ",
      n_markers, ".",
      call. = FALSE
    )
  }
  mat
}

.assignpvalues_coef_block <- function(coef, pheno_idx, n_markers, n_pheno) {
  if (is.null(coef)) {
    return(NULL)
  }
  if (is.vector(coef) && !is.data.frame(coef)) {
    if (n_pheno > 1L) {
      stop(
        "A numeric vector 'coef' requires a single phenotype in 'pheno_name'. ",
        "Use a matrix with one column per phenotype, or call assignPvalues() ",
        "once per trait.",
        call. = FALSE
      )
    }
    if (length(coef) != n_markers) {
      stop(
        "'coef' length (", length(coef), ") must equal nmar(object) (",
        n_markers, ").",
        call. = FALSE
      )
    }
    return(data.frame(Coef.add = as.numeric(coef), check.names = FALSE))
  }

  cm <- as.data.frame(coef, check.names = FALSE)
  if (nrow(cm) != n_markers) {
    stop(
      "'coef' has ", nrow(cm), " rows, but nmar(object) is ", n_markers, ".",
      call. = FALSE
    )
  }
  if (ncol(cm) == n_pheno && n_pheno > 1L) {
    return(cm[, pheno_idx, drop = FALSE])
  }
  if (ncol(cm) == 1L && (is.null(names(cm)) || !nzchar(names(cm)[1L]))) {
    names(cm)[1L] <- "Coef.add"
  }
  cm
}

.assignpvalues_build_scan_mat <- function(p_col,
                                         coef = NULL,
                                         pheno_idx = 1L,
                                         n_markers,
                                         n_pheno,
                                         any_data = NULL) {
  p_col <- as.numeric(p_col)
  if (length(p_col) != n_markers) {
    stop(
      "'p_values' column length (", length(p_col),
      ") must equal nmar(object) (", n_markers, ").",
      call. = FALSE
    )
  }

  out <- data.frame(
    P.model = p_col,
    FDR = p.adjust(p = p_col, method = "fdr"),
    negLog10P = -log10(p_col),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  coef_block <- .assignpvalues_coef_block(
    coef = coef,
    pheno_idx = pheno_idx,
    n_markers = n_markers,
    n_pheno = n_pheno
  )
  if (!is.null(coef_block)) {
    dup <- intersect(names(out), names(coef_block))
    if (length(dup)) {
      stop(
        "Column name overlap between derived scan columns and 'coef': ",
        paste(dup, collapse = ", "),
        call. = FALSE
      )
    }
    out <- cbind(out, coef_block)
  }

  if (!is.null(any_data)) {
    if (!is.data.frame(any_data)) {
      stop("'any_data' must be a data.frame.", call. = FALSE)
    }
    if (nrow(any_data) != n_markers) {
      stop(
        "'any_data' has ", nrow(any_data), " rows, but nmar(object) is ",
        n_markers, ".",
        call. = FALSE
      )
    }
    dup <- intersect(names(out), names(any_data))
    if (length(dup)) {
      stop(
        "Column name overlap between scan columns and 'any_data': ",
        paste(dup, collapse = ", "),
        call. = FALSE
      )
    }
    out <- cbind(out, any_data)
  }

  as.matrix(out)
}

setMethod("assignPvalues",
          "LazyGas",
          function(object, pheno_name, p_values,
                   coef = NULL,
                   any_data = NULL,
                   geno_format = c("genotype", "dosage", "haplotype"),
                   conv_fun = NULL,
                   formula = "phe ~ add",
                   null_formula = NULL,
                   kruskal = NULL,
                   ...) {
            geno_format <- match.arg(
              arg = geno_format,
              choices = c("genotype", "dosage", "haplotype")
            )

            if (is.null(conv_fun)) {
              formula <- formula("phe ~ add")
              null_formula <- NULL
            }

            hit <- pheno_name %in% object@lazydata$pheno_names
            if (any(!hit)) {
              stop(
                "Following phenotype name(s) was not found in the input LazyGas object: ",
                paste(pheno_name[!hit], collapse = ", "),
                call. = FALSE
              )
            }

            p_mat <- .assignpvalues_as_matrix(
              x = p_values,
              n_markers = nmar(object),
              label = "p_values"
            )
            n_pheno <- length(pheno_name)
            if (ncol(p_mat) == 1L && n_pheno > 1L) {
              stop(
                "p_values has one column but pheno_name has ", n_pheno,
                " traits. Supply a matrix with one column per phenotype.",
                call. = FALSE
              )
            }
            if (ncol(p_mat) != n_pheno) {
              stop(
                "p_values has ", ncol(p_mat), " column(s) but pheno_name has ",
                n_pheno, " trait(s).",
                call. = FALSE
              )
            }

            if (!is.null(coef) && (is.matrix(coef) || is.data.frame(coef))) {
              coef_n <- ncol(as.data.frame(coef))
              if (coef_n > 1L && coef_n != n_pheno && n_pheno > 1L) {
                stop(
                  "coef has ", coef_n, " columns; expected 1 (shared block) or ",
                  n_pheno, " (one per phenotype).",
                  call. = FALSE
                )
              }
            }

            .create_scan_folder(object)

            for (i in seq_along(pheno_name)) {
              i_pheno_names <- pheno_name[i]
              p_values_i <- .assignpvalues_build_scan_mat(
                p_col = p_mat[, i],
                coef = coef,
                pheno_idx = i,
                n_markers = nmar(object),
                n_pheno = n_pheno,
                any_data = any_data
              )

              .store_write_scan(
                object = object,
                pheno_name = i_pheno_names,
                mat = p_values_i,
                colnames_stats = colnames(p_values_i)
              )
            }

            .store_additional_info(
              object = object,
              kruskal = kruskal,
              formula = formula,
              null_formula = null_formula,
              fixed_effect = NULL,
              conv_fun = conv_fun,
              geno_format = geno_format
            )
            invisible(object)
          }
)

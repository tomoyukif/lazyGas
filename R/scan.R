################################################################################
#' Extract genotype values for a single marker
#'
#' Returns genotypes for one marker. Among valid markers,
#' \code{marker_index} is a 1-based index (not a variant ID).
#'
#' @param object A \code{LazyGas} object.
#' @param geno_format One of \code{"genotype"}, \code{"corrected"},
#'   \code{"dosage"}, or \code{"haplotype"}.
#' @param marker_index Integer index of the marker among
#'   \code{validMar(object)} markers (default \code{1L}).
#' @param ... Unused. Legacy \code{marker_id=} is accepted with a warning and
#'   treated as \code{marker_index}.
#'
#' @return
#' \describe{
#'   \item{dosage}{Integer vector of length \code{nsam(object)}.}
#'   \item{genotype / haplotype}{Matrix with dimensions
#'     \code{c(ploidy, nsam)} (rows = alleles, columns = samples).}
#' }
#'
#' @export
#'
setGeneric("getGenoPerMarker",
           function(object, geno_format, marker_index = 1L, ...)
             standardGeneric("getGenoPerMarker"))

setMethod("getGenoPerMarker",
          "LazyGas",
          function(object,
                   geno_format = c("genotype", "corrected", "dosage", "haplotype"),
                   marker_index = 1L,
                   ...) {
            geno_format <- match.arg(
              arg = geno_format,
              choices = c("genotype", "corrected", "dosage", "haplotype")
            )

            dots <- list(...)
            if (!is.null(dots$marker_id)) {
              warning(
                "Argument 'marker_id' is deprecated; use 'marker_index' instead.",
                call. = FALSE
              )
              marker_index <- dots$marker_id
            }

            marker_index <- as.integer(marker_index)
            if (length(marker_index) != 1L || is.na(marker_index) ||
                marker_index < 1L) {
              stop("'marker_index' must be a single positive integer.",
                   call. = FALSE)
            }

            n_valid <- sum(validMar(object = object))
            if (marker_index > n_valid) {
              stop(
                "'marker_index' (", marker_index,
                ") exceeds the number of valid markers (", n_valid, ").",
                call. = FALSE
              )
            }

            path <- switch(
              geno_format,
              "genotype" = "genotype/data",
              "corrected" = "annotation/format/CGT/data",
              "dosage" = "annotation/format/EDS/data",
              "haplotype" = "annotation/format/HAP/data"
            )

            if (geno_format == "dosage") {
              selection <- list(
                validSam(object = object),
                validMar(object = object)
              )
              selection[[2]][selection[[2]]][-marker_index] <- FALSE
            } else {
              target_node_index <- index.gdsn(node = object, path = path)
              obj_desp <- objdesp.gdsn(node = target_node_index)
              selection <- list(
                rep(TRUE, obj_desp$dim[1]),
                validSam(object = object),
                validMar(object = object)
              )
              selection[[3]][selection[[3]]][-marker_index] <- FALSE
            }

            out <- .get_data(object = object, node = path, sel = selection)
            if (identical(geno_format, "haplotype")) {
              out[out == 0] <- NA
            } else if (identical(geno_format, "dosage")) {
              out[out == 63] <- NA
            } else {
              out[out == 3] <- NA
            }
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
#' Continuous phenotypes are standardized internally before regression. Stored
#' \code{Coef.*} columns are rescaled to the original phenotype units (change in
#' phenotype per genotype unit). Binary phenotypes are not standardized; their
#' coefficients remain on the GLM scale (log-odds).
#'
#' @param object A LazyGas object
#' @param out_fn Prefix of output file name
#' @param formula The formula of the regression model
#' @param conv_fun The function to convert genotype data to a model matrix for the regression
#' @param fixed_effect A data.frame of fixed effect(s) incorporated in the regression model. The column names must match the terms specified in `formula` and `null_formula.` The number of rows must match the number of samples.
#'
#' @importFrom gdsfmt apply.gdsn objdesp.gdsn
#'
#' @export
#'
setGeneric("scanAssoc", function(object,
                                 formula = "phe ~ add",
                                 null_formula = NULL,
                                 conv_fun = NULL,
                                 fixed_effect = NULL,
                                 geno_format = c("genotype", "corrected", "dosage", "haplotype"),
                                 kruskal = NULL,
                                 method = c("glm", "mlm"),
                                 ...)
  standardGeneric("scanAssoc"))

## Define the scanAssoc method for the LazyGas class
setMethod("scanAssoc",
          "LazyGas",
          function(object,
                   formula = "phe ~ add",
                   null_formula = NULL,
                   conv_fun = NULL,
                   fixed_effect = NULL,
                   geno_format = c("genotype", "corrected", "dosage", "haplotype"),
                   kruskal = NULL,
                   method = c("glm", "mlm")){

            method <- match.arg(arg = method, choices = c("glm", "mlm"))
            geno_format <- match.arg(arg = geno_format, choices = c("genotype", "corrected", "dosage", "haplotype"))

            if(method == "mlm"){
              if(geno_format %in% c("dosage", "haplotype")){
                stop('When method = "mlm", geno_format can only be "genotype" or "corrected".',
                     call. = FALSE)
              }
              conv_fun <- NULL
            }

            if(is.null(conv_fun)){
              formula <- "add"  # Define the formula for regression
              null_formula <- NULL
            }

            if(!is.null(null_formula)){
              formula_terms <- unlist(strsplit(formula, "\\+|\\*|\\:"))
              formula_terms <- gsub("\\s", "", formula_terms)
              null_formula_terms <- unlist(strsplit(null_formula, "\\+|\\*|\\:"))
              null_formula_terms <- gsub("\\s", "", null_formula_terms)
              check <- all(null_formula_terms %in% formula_terms)
              if(!check){
                stop("All terms in null_formula must be included in formula.",
                     call. = FALSE)
              }
              if(is.null(fixed_effect)){
                terms <- names(conv_fun(1))
                check <- all(null_formula_terms %in% terms)
                if(!check){
                  stop("The output of conv_fun does not contain term(s) appeared in null_formula.\nYou may need to specify fixed_effect", call. = FALSE)
                }
              }
            }

            if(!is.null(fixed_effect)){
              check <- all(apply(fixed_effect, 2, is.numeric))
              if(!check){
                stop("All columns in fixed_effect must be numeric.",
                     call. = FALSE)
              }

              terms <- unlist(strsplit(formula, "\\+|\\*|\\:"))
              terms <- gsub("\\s", "", terms)
              if(!is.null(null_formula)){
                null_formula_terms <- unlist(strsplit(null_formula, "\\+|\\*|\\:"))
                null_formula_terms <- gsub("\\s", "", null_formula_terms)
                terms <- c(terms, null_formula_terms)
              }
              check <- all(names(fixed_effect) %in% terms)
              if(!check){
                stop("All terms in fixed_effect must be included in formula and null_formula.",
                     call. = FALSE)
              }
            }

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
              i_pheno_raw <- getPheno(object = object)$pheno[, i_pheno_names]
              pheno_scale <- if (isTRUE(binary)) {
                1
              } else {
                stats::sd(i_pheno_raw, na.rm = TRUE)
              }
              i_pheno <- .standardize(val = i_pheno_raw, binary = binary)

              # Perform regression analysis and store results
              .perform_regression(object = object,
                                  i_pheno = i_pheno,
                                  geno_format = geno_format,
                                  conv_fun = conv_fun,
                                  formula = formula,
                                  null_formula = null_formula,
                                  fixed_effect = fixed_effect,
                                  dokruskal = dokruskal,
                                  i_pheno_names = i_pheno_names,
                                  binary = binary,
                                  method = method,
                                  pheno_scale = pheno_scale)
            }

            .store_flush_scan_meta(object = object)

            # Add additional information to the root node in the GDS object
            .store_additional_info(object = object,
                                   kruskal = kruskal,
                                   formula = formula,
                                   null_formula = null_formula,
                                   fixed_effect = fixed_effect,
                                   conv_fun = conv_fun,
                                   geno_format = geno_format)
          }
)

## Sub-function to create the scan folder in the store
.create_scan_folder <- function(object) {
  .store_ensure_scan_folder(object)
}

.standardize <- function(val, binary){
  if(binary){
    return(val)

  } else {
    return((val - mean(val, na.rm = TRUE)) / sd(val, na.rm = TRUE))
  }
}

.rescale_scan_coefs <- function(mat, pheno_scale) {
  if (is.null(mat)) {
    return(mat)
  }
  if (length(pheno_scale) != 1L || !is.finite(pheno_scale) || pheno_scale <= 0) {
    return(mat)
  }
  if (abs(pheno_scale - 1) < .Machine$double.eps^0.5) {
    return(mat)
  }
  coef_cols <- grep("^Coef\\.", colnames(mat), value = TRUE)
  if (length(coef_cols) == 0L) {
    return(mat)
  }
  mat[, coef_cols] <- mat[, coef_cols, drop = FALSE] * pheno_scale
  mat
}

## Sub-function to perform regression analysis and store results
#' @importFrom gdsfmt apply.gdsn
#' @importFrom gaston as.bed.matrix GRM association.test
.perform_regression <- function(object,
                                i_pheno,
                                geno_format,
                                conv_fun,
                                formula,
                                null_formula,
                                fixed_effect,
                                dokruskal,
                                i_pheno_names,
                                binary,
                                method,
                                pheno_scale = 1) {
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

  } else if(geno_format == "dosage"){
    na_val <- 63

  } else {
    na_val <- 3
  }

  if(method == "glm"){
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
                           null_formula = null_formula,
                           fixed_effect = fixed_effect,
                           dokruskal = dokruskal)
    # Fill NA values
    p_values <- .fill.na(p_values = p_values)

    # Combine results into a single data frame
    p_values <- do.call("rbind", p_values)

  } else {
    sample_id <- getSamID(object)

    if(geno_format == "genotype"){
      gen <- getGenotype(object, "raw")

    } else {
      gen <- getGenotype(object, "cor")
    }
    snp_id <- getMarID(object)
    chr <- as.numeric(factor(getChromosome(object)))
    pos <- getPosition(object)
    allele <- getAllele(object)

    sample_id <- sample_id[!is.na(i_pheno)]
    gen <- gen[!is.na(i_pheno), ]
    i_pheno <- i_pheno[!is.na(i_pheno)]
    fam <- data.frame(famid = sample_id, id = sample_id, father = 0, mother = 0,
                      sex = 0, pheno = i_pheno)
    bim <- data.frame(chr = chr, id = snp_id, dist = 0,
                      pos = pos, A1 = sub(",.+", "", allele),
                      A2 = sub(".+,", "", allele))
    bed_mat <- as.bed.matrix(x = gen, fam = fam, bim = bim)
    k <- GRM(x = bed_mat)
    if(is.null(fixed_effect)){
      fixed_effect <- matrix(1, nrow(bed_mat))
    }
    p_values <- association.test(x = bed_mat, Y = i_pheno, X = fixed_effect,
                                 K = k, p = 2,
                                 eigenK = eigen(k), test = "wald",
                                 method = "lmm", response = "quantitative")
    p_values <- cbind(P.model = p_values$p,
                      P.add = p_values$p,
                      Coef.add = p_values$beta,
                      PVE = p_values$h2)
  }


  # Calculate FDR and -log10(p-values)
  p_values <- cbind(p_values,
                    FDR = p.adjust(p = p_values[, "P.model"], method = "fdr"),
                    negLog10P = -log10(p_values[, "P.model"]))

  p_values <- .rescale_scan_coefs(mat = p_values, pheno_scale = pheno_scale)

  .store_write_scan(
    object = object,
    pheno_name = i_pheno_names,
    mat = p_values,
    colnames_stats = colnames(p_values),
    update_meta = FALSE
  )
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

## Sub-function to store additional scan metadata
.store_additional_info <- function(object,
                                   kruskal,
                                   formula,
                                   null_formula,
                                   fixed_effect,
                                   conv_fun,
                                   geno_format) {
  .store_write_scan_meta(
    object = object,
    kruskal = kruskal,
    formula = formula,
    null_formula = null_formula,
    fixed_effect = fixed_effect,
    conv_fun = conv_fun,
    geno_format = geno_format
  )
}

## Perform regression analysis on genotype data
.regression <- function(g,
                        pheno,
                        binary,
                        na_val,
                        geno_format,
                        conv_fun,
                        formula,
                        null_formula,
                        fixed_effect,
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

    if(!is.null(fixed_effect)){
      df$df <- cbind(df$df, fixed_effect)
    }

    if(!is.null(null_formula)){
      null_df <- .makeDF(g = g,
                         phe = pheno,
                         conv_fun = conv_fun,
                         formula = null_formula)  # Create a data frame for GLM
      if(!is.null(fixed_effect)){
        null_df$df <- cbind(null_df$df, fixed_effect)
      }
    } else {
      null_df <- NULL
    }

    # Set the family for GLM based on whether the phenotype is binary or continuous
    if(binary){
      family <- "binomial"

    } else {
      family <- "gaussian"
    }
    out <- .doGLM(df = df, null_df = null_df, family = family)  # Perform GLM
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
    # Calculate the proportion of variance explained by the model
    tss <- res$null.deviance
    rss <- sum(res$residuals^2)
    pervar <- 1 - rss / tss

  } else {
    # If a null model is provided, fit the null model and calculate p-value
    null_res <- try(glm(formula = null_df$fml,
                        data = null_df$df,
                        family = family))
    p <- pchisq(null_res$deviance - res$deviance,
                null_res$df.residual - res$df.residual,
                lower.tail = FALSE)
    # Calculate the proportion of variance explained by the model
    tss <- res$null.deviance
    rss1 <- sum(res$residuals^2)
    rss2 <- sum(null_res$residuals^2)
    pervar <- (rss2 - rss1) / tss
  }

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
                                     signif = 0.05,
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
                   signif = 0.05){
            # Check if scan data exists in the GDS object
            .check_scan_data_exists(object = object)

            # Determine the phenotype index based on user input
            pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)

            # Retrieve the scan data for the specified phenotype
            pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))

            # Generate the Manhattan plot
            p <- .mhplot_draw(x = pvalues, chr = chr, start = start, end = end, signif = signif)

            return(p)  # Return the plot
          }
)

## Function to check if scan data exists
.check_scan_data_exists <- function(object) {
  if (!.store_check_scan_exists(object)) {
    stop("No scan data in the input LazyGas object.\n",
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
  x_signif <- subset(x, subset = FDR <= signif)
  x <- subset(x, subset = FDR > signif)

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

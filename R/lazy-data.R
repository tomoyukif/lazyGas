#' Extract lazyGas result data from the companion store
#'
#' @param object A \code{LazyGas} object.
#' @param dataset The dataset to retrieve. One of `"scan"`, `"peakcall"`,
#'   `"recalc"`, `"groups"`, `"candidate"`, `"snpeff"`, `"qc"`, `"multitrait"`,
#'   `"conditional"`, `"credible_set"`, or `"pipeline"`.
#' @param pheno The phenotype to retrieve.
#' @param kind Table kind for section datasets (e.g. `"summary"`, `"peak_1"`).
#' @param ... Additional arguments.
#'
#' @return The requested data from the LazyGas object.
#' @export
#'
setGeneric("lazyData", function(object,
                                dataset = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff",
                                            "qc", "multitrait", "conditional", "credible_set", "pipeline"),
                                pheno,
                                kind = NULL,
                                ...)
  standardGeneric("lazyData"))

#'
#' @rdname lazyData
#' @method lazyData LazyGas
#' @export
#'
setMethod("lazyData",
          "LazyGas",
          function(object, dataset, pheno, kind = NULL){
            dataset <- match.arg(
              arg = dataset,
              choices = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff",
                          "qc", "multitrait", "conditional", "credible_set", "pipeline")
            )
            if (dataset == "qc") {
              k <- if (is.null(kind)) "summary" else kind
              return(.store_read_section_df(object, "qc", k, pheno_name = NULL))
            }
            if (dataset == "multitrait") {
              k <- if (is.null(kind)) "clusters" else kind
              return(.store_read_section_df(object, "multitrait", k, pheno_name = NULL))
            }
            if (dataset == "pipeline") {
              return(.store_read_section_df(object, "pipeline", "history", pheno_name = NULL))
            }
            if (dataset %in% c("conditional", "credible_set")) {
              pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
              if (is.null(kind)) {
                stop("kind is required for dataset '", dataset, "' (e.g. 'peak_1').",
                     call. = FALSE)
              }
              return(.store_read_section_df(object, dataset, kind, pheno_name))
            }
            # Determine phenotype index based on input
            pheno_name <- .determine_phenotype_name(object = object,
                                                    pheno = pheno)

            if (dataset == "scan") {
              if (!.store_dataset_exists(object, "scan", pheno_name)) {
                return(NULL)
              }
              out <- .get_scan(object = object, pheno_name = pheno_name)

            } else if (dataset %in% c("candidate", "snpeff", "groups")) {
              if (!.store_dataset_exists(object, dataset, pheno_name)) {
                return(NULL)
              }
              out <- .store_read_matrix_dataset(object, dataset, pheno_name)
              if (is.null(out) || nrow(out) == 0L) {
                return(NULL)
              }

            } else {
              if (!.store_dataset_exists(object, dataset, pheno_name)) {
                return(NULL)
              }
              out <- .get_peakcall(
                object = object,
                pheno_name = pheno_name,
                recalc = (dataset == "recalc")
              )
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

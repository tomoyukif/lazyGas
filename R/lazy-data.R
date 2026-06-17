#' @param dataset The dataset to retrieve. One of "scan", "peakcall", or "recalc".
#' @param pheno The phenotype to retrieve.
#' @param ... Additional arguments.
#'
#' @return The requested data from the LazyGas object.
#' @export
#'
setGeneric("lazyData", function(object,
                                dataset = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff"),
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
                                 choices = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff"))
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

#' Extract lazyGas result data from the companion store
#'
#' @param object A \code{LazyGas} object.
#' @param dataset The dataset to retrieve. One of `"scan"`, `"peakcall"`,
#'   `"recalc"`, `"groups"`, `"candidate"`, `"snpeff"`, `"qc"`, `"multitrait"`,
#'   `"conditional"`, `"credible_set"`, `"fine_mapping"`, `"pipeline"`,
#'   `"phenotype_rank"`, or `"phenotype_query"`.
#' @param pheno The phenotype to retrieve.
#' @param kind Table kind for section datasets (e.g. `"summary"`, `"peak_1"`,
#'   or a phenotype-rank query id such as `"rank_pq_..."`).
#' @param ... Additional arguments.
#'
#' @return The requested data from the LazyGas object.
#' @export
#'
setGeneric("lazyData", function(object,
                                dataset = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff",
                                            "qc", "multitrait", "conditional", "credible_set", "fine_mapping",
                                            "pipeline", "phenotype_rank", "phenotype_query"),
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
                          "qc", "multitrait", "conditional", "credible_set", "fine_mapping",
                          "pipeline", "phenotype_rank", "phenotype_query")
            )
            if (dataset == "phenotype_query") {
              if (is.null(kind)) {
                return(.store_list_phenotype_queries(object))
              }
              q <- .store_read_phenotype_query(object, kind)
              return(q)
            }
            if (dataset == "phenotype_rank") {
              pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
              if (is.null(kind)) {
                meta <- .store_read_meta(object)
                latest <- meta$phenotype_rank_latest
                if (is.null(latest) || is.null(latest[[pheno_name]])) {
                  return(NULL)
                }
                kind <- paste0("rank_", .store_safe_name(latest[[pheno_name]]))
              } else if (!grepl("^rank_", kind)) {
                kind <- paste0("rank_", .store_safe_name(kind))
              }
              return(.store_read_section_df(object, "phenotype_rank", kind, pheno_name))
            }
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
            if (dataset %in% c("conditional", "credible_set", "fine_mapping")) {
              pheno_name <- .determine_phenotype_name(object = object, pheno = pheno)
              if (is.null(kind)) {
                default_kind <- if (dataset == "fine_mapping") "summary" else NULL
                if (is.null(default_kind)) {
                  stop("kind is required for dataset '", dataset, "' (e.g. 'peak_1').",
                       call. = FALSE)
                }
                kind <- default_kind
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

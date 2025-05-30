% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_functions.R
\name{lazyData}
\alias{lazyData}
\alias{lazyData,LazyGas-method}
\title{Retrieve LazyGas Data}
\usage{
lazyData(
  object,
  dataset = c("scan", "peakcall", "recalc", "groups", "candidate", "snpeff"),
  pheno,
  ...
)

\S4method{lazyData}{LazyGas}(object, dataset, pheno)
}
\arguments{
\item{object}{A LazyGas object.}

\item{dataset}{The dataset to retrieve. One of "scan", "peakcall", or "recalc".}

\item{pheno}{The phenotype to retrieve.}

\item{...}{Additional arguments.}
}
\value{
The requested data from the LazyGas object.
}
\description{
This function retrieves data from the LazyGas object based on the specified dataset and phenotype.
}

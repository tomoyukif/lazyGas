% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_functions.R
\name{recalcAssoc}
\alias{recalcAssoc}
\alias{recalcAssoc,LazyGas-method}
\title{Recalculate Associations}
\usage{
recalcAssoc(object, n_threads = NULL, ...)

\S4method{recalcAssoc}{LazyGas}(object, n_threads)
}
\arguments{
\item{object}{A LazyGas object.}
}
\value{
None. The function modifies the LazyGas object in place.
}
\description{
Recalculate associations using the LazyGas object and specified r-square threshold.
}

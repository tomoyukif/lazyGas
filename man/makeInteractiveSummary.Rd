% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_functions.R
\name{makeInteractiveSummary}
\alias{makeInteractiveSummary}
\title{Generate interactive summary}
\usage{
makeInteractiveSummary(
  object,
  pheno,
  what = c("scan_png", "peakcall", "recalc", "groups", "recalc_haplo", "candidate"),
  out_fn
)
}
\arguments{
\item{object}{QTLscan object}

\item{out_fn}{Prefix of output file}

\item{pehno}{Phenotype names to be drawn}
}
\description{
The \code{what} argument accept a character vector of one or more of the following
strings:
\itemize{
\item{"scan"}{Draw an interactive Manhattan plot}
\item{"scan_png"}{Draw an static Manhattan plot}
\item{"peakcall"}{Draw an interactive peackcall plot}
\item{"recalc"}{Draw an interactive recalculated peackll plot}
\item{"groups"}{Draw an interactive peack groping list}
\item{"preakcall_haplo"}{Draw an interactive haplotype-wise phenotype distribution plot for peaks before the recalculation}
\item{"recalc_haplo"}{Draw an interactive haplotype-wise phenotype distribution plot for peaks after the recalculation}
\item{"candidate"}{Draw an interactive candidate list}
}
}

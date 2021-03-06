% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{training_human}
\alias{training_human}
\title{Training Data from the Top 101 Cyclic Genes in 888
  Single-cell Samples}
\format{
A list with the following elements:

\describe{
\item{predict.yy}{Estimated cyclic expression values in the
  training data}

\item{cellcycle_peco_ordered}{Training labels ordered from 0 to
  \code{2*pi}.}

\item{cell_cycle_function}{Nonparametric function of cyclic gene
  expression trend obtained by the \code{trendfilter} function in the
  genlasso package.}

\item{pve}{Proportion of variance explained in each gene by the
  cell cycle phase label.}}
}
\description{
Pre-computed peco results obtained by applying
  \code{fit_cyclic_many} to 888 single-cell samples that have both
  normalized gene expression values and cell cycle labels. These
  training results can be used to predicting cell cycle phase in
  other data.
}
\examples{
data(training_human)

}
\keyword{data}

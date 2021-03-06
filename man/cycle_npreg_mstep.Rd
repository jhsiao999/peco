% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cycle_npreg.R
\name{cycle_npreg_mstep}
\alias{cycle_npreg_mstep}
\title{Estimate Parameters of the Cyclic Trends}
\usage{
cycle_npreg_mstep(
  Y,
  theta,
  method.trend = c("trendfilter", "loess", "bspline"),
  polyorder = 2,
  ncores = 2
)
}
\arguments{
\item{Y}{Gene by sample expression matrix (log2CPM).}

\item{theta}{Observed cell times.}

\item{method.trend}{How to estimate cyclic trend of gene expression
values. We offer three options: \code{method.trend =
"trendfilter"}, uses \code{\link{fit_trendfilter}},
\code{method.trend = "loess"} uses \code{\link{fit_loess}}), and
\code{method.trend = "bsplines"} uses \code{\link{fit_bspline}}.  We
found that \code{"trendfilter"} provided the best fits in our
experiments. However, trend-filtering may require more
computational effort since it uses cross-validation, so for fast
results we recommend using \code{"bspline"}.}

\item{polyorder}{We estimate cyclic trends of gene expression
levels using nonparamtric trend filtering.}

\item{ncores}{How many computing cores to use? We use the
\code{doParallel} package for parallel computing.}
}
\value{
A list with the following elements:

\item{Y}{Input gene expression data.}

\item{theta}{Input angles.}

\item{mu_est}{Estimated expression levels given the cyclic function
      for each gene.}

\item{sigma_est}{Estimated standard error of the cyclic trends for
  each gene.}

\item{funs}{Estimated cyclic functions.}
}
\description{
This is used by \code{\link{cycle_npreg_insample}}
  (model fitting from training data) and
  \code{\link{cycle_npreg_outsample}} (prediction in test data
  estimate cyclic trends of gene expression values. The function
  outputs, for each gene, a standard error of the cyclic trend, a
  cyclic function, and the estimated expression levels from the
  cyclic function.
}
\keyword{internal}

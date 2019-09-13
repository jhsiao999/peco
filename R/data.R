#' Molecule counts of the 101 significant cyclical genes in the 888 samples
#' analyzed in the study.
#'
#' An ExpressionSet object (require Biobase package) including
#' molecule count data after gene and smaple filtering. The `phenotypeData()`
#' slot contains sample phenotype information and the `featureData()` slot
#' contains gene feature information.
#'
#' @format An ExpressionSet object with 888 samples and the 101 significant
#'     cyclic genes,
#' \describe{
#'   \item{theta}{Inferred angles of each cell along
#'     a circle, also known as FUCCI phase.}
#' }
#'
#' @docType data
#'
#' @usage data(eset_sub)
#'
#' @keywords data
"eset_sub"

#' @title Traing model results among samples from 5 individuals.
#'
#' @description Pre-computed results. Applied \emph{cycle_npreg_insample} to
#'   obtain gene-specific cyclic trend parameters using samples from 5
#'   individuals
#'
#' @format A list with the follwing elements
#' \describe{
#'   \item{Y}{a data.frame (gene by sample) of quantile-normailzed gene
#'   expression values}
#'   \item{theta}{a vector of cell cycl phase values (range between 0 to 2pi)}
#'   \item{sigma_est}{a vector of estimated standard errors}
#'   \item{funs_est}{a list of estimated cyclic functions}
#' }
#'
#' @docType data
#'
#' @usage data(fit_train)
#'
#' @keywords data
"fit_train"

#' @title Results of predicting cell cycle phase for samples from NA18511.
#'
#' @description Pre-computed results. Applied \emph{cycle_npreg_outsample} and
#'   results stored in \emph{fit_train} to predict cell cycle phase for
#'   single-cell samples of NA19098.
#'
#' @format A list with the follwing elements
#' \describe{
#'   \item{Y}{a data.frame (gene by sample) of quantile-normailzed gene
#'   expression values}
#'   \item{cell_times_est}{a vector of estimated cell cycl phase values
#'   (range between 0 to 2pi)}
#'   \item{loglik_est}{a vector of estimated log-likelihoods for each cell}
#'   \item{Y_reordered}{Y ordered along columns by \emph{cell_times-est}}
#'   \item{cell_times_reordered}{cell_times_est ordered from 0 to 2pi}
#'   \item{funs_reordered}{A lists of estimated cyclic functions for each gene
#'         using re-ordered gene expression values (rows of Y_reordered) }
#'   \item{mu_reordered}{A data.frame (gene by sample) of estimated expression
#'   for each gene across cells in the test data using \emph{funs_reordered}}
#'   \item{sigma_reordered}{A vector of estimated standard error for each gene
#'   in the test data across cells using \emph{funs_reordered}}
#'   \item{prob_per_cell_by_celltimes}{A data.frame (gene by sample) of
#'   log-likelihoods of of gene expression values in the test data}
#' }
#'
#' \describe{
#'   \item{Y_reordered}{Quntile-normalized expression values}
#'   \item{cell_times_reordered}{Estimated phase}
#'  }
#'
#' @docType data
#'
#' @usage data(fit_predict)
#'
#' @keywords data
"fit_predict"

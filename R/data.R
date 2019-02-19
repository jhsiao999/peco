#' @name eset_final_sub
#'
#' @docType data
#'
#' @title Molecule counts of top 100 cyclical genes in 888 samples analyzed
#'   in the study.
#'
#' @description An ExpressionSet object (require Biobase package) including
#'   molecule count data after gene and smaple filtering. The `phenotypeData()` slot
#'   contains sample phenotype information and the `featureData()` slot contains
#'   gene feature information.
#'
#' @format An ExpressionSet object with 888 samples and top 100 cyclic genes.
#'
#' \describe{
#'   \item{`pData(eset_final_sub)$theta`}{Inferred angles of each cell along
#'     a circle, also known as FUCCI phase.}
#'   \item{`exprs(est_final_sub)`}{Molecule counts of top 100 cyclical genes.}
#' }
#'
#' @keywords data
#'
NULL

#' @name genes_cyclic_list
#'
#' @docType data
#'
#' @title List of top 100 cyclic genes
#'
#' @description Top 100 cyclic genes ordered by their cyclic trend (strong to weak)
#'   across 888 samples analyzed in the study.
#'
#' @format A data frame with the following columns, in which the rows
#' are ordered according to the proportion of variance explained
#' (large to small):
#'
#' \describe{
#'
#'   \item{ensg}{ENSG gene ID.}
#'
#'   \item{pve}{Proportion of variance explained in the expression
#'   values by the estimated cyclic trend.}
#' }
#'
#' @keywords data
#'
NULL

#' @name fit_train
#'
#' @docType data
#'
#' @title Traing model results among samples from 5 individuals.
#'
#' @description Pre-computed results. Applied \code{cycle_npreg_insample} to
#'   obtain gene-specific cyclic trend parameters using samples from 5
#'   individuals
#'
#' @format
#'
#' \describe{
#'   \item{\code{fit_train}}{training results}
#'   }
#'
#' @keywords data
#'
NULL

#' @name fit_pred
#'
#' @docType data
#'
#' @title Results of predicting cell cycle phase for samples from NA19098.
#'
#' @description Pre-computed results. Applied \code{cycle_npreg_outsample} and
#'   results stored in \code{fit_train} to predict cell cycle phase for
#'   single-cell samples of NA19098.
#'
#' @format
#'
#' \describe{
#'   \item{\code{Y_reordered}}{Normalized expression values (log2CPM , )}
#'   \item{\code{cell_times_reordered}}{Estimated phase}
#'  }
#'
#' @keywords data
#'
NULL

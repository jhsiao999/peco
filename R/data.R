#' @name training_human
#'
#' @title Training Data from 888 Single-cell Samples and 101 Top
#'   Cyclic Genes
#'
#' @description Pre-computed results. Applied \code{fit_cyclic_many}
#' to 888 single-cell samples that have both normalized gene
#' expression values and cell cycle labels to obtain training results
#' that can be used as input for predicting cell cycle phase in other
#' data.
#'
#' @format A list with the follwing elements:
#'
#' \describe{
#'   \item{predict.yy}{Estimated cyclic expression values in the
#' training data}
#' 
#'   \item{cellcycle_peco_ordered}{Training labels ordered from 0 to
#'     2*pi}
#' 
#'   \item{cell_cycle function}{Nonparametric function of cyclic gene
#' expression trend obtained by trendfilter function in genlasso}
#' 
#'   \item{pve}{Proportion of variance explained in each gene by the
#'   cell cycl phase label}
#'  }
#'
#' @docType data
#'
#' @usage data(training_human)
#'
#' @keywords data
#' 
NULL

#' @name sce_top101genes
#' 
#' @title Molecule counts of the 101 significant cyclical genes in the
#' 888 samples analyzed in the study.
#'
#' @description A SingleCellExperiment object (require
#' SingleCellExperiment package) including molecule count data after
#' gene and smaple filtering.  The `colData()` slot contains sample
#' phenotype information and the `rowData()` slot contains gene
#' feature information.
#'
#' @format A SingleCellExperiment object with 888 samples and the 101
#'     significant cyclic genes,
#' 
#' \describe{
#'     \item{theta}{Inferred angles of each cell along
#'         a circle, also known as FUCCI phase.}
#' }
#'
#' @docType data
#'
#' @usage data(sce_top101genes)
#'
#' @keywords data
#' 
NULL

#' @name model_5genes_train
#' 
#' @title Traing model results among samples from 5 individuals.
#'
#' @description Pre-computed results. Applied
#'   \code{cycle_npreg_insample} to obtain gene-specific cyclic trend
#'   parameters using samples from 5 individuals
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
#' @usage data(model_5genes_train)
#'
#' @keywords data
#' 
NULL

#' @name model_5genes_predict
#' 
#' @title A SingleCellExperiment object
#'
#' @description Pre-computed results. Applied
#' \emph{cycle_npreg_outsample} and results stored in
#' \emph{model_5genes_train} to predict cell cycle phase for
#' single-cell samples of NA19098. The predicted cell cycle is stored
#' as variable \emph{cellcycle_peco}.
#'
#' @format A list with the follwing elements
#'
#' \describe{
#'   \item{cellcycle_peco}{Predict cell cycle,
#'    the values ranged between 0 to 2pi}
#'  }
#'
#' @docType data
#'
#' @usage data(model_5genes_predict)
#'
#' @keywords data
#'
NULL

#' @name cellcyclegenes_whitfield2002
#' 
#' @title List of Cell Cycle Genes Identified in Whitfield et al 2002.
#'
#' @description List of cell cycle genes and their associated cell
#'   cycle state as reported in Whitfield et al. 2002.
#'
#' @format A list with the follwing elements:
#'
#' \describe{
#'   \item{hgnc}{Gene symbol}
#'   \item{ensembl}{ENSEMBL gene ID}
#'   \item{phase}{Marker phase identified in Whitfield et al. 2002}
#'  }
#'
#' @docType data
#'
#' @usage data(cellcyclegenes_whitfield2002)
#'
#' @keywords data
#'
NULL


#' @name eset_final_sub
#'
#' @docType data
#' 
#' @title Filtered gene expression data
#' 
#' @description Filtered gene expression data, along with summarized FUCCI
#' intensities and sample information. Contains the top 100 cyclic
#' genes from our data.
#'
#'
#' @format An ExpressionSet object with 888 samples and top 100 cyclic
#' genes.
#' 
#' \describe{
#'   \item{theta}{Inferred angles.}
#' 
#'   \item{\code{exprs(est_final)}}{Molecule counts.}
#' }
#'
#' @keywords data
#' 
NULL

#' @name genes_cyclic_list
#'
#' @docType data
#' 
#' @title Genes ordered by their cyclic trend
#'
#' @description Genes ordered by their cyclic trend (strong to weak),
#' and sample information.
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

#' @name fit
#'
#' @docType data
#' 
#' @title Training and prediction results
#' 
#' @description Pre-computed results, including training and
#' prediction and sample information.
#'
#' @format
#' 
#' \describe{
#'   \item{\code{fit_train}}{training results}
#' 
#'   \item{\code{fit_predict}}{prediction results}
#'  }
#'
#' @keywords data
#' 
NULL

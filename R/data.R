#' @name training_human
#'
#' @title Training Data from the Top 101 Cyclic Genes in 888
#'   Single-cell Samples
#'
#' @description Pre-computed peco results obtained by applying
#'   \code{fit_cyclic_many} to 888 single-cell samples that have both
#'   normalized gene expression values and cell cycle labels. These
#'   training results can be used to predicting cell cycle phase in
#'   other data.
#'
#' @format A list with the following elements:
#'
#' \describe{
#' \item{predict.yy}{Estimated cyclic expression values in the
#'   training data}
#' 
#' \item{cellcycle_peco_ordered}{Training labels ordered from 0 to
#'   \code{2*pi}.}
#' 
#' \item{cell_cycle_function}{Nonparametric function of cyclic gene
#'   expression trend obtained by the \code{trendfilter} function in the
#'   genlasso package.}
#' 
#' \item{pve}{Proportion of variance explained in each gene by the
#'   cell cycle phase label.}}
#'
#' @docType data
#'
#' @examples
#' data(training_human)
#'
#' @keywords data
#' 
NULL

#' @name sce_top101genes
#' 
#' @title Molecule Counts of the 101 Significant Cyclical Genes in the
#'   888 samples Analyzed in the Study.
#'
#' @description A \code{SingleCellExperiment} object containing
#'   processed molecule count data (obtained after gene and sample
#'   filtering). The \code{colData} slot contains sample phenotype
#'   information.
#'
#' @format A \code{SingleCellExperiment} object with 888 samples and
#'   the 101 significant cyclic genes. Also, {sce_top101genes$theta}
#'   contains the inferred angles (FUCCI phase) for each cell.
#'
#' @docType data
#'
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#' dim(sce_top101genes@colData)
#'
#' @keywords data
#' 
NULL

#' @name model_5genes_train
#' 
#' @title Model Training Results for 5 Genes
#'
#' @description Pre-computed gene-specific cyclic trend results
#'   obtained by applying \code{cycle_npreg_insample}.
#'
#' @format A list with the follwing elements:
#' \describe{
#' 
#' \item{Y}{A matrix (gene by sample) of quantile-normailzed gene
#'   expression values.}
#' 
#' \item{theta}{A vector of cell cycle phase values (ranging between
#'   0 and \code{2*pi}).}
#' 
#' \item{sigma_est}{A vector of estimated standard errors, one for
#'   each gene.}
#' 
#' \item{funs_est}{A list of estimated cyclic functions.}}
#'
#' @docType data
#'
#' @examples
#' data(model_5genes_train)
#'
#' @keywords data
#' 
NULL

#' @name model_5genes_predict
#' 
#' @title Cell-Cycle Predictions using 5 Genes
#'
#' @description Pre-computed peco results obtained by applying
#'   \emph{cycle_npreg_outsample} with model stored in
#'   \emph{model_5genes_train} to predict cell cycle phase for
#'   single-cell samples of NA19098. 
#'
#' @format A \code{SingleCellExperiment} object. The predicted cell
#' cycle is stored in \emph{model_5genes_predict$cellcycle_peco}.
#'
#' @docType data
#'
#' @examples
#' data(model_5genes_predict)
#'
#' @keywords data
#'
NULL

#' @name cellcyclegenes_whitfield2002
#' 
#' @title List of Cell Cycle Genes Identified in Whitfield et al, 2002.
#'
#' @description List of cell cycle genes and their associated cell
#'   cycle state from Whitfield \emph{et al} (2002).
#'
#' @format A data frame with the three columns:
#'
#' \describe{
#'   \item{hgnc}{HGNC gene symbol.}
#'   \item{ensembl}{Ensembl gene identifier.}
#'   \item{phase}{Marker phase identified in Whitfield \emph{et al} (2002).}
#'  }
#'
#' @docType data
#'
#' @examples
#' data(cellcyclegenes_whitfield2002)
#'
#' @keywords data
#'
NULL


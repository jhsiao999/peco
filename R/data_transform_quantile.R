#' @name data_transform_quantile
#'
#' @title Transform counts by first computing counts-per-million (CPM), then
#'   quantile-normalize CPM for each gene
#'
#' @description
#'   For each gene, transform CPM to a normal distribution.
#' This way the zero-count cells are assigned the lowest qunatiles.
#'
#' @param Y_sce SingleCellExperiment Object.
#' @param ncores We use doParallel package for parallel computing.
#'
#' @return SingleCellExperiment Object with an added slot of cpm_quant,
#'     cpm slot is added if it doesn't exist.
#'
#' @examples
#' # use our data
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # perform CPM normalization using scater, and
#' # quantile-normalize the CPM values of each gene to normal distribution
#' sce_normed <- data_transform_quantile(sce_top101genes, ncores=2)
#'
#' plot(y=assay(sce_normed, "cpm_quantNormed")[1,],
#'      x=assay(sce_normed, "cpm")[1,],
#'     xlab = "CPM bbefore quantile-normalization",
#'     ylab = "CPM after quantile-normalization")
#'
#' @author Joyce Hsiao
#'
#' @import SingleCellExperiment
#' @import methods
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import scater
#' @export
data_transform_quantile <- function(sce, ncores=2) {

  if (is.null(ncores)) {
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  }

  cpm(sce) <- calculateCPM(sce)

  G <- nrow(sce)
  cpm_sce <- cpm(sce)

  cpm_quantNormed <- foreach::foreach(g=seq_len(G)) %dopar% {
    y_g <- cpm_sce[g,]
    is.zero <- which(y_g == 0)
    qq.map <- stats::qqnorm(y_g, plot.it=FALSE)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(y_g= yy.qq)
  }
  parallel::stopCluster(cl)
  cpm_quantNormed <- do.call(rbind, cpm_quantNormed)
  colnames(cpm_quantNormed) <- colnames(cpm_sce)
  rownames(cpm_quantNormed) <- rownames(cpm_sce)

  assays(sce)[[3]] <- cpm_quantNormed
  assayNames(sce)[3] <- "cpm_quantNormed"

  return(sce)
  }


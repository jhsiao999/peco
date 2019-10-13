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
#' # normalize expression counts to counts per million
#' sce_normed <- data_transform_quantile(sce_top101genes, ncores=2)
#'
#' plot(y=assay(sce_normed, "cpm_quant")[1,],
#'      x=assay(sce_normed, "umi_counts")[1,],
#'     xlab = "Before quantile-normalization",
#'     ylab = "After quantile-normalization")
#'
#' @author Joyce Hsiao
#'
#' @import SingleCellExperiment
#' @import methods
#' @import doParallel
#' @import parallel
#' @import foreach
#' @export
data_transform_quantile <- function(Y_sce, ncores=2) {

  if (is.null(ncores)) {
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  }

  cpm <- t((10^6)*t(assay(Y_sce, "umi_counts"))/colData(Y_sce)$molecules)

  G <- nrow(cpm)

  df <- foreach::foreach(g=seq_len(G)) %dopar% {
    y_g <- cpm[g,]
    is.zero <- which(y_g == 0)
    qq.map <- stats::qqnorm(y_g, plot.it=FALSE)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(y_g= yy.qq)
  }
  parallel::stopCluster(cl)
  df <- do.call(rbind, df)
  colnames(df) <- colnames(cpm)
  rownames(df) <- rownames(cpm)

  assays(Y_sce) = list(umi_counts = assay(Y_sce, "umi_counts"),
                       cpm = cpm,
                       cpm_quant = df)

  return(Y_sce)
  }


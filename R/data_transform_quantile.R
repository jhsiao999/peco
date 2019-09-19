#' @name data_transform_quantile
#'
#' @title Quantile-normalize log counts of gene expression values
#'
#' @description
#'   For each gene, transform the log counts to a normal distribution.
#' This way the zero-count cells are assigned the lowest qunatiles.
#'
#' @param Y A gene by sample matrix. Contains library size-normalized
#' molecule counts.
#' @param ncores We use doParallel package for parallel computing.
#'
#' @return A gene by sample expression matrix.
#'
#' @examples
#' # use our data
#' data(sce_sub)
#'
#' # normalize expression counts to counts per million
#' counts_normed<-t((10^6)*t(assay(sce_sub)[1:5,])/colData(sce_sub)$molecules)
#' counts_quant <- data_transform_quantile(counts_normed, ncores=2)
#'
#' plot(x=counts_normed[1,], y=counts_quant[1,],
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
data_transform_quantile <- function(Y, ncores=2) {

  if (is.null(ncores)) {
    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))
  }

  G <- nrow(Y)

  # df <- mclapply(seq_len(G), function(g) {
  #   y_g <- Y[g,]
  #   is.zero <- which(y_g == 0)
  #   qq.map <- stats::qqnorm(y_g, plot.it=FALSE)
  #   yy.qq <- qq.map$x
  #   yy.qq[is.zero] <- sample(qq.map$x[is.zero])
  #   return(y_g= yy.qq)
  # }, mc.cores = ncores)
  df <- foreach::foreach(g=seq_len(G)) %dopar% {
    y_g <- Y[g,]
    is.zero <- which(y_g == 0)
    qq.map <- stats::qqnorm(y_g, plot.it=FALSE)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(y_g= yy.qq)
  }
  parallel::stopCluster(cl)

  df <- do.call(rbind, df)
  colnames(df) <- colnames(Y)
  rownames(df) <- rownames(Y)

  return(df)
  }


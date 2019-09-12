#' @title Quantile-normalize log counts of gene expression values
#'
#' @description
#'   For each gene, transform the log counts to a normal distribution.
#' This way the zero-count cells are assigned the lowest qunatiles.
#'
#' @param Y A gene by sample matrix. Contains library size-normalized
#' molecule counts.
#' @param ncores We use mclapply function for parallel computing.
#'
#' @return A gene by sample expression matrix.
#'
#' @examples
#' # use our data
#' data(eset_sub)
#'
#' # normalize expression counts to counts per million
#' counts_normed<-t((10^6)*t(exprs(eset_sub)[1:5,])/pData(eset_sub)$molecules)
#' counts_quant <- data_transform_quantile(counts_normed, ncores=2)
#'
#' plot(x=counts_normed[1,], y=counts_quant[1,],
#'     xlab = "Before quantile-normalization",
#'     ylab = "After quantile-normalization")
#'
#' @author Joyce Hsiao
#' @export
#'
#' @importFrom parallel mclapply
#' @import methods Biobase MASS Matrix ggplot2
NULL
data_transform_quantile <- function(Y, ncores=2) {
  G <- nrow(Y)

  df <- mclapply(seq_len(G), function(g) {
    y_g <- Y[g,]
    is.zero <- which(y_g == 0)
    qq.map <- stats::qqnorm(y_g)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(y_g= yy.qq)
  }, mc.cores = ncores)
  df <- do.call(rbind, df)
  colnames(df) <- colnames(Y)

  return(df)
  }


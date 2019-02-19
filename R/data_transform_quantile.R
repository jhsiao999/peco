#' @title Quantile-normalize log counts of gene expression values
#'
#' @description For each gene, transform the log counts to a normal distribution.
#'   This way the zero-count cells are assigned the lowest qunatiles.
#'
#' @param Y A gene by sample matrix. Contains library size-normalized molecule counts.
#' @param ncores We use mclapply function for parallel computing.
#'
#' @return A gene by sample expression matrix.
#'
#' @author Joyce Hsiao
#'
#' @examples

data_transform_quantile <- function(Y, ncores=2) {
  G <- nrow(Y)

  df <- mclapply(1:G, function(g) {
    y_g <- Y[g,]
    is.zero <- which(y_g == 0)
    qq.map <- qqnorm(y_g)
    yy.qq <- qq.map$x
    yy.qq[is.zero] <- sample(qq.map$x[is.zero])
    return(y_g= yy.qq)
  }, mc.cores = ncores)
  df <- do.call(rbind, df)
  colnames(df) <- colnames(Y)

  return(df)
  }


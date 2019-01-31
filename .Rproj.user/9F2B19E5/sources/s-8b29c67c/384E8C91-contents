#' @title Get cyclical genes
#'
#' @export
# Y <- data_training$log2cpm.quant.nonvalid[,fold_indices[[1]]$train]
# theta <- data_training$theta.nonvalid[fold_indices[[1]]$train]
# ncores=15
# Y <- Y[1:10,]
get.cyclical <- function(Y, theta, polyorder=2, ncores) {
  library(genlasso)
  library(assertthat)
  library(parallel)

  G <- nrow(Y)
  N <- ncol(Y)

  if (!assert_that(all.equal(names(theta), colnames(Y)))) {
    Y_ordered <- Y[,match(names(theta), colnames(Y))]
    ord <- order(theta)
    theta_ordered <- theta[ord]
    Y_ordered <- Y_ordered[,ord]
  } else {
    ord <- order(theta)
    theta_ordered <- theta[ord]
    Y_ordered <- Y[,ord]
  }

  # for each gene, estimate the cyclical pattern of gene expression
  # conditioned on the given cell times
  fit <- mclapply(1:G, function(g) {
    #    print(g)
    y_g <- Y_ordered[g,]

    fit_g <- fit.trendfilter.generic(yy=y_g, polyorder = polyorder)

    return(fit_g$pve)
  }, mc.cores = ncores)

  out <- do.call(rbind, fit)
  out <- as.data.frame(out)
  colnames(out) <- "pve"
  rownames(out) <- rownames(Y_ordered)

  return(out)
}



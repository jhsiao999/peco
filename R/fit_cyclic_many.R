#' @title Estimate Cyclic Trend For Several Genes Using Trendfiltering
#'
#' @description For each gene, apply trend filtering as implemented in
#'   the genlasso package to estimate cell cycle. For more details, see
#'   \code{link{fit_trendfiltering}}.
#'
#' @param Y A matrix (gene by sample) of gene expression values. The
#'   expression values are assumed to have been normalized and
#'   transformed to the standard normal distribution.
#'
#' @param theta A vector of cell cycle phase (angles) for single-cell
#'   samples.
#'
#' @param polyorder Argument passed to \code{\link{fit_trendfilter}}
#'   specifying the degree of the polynomials used in nonparametric
#'   trend filtering.
#'
#' @param ncores Argument passed to
#'   \code{\link[parallel]{makeCluster}} specifying the number of
#'   threads.
#'
#' @return A list containing the following elements:
#'
#' \item{predict.yy}{A matrix of predicted expression values at observed
#'   cell cycle.}
#' 
#' \item{cellcycle_peco_ordered}{A vector of predicted cell
#'   cycle. Values range between 0 to 2pi}
#' 
#' \item{cellcycle_function}{List of predicted cell cycle functions.}
#' 
#' \item{pve}{Vector of proportion of variance explained in each gene by
#'   the predicted cell cycle.}
#'
#' @seealso \code{\link{fit_trendfilter}}
#' 
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # Select top 10 cyclic genes.
#' sce_top10 <- sce_top101genes[order(rowData(sce_top101genes)$pve_fucci,
#'                                   decreasing=TRUE)[1:10],]
#' coldata <- colData(sce_top10)
#'
#' # Get cell cycle phase based on FUCCI scores.
#' theta <- coldata$theta
#' names(theta) <- rownames(coldata)
#'
#' # Normalize expression counts.
#' sce_top10 <- data_transform_quantile(sce_top10, ncores=2)
#' exprs_quant <- assay(sce_top10, "cpm_quantNormed")
#'
#' # Order FUCCI phase and expression.
#' theta_ordered <- theta[order(theta)]
#' yy_ordered <- exprs_quant[, names(theta_ordered)]
#'
#' fit <- fit_cyclical_many(Y=yy_ordered, theta=theta_ordered)
#'
#' @author Joyce Hsiao
#' 
#' @seealso \code{\link{fit_trendfilter}} for fitting one gene
#'   using trend filtering.
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom stats predict
#'
#' @export
#' 
fit_cyclical_many <- function(Y, theta, polyorder=2, ncores=2) {
    if (is.null(ncores)) {
        cl <- makeCluster(2)
        registerDoParallel(cl)
        message(paste("computing with",ncores,"threads"))
    } else {
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        message(paste("computing with",ncores,"threads"))
    }

    G <- nrow(Y)
    N <- ncol(Y)

    Y_ordered <- Y[,names(theta)]
    ord <- order(theta)
    theta_ordered <- theta[ord]
    Y_ordered <- Y_ordered[,ord]

    # For each gene, estimate the cyclical pattern of gene expression
    # conditioned on the given cell times.
    fit <- foreach(g=seq_len(G),
                   .export = c("trendfilter", "fit_trendfilter"),
                   .packages = c("genlasso", "peco")) %dopar% {
        y_g <- Y_ordered[g,]
        fit_g <- fit_trendfilter(yy=y_g, polyorder=polyorder)
        fun_g <- approxfun(x=as.numeric(theta_ordered),
                           y=as.numeric(fit_g$trend.yy), rule=2)
        mu_g <- fit_g$trend.yy
        sigma_g <- sqrt(sum((y_g - mu_g)^2)/N)
        return(list(trend.yy=fit_g$trend.yy,
                    sigma_g=sigma_g,
                    fun_g=fun_g,
                    pve=fit_g$pve))
    }
    stopCluster(cl)

    predict.yy <- do.call(rbind,(lapply(fit, "[[", 1)))
    colnames(predict.yy) <- colnames(Y_ordered)
    rownames(predict.yy) <- rownames(Y_ordered)

    sigma <- vapply(fit, function(x) x[[2]], double(1))
    sigma <- as.data.frame(sigma)
    colnames(sigma) <- "sigma"
    rownames(sigma) <- rownames(Y_ordered)

    pve <- vapply(fit, function(x) x[[4]], double(1))
    pve <- as.data.frame(pve)
    colnames(pve) <- "pve"
    rownames(pve) <- rownames(Y_ordered)

    cellcycle_function <- lapply(fit, "[[", 3)
    names(cellcycle_function) <- rownames(Y_ordered)

    out <- list(predict.yy=predict.yy,
                cellcycle_peco_ordered=theta_ordered,
                cellcycle_function=cellcycle_function,
                sigma=sigma,
                pve=pve)
    return(out)
}

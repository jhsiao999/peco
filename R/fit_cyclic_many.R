#' @name fit_cyclical_many
#'
#' @title Compute proportation of variance explained by cyclic trends
#' in the gene expression levels for each gene.
#'
#' @description We applied quadratic (second order) trend filtering
#' using the trendfilter function in the genlasso package (Tibshirani,
#' 2014). The trendfilter function implements a nonparametric
#' smoothing method which chooses the smoothing parameter by
#' cross-validation and fits a piecewise polynomial regression. In
#' more specifics: The trendfilter method determines the folds in
#' cross-validation in a nonrandom manner. Every k-th data point in
#' the ordered sample is placed in the k-th fold, so the folds contain
#' ordered subsamples. We applied five-fold cross-validation and chose
#' the smoothing penalty using the option lambda.1se: among all
#' possible values of the penalty term, the largest value such that
#' the cross-validation standard error is within one standard error of
#' the minimum. Furthermore, we desired that the estimated expression
#' trend be cyclical. To encourage this, we concatenated the ordered
#' gene expression data three times, with one added after another. The
#' quadratic trend filtering was applied to the concatenated data
#' series of each gene.
#'
#' @param Y A matrix (gene by sample) of gene expression values.  The
#' expression values are assumed to have been normalized and
#' transformed to standard normal distribution.
#'
#' @param theta A vector of cell cycle phase (angles) for single-cell
#' samples.
#'
#' @param polyorder We estimate cyclic trends of gene expression
#' levels using nonparamtric trend filtering. The default fits second
#' degree polynomials.
#'
#' @param ncores doParallel package is used to perform
#' parallel computing to reduce computational time.
#'
#' @return A list containing the following objects
#'
#' \item{predict.yy}{A matrix of predicted expression values at observed
#'   cell cycle.}
#' \item{cellcycle_peco_ordered}{A vector of predicted cell cycle. The values
#'   range between 0 to 2pi}
#' \item{cellcycle_function}{A list of predicted cell cycle functions.}
#' \item{pve}{A vector of proportion of variance explained in each gene by
#'   the predicted cell cycle.}
#'
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # select top 10 cyclic genes
#' sce_top10 <- sce_top101genes[order(rowData(sce_top101genes)$pve_fucci,
#'                                   decreasing=TRUE)[1:10],]
#' coldata <- colData(sce_top10)
#'
#' # cell cycle phase based on FUCCI scores
#' theta <- coldata$theta
#' names(theta) <- rownames(coldata)
#'
#' # normalize expression counts
#' sce_top10 <- data_transform_quantile(sce_top10, ncores=2)
#' exprs_quant <- assay(sce_top10, "cpm_quantNormed")
#'
#' # order FUCCI phase and expression
#' theta_ordered <- theta[order(theta)]
#' yy_ordered <- exprs_quant[, names(theta_ordered)]
#'
#' fit <- fit_cyclical_many(Y=yy_ordered, theta=theta_ordered)
#'
#' @author Joyce Hsiao
#' @seealso
#'     \code{\link{fit_trendfilter}} for fitting one gene
#'     using trendfilter
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom genlasso trendfilter cv.trendfilter
#' @importFrom stats var predict
#' @import foreach
#'
#' @import parallel
#' @export
fit_cyclical_many <- function(Y, theta, polyorder=2, ncores=2) {

    if (is.null(ncores)) {
        cl <- makeCluster(2)
        registerDoParallel(cl)
        message(paste("computing on",ncores,"cores"))
    } else {
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        message(paste("computing on",ncores,"cores"))
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
        fit_g <- fit_trendfilter(yy=y_g, polyorder = polyorder)
        fun_g <- approxfun(x=as.numeric(theta_ordered),
                        y=as.numeric(fit_g$trend.yy), rule=2)
        mu_g <- fit_g$trend.yy
        sigma_g <- sqrt(sum((y_g-mu_g)^2)/N)
        return(list(trend.yy=fit_g$trend.yy,
                    sigma_g = sigma_g,
                    fun_g = fun_g,
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

    out <- list(predict.yy = predict.yy,
                cellcycle_peco_ordered=theta_ordered,
                cellcycle_function = cellcycle_function,
                sigma = sigma,
                pve = pve)

    return(out)
}

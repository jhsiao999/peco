#' @name fit_trendfilter
#'
#' @title Using trendfiltering to estimate cyclic trend of gene
#' expression
#'
#' @description We applied quadratic (second order) trend filtering
#' using the trendfilter function in the genlasso package (Tibshirani,
#' 2014).  The trendfilter function implements a nonparametric
#' smoothing method which chooses the smoothing parameter by
#' cross-validation and fits a piecewise polynomial regression. In
#' more specifics: The trendfilter method determines the folds in
#' cross-validation in a nonrandom manner.  Every k-th data point in
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
#' @param yy A vector of gene expression values for one gene
#' that are ordered by cell cycle phase. Also, the
#' expression values are normalized and
#' transformed to standard normal distribution.
#'
#' @param polyorder We estimate cyclic trends of gene expression
#' levels using nonparamtric trend filtering. The default fits second
#' degree polynomials.
#'
#' @return A list with two elements:
#'     \item{trend.yy}{The estimated cyclic trend.}
#'     \item{pve}{Proportion of variance explained by the cyclic
#'     trend in the gene expression levels.}
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
#' yy_ordered <- exprs_quant[1, names(theta_ordered)]
#'
#' fit <- fit_trendfilter(yy_ordered)
#'
#' plot(x=theta_ordered, y=yy_ordered, pch=16, cex=.7, axes=FALSE,
#'   ylab="quantile-normalized expression values", xlab="FUCCI phase",
#'   main = "trendfilter fit")
#' points(x=theta_ordered, y=fit$trend.yy, col="blue", pch=16, cex=.7)
#' axis(2)
#' axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#'   labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
#'   expression(2*pi)))
#' abline(h=0, lty=1, col="black", lwd=.7)
#'
#' @author Joyce Hsiao
#' @importFrom genlasso trendfilter cv.trendfilter
#' @importFrom stats var predict
#' @export
fit_trendfilter <- function(yy, polyorder=2) {

    yy.rep <- rep(yy,3)
    include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

    fit.trend <- trendfilter(yy.rep,
                            ord=polyorder, approx=FALSE, maxsteps = 1000)
    cv.trend <- cv.trendfilter(fit.trend)
    which.lambda <- cv.trend$i.1se
    yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                            df=fit.trend$df[which.lambda])$fit
    trend.yy <- yy.trend.pred
    pve <- 1-var(yy-trend.yy)/var(yy)

    return(list(trend.yy=trend.yy[include],
                pve=pve))
}

#' @name fit_bspline
#'
#' @title Use bsplies to cyclic trend of gene expression levels
#'
#' @param yy A vector of gene expression values for one gene. The
#' expression values are assumed to have been normalized and
#' transformed to standard normal distribution.
#'
#' @param time A vector of angels (cell cycle phase).
#'
#' @return A list with one element, \code{pred.yy}, giving the
#' estimated cyclic trend.
#'
#' @author Joyce Hsiao
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
#' yy_ordered <- exprs_quant[1, names(theta_ordered)]
#'
#' fit <- fit_bspline(yy_ordered, time=theta_ordered)
#'
#' plot(x=theta_ordered, y=yy_ordered, pch=16, cex=.7, axes=FALSE,
#'   ylab="quantile-normalized expression values", xlab="FUCCI phase",
#'   main = "bspline fit")
#' points(x=theta_ordered, y=fit$pred.yy, col="blue", pch=16, cex=.7)
#' axis(2)
#' axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#'   labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
#'   expression(2*pi)))
#' abline(h=0, lty=1, col="black", lwd=.7)
#'
#' @importFrom stats smooth.spline var predict
#' @export
fit_bspline <- function(yy, time) {

    yy.rep <- rep(yy,3)
    time.rep <- c(time, time+(2*pi), time+(4*pi))
    include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

    fit <- smooth.spline(x=time.rep, y=yy.rep)
    pred.yy <- predict(fit, time.rep)$y[include]

    pve <- 1-var(yy-pred.yy)/var(yy)

    return(list(pred.yy=pred.yy,
                pve=pve))
}

#' @name fit_loess
#'
#' @title Use loess to estimate cyclic trends of expression values
#'
#' @param yy A vector of gene expression values for one gene. The
#' expression values are assumed to have been normalized and
#' transformed to standard normal distribution.
#'
#' @param time A vector of angles (cell cycle phase).
#'
#' @return A list with one element, \code{pred.yy}, giving the
#' estimated cyclic trend.
#'
#' @author Joyce Hsiao
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
#' yy_ordered <- exprs_quant[1, names(theta_ordered)]
#'
#' fit <- fit_loess(yy_ordered, time=theta_ordered)
#'
#' plot(x=theta_ordered, y=yy_ordered, pch=16, cex=0.7, axes=FALSE,
#'   ylab="quantile-normalized expression values", xlab="FUCCI phase",
#'   main = "loess fit")
#' points(x=theta_ordered, y=fit$pred.yy, col="blue", pch=16, cex=0.7)
#' axis(2)
#' axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#'   labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
#'   expression(2*pi)))
#' abline(h=0, lty=1, col="black", lwd=0.7)
#'
#' @importFrom stats loess
#' @importFrom stats var
#' @importFrom stats predict
#' 
#' @export
#' 
fit_loess <- function(yy, time) {
    yy.rep <- rep(yy,3)
    time.rep <- c(time, time + (2*pi), time + (4*pi))
    include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

    fit <- loess(yy.rep ~ time.rep)
    pred.yy <- fit$fitted[include]

    pve <- 1 - var(yy - pred.yy)/var(yy)

    return(list(pred.yy=pred.yy,
                pve=pve))
}

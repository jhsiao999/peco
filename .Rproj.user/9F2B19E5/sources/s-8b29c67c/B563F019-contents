#' Using trendfiltering to estimate cyclic trend of gene expression
#'
#' We applied quadratic (second order) trend filtering using the
#' trendfilter function in the genlasso package (Tibshirani, 2014).
#' The trendfilter function implements a nonparametric smoothing method
#' which chooses the smoothing parameter by cross-validation and fits
#' a piecewise polynomial regression. In more specifics: The trendfilter
#' method determines the folds in cross-validation in a nonrandom manner.
#' Every k-th data point in the ordered sample is placed in the k-th fold,
#' so the folds contain ordered subsamples. We applied five-fold
#' cross-validation and chose the smoothing penalty using
#' the option lambda.1se: among all possible values of the penalty term,
#' the largest value such that the cross-validation standard error is within
#' one standard error of the minimum. Furthermore, we desired that the estimated
#' expression trend be cyclical. To encourage this, we concatenated the ordered
#' gene expression data three times, with one added after another. The quadratic
#' trend filtering was applied to the concatenated data series of each gene.
#'
#' @param yy A vector of gene expression values for one gene. The expression
#'     values are assumed to have been normalized and transformed to
#'     standard normal distribution.
#' @param polyorder We estimate cyclic trends of gene expression levels using
#'    nonparamtric trend filtering. The default fits second degree polynomials
#'    (polyorder=2).
#'
#' @return
#'     \describe{
#'      \item{\code{trend.yy}}{The estimated cyclic trend.}
#'      \item{\code{pve}}{Proportion of variance explained by the cyclic trend
#'          in the gene expression levels}
#'          }
#'
#' @author Joyce Hsiao
#'
#' @export
fit.trendfilter.generic <- function(yy, polyorder=2) {
  library(genlasso)

  yy.rep <- rep(yy,3)
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  fit.trend <- trendfilter(yy.rep,
                           ord=polyorder, approx = F, maxsteps = 1000)
  cv.trend <- cv.trendfilter(fit.trend)
  which.lambda <- cv.trend$i.1se
  yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                           df=fit.trend$df[which.lambda])$fit
  trend.yy <- yy.trend.pred
  pve <- 1-var(yy-trend.yy)/var(yy)

  return(list(trend.yy=trend.yy[which(include)],
              pve=pve))
}







#' Use bsplies to cyclic trend of gene expression levels
#'
#' @param yy A vector of gene expression values for one gene. The expression
#'     values are assumed to have been normalized and transformed to
#'     standard normal distribution.
#' @param time A vector of angels (cell cycle phase).
#'
#' @return
#'     \describe{
#'      \item{\code{pred.yy}}{The estimated cyclic trend.}
#'          }
#'
#' @author Joyce Hsiao
#'
#' @export
fit.bspline <- function(yy, time) {

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  time.rep <- c(time, time+(2*pi), time+(4*pi))
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  # trendfilter
  fit <- smooth.spline(x=time.rep, y=yy.rep)
  pred.yy <- predict(fit, time.rep)$y[which(include==TRUE)]

  return(list(pred.yy=pred.yy))
}




#' Use loess to estimate cyclic trends of expression values
#'
#' @param yy A vector of gene expression values for one gene. The expression
#'     values are assumed to have been normalized and transformed to
#'     standard normal distribution.
#' @param time A vector of angels (cell cycle phase).
#'
#' @return
#'     \describe{
#'      \item{\code{pred.yy}}{The estimated cyclic trend.}
#'          }
#'
#' @author Joyce Hsiao
#'
#' @export
fit.loess <- function(yy, time) {

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  time.rep <- c(time, time+(2*pi), time+(4*pi))
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  # trendfilter
  fit <- loess(yy.rep~time.rep)
  pred.yy <- fit$fitted[which(include==TRUE)]

  return(list(pred.yy=pred.yy))
}

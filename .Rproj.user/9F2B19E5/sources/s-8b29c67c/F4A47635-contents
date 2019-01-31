
#' @title fitting trendfilter
#'
#' @param yy log2cpm gene expression vector
#'
#' @export
fit.trendfilter.generic <- function(yy, pos.yy=c(1:length(yy)), polyorder=2) {
  library(genlasso)

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  pos.rep <- rep(pos.yy,3)
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  # trendfilter
  fit.trend <- trendfilter(yy.rep,
                           ord=polyorder, approx = F, maxsteps = 1000)
  cv.trend <- cv.trendfilter(fit.trend)
  which.lambda <- cv.trend$i.1se
  yy.trend.pred <- predict(fit.trend, lambda=cv.trend$lambda.1se,
                           df=fit.trend$df[which.lambda])$fit
  trend.yy <- yy.trend.pred
  pve <- 1-var(yy-trend.yy[which(include)])/var(yy)

  # trend.mad.pred <- mad(yy.nonzero.rep[include]-yy.trend.pred[include], constant = 1)
  # trend.mad.constant <- mad(yy.nonzero.rep[include]-mean(yy.nonzero.rep[include]), constant = 1)
  # trend.mad.ratio <- trend.mad.pred/trend.mad.constant

  return(list(trend.yy=trend.yy[which(include)],
              pve=pve))
}



#' @title fitting bsplines
#'
#' @param yy log2cpm gene expression vector
#'
#' @export
fit.bspline <- function(yy, time, polyorder=2) {

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  time.rep <- c(time, time+(2*pi), time+(2*pi))
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  # trendfilter
  fit.bspline <- smooth.spline(x=time.rep, y=yy.rep)
  pred.yy <- fit.bspline$y[which(include)]

  return(list(pred.yy=pred.yy))
}




#' @title fitting loess
#'
#' @param yy log2cpm gene expression vector
#'
#' @export
fit.loess <- function(yy, time, polyorder=2) {

  yy.rep <- rep(yy,3)
  #  theta.nonzero.rep <- rep(theta.nonzero,3)
  time.rep <- c(time, time+(2*pi), time+(2*pi))
  include <- rep(c(FALSE, TRUE, FALSE), each = length(yy))

  # trendfilter
  fit.loess <- loess(yy.rep~time.rep)
  pred.yy <- fit.loess$fitted[which(include)]

  return(list(pred.yy=pred.yy))
}


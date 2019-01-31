#' Infer angles for each single-cell samples using fluorescence intensities
#'
#' We use FUCCI intensities to infer the position of the cells in cell
#' cycle progression. The result is a vector of angles on a unit circle
#' corresponding to the positions of the cells in cell cycle progression.
#'
#' @param mat A matrix of two columns of summarized fluorescence intensity.
#'     Rows correspond to samples.
#' @param method The method used to fit the circle. \code{trig} uses trignometry to
#'     transform intensity measurements from cartesian coordinates to
#'     polar coordinates. \code{algebraic} uses an algebraic approach for
#'     circle fitting, using the \code{conicfit} package.
#' @param plot.it T/F. Plot the fitted results.
#'
#' @return The inferred angles on unit circle based on the input intensity
#'     measurements.
#'
#' @author Joyce Hsiao
#'
#' @import stats, circular, conicfit
#'
#' @export
intensity2circle <- function(mat, plot.it=F, method=c("trig",
                                                      "algebraic")) {

  if (is.matrix(mat)==F) { mat <- as.matrix(mat) }

  if (method=="trig") {
#    library(circular)
    pca <- prcomp(mat, scale. = T)
    pca_scores <- pca$x
    theta <- as.numeric(circular::coord2rad(pca_scores[,1:2]))

    if (plot.it==TRUE) {
      rng <- c(-1,1)*max(abs(range(pca_scores[,1:2])))
      plot(pca_scores[,1:2], pch=16,col="gray80", cex=.7,
           xlim=rng, ylim=rng)  #;par(new=TRUE)
      points(x=cos(theta),y=sin(theta), pch=16,col="black", cex=.7)
    }
    # convert to angle
    theta <- theta%%(2*pi)
    return(theta)
  }

  if (method=="algebraic") {
#    library(conicfit)

    pca <- prcomp(mat, scale. = F)
    pca_scores <- pca$x
    ellipDirect <- conicfit::EllipseDirectFit(pca_scores[,1:2])
    ellipDirectG <- conicfit::AtoG(ellipDirect)$ParG
    xyDirect<- conicfit::calculateEllipse(ellipDirectG[1], ellipDirectG[2], ellipDirectG[3],
                                          ellipDirectG[4], 180/pi*ellipDirectG[5])
    ellipProj <- conicfit::Residuals.ellipse(pca_scores[,1:2], ellipDirectG)

    if (plot.it==TRUE) {
      rng <- c(-1,1)*max(abs(range(xyDirect)))
      plot(pca_scores[,1:2], pch=16,col="gray80", cex=.7,
           xlim=rng, ylim=rng)  #;par(new=TRUE)
      points(ellipProj$XYproj, pch=16,col="black", cex=.7)
    }
    # convert to angle
    ang <- atan2(ellipDirectG[3]/ellipDirectG[4]*ellipProj$XYproj[,2], ellipProj$XYproj[,1])
    ang <- ang%%(2*pi)

    return(ang)
  }

}

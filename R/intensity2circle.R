#' @title Infer angles for each single-cell samples using fluorescence
#' intensities
#'
#' @description We use FUCCI intensities to infer the position of the
#' cells in cell cycle progression. The result is a vector of angles
#' on a unit circle corresponding to the positions of the cells in
#' cell cycle progression.
#'
#' @param mat A matrix of two columns of summarized fluorescence
#' intensity. Rows correspond to samples.
#' @param method The method used to fit the circle. \code{trig} uses
#' trignometry to transform intensity measurements from cartesian
#' coordinates to polar coordinates. \code{algebraic} uses an
#' algebraic approach for circle fitting, using the conicfit package.
#' @param plot.it TRUE or FALSE. Plot the fitted results.
#'
#' @return The inferred angles on unit circle based on the input
#' intensity measurements.
#'
#' @examples
#' # use our data
#' library(Biobase)
#'
#' data(eset_sub)
#'
#' # FUCCI scores - log10 transformed sum of intensities that were
#' # corrected for background noise
#' ints <- pData(eset_sub)[,c("rfp.median.log10sum.adjust",
#'                           "gfp.median.log10sum.adjust")]
#' intensity2circle(ints, plot.it=TRUE, method = "trig")
#'
#' @author Joyce Hsiao
#'
#' @importFrom stats prcomp
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom circular coord2rad
#' @importFrom conicfit AtoG
#' @importFrom conicfit EllipseDirectFit
#' @importFrom conicfit calculateEllipse
#' @importFrom conicfit Residuals.ellipse
#'
#' @export
#'
intensity2circle <- function(mat, plot.it=FALSE,
                             method=c("trig","algebraic")) {

  if (is.matrix(mat)==FALSE)
    mat <- as.matrix(mat)

  if (method=="trig") {
    pca <- prcomp(mat, scale. = TRUE)
    pca_scores <- pca$x
    theta <- as.numeric(coord2rad(pca_scores[,c(1,2)]))

    if (plot.it==TRUE) {
      rng <- c(-1,1)*max(abs(range(pca_scores[,c(1,2)])))
      plot(pca_scores[,c(1,2)], pch=16,col="gray80", cex=.7,
           xlim=rng, ylim=rng)
      points(x=cos(theta),y=sin(theta), pch=16,col="black", cex=.7)
    }

    # convert to angle
    theta <- theta%%(2*pi)
    return(theta)
  }
  else if (method=="algebraic") {
    pca <- prcomp(mat, scale. = FALSE)
    pca_scores <- pca$x
    ellipDirect <- EllipseDirectFit(pca_scores[,c(1,2)])
    ellipDirectG <- AtoG(ellipDirect)$ParG
    xyDirect <- calculateEllipse(ellipDirectG[1], ellipDirectG[2],
                                 ellipDirectG[3],
                                 ellipDirectG[4], 180/pi*ellipDirectG[5])
    ellipProj <- Residuals.ellipse(pca_scores[,c(1,2)], ellipDirectG)

    if (plot.it==TRUE) {
      rng <- c(-1,1)*max(abs(range(xyDirect)))
      plot(pca_scores[,c(1,2)], pch=16,col="gray80", cex=.7,
           xlim=rng, ylim=rng)
      points(ellipProj$XYproj, pch=16,col="black", cex=.7)
  }

    # convert to angle
    ang <- atan2(ellipDirectG[3]/ellipDirectG[4]*ellipProj$XYproj[,2],
                 ellipProj$XYproj[,1])
    ang <- ang%%(2*pi)

    return(ang)
  } else
    stop("Invalid value for argument \"method\"")
}

#' @title Infer Angles for Single-cell Samples Using Fluorescence
#'   Intensities
#'
#' @description We use FUCCI intensities to infer the position of the
#'   cells in cell cycle progression. The result is a vector of angles
#'   on a unit circle corresponding to the positions of the cells in
#'   cell cycle progression.
#'
#' @param mat A matrix with two columns giving summarized fluorescence
#'   intensities. Rows correspond to samples.
#' 
#' @param method The method used to fit the circle. \code{method =
#'   "trig"} uses trignometry to transform intensity measurements from
#'   cartesian coordinates to polar coordinates; \code{method =
#'   "algebraic"} uses an algebraic approach for circle fitting,
#'   implemented in the \code{conicfit} package.
#' 
#' @param plot.it \code{TRUE} or \code{FALSE}. If \code{plot.it =
#'   TRUE}, plot the fitted results.
#'
#' @return The inferred angles on a unit circle based on the input
#'   intensity measurements.
#'
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # Compute FUCCI scores---the log10-transformed sum of intensities
#' # corrected for background noise.
#' ints <- colData(sce_top101genes)[,c("rfp.median.log10sum.adjust",
#'                                     "gfp.median.log10sum.adjust")]
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
    method <- match.arg(method)
    if (!is.matrix(mat))
        mat <- as.matrix(mat)

    if (method == "trig") {
        pca <- prcomp(mat, scale. = TRUE)
        pca_scores <- pca$x
        theta <- as.numeric(coord2rad(pca_scores[,1:2]))
        if (plot.it) {
            rng <- c(-1,1) * max(abs(range(pca_scores[,1:2])))
            plot(pca_scores[,1:2], pch=16, col="gray", cex=0.7,
                 xlim=rng, ylim=rng)
            points(x=cos(theta), y=sin(theta), pch=16,col="darkblue", cex=0.5)
        }
        theta <- theta %% (2*pi)
        return(theta)
    } else if (method == "algebraic") {
        pca <- prcomp(mat, scale. = FALSE)
        pca_scores <- pca$x
        ellipDirect <- EllipseDirectFit(pca_scores[,1:2])
        ellipDirectG <- AtoG(ellipDirect)$ParG
        xyDirect <- calculateEllipse(ellipDirectG[1], ellipDirectG[2],
                                     ellipDirectG[3], ellipDirectG[4],
                                     180/pi*ellipDirectG[5])
        ellipProj <- Residuals.ellipse(pca_scores[,1:2], ellipDirectG)
        if (plot.it) {
            rng <- c(-1,1) * max(abs(range(xyDirect)))
            plot(pca_scores[,1:2], pch=16, col="gray", cex=0.7,
                 xlim=rng, ylim=rng)
            points(ellipProj$XYproj, pch=16, col="darkblue", cex=0.5)
        }
        ang <- atan2(ellipDirectG[3] / ellipDirectG[4] * ellipProj$XYproj[,2],
                     ellipProj$XYproj[,1])
        ang <- ang %% (2*pi)

        return(ang)
    } else
        stop("Invalid value for argument \"method\"")
}

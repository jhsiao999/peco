#' @title Quantile-normalize CPM for Each Gene
#' 
#' @title Transform counts by first computing counts-per-million
#'   (CPM), then quantile-normalize CPM for each gene.
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object.
#' 
#' @param ncores Argument passed to
#'   \code{\link[parallel]{makeCluster}} specifying the number of
#'   threads.
#'
#' @return SingleCellExperiment object with an additional
#'   \dQuote{cpm_quant} slot; this is added if it doesn't already exist.
#'
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # Perform CPM normalization using scater and quantile-normalize
#' # the CPM values of each gene to normal distribution.
#' sce_top101genes <- data_transform_quantile(sce_top101genes, ncores=2)
#'
#' plot(y=assay(sce_top101genes, "cpm_quantNormed")[1,],
#'      x=assay(sce_top101genes, "cpm")[1,],
#'      xlab = "CPM before quantile-normalization",
#'      ylab = "CPM after quantile-normalization")
#'
#' @author Joyce Hsiao
#'
#' @importFrom SingleCellExperiment cpm
#' @importFrom SingleCellExperiment cpm<-
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment assays<-
#' @importFrom SummarizedExperiment assayNames<-
#' @importFrom assertthat has_name
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom scater calculateCPM
#' @importFrom stats qqnorm
#'
#' @export
#' 
data_transform_quantile <- function(sce, ncores=2) {

    if (is.null(ncores)) {
        cl <- makeCluster(2)
        registerDoParallel(cl)
        message(paste("computing with",ncores,"threads"))
    } else {
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        message(paste("computing with",ncores,"threads"))
    }

    # Check if there's already cpm normalied data, if yes, skip
    # this step.
    if(inherits(sce, "SingleCellExperiment")) {
        if (!has_name(assays(sce), "cpm"))
            cpm(sce) <- calculateCPM(sce)
        G <- nrow(sce)
        cpm_sce <- assay(sce, "cpm")

        cpm_quantNormed <- foreach(g = seq_len(G)) %dopar% {
            y_g <- cpm_sce[g,]
            is.zero <- which(y_g == 0)
            qq.map <- qqnorm(y_g, plot.it=FALSE)
            if (is.null(is.zero))
              yy.qq <- qq.map$x
            else {
              yy.qq <- qq.map$x
              yy.qq[is.zero] <- sample(qq.map$x[is.zero])
            }
            return(yy.qq)
        }
        stopCluster(cl)
        cpm_quantNormed <- do.call(rbind, cpm_quantNormed)
        colnames(cpm_quantNormed) <- colnames(cpm_sce)
        rownames(cpm_quantNormed) <- rownames(cpm_sce)

        assays(sce)[[3]] <- cpm_quantNormed
        assayNames(sce)[3] <- "cpm_quantNormed"
        return(sce)
    } else {
        cpm_quantNormed <- calculateCPM(sce)
        return(cpm_quantNormed)
    }
}

#' @name data_transform_quantile
#'
#' @title Transform counts by first computing counts-per-million (CPM), then
#'   quantile-normalize CPM for each gene
#'
#' @description
#'   For each gene, transform counts to CPM and then to a normal distribution.
#'
#' @param sce SingleCellExperiment Object.
#' @param ncores We use doParallel package for parallel computing.
#'
#' @return SingleCellExperiment Object with an added slot of cpm_quant,
#'     cpm slot is added if it doesn't exist.
#'
#' @examples
#' # use our data
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # perform CPM normalization using scater, and
#' # quantile-normalize the CPM values of each gene to normal distribution
#' sce_top101genes <- data_transform_quantile(sce_top101genes, ncores=2)
#'
#' plot(y=assay(sce_top101genes, "cpm_quantNormed")[1,],
#'      x=assay(sce_top101genes, "cpm")[1,],
#'     xlab = "CPM bbefore quantile-normalization",
#'     ylab = "CPM after quantile-normalization")
#'
#' @author Joyce Hsiao
#'
#' @importFrom SingleCellExperiment SingleCellExperiment cpm
#' @importFrom SummarizedExperiment assay assays assayNames
#' @importFrom SummarizedExperiment assays<- assayNames<-
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom scater calculateCPM
#' @importFrom stats qqnorm
#' @import foreach
#' @import methods
#'
#' @export
data_transform_quantile <- function(sce, ncores=2) {

    if (is.null(ncores)) {
        cl <- makeCluster(2)
        registerDoParallel(cl)
        message(paste("computing on",ncores,"cores"))
    } else {
        cl <- makeCluster(ncores)
        registerDoParallel(cl)
        message(paste("computing on",ncores,"cores"))
    }

    # check if there's already cpm normalied data,
    # if yes, then skip this step
    if (has_name(assays(sce), "cpm")) {
        sce <- sce
    } else {
        cpm(sce) <- calculateCPM(sce)
    }

    G <- nrow(sce)
    cpm_sce <- cpm(sce)

    cpm_quantNormed <- foreach(g=seq_len(G)) %dopar% {
        y_g <- cpm_sce[g,]
        is.zero <- which(y_g == 0)
        qq.map <- qqnorm(y_g, plot.it=FALSE)
        if (is.null(is.zero)) {
          yy.qq <- qq.map$x
        } else {
          yy.qq <- qq.map$x
          yy.qq[is.zero] <- sample(qq.map$x[is.zero])
        }
        return(y_g= yy.qq)
    }
    stopCluster(cl)
    cpm_quantNormed <- do.call(rbind, cpm_quantNormed)
    colnames(cpm_quantNormed) <- colnames(cpm_sce)
    rownames(cpm_quantNormed) <- rownames(cpm_sce)

    assays(sce)[[3]] <- cpm_quantNormed
    assayNames(sce)[3] <- "cpm_quantNormed"

    return(sce)
}


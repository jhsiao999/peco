# cycle_npreg_insample and cycle_npreg_outsample are the two main
# functions in peco. cycle_npreg_insample generates cyclic trend
# estimates of gene expression levels using training data, and
# cycle_npreg_outsample applies the estimates of cycle_npreg_insample
# to another gene expression dataset to infer an angle or cell cycle
# phase for each cell.


#' @name cycle_npreg_outsample
#'
#' @title Predict test-sample ordering using training labels (no update)
#'
#' @description Apply the estimates of cycle_npreg_insample to another
#' gene expression dataset to infer an angle or cell cycle phase for
#' each cell.
#'
#' @param Y_test A SingleCellExperiment object.
#'
#' @param sigma_est Input from training data. A vector of gene-specific
#'   standard error of the cyclic trends.
#'
#' @param funs_est Input fron training data. A vector of cyclic functions
#'   estimating cyclic trends.
#'
#' @param polyorder We estimate cyclic trends of gene expression levels using
#'    nonparamtric trend filtering. The default fits second degree polynomials.
#' @param normed Is the data already normalized? TRUE or FALSE.
#' @param method.grid Method for defining bins along the circle.
#' @param method.trend Varous methods that can be applied to estimate
#' cyclic trend of gene expression levels.
#'
#' @param grids number of bins to be selected along 0 to 2pi.
#'
#' @param get_trend_estimates To re-estimate the cylic trend based on
#'   the predicted cell cycle phase or not (T or F). Default FALSE. This step
#'   calls trendfilter and is computationally intensive.
#' @param ncores We use doParallel package for parallel computing.
#'
#' @inheritParams cycle_npreg_mstep
#' @inheritParams cycle_npreg_loglik
#'
#' @return A list with the following elements:
#'
#' \item{Y}{The input gene expression marix.}
#' \item{cell_times_est}{Inferred angles or cell cycle phases, NOT
#' ordered.}
#' \item{loglik_est}{Log-likelihood estimates for each gene.}
#' \item{cell_times_reordered}{The inferred angles reordered (in
#' ascending order).}
#' \item{Y_reorded}{The input gene expression matrix reordered by
#' cell_times_reordered.}
#' \item{sigma_reordered}{Estimated standard error of the cyclic trend
#' for each gene, reordered by cell_times_reordered.}
#' \item{funs_reordered}{A list of functions for approximating the
#' cyclic trends of gene express levels for each gene, reordered by
#' cell_times_reordered.}
#' \item{mu_reordered}{Estimated cyclic trend of gene expression
#' values for each gene, reordered by cell_times_reordered.}
#' \item{prob_per_cell_by_celltimes}{Probabilities of each cell belong
#' to each bin.}
#'
#' @examples
#' # import data
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # select top 5 cyclic genes
#' sce_top5 <- sce_top101genes[order(rowData(sce_top101genes)$pve_fucci,
#'                                   decreasing=TRUE)[1:5],]
#'
#' # Select samples from NA18511 for our prediction example
#' coldata <- colData(sce_top5)
#' which_samples_train <- rownames(coldata)[coldata$chip_id != "NA18511"]
#' which_samples_predict <- rownames(coldata)[coldata$chip_id == "NA18511"]
#'
#' # learning cyclic functions of the genes using our training data
#' sce_top5 <- data_transform_quantile(sce_top5)
#' expr_quant <- assay(sce_top5, "cpm_quantNormed")
#' Y_train <- expr_quant[, colnames(expr_quant) %in% which_samples_train]
#' theta_train <-
#'     coldata$theta_shifted[rownames(coldata) %in% which_samples_train]
#' names(theta_train) <-
#'     rownames(coldata)[rownames(coldata) %in% which_samples_train]
#'
#' # obtain cyclic function estimates
#' model_5genes_train <- cycle_npreg_insample(Y = Y_train,
#'                                            theta = theta_train,
#'                                            polyorder=2,
#'                                            ncores=2,
#'                                            method.trend="trendfilter")
#'
#' # predict cell cycle
#' model_5genes_predict <- cycle_npreg_outsample(
#'   Y_test=sce_top5[,colnames(sce_top5) %in% which_samples_predict],
#'   sigma_est=model_5genes_train$sigma_est,
#'   funs_est=model_5genes_train$funs_est,
#'   method.trend="trendfilter",
#'   ncores=2,
#'   get_trend_estimates=FALSE)
#'
#' # estimate cyclic gene expression levels given cell cycle for each gene
#' predict_cyclic <-
#'     fit_cyclical_many(Y=assay(model_5genes_predict,"cpm_quantNormed"),
#'                       theta=colData(model_5genes_predict)$cellcycle_peco)
#' all.equal(rownames(predict_cyclic[[2]]), rownames(predict_cyclic[[1]]))
#'
#' par(mfrow=c(2,3), mar=c(4,4,3,1))
#' for (g in seq_along(rownames(model_5genes_predict))) {
#'   plot(assay(model_5genes_predict,"cpm_quantNormed")[
#'       rownames(model_5genes_predict)[g],],
#'        x=colData(model_5genes_predict)$cellcycle_peco, axes=FALSE,
#'        xlab="FUCCI phase",
#'        ylab="Predicted phase")
#'   points(y=predict_cyclic$cellcycle_function[[
#'            rownames(model_5genes_predict)[g]]](
#'     seq(0, 2*pi, length.out = 100)),
#'     x=seq(0, 2*pi, length.out = 100),
#'     pch=16, col="royalblue")
#'   axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#'                 labels=c(0,expression(pi/2), expression(pi),
#'                          expression(3*pi/2), expression(2*pi)))
#'   abline(h=0, lty=1, col="black", lwd=.7)
#'   title(rownames(model_5genes_predict$Y_reordered)[g])
#' }
#' title("Predicting cell cycle phase for NA18511", outer=TRUE)
#'
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays assayNames colData colData<-
#' @importFrom assertthat has_name
#' @import methods
#'
#' @family peco classifier functions
#' @seealso \code{\link{cycle_npreg_insample}} for obtaining parameteres for
#'     cyclic functions from training data,
#'     \code{\link{cycle_npreg_loglik}} for log-likehood at
#'     angles between 0 to 2pi,
#'     \code{\link{initialize_grids}} for selecting
#'     angles in \code{\link{cycle_npreg_loglik}},
#'     \code{\link{cycle_npreg_mstep}} for estimating cyclic functions given
#'     inferred phases from  \code{\link{cycle_npreg_loglik}}
#' @export
cycle_npreg_outsample <- function(Y_test,
                                    sigma_est,
                                    funs_est,
                                    method.trend=c("trendfilter",
                                                    "loess", "bspline"),
                                    normed = TRUE,
                                    polyorder=2,
                                    method.grid="uniform",
                                    ncores=2,
                                    grids=100,
                                    get_trend_estimates=FALSE) {

    if (!normed & is(Y_test, "SingleCellExperiment")) {
      Y_test <- data_transform_quantile(Y_test)
      exprs_test <- assay(Y_test, "cpm_quantNormed")
    } else if (normed & is(Y_test, "SingleCellExperiment")) {
      exprs_test <- assay(Y_test, "cpm_quantNormed")
    } else if (!normed & !is(Y_test, "SingleCellExperiment")) {
      exprs_test <- data_transform_quantile(Y_test)
    } else {
      exprs_test <- Y_test
    }

    # compute expected cell time for the test samples
    # under mu and sigma estimated from the training samples
    initial_loglik <- cycle_npreg_loglik(Y = exprs_test,
                                        sigma_est = sigma_est,
                                        method.grid = method.grid,
                                        funs_est = funs_est,
                                        grids = grids)

    if (is(Y_test, "SingleCellExperiment")) {
      colData(Y_test)$cellcycle_peco <- initial_loglik$cell_times_est
      colData(Y_test)@elementMetadata$labelDescription[ncol(colData(Y_test))] <-
        "peco predicted continuous cell cycle, between 0 to 2pi"
    }

    if (get_trend_estimates) {
        updated_estimates <- cycle_npreg_mstep(Y = Y_test,
                                        theta = initial_loglik$cell_times_est,
                                        method.trend = method.trend,
                                        polyorder = polyorder,
                                        ncores = ncores)
        out <- list(
            Y=Y_test,
            cell_times_est=initial_loglik$cell_times_est,
            loglik_est=initial_loglik$loglik_est,
            Y_reordered=updated_estimates$Y,
            cell_times_reordered=updated_estimates$theta,
            mu_reordered=updated_estimates$mu_est,
            sigma_reordered=updated_estimates$sigma_est,
            funs_reordered=updated_estimates$funs,
            prob_per_cell_by_celltimes=initial_loglik$prob_per_cell_by_celltimes)
        } else {
            out <- Y_test
        }
    return(out)
}

#' @name cycle_npreg_insample
#'
#' @title Obtain cyclic trend estimates from the training data
#'
#' @description Estimates cyclic trends of gene expression levels
#' using training data.
#'
#' @param Y A matrix of normalized and transformed gene expression
#' values.  Gene by sample.
#'
#' @param theta A vector of angles.
#'
#' @param ncores We use doParallel package for parallel computing.
#'
#' @param polyorder We estimate cyclic trends of gene expression
#' levels using nonparamtric trend filtering. The default fits second
#' degree polynomials.
#'
#' @param method.trend Varous methods that can be applied to estimate
#' cyclic trend of gene expression levels.
#'
#' @inheritParams cycle_npreg_mstep
#'
#' @return A list with four elements:
#'
#' \item{Y}{Gene expression marix.}
#'
#' \item{theta}{Vector of angles or cell cycle phases.}
#'
#' \item{sigma_est}{Estimated standard error of the cyclic trend for
#' each gene.}
#'
#' \item{funs_est}{A list of functions for approximating the cyclic
#' trends of gene express levels for each gene.}
#'
#' @examples
#' # see \code{\link{cycle_npreg_insample}}
#'
#' @author Joyce Hsiao
#'
#' @family peco classifier functions
#' @seealso
#'        \code{\link{cycle_npreg_mstep}} for estimating cyclic functions given
#'        inferred phases from  \code{\link{cycle_npreg_loglik}},
#'        \code{\link{cycle_npreg_outsample}} for predicting cell cycle phase
#'        using parameters learned from \code{\link{cycle_npreg_insample}}
#'
#' @export
cycle_npreg_insample <- function(Y, theta,
                                ncores=2,
                                polyorder=2,
                                method.trend=c("trendfilter",
                                                "loess", "bspline")) {

    # order data by initial cell times
    G <- nrow(Y)
    N <- ncol(Y)
    theta_ordered_initial <- theta[order(theta)]
    Y_ordered <- Y[,order(theta)]

    # initialize mu and sigma
    initial_mstep <- cycle_npreg_mstep(Y = Y_ordered,
                                    theta = theta_ordered_initial,
                                    polyorder=polyorder,
                                    method.trend=method.trend,
                                    ncores = ncores)

    return(list(Y=Y,
                theta=theta,
                sigma_est=initial_mstep$sigma_est,
                funs_est=initial_mstep$funs))
}



#------ Supporting functions

#' @name initialize_grids
#'
#' @title For prediction, initialize grid points for cell cycle phase
#' on a circle.
#'
#' @param Y Gene expression matrix. Gene by sample.
#'
#' @param grids number of bins to be selected along 0 to 2pi.
#'
#' @param method.grid The approach to initialize angles in the
#' computation. \code{uniform} creates k equally-spaced bins
#' (grids). \code{pca} uses gene expression values to infer angles,
#' and then use these pca-based angles to move the cells to the
#' closest bin (as defined by \code{uniform}).
#'
#' @return A vector of initialized angles to be used in
#' \code{cycle_npreg_loglik} to infer angles.
#'
#' @family peco classifier functions
#' @seealso
#'     \code{\link{cycle_npreg_loglik}} for log-likehood at
#'     angles between 0 to 2pi,
#'     \code{\link{cycle_npreg_mstep}} for estimating cyclic functions given
#'     inferred phases from  \code{\link{cycle_npreg_loglik}},
#'     \code{\link{cycle_npreg_outsample}} for predicting cell cycle phase
#'      using parameters learned from \code{\link{cycle_npreg_insample}}
#' @author Joyce Hsiao
#'
initialize_grids <- function(Y, grids=100,
                                method.grid=c("pca", "uniform")) {

    len <- (2*pi)/(2*grids)
    theta_grids <- seq(len, (2*pi)-(len), length.out=grids)

    if (method.grid=="pca") {
        pc_res <- prcomp(t(Y), scale = TRUE)
        grid_approx <- coord2rad(cbind(pc_res$x[,1], pc_res$x[,2]))
        grid_approx <- as.numeric(grid_approx)
        names(grid_approx) <- colnames(Y)

        theta_initial <- rep(0, length(grid_approx))
        theta_initial_ind <- rep(0, length(grid_approx))
        names(theta_initial) <- names(grid_approx)
        names(theta_initial_ind) <- names(grid_approx)

    for (i in seq_along(grid_approx)) {
        theta_initial_ind[i] <-
            which.min(pmin(abs(theta_grids-grid_approx[i]),
                        abs(theta_grids-(2*pi-grid_approx[i]))))
        theta_initial <- theta_grids[theta_initial_ind]
    }
    }

    if (method.grid=="uniform") {
        theta_initial <- theta_grids
    }

    return(theta_initial)
}

#' @name cycle_npreg_loglik
#'
#' @title Infer angles or cell cycle phase based on gene expression data
#'
#' @param Y Gene by sample expression matrix.
#' @param sigma_est A vector of standard errors for each gene from the
#' training data.
#' @param funs_est A vector of cyclic functions estimated for each
#' gene from the training data.
#'
#' @inheritParams initialize_grids
#'
#' @return A list with the following three elements:
#'
#'     \item{cell_times_est}{Inferred angles or cell cycle phases, NOT
#' ordered.}
#'     \item{loglik_est}{Log-likelihood estimates for each gene.}
#'     \item{prob_per_cell_by_celltimes}{Probabilities of each cell belong
#' to each bin.}
#'
#' @importFrom stats dnorm
#' @import methods
#'
#' @family peco classifier functions
#' @seealso \code{\link{initialize_grids}} for selecting
#'     angles in \code{\link{cycle_npreg_loglik}},
#'     \code{\link{cycle_npreg_mstep}} for estimating cyclic functions given
#'     inferred phases from  \code{\link{cycle_npreg_loglik}},
#'     \code{\link{cycle_npreg_outsample}} for predicting cell cycle phase
#'      using parameters learned from \code{\link{cycle_npreg_insample}}
#'
#' @author Joyce Hsiao
cycle_npreg_loglik <- function(Y, sigma_est, funs_est,
                                grids=100,
                                method.grid=c("pca", "uniform")) {

    N <- ncol(Y); G <- nrow(Y)
    theta_choose <- initialize_grids(Y, grids=grids, method.grid=method.grid)

    # for each cell, sum up the loglikelihood for each gene
    # at the observed cell times
    loglik_per_cell_by_celltimes <- matrix(0, N, grids)
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    for (n in seq_len(N)) {
        loglik_per_cell <- do.call(rbind, lapply(seq_len(G), function(g) {
        dnorm(Y[g,n], funs_est[[g]](theta_choose), sigma_est[g], log = TRUE)
    }))
    loglik_per_cell <- colSums(loglik_per_cell)
    loglik_per_cell_by_celltimes[n,] <- loglik_per_cell }

    # use max likelihood to assign samples
    prob_per_cell_by_celltimes <- t(apply(loglik_per_cell_by_celltimes, 1,
                                        function(x) {
            sumll <- sum(exp(x), na.rm=TRUE)
            if (sumll == 0) { return(rep(0, grids))
                } else { return(exp(x)/sumll) } }) )
    colnames(prob_per_cell_by_celltimes) <- theta_choose

    cell_times_samp_ind <- apply(prob_per_cell_by_celltimes, 1, function(x) {
        if (max(x, na.rm=TRUE)==0) { sample(seq_len(grids), 1, replace=FALSE)
            } else { which.max(x) } })
    names(cell_times_samp_ind) <-
        colnames(prob_per_cell_by_celltimes)[cell_times_samp_ind]

    cell_times_est <- do.call(c, lapply(seq_len(N), function(n) {
        theta_choose[cell_times_samp_ind[n]] }) )
    names(cell_times_est) <- colnames(Y)

    # compute likelihood based on the selected cell times
    loglik_max_per_cell <- do.call(c, lapply(seq_len(N), function(n) {
        ll <- loglik_per_cell_by_celltimes[n,]
        ll[cell_times_samp_ind[n]] }) )
    loglik_est <- sum(loglik_max_per_cell)

    return(list(loglik_est=loglik_est,
                cell_times_est=cell_times_est,
                prob_per_cell_by_celltimes=prob_per_cell_by_celltimes))
}


#' @name cycle_npreg_mstep
#'
#' @title Estimate parameters of the cyclic trends
#'
#' @description This is used in both cycle_npreg_insample (training
#' data fitting) and cycle_npreg_outsample (testing data prediction)
#' to estimate cyclic trends of gene expression values. The function
#' outputs for each gene standard error of the cyclic trend, cyclic
#' function, and the estimated expression levels given the cyclic
#' function.
#'
#' @param Y Gene by sample expression matrix (log2CPM).
#' @param theta Observed cell times.
#' @param method.trend How to estimate cyclic trend of gene expression values?
#'     We offer three options: 'trendfilter' (\code{fit_trendfilter_generic()}),
#'     'loess' (\code{fit_loess()}) and 'bsplines' (\code{fit_bspline()}).
#'     'trendfilter' provided the best fit in our study. But 'trendfilter`
#'     uses cross-validation and takes some time. Therefore, we recommend
#'     using bspline for quick results.
#' @param ncores How many computing cores to use? We use doParallel package for
#'     parallel computing.
#'
#' @inheritParams cycle_npreg_insample
#'
#' @return A list with the following elements:
#'
#' \item{Y}{Input gene expression data.}
#' \item{theta}{Input angles.}
#' \item{mu_est}{Estimated expression levels given the cyclic function
#'       for each gene.}
#' \item{sigma_est}{Estimated standard error of the cyclic trends for
#'       each gene}
#' \item{funs}{Estimated cyclic functions}
#'
#' @importFrom stats approxfun
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @import foreach
#'
#' @family peco classifier functions
#' @seealso
#'     \code{\link{cycle_npreg_insample}} for estimating cyclic functions
#'     given known phasesfrom training data,
#'     \code{\link{cycle_npreg_outsample}} for predicting cell cycle phase
#'      using parameters learned from \code{\link{cycle_npreg_insample}}
#'
#' @author Joyce Hsiao
cycle_npreg_mstep <- function(Y, theta, method.trend=c("trendfilter",
                                                        "loess", "bspline"),
                                polyorder=2, ncores=2) {

    if (is.null(ncores)) {
        cl <- makeCluster(2); registerDoParallel(cl)
        message(paste("computing on",ncores,"cores"))
    } else {
        cl <- makeCluster(ncores); registerDoParallel(cl)
        message(paste("computing on",ncores,"cores")) }

    G <- nrow(Y); N <- ncol(Y)

    Y_ordered <- Y[,names(theta)]
    ord <- order(theta)
    theta_ordered <- theta[ord]
    Y_ordered <- Y_ordered[,ord]

    fit <- foreach(g=seq_len(G)) %dopar% {
        y_g <- Y_ordered[g,]

        if (method.trend=="trendfilter") {
            fit_g <- fit_trendfilter_generic(yy=y_g, polyorder = polyorder)
            fun_g <- approxfun(x=as.numeric(theta_ordered),
                            y=as.numeric(fit_g$trend.yy), rule=2)
            mu_g <- fit_g$trend.yy
        }
        if (method.trend=="bspline") {
            fit_g <- fit_bspline(yy=y_g, time = theta_ordered)
            fun_g <- approxfun(x=as.numeric(theta_ordered),
                            y=as.numeric(fit_g$pred.yy), rule=2)
            mu_g <- fit_g$pred.yy
        }

        if (method.trend=="loess") {
            fit_g <- fit_loess(yy=y_g, time = theta_ordered)
            fun_g <- approxfun(x=as.numeric(theta_ordered),
                            y=as.numeric(fit_g$pred.yy), rule=2)
            mu_g <- fit_g$pred.yy
        }

        sigma_g <- sqrt(sum((y_g-mu_g)^2)/N)

        list(y_g =y_g,
            mu_g=mu_g,
            sigma_g=sigma_g,
            fun_g=fun_g)
    }
    stopCluster(cl)

    sigma_est <- do.call(c, lapply(fit, "[[", "sigma_g"))
    names(sigma_est) <- rownames(Y_ordered)

    mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
    colnames(mu_est) <- colnames(Y_ordered)
    rownames(mu_est) <- rownames(Y_ordered)

    funs <- do.call(c, lapply(fit, "[[", "fun_g"))
    names(funs) <- rownames(Y_ordered)

    return(list(Y = Y_ordered,
                theta = theta_ordered,
                mu_est = mu_est,
                sigma_est = sigma_est,
                funs = funs))
}

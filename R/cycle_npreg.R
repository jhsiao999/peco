#' @title Predict Test-sample Ordering Using Training Labels
#'
#' @description Apply the estimates from \code{cycle_npreg_insample}
#'   to another gene expression dataset to infer an angle, or cell cycle
#'   phase, for each cell.
#'
#' @param Y_test A \code{SingleCellExperiment} object.
#'
#' @param sigma_est A vector of gene-specific standard errors of the
#'   cyclic trends.
#'
#' @param funs_est A vector of cyclic functions used to estimate
#'   cyclic trends.
#'
#' @param polyorder We estimate cyclic trends of gene expression levels using
#'   nonparamtric trend filtering.
#' 
#' @param normed Is the data already normalized? \code{TRUE} or
#'   \code{FALSE}.
#' 
#' @param method.grid The approach used to initialize angles in the
#'   computation. \code{method.grid = "uniform"} creates k
#'   equally-spaced bins ("grids"). \code{method.grid = "pca"} uses gene
#'   expression values to infer angles, and then these angles are used
#'   to project the cells to the closest bin.
#'
#' @param method.trend How to estimate cyclic trend of gene expression
#'   values. We offer three options: \code{method.trend =
#'   "trendfilter"}, uses \code{\link{fit_trendfilter}},
#'   \code{method.trend = "loess"} uses \code{\link{fit_loess}}), and
#'   \code{method.trend = "bsplines"} uses \code{\link{fit_bspline}}.  We
#'   found that \code{"trendfilter"} provided the best fits in our
#'   experiments. However, trend-filtering may require more
#'   computational effort since it uses cross-validation, so for fast
#'   results we recommend using \code{"bspline"}.
#'
#' @param grids number of bins to be selected along 0 to 2pi.
#'
#' @param get_trend_estimates Re-estimate the cylic trend based on the
#'   predicted cell cycle phase, or not (\code{TRUE} or
#'   \code{FALSE}). This step calls trendfilter and is
#'   computationally intensive.
#' 
#' @param ncores We use doParallel package for parallel computing.
#'
#' @return A list with the following elements:
#'
#' \item{Y}{The inputted gene expression marix.}
#' 
#' \item{cell_times_est}{Inferred angles or cell cycle phases
#'   (not ordered).}
#' 
#' \item{loglik_est}{Log-likelihoods for each gene.}
#' 
#' \item{cell_times_reordered}{The inferred angles, in ascending
#'   order.}
#' 
#' \item{Y_reorded}{The input gene expression matrix reordered by
#'   \code{cell_times_reordered}.}
#' 
#' \item{sigma_reordered}{Standard error of the cyclic trend for each
#'   gene, reordered by \code{cell_times_reordered}.}
#' 
#' \item{funs_reordered}{A list of functions for approximating the
#'   cyclic trends of gene express levels for each gene, reordered by
#'   cell_times_reordered.}
#' 
#' \item{mu_reordered}{Estimated cyclic trend of gene expression
#'   values for each gene, reordered by cell_times_reordered.}
#' 
#' \item{prob_per_cell_by_celltimes}{Probabilities of each cell belong
#'   to each bin.}
#'
#' @examples
#' library(SingleCellExperiment)
#' data(sce_top101genes)
#'
#' # Select top 5 cyclic genes.
#' sce_top5 <- sce_top101genes[order(rowData(sce_top101genes)$pve_fucci,
#'                                   decreasing = TRUE)[1:5],]
#'
#' # Select NA18511 samples for training and test sets.
#' coldata               <- colData(sce_top5)
#' which_samples_train   <- rownames(coldata)[coldata$chip_id != "NA18511"]
#' which_samples_predict <- rownames(coldata)[coldata$chip_id == "NA18511"]
#' 
#' # Learn gene cyclic functions from the training data using the
#' # trend filtering method.
#' sce_top5    <- data_transform_quantile(sce_top5)
#' expr_quant  <- assay(sce_top5, "cpm_quantNormed")
#' Y_train     <- expr_quant[, colnames(expr_quant) %in% which_samples_train]
#' theta_train <- coldata$theta_shifted[rownames(coldata) %in%
#'                  which_samples_train]
#' names(theta_train) <- rownames(coldata)[rownames(coldata) %in%
#'                                         which_samples_train]
#' model_5genes_train <- cycle_npreg_insample(Y = Y_train,theta = theta_train,
#'                                            polyorder = 2,ncores = 2,
#'                                            method.trend = "trendfilter")
#'
#' # Predict cell cycle in the test samples.
#' Y_test <- sce_top5[,colnames(sce_top5) %in% which_samples_predict]
#' model_5genes_predict <-
#'   cycle_npreg_outsample(Y_test = Y_test,
#'                         sigma_est = model_5genes_train$sigma_est,
#'                         funs_est = model_5genes_train$funs_est,
#'                         method.trend = "trendfilter",
#'                         ncores = 2,get_trend_estimates = FALSE)
#'
#' # Estimate cyclic gene expression levels given the cell cycle for each
#' # gene.
#' predict_cyclic <-
#'   fit_cyclical_many(Y = assay(model_5genes_predict$Y,"cpm_quantNormed"),
#'                     theta = colData(model_5genes_predict$Y)$cellcycle_peco)
#'
#' # Plot the cell cycle predictions vs. FUCCI phase.
#' par(mfrow = c(2,3),mar = c(4,4,3,1))
#' for (g in seq_along(rownames(model_5genes_predict$Y))) {
#'   plot(assay(model_5genes_predict$Y,"cpm_quantNormed")[
#'         rownames(model_5genes_predict$Y)[g],],
#'        x = colData(model_5genes_predict$Y)$cellcycle_peco,axes = FALSE,
#'        xlab = "FUCCI phase",ylab = "predicted phase")
#'   x <- seq(0,2*pi,length.out = 100)
#'   lines(x = x,y = predict_cyclic$cellcycle_function[[
#'                      rownames(model_5genes_predict$Y)[g]]](x),
#'         lwd = 1.5,col = "tomato")
#'   axis(2)
#'   axis(1,at = seq(0,2*pi,length.out = 5),
#'        labels = c(0,expression(pi/2),expression(pi),expression(3*pi/2),
#'                   expression(2*pi)))
#'   abline(h = 0,lwd = 1,col = "skyblue",lty = "dotted")
#'   title(names(model_5genes_predict$Y)[g])
#' }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom assertthat has_name
#'
#' @seealso
#'
#' \code{\link{cycle_npreg_insample}} for estimating parameters for
#' cyclic functions from training data,
#' \code{\link{cycle_npreg_loglik}} for log-likehood at angles between
#' 0 to 2pi, and \code{\link{cycle_npreg_mstep}} for estimating cyclic
#' functions given inferred phases from
#' \code{\link{cycle_npreg_loglik}}.
#' 
#' @export
#' 
cycle_npreg_outsample <- function(Y_test, sigma_est, funs_est,
                                  method.trend = c("trendfilter", "loess",
                                                   "bspline"),
                                  normed = TRUE, polyorder = 2,
                                  method.grid = "uniform", ncores = 2,
                                  grids = 100, get_trend_estimates = FALSE) {

    if (!normed & inherits(Y_test, "SingleCellExperiment")) {
        Y_test <- data_transform_quantile(Y_test)
        exprs_test <- assay(Y_test, "cpm_quantNormed")
    } else if (normed & inherits(Y_test, "SingleCellExperiment"))
      exprs_test <- assay(Y_test, "cpm_quantNormed")
    else if (!normed & !inherits(Y_test, "SingleCellExperiment"))
      exprs_test <- data_transform_quantile(Y_test)
    else
      exprs_test <- Y_test

    # Compute expected cell time for the test samples using mu and
    # sigma estimated from the training samples.
    initial_loglik <- cycle_npreg_loglik(Y = exprs_test,
                                         sigma_est = sigma_est,
                                         method.grid = method.grid,
                                         funs_est = funs_est,
                                         grids = grids)

    if (inherits(Y_test, "SingleCellExperiment")) {
      colData(Y_test)$cellcycle_peco <- initial_loglik$cell_times_est
      if (!is.null(colData(Y_test)@elementMetadata$labelDescription)) 
        colData(Y_test)@elementMetadata$labelDescription[
          ncol(colData(Y_test))] <-
            "peco predicted continuous cell cycle, between 0 to 2pi"
    }

    if (get_trend_estimates) {
      updated_estimates <-
        cycle_npreg_mstep(Y = Y_test,
                          theta = initial_loglik$cell_times_est,
                          method.trend = method.trend,
                          polyorder = polyorder,
                          ncores = ncores)
      out <- list(Y               = Y_test,
                  cell_times_est  = initial_loglik$cell_times_est,
                  loglik_est      = initial_loglik$loglik_est,
                  Y_reordered     = updated_estimates$Y,
                  cell_times_reordered = updated_estimates$theta,
                  mu_reordered    = updated_estimates$mu_est,
                  sigma_reordered = updated_estimates$sigma_est,
                  funs_reordered  = updated_estimates$funs,
                  prob_per_cell_by_celltimes =
                    initial_loglik$prob_per_cell_by_celltimes)
    } else
      out <- list(Y              = Y_test,
                  cell_times_est = initial_loglik$cell_times_est,
                  loglik_est     = initial_loglik$loglik_est,
                  prob_per_cell_by_celltimes =
                    initial_loglik$prob_per_cell_by_celltimes)
    return(out)
}

#' @title Obtain Cyclic Trend Estimates from Training Data
#'
#' @description Estimate cyclic trends of gene expression levels
#'   using training data.
#'
#' @param Y A matrix (gene by sample) of normalized and transformed
#'   gene expression values.
#'
#' @param theta A vector of angles.
#'
#' @param ncores We use the doParallel package for parallel computing.
#'
#' @param polyorder We estimate cyclic trends of gene expression
#'   levels using nonparamtric trend filtering.
#'
#' @param method.trend How to estimate cyclic trend of gene expression
#'   values. We offer three options: \code{method.trend =
#'   "trendfilter"}, uses \code{\link{fit_trendfilter}},
#'   \code{method.trend = "loess"} uses \code{\link{fit_loess}}), and
#'   \code{method.trend = "bsplines"} uses \code{\link{fit_bspline}}.  We
#'   found that \code{"trendfilter"} provided the best fits in our
#'   experiments. However, trend-filtering may require more
#'   computational effort since it uses cross-validation, so for fast
#'   results we recommend using \code{"bspline"}.
#'
#' @return A list with four elements:
#'
#' \item{Y}{The gene expression marix.}
#'
#' \item{theta}{Vector of angles or cell cycle phases.}
#'
#' \item{sigma_est}{Estimated standard error of the cyclic trend for
#'   each gene.}
#'
#' \item{funs_est}{A list of functions for approximating the cyclic
#' trends of gene express levels for each gene.}
#'
#' @author Joyce Hsiao
#'
#' @examples
#' # See \code{help(cycle_npreg_outsample)} for an example.
#' 
#' @seealso
#'   \code{\link{cycle_npreg_mstep}} for estimating cyclic functions
#'   given inferred phases from \code{\link{cycle_npreg_loglik}}, and
#'   \code{\link{cycle_npreg_outsample}} for predicting cell cycle phase
#'   using parameters learned from \code{\link{cycle_npreg_insample}}.
#'
#' @export
#' 
cycle_npreg_insample <- function(Y, theta, ncores = 2, polyorder = 2,
                                 method.trend = c("trendfilter", "loess",
                                                  "bspline")) {

    # Order data by initial cell times.
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

    return(list(Y         = Y,
                theta     = theta,
                sigma_est = initial_mstep$sigma_est,
                funs_est  = initial_mstep$funs))
}

#' @title Initialize Grid for cycle_npreg_loglik
#' 
#' @description For prediction, initialize grid points for cell cycle phase
#'   on a circle.
#'
#' @param Y Gene expression matrix (gene by sample).
#'
#' @param grids Number of bins used over interval 0 to \code{2*pi}.
#'
#' @param method.grid The approach used to initialize angles in the
#'   computation. \code{method.grid = "uniform"} creates k
#'   equally-spaced bins ("grids"). \code{method.grid = "pca"} uses gene
#'   expression values to infer angles, and then these angles are used
#'   to project the cells to the closest bin.
#'
#' @return A vector of initialized angles to be used in
#'   \code{\link{cycle_npreg_loglik}}.
#'
#' @keywords internal
#'
initialize_grids <- function(Y, grids = 100, method.grid = c("pca","uniform")) {

    len <- (2*pi)/(2*grids)
    theta_grids <- seq(len, (2*pi)-(len), length.out = grids)

    if (method.grid == "pca") {
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
                which.min(pmin(abs(theta_grids - grid_approx[i]),
                               abs(theta_grids - (2*pi - grid_approx[i]))))
            theta_initial <- theta_grids[theta_initial_ind]
        }
    }

    if (method.grid == "uniform")
        theta_initial <- theta_grids
    return(theta_initial)
}

#' @title Infer Angles or Cell Cycle Phase Using Gene Expression Data
#'
#' @param Y Gene by sample expression matrix.
#' 
#' @param sigma_est A vector of standard errors for each gene from the
#'   training data.
#' 
#' @param funs_est A vector of cyclic functions estimated for each
#'   gene from the training data.
#'
#' @param grids Number of bins used over interval 0 to \code{2*pi}.
#'
#' @param method.grid The approach used to initialize angles in the
#'   computation. \code{method.grid = "uniform"} creates k
#'   equally-spaced bins ("grids"). \code{method.grid = "pca"} uses gene
#'   expression values to infer angles, and then these angles are used
#'   to project the cells to the closest bin.
#'
#' @return A list with the following three elements:
#'
#' \item{cell_times_est}{Inferred angles or cell cycle phases (not
#'   ordered).}
#' 
#' \item{loglik_est}{Log-likelihood estimates for each gene.}
#' 
#' \item{prob_per_cell_by_celltimes}{Probabilities of each cell belonging
#'   to each bin.}
#'
#' @importFrom stats dnorm
#'
#' @keywords internal
#' 
cycle_npreg_loglik <- function(Y, sigma_est, funs_est, grids = 100,
                               method.grid = c("pca","uniform")) {

    N <- ncol(Y)
    G <- nrow(Y)
    theta_choose <- initialize_grids(Y, grids = grids,
                                     method.grid = method.grid)

    # For each cell, sum up the loglikelihood for each gene
    # at the observed cell times.
    loglik_per_cell_by_celltimes <- matrix(0, N, grids)
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    for (n in seq_len(N)) {
        loglik_per_cell <- do.call(rbind, lapply(seq_len(G), function(g) 
            dnorm(Y[g,n], funs_est[[g]](theta_choose), sigma_est[g], log = TRUE)
        ))
        loglik_per_cell <- colSums(loglik_per_cell)
        loglik_per_cell_by_celltimes[n,] <- loglik_per_cell
    }

    # Use maximum-likelihood to assign samples.
    prob_per_cell_by_celltimes <-
        t(apply(loglik_per_cell_by_celltimes, 1,
                function(x) {
                  sumll <- sum(exp(x), na.rm=TRUE)
                  if (sumll == 0)
                    return(rep(0, grids))
                  else
                    return(exp(x)/sumll)
                 }))
    colnames(prob_per_cell_by_celltimes) <- theta_choose

    cell_times_samp_ind <-
        apply(prob_per_cell_by_celltimes, 1, function(x)
          if (max(x, na.rm = TRUE) == 0)
            sample(seq_len(grids), 1, replace = FALSE)
          else
            which.max(x))
    names(cell_times_samp_ind) <-
        colnames(prob_per_cell_by_celltimes)[cell_times_samp_ind]

    cell_times_est <-
        do.call(c, lapply(seq_len(N), function(n)
          theta_choose[cell_times_samp_ind[n]]))
    names(cell_times_est) <- colnames(Y)

    # Compute likelihood based on the selected cell times.
    loglik_max_per_cell <-
        do.call(c, lapply(seq_len(N), function(n) {
          ll <- loglik_per_cell_by_celltimes[n,]
          return(ll[cell_times_samp_ind[n]])
        }))
    loglik_est <- sum(loglik_max_per_cell)

    return(list(loglik_est = loglik_est,
                cell_times_est = cell_times_est,
                prob_per_cell_by_celltimes = prob_per_cell_by_celltimes))
}

#' @title Estimate Parameters of the Cyclic Trends
#'
#' @description This is used by \code{\link{cycle_npreg_insample}}
#'   (model fitting from training data) and
#'   \code{\link{cycle_npreg_outsample}} (prediction in test data
#'   estimate cyclic trends of gene expression values. The function
#'   outputs, for each gene, a standard error of the cyclic trend, a
#'   cyclic function, and the estimated expression levels from the
#'   cyclic function.
#'
#' @param Y Gene by sample expression matrix (log2CPM).
#' 
#' @param theta Observed cell times.
#' 
#' @param method.trend How to estimate cyclic trend of gene expression
#'   values. We offer three options: \code{method.trend =
#'   "trendfilter"}, uses \code{\link{fit_trendfilter}},
#'   \code{method.trend = "loess"} uses \code{\link{fit_loess}}), and
#'   \code{method.trend = "bsplines"} uses \code{\link{fit_bspline}}.  We
#'   found that \code{"trendfilter"} provided the best fits in our
#'   experiments. However, trend-filtering may require more
#'   computational effort since it uses cross-validation, so for fast
#'   results we recommend using \code{"bspline"}.
#' 
#' @param polyorder We estimate cyclic trends of gene expression
#'   levels using nonparamtric trend filtering.
#'
#' @param ncores How many computing cores to use? We use the
#'   \code{doParallel} package for parallel computing.
#'
#' @return A list with the following elements:
#'
#' \item{Y}{Input gene expression data.}
#' 
#' \item{theta}{Input angles.}
#' 
#' \item{mu_est}{Estimated expression levels given the cyclic function
#'       for each gene.}
#' 
#' \item{sigma_est}{Estimated standard error of the cyclic trends for
#'   each gene.}
#' 
#' \item{funs}{Estimated cyclic functions.}
#'
#' @importFrom stats approxfun
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @keywords internal
#' 
cycle_npreg_mstep <- function(Y, theta,
                              method.trend = c("trendfilter","loess","bspline"),
                              polyorder = 2, ncores = 2) {

    if (inherits(Y, "SingleCellExperiment"))
      exprs_test <- assay(Y, "cpm_quantNormed")
    else
      exprs_test <- Y

    if (is.null(ncores))
      ncores <- 2
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    message(paste("computing on",ncores,"cores"))

    G <- nrow(exprs_test)
    N <- ncol(exprs_test)

    exprs_test_ordered <- exprs_test[,names(theta)]
    ord <- order(theta)
    theta_ordered <- theta[ord]
    exprs_test_ordered <- exprs_test_ordered[,ord]

    fit <- foreach(g = seq_len(G)) %dopar% {
        y_g <- exprs_test_ordered[g,]
        if (method.trend == "trendfilter") {
            fit_g <- fit_trendfilter(yy = y_g, polyorder = polyorder)
            fun_g <- approxfun(x = as.numeric(theta_ordered),
                               y = as.numeric(fit_g$trend.yy),
                               rule = 2)
            mu_g  <- fit_g$trend.yy
        }
        if (method.trend == "bspline") {
            fit_g <- fit_bspline(yy = y_g, time = theta_ordered)
            fun_g <- approxfun(x = as.numeric(theta_ordered),
                               y = as.numeric(fit_g$pred.yy),
                               rule = 2)
            mu_g  <- fit_g$pred.yy
        }

        if (method.trend == "loess") {
            fit_g <- fit_loess(yy = y_g, time = theta_ordered)
            fun_g <- approxfun(x = as.numeric(theta_ordered),
                               y = as.numeric(fit_g$pred.yy),
                               rule = 2)
            mu_g  <- fit_g$pred.yy
        }

        sigma_g <- sqrt(sum((y_g - mu_g)^2)/N)
        return(list(y_g     = y_g,
                    mu_g    = mu_g,
                    sigma_g = sigma_g,
                    fun_g   = fun_g))
    }
    stopCluster(cl)

    sigma_est <- do.call(c, lapply(fit, "[[", "sigma_g"))
    names(sigma_est) <- rownames(exprs_test_ordered)

    mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
    colnames(mu_est) <- colnames(exprs_test_ordered)
    rownames(mu_est) <- rownames(exprs_test_ordered)

    funs <- do.call(c, lapply(fit, "[[", "fun_g"))
    names(funs) <- rownames(exprs_test_ordered)

    return(list(Y         = exprs_test_ordered,
                theta     = theta_ordered,
                mu_est    = mu_est,
                sigma_est = sigma_est,
                funs      = funs))
}

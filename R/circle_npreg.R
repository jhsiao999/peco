# cycle_npreg_insample and cycle_npreg_outsample are the two main
# functions in peco. cycle_npreg_insample generates cyclic trend
# estimates of gene expression levels using training data, and
# cycle_npreg_outsample applies the estimates of cycle_npreg_insample
# to another gene expression dataset to infer an angle or cell cycle
# phase for each cell.

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
#' @param ncores We use mclapply function for parallel computing.
#' 
#' @param polyorder We estimate cyclic trends of gene expression
#' levels using nonparamtric trend filtering. The default fits second
#' degree polynomials.
#' 
#' @param method.trend Varous methods that can be applied to estimate
#' cyclic trend of gene expression levels.
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
#' @author Joyce Hsiao
#'
#' @export
#' 
cycle_npreg_insample <- function(Y, theta,
                                 ncores=4,
                                 polyorder=2,
                                 method.trend=c("trendfilter",
                                                "loess", "bspline"),...) {

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

#' @title Predict test-sample ordering using training labels (no update)
#'
#' @description Apply the estimates of cycle_npreg_insample to another
#' gene expression dataset to infer an angle or cell cycle phase for
#' each cell.
#'
#' @param Y_test Gene expression data to be used for prediction.
#'     Gene by sample.
#' 
#' @param sigma_est A vector of gene-specific standard error of the
#' cyclic trends.
#' 
#' @param funs_est A vector of cyclic functions estimating cyclic trends.
#' 
#' @param polyorder We estimate cyclic trends of gene expression levels using
#'    nonparamtric trend filtering. The default fits second degree polynomials.
#' 
#' @param method.grid Method for defining bins along the circle.
#' 
#' @param method.trend Varous methods that can be applied to estimate
#' cyclic trend of gene expression levels.
#' 
#' @param ncores We use mclapply function for parallel computing.
#'
#' @return A list with the following elements:
#' 
#' \item{Y}{The input gene expression marix.}
#' 
#' \item{cell_times_est}{Inferred angles or cell cycle phases, NOT
#' ordered.}
#' 
#' \item{loglik_est}{Log-likelihood estimates for each gene.}
#' 
#' \item{cell_times_reordered}{The inferred angles reordered (in
#' ascending order).}
#' 
#' \item{Y_reorded}{The input gene expression matrix reordered by
#' cell_times_reordered.}
#' 
#' \item{sigma_reordered}{Estimated standard error of the cyclic trend
#' for each gene, reordered by cell_times_reordered.}
#' 
#' \item{funs_reordered}{A list of functions for approximating the
#' cyclic trends of gene express levels for each gene, reordered by
#' cell_times_reordered.}
#' 
#' \item{mu_reordered}{Estimated cyclic trend of gene expression
#' values for each gene, reordered by cell_times_reordered.}
#' 
#' \item{prob_per_cell_by_celltimes}{Probabilities of each cell belong
#' to each bin.}
#'
#' @export
#' 
cycle_npreg_outsample <- function(Y_test,
                                  sigma_est,
                                  funs_est,
                                  method.trend=c("trendfilter",
                                                 "loess", "bspline"),
                                  polyorder=2,
                                  method.grid=c("pca", "uniform"),
                                  ncores=4,...) {

  # compute expected cell time for the test samples
  # under mu and sigma estimated from the training samples
  initial_loglik <- cycle_npreg_loglik(Y = Y_test,
                                       sigma_est = sigma_est,
                                       method.type="supervised",
                                       method.grid=method.grid,
                                       funs_est=funs_est)
  updated_estimates <- cycle_npreg_mstep(Y = Y_test,
                                         theta = initial_loglik$cell_times_est,
                                         method.trend = method.trend,
                                         polyorder=polyorder,
                                         ncores = ncores)

  out <- list(Y=Y_test,
              cell_times_est=initial_loglik$cell_times_est,
              loglik_est=initial_loglik$loglik_est,
              Y_reordered=updated_estimates$Y,
              cell_times_reordered=updated_estimates$theta,
              mu_reordered=updated_estimates$mu_est,
              sigma_reordered=updated_estimates$sigma_est,
              funs_reordered=updated_estimates$funs,
              prob_per_cell_by_celltimes=initial_loglik$prob_per_cell_by_celltimes)
  return(out)
}

#------ Supporting functions

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
#' @author Joyce Hsiao
#'
#' @importFrom stats prcomp
#' @importFrom circular coord2rad
#' 
#' @export
#' 
initialize_grids <- function(Y, grids=100,
                             method.grid=c("pca", "uniform"),...) {
    
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

    for (i in 1:length(grid_approx)) {
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

#' @title Infer angles or cell cycle phase based on gene expression data
#'
#' @param Y Gene by sample expression matrix.
#' 
#' @param sigma_est A vector of standard errors for each gene from the
#' training data.
#' 
#' @param funs_est A vector of cyclic functions estimated for each
#' gene from the training data.
#'
#' @return A list with the following three elements:
#' 
#' \item{cell_times_est}{Inferred angles or cell cycle phases, NOT
#' ordered.}
#' 
#' \item{loglik_est}{Log-likelihood estimates for each gene.}
#' 
#' \item{prob_per_cell_by_celltimes}{Probabilities of each cell belong
#' to each bin.}
#'
#' @author Joyce Hsiao
#'
#' @importFrom stats dnorm
#' 
#' @export
#' 
cycle_npreg_loglik <- function(Y, sigma_est, funs_est,
                               grids=100,
                               method.type=c("supervised", "unsupervised"),
                               method.grid=c("pca", "uniform"),...) {
# import circular
    
  N <- ncol(Y)
  G <- nrow(Y)

  if (method.type=="unsupervised") {
    theta_choose <- initialize_grids(Y, grids=grids, method.grid="pca")
    loglik_per_cell_by_celltimes <- matrix(0, N, length(theta_choose))
    prob_per_cell_by_celltimes <- matrix(0, N, length(theta_choose))
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    colnames(prob_per_cell_by_celltimes) <- theta_choose
  }
  if (method.type=="supervised") {
    theta_choose <- initialize_grids(Y, grids=grids, method.grid="uniform")
    loglik_per_cell_by_celltimes <- matrix(0, N, grids)
    prob_per_cell_by_celltimes <- matrix(0, N, grids)
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    colnames(prob_per_cell_by_celltimes) <- theta_choose
  }

  for (n in 1:N) {
      
    # for each cell, sum up the loglikelihood for each gene
    # at the observed cell times
    loglik_per_cell <- do.call(rbind, lapply(1:G, function(g) {
      dnorm(Y[g,n], funs_est[[g]](theta_choose), sigma_est[g], log = TRUE)
    }))
    loglik_per_cell <- colSums(loglik_per_cell)

    loglik_per_cell_by_celltimes[n,] <- loglik_per_cell
  }

  # use max likelihood to assign samples
  for (n in 1:N) {
    sumll <- sum(exp(loglik_per_cell_by_celltimes)[n,], na.rm=T)
    if (sumll == 0) {
      prob_per_cell_by_celltimes[n,] <- rep(0, grids)
    } else {
      prob_per_cell_by_celltimes[n,] <- exp(loglik_per_cell_by_celltimes)[n,]/sumll
    }
  }
  cell_times_samp_ind <- sapply(1:N, function(n) {
    if (max(prob_per_cell_by_celltimes[n,], na.rm=T)==0) {
      sample(1:grids, 1, replace=F)
    } else {
      which.max(prob_per_cell_by_celltimes[n,])
    }
  })
  cell_times_est <- sapply(1:N, function(n) {
    theta_choose[cell_times_samp_ind[n]]
  })
  names(cell_times_est) <- colnames(Y)

  # compute likelihood based on the selected cell times
  loglik_max_per_cell <- sapply(1:N, function(n) {
    ll <- loglik_per_cell_by_celltimes[n,]
    ll[cell_times_samp_ind[n]]
  })
  loglik_est <- sum(loglik_max_per_cell)

  return(list(loglik_est=loglik_est,
              cell_times_est=cell_times_est,
              prob_per_cell_by_celltimes=prob_per_cell_by_celltimes))
}

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
#' 
#' @param theta Observed cell times.
#'
#' @return A list with the following elements: 
#'     
#' \item{Y}{Input gene expression data.}
#' 
#' \item{theta}{Input angles.}
#' 
#' \item{mu_est}{Estimated expression levels given the cyclic function
#' for each gene.}
#' 
#' \item{sigma_est}{Estimated standard error of the cyclic trends for
#' each gene.}
#' 
#' \item{funs_est}{Estimated cyclic functions.}
#'
#' @author Joyce Hsiao
#'
#' @importFrom assertthat assert_that
#' @importFrom parallel mclapply
#' @importFrom stats approxfun
#' 
#' @export
#' 
cycle_npreg_mstep <- function(Y, theta, method.trend=c("trendfilter",
                                                       "loess", "bspline"),
                              polyorder=2,
                              ncores=4,...) {

      G <- nrow(Y)
      N <- ncol(Y)

      if (!assert_that(all.equal(names(theta), colnames(Y)))) {
        Y_ordered <- Y[,match(names(theta), colnames(Y))]
        ord <- order(theta)
        theta_ordered <- theta[ord]
        Y_ordered <- Y_ordered[,ord]
      } else {
        ord <- order(theta)
        theta_ordered <- theta[ord]
        Y_ordered <- Y[,ord]
      }

      # for each gene, estimate the cyclical pattern of gene expression
      # conditioned on the given cell times
      fit <- mclapply(1:G, function(g) {
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
      }, mc.cores = ncores)

      sigma_est <- sapply(fit, "[[", "sigma_g")
      names(sigma_est) <- rownames(Y_ordered)

      mu_est <- do.call(rbind, lapply(fit, "[[", "mu_g"))
      colnames(mu_est) <- colnames(Y_ordered)
      rownames(mu_est) <- rownames(Y_ordered)

      funs <- sapply(fit, "[[", "fun_g")
      names(funs) <- rownames(Y_ordered)

      return(list(Y = Y_ordered,
                  theta = theta_ordered,
                  mu_est = mu_est,
                  sigma_est = sigma_est,
                  funs = funs))
    }

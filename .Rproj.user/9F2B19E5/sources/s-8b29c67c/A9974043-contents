library(parallel)

#' @title Initialize grid points for estimation
#'
#' @export
initialize_grids <- function(Y, grids=100,
                             method.grid=c("pca", "uniform"), ...) {

  len <- (2*pi)/(2*grids)
  theta_grids <- seq(len, (2*pi)-(len), length.out=grids)

  library(circular)
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
      theta_initial_ind[i] <- which.min(pmin(abs(theta_grids-grid_approx[i]),
                                             abs(theta_grids-(2*pi-grid_approx[i]))))
      theta_initial <- theta_grids[theta_initial_ind]
    }
  }

  if (method.grid=="uniform") {
    theta_initial <- theta_grids
  }

  return(theta_initial)
}





#' @title Estimate parameters for the cyclial ordering using nonparametric smoothing
#'
#' @description Conditioned on observed cell times, estimate the cyclical trend
#'   for each gene, and compute the likelihood of the observed cell times.
#'
#' @param Y gene by sample expression matrix (log2CPM).
#' @param theta observed cellt times
cycle.npreg.mstep <- function(Y, theta, method.trend=c("trendfilter",
                                                       "npcirc.nw",
                                                       "npcirc.ll",
                                                       "loess", "bspline"),
                              polyorder=3,
                              ncores=12, ...) {
  #      library(NPCirc)
  library(genlasso)
  library(assertthat)

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
    #    print(g)
    y_g <- Y_ordered[g,]

    if (method.trend=="npcirc.nw") {
      fit_g <- kern.reg.circ.lin(theta_ordered, y_g, method = "NW",
                                 t=theta_ordered)
      fun_g <- approxfun(x=as.numeric(fit_g$x), y=fit_g$y, rule=2)
      mu_g <- fun_g(theta_ordered)
    }
    if (method.trend=="npcirc.ll") {
      fit_g <- kern.reg.circ.lin(theta_ordered, y_g, method = "LL",
                                 t=theta_ordered)
      fun_g <- approxfun(x=as.numeric(fit_g$x), y=fit_g$y, rule=2)
      mu_g <- fun_g(theta_ordered)
    }
    if (method.trend=="trendfilter") {
      fit_g <- fit.trendfilter.generic(yy=y_g, polyorder = polyorder)
      fun_g <- approxfun(x=as.numeric(theta_ordered),
                         y=as.numeric(fit_g$trend.yy), rule=2)
      mu_g <- fit_g$trend.yy
    }
    if (method.trend=="bspline") {
      fit_g <- fit.bspline(yy=y_g, time = theta_ordered)
      fun_g <- approxfun(x=as.numeric(theta_ordered),
                         y=as.numeric(fit_g$pred.yy), rule=2)
      mu_g <- fit_g$pred.yy
    }

    if (method.trend=="loess") {
      fit_g <- fit.loess(yy=y_g, time = theta_ordered)
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



#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
cycle.npreg.loglik <- function(Y, sigma_est, funs_est,
                               grids=100,
                               method.type=c("supervised", "unsupervised"),
                               method.grid=c("pca", "uniform"), ...) {

  N <- ncol(Y)
  G <- nrow(Y)

  if (method.type=="unsupervised") {
    theta_choose <- initialize_grids(Y, grids=grids, method.grid="pca")
    loglik_per_cell_by_celltimes <- matrix(0, N, length(theta_choose))
    prob_per_cell_by_celltimes <- matrix(0, N, length(theta_choose))
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    colnames(prob_per_cell_by_celltimes) <- theta_choose
    rownames(loglik_per_cell_by_celltimes) <- colnames(Y)
    rownames(prob_per_cell_by_celltimes) <- colnames(Y)
  }
  if (method.type=="supervised") {
    theta_choose <- initialize_grids(Y, grids=grids, method.grid="uniform")
    loglik_per_cell_by_celltimes <- matrix(0, N, grids)
    prob_per_cell_by_celltimes <- matrix(0, N, grids)
    colnames(loglik_per_cell_by_celltimes) <- theta_choose
    colnames(prob_per_cell_by_celltimes) <- theta_choose
    rownames(loglik_per_cell_by_celltimes) <- colnames(Y)
    rownames(prob_per_cell_by_celltimes) <- colnames(Y)
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
    # print(n)
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
              loglik_per_cell_by_celltimes=loglik_per_cell_by_celltimes,
              prob_per_cell_by_celltimes=prob_per_cell_by_celltimes))
}




#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
#' @title log-likelihood of nonparametric smoothing
#'
#' @param Y gene by sample expression matrix
#' @param mu_est gene by sample matrix of expected mean
#' @param sigma_est vector of standard errors for each gene
#'
#' @export
cycle.npreg.loglik.post <- function(loglik,
                                    pi_bins, ...) {


  #theta_choose <- initialize_grids(Y, grids=grids, method.grid="uniform")
  post_loglik_per_cell_by_celltime <- matrix(0, nrow=nrow(loglik), ncol=ncol(loglik))
  post_prob_per_cell_by_celltime <- matrix(0, nrow=nrow(loglik), ncol=ncol(loglik))
  theta_choose <- as.numeric(colnames(loglik))
  cell_times_est <- matrix(0,nrow=nrow(loglik),
                           dimnames=list(rownames(loglik)))
  names(cell_times_est) <- rownames(loglik)

  post_loglik_est <- 0
  for (n in 1:ncol(loglik)) {
    post_loglik <- loglik[n,] + log(pi_bins)
    post_loglik_per_cell_by_celltime[n,] <- post_loglik

    post_prob <- exp(post_loglik)/sum(exp(post_loglik))
    which_grid <- which.max(post_prob)
    cell_times_est[n] <- theta_choose[which_grid]
    post_loglik_est <- post_loglik_est + sum(post_loglik)
  }

  return(list(loglik_est=post_loglik_est,
              cell_times_est=cell_times_est,
              loglik_per_cell_by_celltimes=post_loglik_per_cell_by_celltime,
              prob_per_cell_by_celltimes=post_prob_per_cell_by_celltime))
}




#' @title Prior probabilty distribution for cell times
#'
#' @export
cycle.npreg.prior <- function(theta, grids=100, plot.it=F,...) {

  breaks <- seq(0,2*pi, by=2*pi/grids)
  bins <- cut(theta, breaks=breaks, include.lowest = T)
  n_eachbin <- table(bins)
  pi_bins <- (n_eachbin+1)/(length(theta)+1)
  pi_bins <- (pi_bins)/sum(pi_bins)

  if (plot.it==TRUE) {
    par(mfrow=c(1,2))
    hist(theta, nclass=100, main = "Theta observed",
         xlab="theta")
    plot(pi_bins, xlab="Theta bins", ylab="Prior probability")
  }
  return(pi_bins)
}




#' @title Estimate cell cycle ordering in the current sample
#'
#' @param update T/F to update cell times
#'
#' @export
cycle.npreg.insample <- function(Y, theta,
                                 ncores=12,
                                 polyorder=3,
                                 method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter",
                                                "loess", "bspline"),
                                 ...) {

  # order data by initial cell times
  G <- nrow(Y)
  N <- ncol(Y)
  theta_ordered_initial <- theta[order(theta)]
  Y_ordered <- Y[,order(theta)]

  # initialize mu and sigma
  initial_mstep <- cycle.npreg.mstep(Y = Y_ordered,
                                     theta = theta_ordered_initial,
                                     polyorder=polyorder,
                                     method.trend=method.trend,
                                     ncores = ncores)

  out <- list(Y=Y,
              theta=theta,
              sigma_est=initial_mstep$sigma_est,
              funs_est=initial_mstep$funs)
  return(out)
}





#' @title Predict test-sample ordering using training lables (no update)
#'
#' @export
cycle.npreg.outsample <- function(Y_test,
                                  sigma_est,
                                  funs_est,
                                  theta_prior=NULL,
                                  method.trend=c("npcirc.nw", "npcirc.ll", "trendfilter",
                                                 "loess", "bspline"),
                                  polyorder=3,
                                  method.grid=c("pca", "uniform"),
                                  ncores=12,...) {


  # compute expected cell time for the test samples
  # under mu and sigma estimated from the training samples
  initial_loglik <- cycle.npreg.loglik(Y = Y_test,
                                       sigma_est = sigma_est,
                                       method.type="supervised",
                                       method.grid=method.grid,
                                       funs_est=funs_est)

  # pi_bins <- cycle.npreg.prior(theta_prior, grids=100)
  # initial_postlik <- cycle.npreg.loglik.post(loglik=initial_loglik$loglik_per_cell_by_celltimes,
  #                                            pi_bins = pi_bins)

  updated_estimates <- cycle.npreg.mstep(Y = Y_test,
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
              #loglik_per_cell_by_celltimes=initial_loglik$loglik_per_cell_by_celltimes,
              prob_per_cell_by_celltimes=initial_loglik$prob_per_cell_by_celltimes)
  return(out)
}







dir <-"/project2/gilad/joycehsiao/fucci-seq"

source(file.path(dir,"/peco/R/fit.cyclical.R"))
source(file.path(dir, "/peco/R/cycle.npreg.R"))
source(file.path(dir, "/peco/R/run_seurat.R"))
source(file.path(dir, "/peco/R/run_seurat.R"))
source(file.path(dir, "/peco/R/unsupervised.R"))
source(file.path(dir, "/peco/R/cycle.corr.R"))


#' @title  run all methods for withheld samples
#'
#' @export
run_methods <- function(Y_test, Y_test_normed,
                             theta_test, training_model,
                             seurat.genes,
                             pdata_test,
                             fdata,
                             ncores=15, maxiter=30, ...) {

  cycle.genes <- rownames(training_model$Y)
  Y_test_normed.cycle <- Y_test_normed[which(rownames(Y_test_normed) %in% cycle.genes),]
  Y_test.cycle <- Y_test[which(rownames(Y_test) %in% cycle.genes),]

  ### supervised methods
  message("Begin PCA method...")
  library(circular)
  prcomp_raw <- prcomp(t(Y_test.cycle), scale.=TRUE)
  coord_raw <- coord2rad(prcomp_raw$x[,1:2])
  # prcomp_normed <- prcomp(t(Y_test_normed.cycle), scale.=TRUE)
  # coord_normed <- coord2rad(prcomp_normed$x[,1:2])
  fit.pca <- list(cell_times_est=as.numeric(coord_raw))
  names(fit.pca$cell_times_est) <- colnames(Y_test.cycle)

  ### supervised methods
  message("Begin supervised method...")
  fit.supervised <- cycle.npreg.outsample(Y_test=Y_test_normed.cycle,
                                          sigma_est=training_model$sigma_est,
                                          funs_est=training_model$funs_est,
                                          theta_prior=training_model$theta,
                                          method.grid = "uniform",
                                          method.trend="trendfilter",
                                          polyorder=2,
                                          ncores=ncores)

  ### unsuperviesd methods
  theta_initial <- initialize_grids(Y_test_normed.cycle, method.grid="pca")
  names(theta_initial) <- colnames(Y_test_normed.cycle)

  # message("Begin unsupervised bspline...")
  # fit.bspline.unsup <- cycle.npreg.unsupervised(Y=Y_test_normed.cycle, theta=theta_initial,
  #                                               ncores=ncores,
  #                                               method.trend="bspline",
  #                                               maxiter=maxiter, verbose=TRUE, tol=1)
  #
  # message("Begin unsupervised loess...")
  # fit.loess.unsup <- cycle.npreg.unsupervised(Y=Y_test_normed.cycle, theta=theta_initial,
  #                                             ncores=ncores,
  #                                             method.trend="loess",
  #                                             maxiter=maxiter, verbose=TRUE, tol=1)

  # message("Begin unsupervised trendfilter...")
  # fit.trend2.unsup <- cycle.npreg.unsupervised(Y=Y_test_normed.cycle, theta=theta_initial,
  #                                              ncores=ncores,
  #                                              method.trend="trendfilter",
  #                                              polyorder=2,
  #                                              maxiter=maxiter, verbose=TRUE, tol=1)

  message("Begin Seurat...")
  # replace Y_test rownames with symbols
  symbols <- fdata$name[match(rownames(Y_test), rownames(fdata))]
  Y_test_seurat <- Y_test
  rownames(Y_test_seurat) <- symbols

  fit.seurat <- run_seurat(Y=Y_test_seurat,
                           s.genes=seurat.genes$s.genes,
                           g2m.genes=seurat.genes$g2m.genes,
                           n.bin=25,
                           seed.use=1, random.seed=1)

  fit.seurat <- as.list(fit.seurat)

  seurat.pca <- prcomp(cbind(fit.seurat$G2M, fit.seurat$S), scale=TRUE)
  seurat.cell_times_est <- as.numeric(coord2rad(cbind(seurat.pca$x[,1],seurat.pca$x[,2])))
  names(seurat.cell_times_est) <- colnames(Y_test_seurat)
  fit.seurat$cell_times_est <- seurat.cell_times_est


  message("Tidying up results...")

  out <- list(fit.supervised=fit.supervised,
              #fit.trend2.unsup=fit.trend2.unsup,
              fit.pca=fit.pca,
              #fit.bspline.unsup=fit.bspline.unsup,
              #fit.loess.unsup=fit.loess.unsup,
              fit.seurat=fit.seurat)

  # set.seed(111)
  # theta_test_null <- sample(theta_test)

  for (i in 1:length(out)) {
    print(i)
    out[[i]]$ref_time <- theta_test
    out[[i]]$pred_time <- with(out[[i]], cell_times_est[match(names(ref_time),
                                                        names(cell_times_est))])
    out[[i]]$pred_time_shift <- with(out[[i]], rotation(ref_time, pred_time))
    out[[i]]$diff_time <- with(out[[i]], circ_dist(pred_time_shift, ref_time))
    out[[i]]$dapi <- pdata_test$dapi.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
    out[[i]]$gfp <- pdata_test$gfp.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
    out[[i]]$rfp <- pdata_test$rfp.median.log10sum.adjust[
      match(names(theta_test), rownames(pdata_test))]
  }

  return(out)
}



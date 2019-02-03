
#' @title Perform Principal Components Analysis (PCA).
#
#' @param x gene-by-sample matrix
#' @param retx center, scale - see ?prcomp
#
#' @return a list with the following elements:
#   PCs - sample-by-PC matrix of principal components
#   explained - proportion of variance explained by each PC
#
#' @references Zhang et al. 2009 (http://www.ncbi.nlm.nih.gov/pubmed/19763933)
#'
#' @export
run_pca <- function(x, retx = TRUE, center = TRUE, scale = TRUE) {
  library("testit")

  pca <- prcomp(t(x), retx = retx, center = center, scale. = scale)
  variances <- pca$sdev^2
  explained <- variances / sum(variances)
  assert("Variance explained is calculated correctly.",
         explained[1:2] - summary(pca)$importance[2, 1:2] < 0.0001)
  return(list(PCs = pca$x, explained = explained))
}

#' @title Plot PCA results.
#
#' @param x numeric matrix of PCs
#' @param pcx PC to plot on x-axis (default: 1)
#' @param pcy PC to plot on y-axis (default: 2)
#' @param explained numeric vector of fractions of variance explained by each PC
#' @param metadata data frame or matrix that contains the metadata used to annotate
#             the plot
#' @param color, shape, size: column name of metadata used to pass column to ggplot
#                       aesthetic
#' @param factors character vector which contains the column names of metadata that
#            need to be explicitly converted to a factor
#' @param ...  Additional arguments passed to geom_point
#'
#' @export
plot_pca <- function(x, pcx = 1, pcy = 2, explained = NULL, metadata = NULL,
                     color = NULL, shape = NULL, factors = NULL,
                     ...) {
  library("ggplot2")
  library("testit")

  # Prepare data frame to pass to ggplot
  if (!is.null(metadata)) {
    assert("PC and metadata have same number of rows.",
           nrow(x) == nrow(metadata))
    plot_data <- cbind(x, metadata)
    plot_data <- as.data.frame(plot_data)
    # Convert numeric factors to class "factor"
    for (f in factors) {
      plot_data[, f] <- as.factor(plot_data[, f])
    }
  } else {
    plot_data <- as.data.frame(x)
  }
  # Prepare axis labels
  if (!is.null(explained)) {
    assert("Number of PCs differs between x and explained.",
           length(explained) == ncol(x))
    xaxis <- sprintf("PC%d (%.2f%%)", pcx, round(explained[pcx] * 100, 2))
    yaxis <- sprintf("PC%d (%.2f%%)", pcy, round(explained[pcy] * 100, 2))
  } else {
    xaxis <- paste0("PC", pcx)
    yaxis <- paste0("PC", pcy)
  }
  # Plot
  p <- ggplot(plot_data, aes_string(x = paste0("PC", pcx),
                                    y = paste0("PC", pcy),
                                    color = color,
                                    shape = shape)) +
    geom_point(...) +
    labs(x = xaxis, y = yaxis)
  p
}

#' @title Print prime numbers
#'
#' @param nprimes number of prime numbers needed (default starting from 1)
#' @param prime.start start search for prime numbers from this value
#'
#' @export
primes <- function(nprimes, prime.start=1) {
  is.prime <- function(n) n == 2L || all(n %% 2L:max(2,floor(sqrt(n))) != 0)
  out <- c()
  num <- prime.start

  while(length(out)<nprimes) {
    if (is.prime(num) ) {
      out <- c(out, num) } else {
      out <- out
      }
    num <- num+1
  }
  return(out)
}

#' @title partition samples by specified sizes of training and testing data
#'
#' @param y data vector
#' @param runs number of times that data vector will be partitioned
#' @param nsize_each number of samples in each partition
#'
#' @export
partitionSamples <- function(y, runs, nsize.each = NULL) {
  dat <- data.frame(y)
  dat$index <- c(1:length(y))
  nseeds <- runs*length(nsize.each)
  seeds <- matrix(primes(nseeds), nrow=runs, ncol=length(nsize.each))

  out <- lapply(1:runs, function(r) {
    lapply(1:length(nsize.each), function(i) {
      set.seed(seeds[r,i])
      indices <- sample(dat$index, size=nsize.each[i])
      indices <- sort(indices)
      dat$y[indices]
    })
  })
  return(list(partitions=out,
              seeds=seeds))
}






#' @title Estimate gene weights for cell time
#'
#' @param Y gene by sample expression matrix
#' @param theta sample cell time vector
#'
#' @export
cycle.spml.trainmodel <- function(Y, theta) {

  library(Rfast)
  library(assertthat)
  fit <- spml.reg(theta, t(Y), seb=TRUE)
  return(fit)
}


#' @title Estimate gene weights for cell time
#'
#' @param Y_test gene by testing samples
#' @param theta_test gene by training samples
#' @param theta_train cell times for training samples
#' @param theta_test cell times for test samples
#'
#' @export
cycle.spml.testmodel <- function(Y_test, Y_train, theta_test, theta_train) {

  library(Rfast)
  library(assertthat)
  assert_that(is.matrix(Y_test))
  assert_that(dim(Y_test)[2]==length(theta_test),
              msg = "dimension of testing expression matrix doesn't match length of cell time vector")
  assert_that(is.matrix(Y_train))
  assert_that(dim(Y_train)[2]==length(theta_train),
              msg = "dimension of training expression matrix doesn't match length of cell time vector")

  fit_train <- cycle.spml.trainmodel(Y_train, theta_train)

  pred_cart <- cbind(1,t(Y_test))%*%fit_train$be
  pred_polar <- atan( pred_cart[, 2] / pred_cart[, 1] ) + pi * I(pred_cart[, 1] < 0)

  rho_test <- rFLIndTestRand(pred_polar, theta_test, 9999)
  boot_ci <- rhoFLCIBoot(pred_polar, theta_test, 95, 9999)

  return(list(betahat=fit_train$be,
              theta_pred=pred_polar,
              theta_test=theta_test,
              rho=rho_test[1],
              boot_95ci_low=boot_ci[1],
              boot_95ci_high=boot_ci[2],
              pval=rho_test[2]))
}





#' @title compute Fisher-Lee correlation coefficient
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rFLCorrCoeff <- function(lcdat1, lcdat2) {
  A <- sum(cos(lcdat1)*cos(lcdat2))
  B <- sum(sin(lcdat1)*sin(lcdat2))
  C <- sum(cos(lcdat1)*sin(lcdat2))
  D <- sum(sin(lcdat1)*cos(lcdat2))
  E <- sum(cos(2*lcdat1)); F <- sum(sin(2*lcdat1))
  G <- sum(cos(2*lcdat2)); H <- sum(sin(2*lcdat2))
  n <- length(lcdat1)
  denom <- sqrt(((n*n)-(E*E)-(F*F))*((n*n)-(G*G)-(H*H)))
  rFLVal <- 4*((A*B)-(C*D))/denom
  return(rFLVal)
}

#' @title compute statistical significance of Fisher-Lee correlation coefficient
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param NR number of permutations
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rFLIndTestRand <- function(lcdat1, lcdat2, NR) {
  rFLObs <- rFLCorrCoeff(lcdat1, lcdat2)
  nxtrm <- 1
  for (r in 1:NR) {
    lcdat1Rand <- sample(lcdat1)
    rFLRand <- rFLCorrCoeff(lcdat1Rand, lcdat2)
    if (abs(rFLRand) >= abs(rFLObs)) {nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1); return(c(rFLObs,pval))
}


#' @title compute boostrap confidence interval for correlation coefficient
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param ConfLevel percent confidence interval
#' @param B number of bootstrap sample
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
rhoFLCIBoot <- function(lcdat1, lcdat2, ConfLevel, B) {
  alpha <- (100-ConfLevel)/100; n <- length(lcdat1)
  rFL <- 0
  for (b in 1:B) {
    randind <- sample(1:n, n, replace = TRUE)
    boot1 <- lcdat1[randind]
    boot2 <- lcdat2[randind]
    rFL[b] <- rFLCorrCoeff(boot1, boot2) }

  rFL[B+1] <- rFLCorrCoeff(lcdat1, lcdat2); rFLsort <- sort(rFL)
  return(c(rFLsort[(alpha/2)*(B+1)], rFLsort[(1-alpha/2)*(B+1)]))
}



#' @title JS test for rotational dependence between circular variables
#'
#' @param lcdat1 length n vector of radians
#' @param lcdat2 length n vector of radians
#' @param ConfLevel percent confidence interval
#' @param B number of bootstrap sample
#'
#' @references Pewsey et al. Circular statistics in R
#'
#' @export
JSTestRand <- function(cdat1, cdat2, NR) {
  library(circular)
  CorrJSObs <- cor.circular(cdat1, cdat2); nxtrm <- 1
  for (r in 1:NR) {
    cdat1Rand <- sample(cdat1); CorrJSRand <- cor.circular(cdat1Rand, cdat2)
    if (abs(CorrJSRand) >= abs(CorrJSObs)) {nxtrm <- nxtrm + 1} }
  pval <- nxtrm/(NR+1); return(c(CorrJSObs, pval))
}



#' @title enrichment of neighborhoold samples
circ.dist.neighbors <- function(labels, k) {
  mat_neighbors <- matrix(0, ncol=length(labels), nrow=length(labels))
  colnames(mat_neighbors) <- labels
  rownames(mat_neighbors) <- labels
  N <- length(labels)
  band <- round(k/2)

  for (i in 1:ncol(mat_neighbors)) {
    if (i == 1) {
      neighbors <- c(labels[c((N-band+i):N)], labels[c((i+1):(i+band))])
      mat_neighbors[,i] <- rownames(mat_neighbors) %in% neighbors
    }
    if (i > 1 & i <= band) {
      neighbors <- c(labels[c(1:(i-1), c((N-band+i):N))], labels[c((i+1):(i+band))])
      mat_neighbors[,i] <- rownames(mat_neighbors) %in% neighbors
    }
    if (i > band & i <= (N-band)) {
      neighbors <- c(labels[c((i-band):(i-1))], labels[c((i+1):(i+band))])
      mat_neighbors[,i] <- rownames(mat_neighbors) %in% neighbors
    }
    if (i > (N-band)) {
      neighbors <- c(labels[c((i-band):(i-1))], labels[c((i+1):N, (band -(N-i)):1)])
      mat_neighbors[,i] <- rownames(mat_neighbors) %in% neighbors
    }
    if (i == N) {
      neighbors <- c(labels[c((N-band):(N-1))], labels[c(1:band)])
      mat_neighbors[,i] <- rownames(mat_neighbors) %in% neighbors
    }

  }
  return(mat_neighbors)
}


#' @title Print prime numbers
#'
#' @param nprimes number of prime numbers needed (default starting from 1)
#' @param prime.start start search for prime numbers from this value
#'
#' @export
primes <- function(nprimes, prime.start=1) {
  is.prime <- function(n) n == 2L || all(n %% 2L:max(2,floor(sqrt(n))) != 0)
  out <- c()
  num <- prime.start

  while(length(out)<nprimes) {
    if (is.prime(num) ) {
      out <- c(out, num) } else {
        out <- out
      }
    num <- num+1
  }
  return(out)
}



#' @title partition samples by specified sizes of training and testing data
#'
#' @param y data vector
#' @param runs number of times that data vector will be partitioned
#' @param nsize_each number of samples in each partition
#'
#' @export
partitionSamples <- function(y, runs, nsize.each = NULL) {
  dat <- data.frame(y)
  dat$index <- c(1:length(y))
  nseeds <- runs*length(nsize.each)
  seeds <- matrix(primes(nseeds), nrow=runs, ncol=length(nsize.each))

  out <- lapply(1:runs, function(r) {
    lapply(1:length(nsize.each), function(i) {
      set.seed(seeds[r,i])
      indices <- sample(dat$index, size=nsize.each[i])
      indices <- sort(indices)
      dat$y[indices]
    })
  })

  part_indices <- vector("list", runs)
  for (i in 1:length(out)) {
    part_indices[[i]] <- list(test=c(out[[i]][[i]]))
    part_indices[[i]]$train <- setdiff(c(1:length(y)), part_indices[[i]]$test)
  }

  return(list(partitions=part_indices,
              seeds=seeds))
}



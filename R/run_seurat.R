#' @title Seurat code for cell cycle scoring
#' 
#' @description Compute phase specific cell cycle scores using the
#' code outlined in Seurat.
#'
#' @param Y log2-normalized gene expression matrix.
#'
#' @references Butler, A., Hoffman, P., Smibert, P., Papalexi, E, and
#' Satija, R. (2018). Integrating single-cell transcriptomic data
#' across different conditions, technologies, and species. \emph{Nature
#' Biotechnology} \bold{36}, 411-420. doi:10.1038/nbt.4096
#'
run_seurat <- function(Y, s.genes, g2m.genes, n.bin=25,
                       seed.use=1, random.seed=1) {
    
  set.seed(random.seed)
  genes.list <- list(S.Score = s.genes, G2M.Score = g2m.genes)
  genes.list <- lapply(X = genes.list, FUN = function(x) {
    return(intersect(x = x, y = rownames(Y)))
  })
  cluster.length <- length(x = genes.list)

  ctrl.size = min(vapply(X = genes.list,
                         FUN = length, FUN.VALUE = numeric(1)))

  # order genes by gene-specific averages
  data.avg <- Matrix::rowMeans(Y)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(cut(x = data.avg,
                             breaks=n.bin,
                             include.lowest = TRUE))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    for (j in 1:length(genes.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(sample(data.cut[which(data.cut ==
                                  data.cut[genes.use[j]])],
                                  size = ctrl.size, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(x = ctrl.use),
                        ncol = ncol(Y))
  for (i in 1:length(ctrl.use)) {
    genes.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = Y[which(rownames(Y) %in% genes.use),])
  }

  # compute phase-specific scores
  genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length,
                         ncol = ncol(Y))
  for (i in 1:cluster.length) {
    genes.use <- genes.list[[i]]
    data.use <- Y[which(rownames(Y) %in% genes.use), ,drop = FALSE]
    genes.scores[i, ] <- Matrix::colMeans(data.use)
  }
  genes.scores.use <- genes.scores - ctrl.scores
  rownames(genes.scores.use) <- c("S", "G2M")
  genes.scores.use <- data.frame(t(genes.scores.use))

  assignments <- apply(X = genes.scores.use,
                       MARGIN = 1,
                       FUN = function(scores, first = "S", second = "G2M",
                           null = "G1") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  })

  out <- data.frame(genes.scores.use, assignments)
  rownames(out) <- colnames(Y)

  out$assignments <- factor(out$assignments, levels=c("G1", "S", "G2M"))

  return(out)
}

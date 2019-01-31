#!/usr/bin/env Rscript

# Description:
#   This script was written for extracting summary information
#  from the intensity values.

args <- commandArgs(trailingOnly = TRUE)

for (i in seq_along(args)) {
  cat(sprintf("Argument %d: %s", i, args[i]), sep = "\n")
}


ncores <- as.numeric(args[1])
ngenes <- as.numeric(args[2])
fold <- as.numeric(args[3])

dir <-"/project2/gilad/joycehsiao/fucci-seq"
source(file.path(dir,"code/working/run_methods.R"))
source(file.path(dir,"peco/R/fit.trendfilter.generic.R"))

fold_indices <- readRDS(file=file.path(dir,"data/results/mixed_overall_fold_indices.rds"))

log2cpm.quant <- readRDS(file.path(dir,"output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds"))

data_training <- log2cpm.quant[,fold_indices[[fold]]$train]
data_withheld <- log2cpm.quant[,fold_indices[[fold]]$test]

seurat.genes <- readLines(con = file.path(dir,
                  "data/cellcycle-genes-previous-studies/seurat_cellcycle/regev_lab_cell_cycle_genes.txt"))
seurat.genes <- list(s.genes=seurat.genes[1:43],
                     g2m.genes=seurat.genes[44:97])


fits_all <- readRDS(file.path(dir,"output/npreg-trendfilter-quantile.Rmd/fit.quant.rds"))
genes_all <- names(fits_all)[order(sapply(fits_all,"[[",3), decreasing=T)]
which_genes <- genes_all[1:ngenes]

# cyclical_genes <- readRDS(file=file.path(dir, paste0(
#                            "data/results/data_training_cyclical_genes.fold.",fold,".rds")))

# which_genes <- rownames(cyclical_genes)[order(cyclical_genes$pve,
#                                               decreasing = T)[1:ngenes]]

print(fold)

eset <- readRDS(file=file.path(dir,"data/eset-final.rds"))
library(Biobase)
pdata <- pData(eset)
theta <- pdata$theta
names(theta) <- rownames(pdata)

log2cpm <- t(log2(1+(10^6)*(t(exprs(eset))/pdata$molecules)))
log2cpm <- log2cpm[,match(colnames(log2cpm.quant),colnames(log2cpm))]
theta_ordered <- theta[match(colnames(log2cpm.quant),names(theta))]
pdata_ordered <- pdata[match(colnames(log2cpm.quant),rownames(pdata)),]

Y_train_normed_fold_topX <- data_training[which(rownames(data_training) %in% which_genes),]
theta_train_fold <- theta_ordered[fold_indices[[fold]]$train]


fit.train <- cycle.npreg.insample(Y = Y_train_normed_fold_topX,
                                  theta = theta_train_fold,
                                  polyorder=2,
                                  ncores=ncores,
                                  method.trend="trendfilter")

Y_test_normed_fold_topX <- data_withheld[which(rownames(data_withheld) %in% which_genes),]
Y_test_fold <- log2cpm[,fold_indices[[fold]]$test]
theta_test_fold <- theta_ordered[fold_indices[[fold]]$test]
pdata_test_fold <- pdata_ordered[fold_indices[[fold]]$test,]

fit.test <- run_methods(Y_test=Y_test_fold,
                     Y_test_normed=Y_test_normed_fold_topX,
                     theta_test=theta_test_fold,
                     training_model=fit.train,
                     seurat.genes=seurat.genes,
                     pdata_test=pdata_test_fold,
                     fdata=fData(eset),
                     ncores=ncores, maxiter=30)

out <- list(fit.train=fit.train,
            fit.test=fit.test)
#}
#names(fits) <- paste0("fold.", 1:length(fold_indices))
saveRDS(out,
        file=file.path(dir, paste0("data/results/finalizing/mixed_results.overallcyclical.fold.",fold,".top",ngenes,".rds")))






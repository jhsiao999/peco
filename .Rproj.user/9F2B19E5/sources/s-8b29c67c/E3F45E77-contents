# eset_raw <- readRDS("data/eset-raw.rds")
# pData(eset_raw)$cell_number

# Reproduce Linear discriminant analysis ---------------------------------------------------------

eset <- readRDS("data/eset-raw.rds")

counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

group_3 <- rep(2,dim(pdata)[1])
group_3[grep("0", pdata$cell_number)] <- 0
group_3[grep("1", pdata$cell_number)] <- 1

## create data frame
data <- pdata %>% dplyr::select(experiment:concentration, mapped, molecules)
data <- data.frame(data, group = group_3)

library("cowplot")
library("dplyr")
library("edgeR")
library("ggplot2")
library("MASS")
library("tibble")
library("Biobase")

## perform lda
data_lda <- lda(group ~ concentration + molecules, data = data)
data_lda_p <- predict(data_lda, newdata = data[,c("concentration", "molecules")])$class

## determine how well the model fix
table(data_lda_p, data[, "group"])

data$data_lda_p <- data_lda_p

## identify the outlier
outliers_lda <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_lda_p == "two")
outliers_lda

anno <- pdata
anno$molecule_outlier <- row.names(anno) %in% outliers_lda$sample_id

## plot before and after
# data$group[data$group=="no"] <- 0
# data$group[data$group=="one"] <- 1
# data$group[data$group=="two"] <- 2
# data$data_lda_p[data$data_lda_p=="no"] <- 0
# data$data_lda_p[data$data_lda_p=="one"] <- 1
# data$data_lda_p[data$data_lda_p=="two"] <- 2
plot_before <- ggplot(data, aes(x = concentration, y = molecules / 10^3,
                                color = as.factor(group))) +
  geom_text(aes(label = cell_number), alpha = 0.5) + labs(col="No. cells") +
  labs(x = "Concentration", y = "Gene molecules (thousands)", title = "Observed data") +
#  scale_color_brewer(palette = "Dark2") +
#  theme(legend.position = "none") +
  scale_color_manual(values= brewer.pal(8,"Dark2")[1:3])

plot_after <- ggplot(data, aes(x = concentration, y = molecules / 10^3,
                               color = as.factor(data_lda_p))) +
  geom_text(aes(label = cell_number), alpha = 0.5) + labs(col="No. cells") +
  labs(x = "Concentration", y = "Gene molecules (thousands)", title = "Inferred number of cells") +
#  scale_color_brewer(palette = "Dark2") +
#  theme(legend.position = "none") +
  scale_color_manual(breaks=c("0", "1", "2"),
                     values= brewer.pal(8,"Dark2")[c(2,1,3)])


plot_grid(plot_before + theme(legend.position=c(.8,.85)),
          plot_after + theme(legend.position=c(.8,.85)),
          labels = LETTERS[1:2])


## calculate convertion
anno$ercc_conversion <- anno$mol_ercc / anno$reads_ercc

anno$conversion <- anno$mol_hs / anno$reads_hs

## try lda
data$conversion <- anno$conversion
data$ercc_conversion <- anno$ercc_conversion

data_ercc_lda <- lda(group ~ ercc_conversion + conversion, data = data)

data_ercc_lda_p <- predict(data_ercc_lda,  newdata = data[,c("ercc_conversion", "conversion")])$class

## determine how well the model fix
table(data_ercc_lda_p, data[, "group"])


data$data_ercc_lda_p <- data_ercc_lda_p

## identify the outlier
outliers_conversion <- data %>% rownames_to_column("sample_id") %>% filter(cell_number == 1, data_ercc_lda_p == "two")
outliers_conversion


## create filter
anno$conversion_outlier <- row.names(anno) %in% outliers_conversion$sample_id

## plot before and after
plot_ercc_before <- ggplot(data, aes(x = ercc_conversion, y = conversion,
                                     color = as.factor(group))) +
  geom_text(aes(label = cell_number), alpha = 0.5) + labs(col="No. cells")+
  labs(x = "Convertion of ERCC spike-ins", y = "Conversion of genes", title = "Observed data") +
  scale_color_brewer(palette = "Dark2")
#  theme(legend.position = "none")

plot_ercc_after <- ggplot(data, aes(x = ercc_conversion, y = conversion,
                                    color = as.factor(data_ercc_lda_p))) +
  geom_text(aes(label = cell_number), alpha = 0.5) + labs(col="No. cells")+
  labs(x = "Convertion of ERCC spike-ins", y = "Conversion of genes", title = "Inferred number of cells") +
  scale_color_brewer(palette = "Dark2")
#  theme(legend.position = "none")

plot_grid(plot_ercc_before + theme(legend.position=c(.8,.85)),
          plot_ercc_after + theme(legend.position=c(.8,.85)),
          labels = LETTERS[3:4])

# RNA-seq data quality, PCA, etc ---------------------------------------------------------

eset <- readRDS("data/eset-final.rds")

counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))



# drop-out rate, number of singletons
mean(rowMeans(log2cpm.all>0))
sd(rowMeans(log2cpm.all>0))

inds <- unique(pdata$chip_id)
cbind(inds, sapply(1:length(inds), function(i) {
  tmp <- log2cpm.all[,pdata$chip_id == inds[i]]
  mean(rowMeans(tmp>0)) }))

library(ggplot2)
library(cowplot)
pdata$dropout <- colMeans(log2cpm.all >0)
pdata$experiment <- factor(pdata$experiment)
pdata$chip_id <- factor(pdata$chip_id)
plot_grid(ggplot(pdata, aes(x=experiment, y=log10(molecules),
                            group=experiment, col=experiment)) +
            geom_boxplot() + coord_flip() +
            ylab("log10 molecule count") + xlab("C1 plate") +
            ggtitle("Library size and C1 plate") + theme(legend.position="none"),
          ggplot(pdata, aes(x=experiment, y=dropout,
                            group=experiment, col=experiment)) +
            geom_boxplot() + coord_flip() +
            ylab("Sample detection rate") + xlab("C1 plate") +
            ggtitle("Detection rate and C1 plate") + theme(legend.position="none"),
          ggplot(pdata, aes(x=experiment,
                            group=experiment, fill=experiment)) +
            geom_bar() + coord_flip() +
            ylab("Singletons") + xlab("C1 plate") +
            ggtitle("Singletons and C1 plate") + theme(legend.position="none"),
          ggplot(pdata, aes(x=chip_id, y=log10(molecules),
                            group=chip_id, col=chip_id)) +
            geom_boxplot() + coord_flip() +
            ylab("log10 molecule count") + xlab("Individual") +
            scale_color_brewer(palette="Dark2") +
            ggtitle("Library size and individual") + theme(legend.position="none"),
          ggplot(pdata, aes(x=chip_id, y=dropout,
                            group=chip_id, col=chip_id)) +
            geom_boxplot() + coord_flip() +
            ylab("Sample detection rate") + xlab("Individual") +
            scale_color_brewer(palette="Dark2") +
            ggtitle("Detection rate and individual") + theme(legend.position="none"),
          ggplot(pdata, aes(x=chip_id, fill=chip_id)) +
            geom_bar() + coord_flip() +
            ylab("Singletons") + xlab("Individual") +
#            scale_fill_manual("legend", values = c("A" = "black", "B" = "orange", "C" = "blue"))
            scale_fill_brewer("Individual",palette="Dark2") +
            ggtitle("Singletons and individual") + theme(legend.position="none"),
          nrow=2, labels=LETTERS[1:6])


# significance tests
kruskal.test(pdata$dropout ~ pdata$chip_id)
kruskal.test(log10(pdata$molecules) ~ pdata$chip_id)

kruskal.test(pdata$dropout ~ pdata$experiment)
kruskal.test(log10(pdata$molecules) ~ pdata$experiment)

dd <- table(pdata$experiment)/96
dd <- as.numeric(dd)
mean(dd)
sd(dd)
range(dd)

table(pdata$chip_id)

# Result section on sequencing reads
eset_raw <- readRDS("data/eset-raw.rds")
mean(pData(eset_raw)$mapped)
sd(pData(eset_raw)$mapped)
range(pData(eset_raw)$mapped)


# Supplemental Table listing sequencing results by individual
# contains all data before filtering of samples or genes
eset <- readRDS("data/eset-final.rds")

pdata <- pData(eset)

inds <- unique(pdata$chip_id)

# number of reads with a valid UMI and also mapped to the index
do.call(rbind, lapply(1:length(inds), function(i) {
  tmp <- pdata$mapped[pdata$chip_id == inds[i]]
  data.frame(ind=inds[i],
             mn=mean(tmp),
             std=sd(tmp), stringsAsFactors = F) }))
cbind(mean(pdata$mapped), sd(pdata$mapped))

# proportin of unmapped reads with a valid UMI
percents_unmapped <- pData(eset)$unmapped/pData(eset)$umi
do.call(rbind, lapply(1:length(inds), function(i) {
  tmp <- percents_unmapped[pdata$chip_id == inds[i]]
  data.frame(ind=inds[i],
             mn=mean(tmp)*100,
             std=sd(tmp)*100, stringsAsFactors = F) }))
cbind(mean(percents_unmapped)*100, sd(percents_unmapped)*100)


# number of singletons
table(pdata$chip_id)
length(pdata$chip_id)


# fraction of ERCC reads: reads_ercc/mapped
percentage_ercc <- pData(eset)$ercc_percentage
do.call(rbind, lapply(1:length(inds), function(i) {
  tmp <- percentage_ercc[pdata$chip_id == inds[i]]
  data.frame(ind=inds[i],
             mn=mean(tmp)*100,
             std=sd(tmp)*100, stringsAsFactors = F) }))
cbind(mean(percentage_ercc)*100, sd(percentage_ercc)*100)


# total molecular count
molecules <- pData(eset)$molecules
do.call(rbind, lapply(1:length(inds), function(i) {
  tmp <- molecules[pdata$chip_id == inds[i]]
  data.frame(ind=inds[i],
             mn=mean(tmp)/(10^6),
             std=sd(tmp)/(10^6), stringsAsFactors = F) }))
cbind(mean(molecules)/(10^6), sd(molecules)/(10^6))




# PCA variation -------------------------------------------------------------------
# selection of technical factor
library(dplyr)
covariates <- pData(eset) %>% dplyr::select(experiment, well, chip_id,
                                                   concentration, raw:unmapped,
                                                   starts_with("detect"),  molecules)
# look at the first 6 PCs
pca_log2cpm <- prcomp(t(log2cpm.all), scale. = TRUE, center = TRUE)
pcs <- pca_log2cpm$x[, 1:6]


# R-square between PCs and the covariates
get_r2 <- function(x, y) {
  stopifnot(length(x) == length(y))
  model <- lm(y ~ x)
  stats <- summary(model)
  return(stats$adj.r.squared)
}
r2 <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
             dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

# plot heatmap
library(heatmap3)
heatmap3(r2, cexRow=1, cexCol=1, margins=c(4,15), scale="none",
         ylab="", main = "", symm=F,
         Colv=NA, showColDendro = F,
         labRow = c("C1: C1 batch",
                    "Well: C1 capture site (well) ID",
                    "Individual",
                    "cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"
                    ))

# r-squared between PC1 and individual
f1 <- lm(pcs[,1] ~ covariates$chip_id)
f2 <- lm(pcs[,1] ~ covariates$molecules)
summary(f1)

summary(f1)$adj.r.squared
summary(f2)$adj.r.squared


# PC proportions
100*((pca_log2cpm$sdev^2)/sum(pca_log2cpm$sdev^2))[1:6]


# Correlation between Technical factors
cor_tech <- cor(as.matrix(covariates[,4:11]),use="pairwise.complete.obs")
library(RColorBrewer)
heatmap3(cor_tech, symm = TRUE, margins=c(7,15),
         col=brewer.pal (9, "Blues" ), cexRow=1, cexCol=1, scale="none",
         labRow = c("cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"),
         labCol = c("cDNA",
                    "Reads raw",
                    "Reads UMI",
                    "Reads mapped",
                    "Reads unmapped",
                    "ERCC",
                    "ENSG",
                    "Molecules"))

summary(lm(covariates$detect_hs ~ covariates$experiment))

summary(lm(covariates$detect_hs ~ covariates$mapped))

summary(lm(covariates$detect_hs ~ covariates$chip_id))

summary(lm(covariates$detect_hs ~ covariates$molecules))

summary(lm(covariates$detect_hs ~ covariates$concentration))

# pca of top 10% experssing genes
log2cpm_mean <- rowMeans(log2cpm.all)
log2cpm_top <- log2cpm.all[rank(log2cpm_mean) / length(log2cpm_mean) > 1 - 0.1, ]
dim(log2cpm_top)

pca_top <- prcomp(t(log2cpm_top), scale. = T, center = T)

## look at the first 6 PCs
pcs <- pca_top$x[, 1:6]

## generate the data
r2_top <- matrix(NA, nrow = ncol(covariates), ncol = ncol(pcs),
                 dimnames = list(colnames(covariates), colnames(pcs)))
for (cov in colnames(covariates)) {
  for (pc in colnames(pcs)) {
    r2_top[cov, pc] <- get_r2(covariates[, cov], pcs[, pc])
  }
}

## plot heatmap
heatmap3(r2_top, cexRow=1, cexCol=1, margins=c(8,8),
         Colv=F, showColDendro = F,
         labRow = c("C1: C1 batch",
                    "Well: C1 capture site (well) ID",
                    "Individual",
                    "cDNA: cDNA concentration",
                    "Reads raw: Raw reads",
                    "Reads UMI: Reads with valid UMI",
                    "Reads mapped: Reads with valid UMI mapped to genome",
                    "Reads unmapped: Reads with valid UMI not mapped to genome",
                    "ERCC: ERCC gene proportion detected",
                    "ENSG: Endogeneous gene proportion detected",
                    "Molecules: Total sample molecule count"
         ))


# identify variable genes ------------------------------------------
eset <- readRDS("data/eset-final.rds")
counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)
log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

library(matrixStats)
ss <- rownames(log2cpm.all)[order(rowVars(log2cpm.all), decreasing = T)]
pca_log2cpm <- prcomp(t(log2cpm.all), scale. = TRUE, center = TRUE)

cycle <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")
counts.pseudo <- counts
counts.pseudo[which(counts==0, arr.ind = T)] <- NA
pp <- data.frame(mn.exp=rowMeans(counts.pseudo, na.rm = T),
                 mn=rowMeans(counts),
                 cv=rowSds(counts)/rowMeans(counts),
                 ss=rowSds(counts))


# given mean the genes with high CV
pp$mn_bin <- cut(pp$mn, breaks=20, include.lowest=T)
for (i in 1:20) {
  uu <- levels(pp$mn_bin)
  pp$cv_z[pp$mn_bin == uu[i]] <- scale(pp$cv[pp$mn_bin == uu[i]])
}
pp$cv_z <- scale(pp$cv)
pp$out <- (abs(pp$cv_z) > 2)
rownames(pp) <- rownames(counts)

plot(x=log10(pp$mn), y=pp$cv_z, cex=.6, pch=16, col="gray80")
points(x=log10(pp$mn), y=pp$cv_z,
       col=c("gray20", "red")[pp$out+1])

sum(rownames(pp)[pp$out==T] %in% cycle$ensembl)
which(rownames(pp)[pp$out==T] %in% cycle$ensembl)
rownames(pp)[192]


# plot(pca_log2cpm$rotation[,1],
#      col=c("gray20", "red")[pp$out+1])
# plot(pca_log2cpm$rotation[,2],
#      col=c("gray20", "red")[pp$out+1])
# plot(pca_log2cpm$rotation[,3],
#      col=c("gray20", "red")[pp$out+1])
# plot(pca_log2cpm$rotation[,4],
#      col=c("gray20", "red")[pp$out+1])





# Intensity batch effect ---------------------------------------------------------
lm.rfp <- lm(rfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.gfp <- lm(gfp.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
             data = pdata)
lm.dapi <- lm(dapi.median.log10sum~factor(chip_id)+factor(experiment) + factor(image_label),
              data = pdata)

library(ibd)
aov.lm.rfp <- Anova(lm.rfp, type = "III")
aov.lm.gfp <- Anova(lm.gfp, type = "III")
aov.lm.dapi <- Anova(lm.dapi, type = "III")

aov.lm.rfp
aov.lm.gfp

library(ggplot2)
# ggplot(pdata, aes(gfp.median.log10sum, col=factor(experiment))) +
#   geom_density()
# ggplot(pdata, aes(rfp.median.log10sum, col=factor(experiment))) +
#   geom_density()
# ggplot(pdata, aes(gfp.median.log10sum.adjust, col=factor(experiment))) +
#   geom_density()

library(ggplot2)
library(RColorBrewer)
library(cowplot)
plot_grid(
  ggplot(pdata,
       aes(x=experiment, y=gfp.median.log10sum, fill=experiment)) +
  geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
  ylab("FUCC score") + xlab("Individual") +
  ggtitle("GFP before batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)),
  ggplot(pdata,
         aes(x=experiment, y=rfp.median.log10sum, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("RFP before batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)),
  ggplot(pdata,
         aes(x=experiment, y=gfp.median.log10sum.adjust, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("GFP after batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)),
  ggplot(pdata,
         aes(x=experiment, y=rfp.median.log10sum.adjust, fill=experiment)) +
    geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
    ylab("FUCC score") + xlab("Individual") +
    ggtitle("RFP after batch-correction") + labs(fill="experiment") +
    theme(axis.text.x  = element_text(angle=90, vjust=1, hjust=1,size=10),
          axis.text.y  = element_text(angle=0, vjust=0.2, size=10)),
  labels=LETTERS[1:4])






# Fucci phase making and properties --------------------------------------------------
df <- readRDS("data/eset-final.rds")

par(mfrow=c(1,1))
plot(x=pData(df)$rfp.median.log10sum.adjust, xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),
     y=pData(df)$gfp.median.log10sum.adjust, pch=16, cex=.5, col="gray50",
     xlab="RFP log10 sum intensity",
     ylab="GFP log10 sum intensity", axes=F)
axis(1); axis(2)

pca <- prcomp(cbind(pData(df)$rfp.median.log10sum.adjust,
                    pData(df)$gfp.median.log10sum.adjust))
(pca$sdev^2)/sum(pca$sdev^2)
plot(pca$x[,1], pca$x[,2], pch=16, cex=.5, xlim=c(-4, 4), ylim=c(-4,4),
     xlab="PC1 (67%)", ylab="PC2 (33%)",
     main = "fucci intensities PC1 vs PC2", col="gray50", axes=F)
axis(1);axis(2)
abline(h=0,v=0, col="gray50", lty=2)
par(new=TRUE)

library(circular)
theta <- coord2rad(pca$x)
plot(circular(theta), stack=T, shrink=1.3, cex=.5, bins=200)


# fucci intensity and expression by phase
fits_all <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)
pve_all_ord <- pve_all[order(pve_all, decreasing = T)]
genes <- names(pve_all_ord)[1:5]
data <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")

source("peco/R/fit.trendfilter.generic.R")
data_quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
theta_final <- (2*pi-as.numeric(theta))%%(2*pi)
sample_ord <- rownames(pData(df))[order(theta_final)]
data_quant_ord <- data_quant[,match(sample_ord,colnames(data_quant))]
fits_tmp <- lapply(1:5, function(i) {
  ii <- which(rownames(data_quant_ord)==genes[i])
  fit_g <- fit.trendfilter.generic(data_quant_ord[ii,])
  fun_g <- approxfun(x=as.numeric(theta_final[order(theta_final)]),
                     y=as.numeric(fit_g$trend.yy), rule=2)
  pred_g <- fun_g(as.numeric(theta_final[order(theta_final)]))
  return(pred_g)
})
names(fits_tmp) <- genes

par(mfrow=c(1,1))
plot(x=(2*pi-as.numeric(theta))%%(2*pi),
     y=pData(df)$gfp.median.log10sum.adjust, col="forestgreen",
     ylim=c(-1.5, 1.5), pch=16, cex=.5,
     xlab="Fucci phase", ylab="Fucci intensities",
     main="Fucci intensities", axes=F)
axis(2); axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
              labels=c(0,expression(pi/2), expression(pi), expression(3*pi/2),
                       expression(2*pi)))
abline(h=0, lty=1, col="black", lwd=.7)
points(x=(2*pi-as.numeric(theta))%%(2*pi),
       y=pData(df)$rfp.median.log10sum.adjust, col="firebrick",
       ylim=c(-1.5, 1.5), pch=16, cex=.5)
theta_tmp <- (2*pi-as.numeric(theta))%%(2*pi)
theta_tmp <- theta_tmp[order(theta_tmp)]
# labs_at <- sapply(c(0,pi/2, pi, 3*pi/2, 2*pi),
#                   function(x) {which.min(abs(theta_tmp-x))})
par(mfrow=c(1,5), mar=c(3,4,3,1))
for (i in 1:5) {
  plot(x=theta_final[order(theta_final)],
       y=data_quant_ord[rownames(data_quant_ord)==genes[i],], col="gray50",
       xlab="Fucci phase",
       ylab="Normalized log2CPM", axes=F, cex=.5, pch=16,
       main = labs[i])
  points(x=theta_final[order(theta_final)],
         y=fits_tmp[[i]], col="blue", pch=16, cex=.7)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                              expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
}


# distribution of Fucci phase
par(mfrow=c(1,1))
hist((2*pi-as.numeric(theta)), nclass=25,
     main="", xlab="Fucci phase",
     axes=F)
axis(2);
axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                             expression(3*pi/2), expression(2*pi)))


# Top cyclical genes before training ----------------------------------------------
fits_all <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve_all <- sapply(fits_all, "[[", 3)

pve_all_ord <- pve_all[order(pve_all, decreasing = T)]

head(pve_all_ord)

genes <- names(pve_all_ord)[1:5]
data <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")
par(mfrow=c(2,3))
for (i in 1:5) {
  ii <- which(names(fits_all)==genes[i])
  plot(data[rownames(data)==genes[i],], col="gray50",
       xlab="Fucci ordering",
       ylab="Normalized log2CPM", axes=F,
       main = labs[i])
  axis(1); axis(2)
  points(fits_all[[ii]]$trend.yy, col="red", pch=16)
}


# Top cyclical genes enriched for GO cell-cycle term
library('org.Hs.eg.db')
key = 'GO:0007049'
columns = c('ENSEMBL', 'GENENAME')
s = select(org.Hs.eg.db, key, columns, keytype='GO')

top100 <- names(pve_all_ord)[1:100]

genesInGo <- names(pve_all_ord) %in% s$ENSEMBL
genesInTop <- names(pve_all_ord) %in% top100

mat <- table(genesInTop, genesInGo)
mat <- mat[2:1, 2:1]

# odds of top100 in GO / odds of non-top100 in GO
# test if this OR is greater than 1
fisher.test(mat,alternative="greater")




# Top cyclical genes enriched for Macosko cell-cycle genes
cc <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")

for (pp in c(5, 50, 100)) {
  top100 <- names(pve_all_ord)[1:pp]
  genesInCC <- names(pve_all_ord) %in% cc$ensembl
  genesInTop <- names(pve_all_ord) %in% top100
  mat <- table(genesInTop, genesInCC)
  mat <- mat[2:1, 2:1]
  # odds of top100 in GO / odds of non-top100 in GO
  # test if this OR is greater than 1
  print(fisher.test(mat,alternative="greater"))
}


# Proportion of variance explained vs dropout rate ------------------------
eset <- readRDS("data/eset-final.rds")
counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

dropout <- rowMeans(log2cpm.all == 0)
dropout <- dropout[match(names(pve_all_ord),names(dropout))]
head(names(dropout))
head(names(pve_all_ord))

mean(dropout[1:100])
mean(dropout[1:50])
mean(dropout[1:5])



# proportion of variance: permutation ---------------------------------------------------
fit.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.quant.rds")
pve <- sapply(fit.quant, "[[", "trend.pve")
log2cpm.quant <- readRDS("output/npreg-trendfilter-quantile.Rmd/log2cpm.quant.rds")
macosko <- readRDS("data/cellcycle-genes-previous-studies/rds/macosko-2015.rds")


perm.lowmiss <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.trend.perm.lowmiss.rds")
perm.highmiss <- readRDS("output/npreg-trendfilter-quantile.Rmd/fit.trend.perm.highmiss.rds")

pve.perm.lowmiss <- sapply(perm.lowmiss, "[[", "trend.pve")
pve.perm.highmiss <- sapply(perm.highmiss, "[[", "trend.pve")

summary(pve.perm.lowmiss)

#Compute p-value based on two different distributions. High consistency between the two.
B <- length(pve.perm.lowmiss)
pval.perm.low <- sapply(fit.quant, function(x) (1+sum(pve.perm.lowmiss > x$trend.pve))/(1+B))
pval.perm.high <- sapply(fit.quant, function(x) (1+sum(pve.perm.highmiss > x$trend.pve))/(1+B))

#quickily check top 100 enrichment for cell cycle genes.
enrich.order <- function(cutoffs, metrics, cyclegenes, allgenes) {
  #  out <- order(mad.ratio$smash.mad.ratio)
  # cutoffs <- c(100, 200, 300)
  cycle.rich <- lapply(cutoffs, function(x) {
    which_top <- order(metrics, decreasing = T)[1:x]
    sig.cycle <- sum(allgenes[which_top] %in% cyclegenes)/x
    non.cycle <- sum(allgenes[-which_top] %in% cyclegenes)/(length(allgenes)-x)
    mat <- table(allgenes %in% cyclegenes,
                 allgenes %in% allgenes[which_top])
    mat <- mat[2:1,2:1]
    res <- fisher.test(mat, alternative = "greater")
    # cbind(as.numeric(sum(allgenes[which_top] %in% cyclegenes)),
    #       sig.cycle/non.cycle)
  })
  # colnames(cycle.rich) <- cutoffs
  # rownames(cycle.rich) <- c("nsig.genes.cycle", "fold.sig.vs.nonsig.cycle")
  return(cycle.rich)
}

enrich.order(cutoffs = c(5, 50, 100),
             metrics = 1-pval.perm.high, cyclegenes = macosko$ensembl,
             allgenes = rownames(log2cpm.quant))

enrich.order(cutoffs = c(5, 50, 100),
             metrics = 1-pval.perm.low, cyclegenes = macosko$ensembl,
             allgenes = rownames(log2cpm.quant))


# dropout
eset <- readRDS("data/eset-final.rds")
counts <- exprs(eset)
counts <- counts[grep("ENSG", rownames(counts)), ]
pdata <- pData(eset)

log2cpm.all <- t(log2(1+(10^6)*(t(counts)/pdata$molecules)))

dropout <- rowMeans(log2cpm.all == 0)
dropout <- dropout[match(names(1-pval.perm.low),names(dropout))]
head(names(dropout))
head(names(pval.perm.low))

mean(dropout[1:100])
mean(dropout[1:50])
mean(dropout[1:5])


nn <- names(pve_all_ord)[1:5]
pval.perm.low[names(pval.perm.low) %in% nn]

# Top cyclical genes in other datasets ----------------------------------------

# import leng data

#log2cpm_quant <- readRDS("data/rnaseq-previous-studies/leng/log2cpm_quant.rds")
df <- readRDS("data/rnaseq-previous-studies/leng/HumanLengESC.rds")
log2cpm <- log2(10^6*t(t(exprs(df)+1)/colSums(exprs(df))))
log2cpm_quant <- t(scale(t(log2cpm)))
pdata <- readRDS("data/rnaseq-previous-studies/leng/pdata_filtered.rds")


log2cpm_quant_sub <- log2cpm_quant[,match(rownames(pdata),
                                          colnames(log2cpm_quant))]

genes <- c("CDK1", "TOP2A", "UBE2C", "HIST1H4E", "HIST1H4C")
par(mfrow=c(2,3))
for (g in 1:length(genes)) {
  boxplot(log2cpm_quant_sub[which(rownames(log2cpm_quant_sub)==genes[g]),]~pdata$cell_state)
}







# prediction error -----------------------------------

diff_time_wrapper <- function(results_list) {

  methods_list <- sapply(names(results_list),
                         function(x) strsplit(x, split=".", fixed=TRUE)[[1]][2])

  diff_time_list <- do.call(rbind, lapply(1:length(results_list), function(i) {
    diff_time <- results_list[[i]]$diff_time
    diff_mean <- mean(diff_time/2/pi)
    diff_sd <- sd(diff_time/2/pi) #/sqrt(ncol(results_list[[1]]$Y))
    diff_se <- sd(diff_time/2/pi)/sqrt(length(diff_time))
    diff_median <- median(diff_time/2/pi)
    diff_mad <- mad(diff_time/2/pi)

    return(data.frame(diff_mean=diff_mean,
                      diff_sd=diff_sd,
                      diff_se=diff_se,
                      diff_median=diff_median,
                      diff_mad=diff_mad,
                      methods=methods_list[i]))
  }) )
  return(diff_time_list)
}

# code saved in code/working/finalizing/job_run_methods.train.ind.R

ngenes <- c(5, seq(10,500, by=10))
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")

eval_res <- do.call(rbind, lapply(1:length(inds), function(j) {
  ind <- inds[j]
  foo <- do.call(rbind, lapply(1:length(ngenes), function(i) {
    ngene <- ngenes[i]
    #  train_topX <- do.call(rbind, lapply(1:5, function(fold) {
    fl_name <- list.files("data/results/finalizing",
                          pattern=paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"),
                          full.names = TRUE)
    df <- readRDS(fl_name)
    out <- diff_time_wrapper(df$fit.test)
    out$ngenes <- ngene
    return(out)
  }) )
  foo$ind <- ind
  return(foo)
}) )


eval_res$ind <- factor(eval_res$ind,
                       levels=c("NA19098", "NA18511", "NA18855", "NA19101", "NA18870", "NA19160"))
library(ggplot2)
ggplot(subset(eval_res,methods=="supervised"),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind), alpha=.5) + #geom_line(lty=3) +
  ylab("Prediction error (% circle)") + xlab("Top X cyclical genes")  +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.30) +
  geom_hline(yintercept=.25, col="red") +
  theme_light() + labs(color="individual") +
  scale_color_brewer(palette="Dark2")

ggplot(subset(eval_res,methods=="supervised" & ngenes <=50),
       aes(x=ngenes, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind), alpha=.5) + #geom_line(lty=3) +
  ylab("Prediction error (% circle)") + xlab("Top X cyclical genes")  +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2) +
  geom_line(aes(color=ind), alpha=.5) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1, alpha=.5) +
  ylim(0,.30) +
  geom_hline(yintercept=.25, col="red") +
  theme_light() + labs(col="individual") +
  scale_color_brewer(palette="Dark2")


# mean prediction error betewen 5 to 100
eval_res_sub <- subset(eval_res,methods=="supervised")
aggregate(diff_mean ~ ngenes, data=eval_res_sub, mean)

# prediction error range of 5 genes predictor
subset(eval_res,ngenes==5 & ind=="NA18511")


# plot out the fit using a simple predictor of 5 genes
# include these in the supplemental figure
source("peco/R/cycle.corr.R")
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")
ngene <- 5
eval_fit <- lapply(1:length(inds), function(i) {
  ind <- inds[i]
  fl_name <- list.files("data/results/finalizing",
                        pattern=paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"),
                        full.names = TRUE)
  df <- readRDS(fl_name)
  return(df$fit.test)
})
names(eval_fit) <- inds



par(mfrow=c(2,5), mar=c(3,4,3,1))
#labs <- rownames(eval_seurat[[1]]$`5`$fit.supervised$Y_reordered)
genes <- c("ENSG00000170312","ENSG00000175063",
           "ENSG00000131747", "ENSG00000198518", "ENSG00000197061")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")
for (i in 1:5) {
  if (i>1) { ylab <- ""} else { ylab <- "Normalized log2CPM"}
  ii <- which(rownames(eval_fit$NA19098$fit.supervised$Y_reordered) == genes[i])
  plot(x=eval_fit$NA19098$fit.supervised$cell_times_reordered,
       y=eval_fit$NA19098$fit.supervised$Y_reordered[ii,], col="gray50",
       xlab="Fucci phase",
       ylab=ylab, axes=F, cex=.5, pch=16,
       main = labs[i])
  points(x=eval_fit$NA19098$fit.supervised$cell_times_reordered,
         y=eval_fit$NA19098$fit.supervised$mu_reordered[ii,], col="blue", pch=16, cex=.7)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
}
for (i in 1:5) {
  if (i>1) { ylab <- ""} else { ylab <- "Normalized log2CPM"}
  ii <- which(rownames(eval_fit$NA19098$fit.supervised$Y_reordered) == genes[i])
  plot(x=eval_fit$NA19160$fit.supervised$cell_times_reordered,
       y=eval_fit$NA19160$fit.supervised$Y_reordered[ii,], col="gray50",
       xlab="Fucci phase",
       ylab=ylab, axes=F, cex=.5, pch=16,
       main = labs[i])
  points(x=eval_fit$NA19160$fit.supervised$cell_times_reordered,
         y=eval_fit$NA19160$fit.supervised$mu_reordered[ii,], col="blue", pch=16, cex=.7)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
}





df_tmp <- data.frame(diff_time=c(eval_fit$NA19098$fit.supervised$diff_time,
                                 eval_fit$NA18511$fit.supervised$diff_time,
                                 eval_fit$NA18870$fit.supervised$diff_time,
                                 eval_fit$NA19101$fit.supervised$diff_time,
                                 eval_fit$NA18855$fit.supervised$diff_time,
                                 eval_fit$NA19160$fit.supervised$diff_time),
                     ind=c(rep("NA19098", length(eval_fit$NA19098$fit.supervised$diff_time)),
                           rep("NA18511", length(eval_fit$NA18511$fit.supervised$diff_time)),
                           rep("NA18870", length(eval_fit$NA18870$fit.supervised$diff_time)),
                           rep("NA19101", length(eval_fit$NA19101$fit.supervised$diff_time)),
                           rep("NA18855", length(eval_fit$NA18855$fit.supervised$diff_time)),
                           rep("NA19160", length(eval_fit$NA19160$fit.supervised$diff_time))),
                     stringsAsFactors = F)
df_tmp$ind <- factor(df_tmp$ind,
                     levels=rev(c("NA19098", "NA18511", "NA18855", "NA19101", "NA18870", "NA19160")))
library(ggplot2)
ggplot(df_tmp,
       aes(x=ind, y=diff_time/2/pi, fill=ind)) +
  geom_violin(alpha=.5) + geom_boxplot(width=.1, fill="gray90") +
  scale_fill_brewer(palette="Dark2") +
  ylab("Prediction error") + xlab("Individual") +
  geom_hline(yintercept=.25, col="red")



cbind(mean(eval_fit$NA19098$fit.supervised$diff_time/2/pi),
      sd(eval_fit$NA19098$fit.supervised$diff_time/2/pi)/sqrt(length(eval_fit$NA19098$fit.supervised$diff_time)))
cbind(mean(eval_fit$NA19160$fit.supervised$diff_time/2/pi),
      sd(eval_fit$NA19160$fit.supervised$diff_time/2/pi)/sqrt(length(eval_fit$NA19160$fit.supervised$diff_time)))

inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")
ngene <- 5
eval_fit <- lapply(1:length(inds), function(i) {
  ind <- inds[i]
  fl_name <- list.files("data/results/finalizing",
                        pattern=paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"),
                        full.names = TRUE)
  df <- readRDS(fl_name)
  return(df$fit.test)
})
names(eval_fit) <- inds


eval_fit_mixed <- lapply(1:5, function(i) {
  fl_name <- file.path("data/results",
                        paste0("results_train.fold.",i,".top",ngene,".rds"))
  df <- readRDS(fl_name)
  return(df$fit.test)
})
names(eval_fit_mixed) <- c(1:5)

for (i in 1:length(eval_fit_mixed)) {
  print(mean(eval_fit_mixed[[i]]$fit.supervised$diff_time)/2/pi)
}

mat_pval <- matrix(0, ncol=length(inds), nrow=length(inds))
rownames(mat_pval) <- names(eval_fit)
colnames(mat_pval) <- names(eval_fit)

for (i in 1:length(inds)) {
  for (j in 1:length(inds)) {
    res <- wilcox.test(eval_fit[[i]]$fit.supervised$diff_time/2/pi,
                       eval_fit[[j]]$fit.supervised$diff_time/2/pi)
    mat_pval[i,j] <- res$p.value
  }
}


ind_vs_mixed <- sapply(1:length(inds), function(i) {
  res <- wilcox.test(eval_fit[[i]]$fit.supervised$diff_time/2/pi,
                            eval_fit_mixed$`1`$fit.supervised$diff_time/2/pi)
  return(res$p.value)
})

mat_pval <- rbind(mat_pval, ind_vs_mixed)
which(mat_pval < .01, arr.ind = T)

ind_vs_seurat <- sapply(1:length(inds), function(i) {
  res <- wilcox.test(eval_fit[[i]]$fit.supervised$diff_time/2/pi,
                     eval_fit[[i]]$fit.seurat$diff_time/2/pi)
  return(res$p.value)
})
mat_pval <- rbind(mat_pval, ind_vs_seurat)
which(mat_pval < .01, arr.ind = T)

df_plot <- rbind(data.frame(diff_time=eval_fit$NA19098$fit.supervised$diff_time,
                      ind="NA19098",
                      method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA19098$fit.seurat$diff_time,
                            ind="NA19098",
                            method="seurat",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18511$fit.supervised$diff_time,
                            ind="NA18511",
                            method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18511$fit.seurat$diff_time,
                            ind="NA18511",
                            method="seurat",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18870$fit.supervised$diff_time,
                            ind="NA18870",
                            method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18870$fit.seurat$diff_time,
                            ind="NA18870",
                            method="seurat",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA19101$fit.supervised$diff_time,
                            ind="NA19101",
                            method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA19101$fit.seurat$diff_time,
                            ind="NA19101",
                            method="seurat",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18855$fit.supervised$diff_time,
                            ind="NA18855",
                            method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA18855$fit.seurat$diff_time,
                            ind="NA18855",
                            method="seurat",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA19160$fit.supervised$diff_time,
                            ind="NA19160",
                            method="supervised",stringsAsFactors = F),
                 data.frame(diff_time=eval_fit$NA19160$fit.seurat$diff_time,
                            ind="NA19160",
                            method="seurat",stringsAsFactors = F) )

library(ggplot2)
df_plot$ind <- factor(df_plot$ind, levels=c("NA19098", "NA18511",
                                            "NA18855", "NA19101", "NA18870", "NA19160"))
df_plot$method <- factor(df_plot$method, levels=c("supervised", "seurat"))
ggplot(df_plot, aes(x=method,y=diff_time/2/pi, fill=ind)) +
#  geom_violin(alpha=.5) +
  geom_boxplot(col="gray30") + labs(fill="individual") +
  scale_fill_brewer(palette="Dark2") +
  ylab("Prediction error") + xlab("") +
  geom_hline(yintercept=.25, col="red")




par(mfcol=c(1,2), mar=c(5,4,2,1))
plot(x=eval_fit$NA19098$fit.supervised$ref_time,
     y=eval_fit$NA19098$fit.supervised$pred_time_shift, col="gray50",
     ylab="Predicted phase",
     xlab="FUCCI phase", cex=.5, pch=16, axes=F)
axis(2,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                             expression(3*pi/2), expression(2*pi)))
axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                             expression(3*pi/2), expression(2*pi)))
abline(a=0,b=1, col="black", lwd=.5)
plot(x=eval_fit$NA19160$fit.supervised$ref_time,
     y=eval_fit$NA19160$fit.supervised$pred_time_shift, col="gray50",
     ylab="Predicted phase",
     xlab="FUCCI phase", cex=.5, pch=16, axes=F)
axis(2,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                             expression(3*pi/2), expression(2*pi)))
axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                             expression(3*pi/2), expression(2*pi)))
abline(a=0,b=1, col="black", lwd=.5)


# get results from mixed indivdiuals ----------------------------------------
#mixed <- readRDS("data/results/results_train.fold.1.top5.rds")

ngenes <- c(5, seq(10,500, by=10))

eval_res_mixed <- do.call(rbind, lapply(c(1:6), function(j) {

  foo <- do.call(rbind, lapply(1:length(ngenes), function(i) {
    ngene <- ngenes[i]
    #  train_topX <- do.call(rbind, lapply(1:5, function(fold) {
    fl_name <- file.path(paste0("data/results/finalizing/mixed_results.overallcyclical.fold.",
                                j,".top",ngene,".rds"))
    df <- readRDS(fl_name)
    out <- diff_time_wrapper(df$fit.test)
    out$ngenes <- ngene
    return(out)
  }) )
  foo$fold <- j
  return(foo)
}) )


library(ggplot2)
eval_res_mixed$fold <- factor(eval_res$fold, levels=c(1:6))
ggplot(subset(eval_res_mixed,methods=="supervised"),
       aes(x=ngenes, y=diff_mean, group=fold)) +
  geom_point(aes(color=fold), alpha=.5) + #geom_line(lty=3) +
  ylab("Prediction error (% circle)") + xlab("Top X cyclical genes")  +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=fold), width=.2) +
  stat_summary(fun.y=mean,geom="line",lwd=.5, group=1) +
  ylim(0,.30) +
  geom_hline(yintercept=.25, col="red") +
  theme_light() +
  scale_color_brewer(palette="Spectral")



eval_res_mixed[eval_res_mixed$methods=="supervised",]




# sample size and prediction error
ns_error <- do.call(rbind, lapply(1:length(eval_fit), function(i) {
  data.frame(ind=names(eval_fit)[i],
        error=mean(eval_fit[[i]]$fit.supervised$diff_time),
        ns=ncol(eval_fit[[i]]$fit.supervised$Y_reordered),
        stringsAsFactors = F)
}))
cols <- RColorBrewer::brewer.pal(9,"Dark2")
par(mfrow=c(1,1))
plot(y=(ns_error$error/2/pi),x=ns_error$ns, pch=16,
     ylab="Prediction error",
     xlab="Sample size in test sample", axes=F,
     xlim=c(80, 220), col=cols[1:6])
axis(1); axis(2)
text(labels=ns_error$ind[1],
     y=(ns_error$error[1]/2/pi),
     x=ns_error$ns[1], pos=4, offset=1, col=cols[1])
text(labels=ns_error$ind[2],
     y=(ns_error$error[2]/2/pi),
     x=ns_error$ns[2], pos=3, offset=1, col=cols[2])
text(labels=ns_error$ind[3],
     y=(ns_error$error[3]/2/pi),
     x=ns_error$ns[3], pos=3, offset=1, col=cols[3])
text(labels=ns_error$ind[4],
     y=(ns_error$error[4]/2/pi),
     x=ns_error$ns[4], pos=1, offset=1, col=cols[4])
text(labels=ns_error$ind[5],
     y=(ns_error$error[5]/2/pi),
     x=ns_error$ns[5], pos=3, offset=1, col=cols[5])
text(labels=ns_error$ind[6],
     y=(ns_error$error[6]/2/pi),
     x=ns_error$ns[6], pos=1, offset=1, col=cols[6])


#  compute proportion of variance explained
source("peco/R/utility.R")
for (i in 1:length(inds)) {
  yy <- eval_fit[[i]]$fit.supervised$Y_reordered[1,]
  trend.yy <- eval_fit[[i]]$fit.supervised$mu_reordered[1,]
  ss_total <- (length(yy)-1)*var(yy)
  ss_error <- (length(yy)-1)*var(yy-trend.yy)
  pve <- 1-ss_error/ss_total
  print(cbind(inds[i],ncol(eval_fit[[i]]$fit.supervised$Y_reordered),
              pve, pve/ncol(eval_fit[[i]]$fit.supervised$Y_reordered)))
}
# PVE: NA18511>NA18855>NA19098>NA19160>NA18870>NA19101
# Prediction error (best to worst): NA19098>NA18855>NA18511>NA19101>NA18870>NA19160
# PVE/n: NA18511>NA19160>NA19101>NA18855>NA19098>NA18870













# make figure 2
eval_res_sub2 <- subset(eval_res,
                        (ngenes <= 10) &( methods == "supervised"|methods=="seurat"))
eval_res_sub2 <- subset(eval_res_sub2,
                        !(methods == "seurat" & (ngenes > 5)))
eval_res_sub2$predictor <- NULL

eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==5] <- "peco 5 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="supervised" & eval_res_sub2$ngenes ==10] <- "peco 10 genes"
eval_res_sub2$predictor[eval_res_sub2$methods=="seurat" & eval_res_sub2$ngenes ==5] <- "seurat 97 genes"

eval_res_sub2$predictor <- factor(eval_res_sub2$predictor)
eval_res_sub2$predictor <- factor(eval_res_sub2$predictor,
                                  levels(eval_res_sub2$predictor)[c(2,1,3)])

aggregate(diff_mean ~ predictor, data=eval_res_sub2, mean)

ggplot(eval_res_sub2,
       aes(x=predictor, y=diff_mean, group=ind)) +
  geom_point(aes(color=ind), alpha=.7) +
  geom_errorbar(aes(ymin=diff_mean-diff_se,
                    ymax=diff_mean+diff_se, col=ind), width=.2) +
#  stat_summary(fun.y=mean,geom="point",shape="_",stroke=10, group=1) +
  ylim(0,.30) +
  geom_hline(yintercept=.25, col="blue") +
  ylab("Predictor error (percent circle)") + xlab("Predictors") +
  theme_light() +
  theme(axis.text.x = element_text(angle=35, vjust=.5, hjust=.3))



# Seurat to assign cell-cycle phases ---------------------------------------------------

#ngenes <- c(5, seq(10,500, by=10))
#genes_list <- readRDS("data/results/finalizing/ind_NA18511_cyclical_genes.rds")
inds <- c("NA19098","NA18511","NA18870","NA19101","NA18855","NA19160")
ngenes <- c(5,10)
eval_seurat <- lapply(1:length(inds), function(j) {
  ind <- inds[j]
  foo <- lapply(1:length(ngenes), function(i) {
    ngene <- ngenes[i]
    #  train_topX <- do.call(rbind, lapply(1:5, function(fold) {
    fl_name <- list.files("data/results/finalizing",
                          pattern=paste0("ind_",ind,"_results_overallcyclical.top",ngene,".rds"),
                          full.names = TRUE)
    df <- readRDS(fl_name)
    return(df$fit.test)
  })
  names(foo) <- ngenes
  return(foo)
})
names(eval_seurat) <- inds


par(mfrow=c(1,1), mar=c(5,4,3,1))
plot(eval_seurat[[1]]$`5`$fit.seurat$S,
     eval_seurat[[1]]$`5`$fit.seurat$G2M,
     xlab="S score",
     ylab="G2M score",
     pch=16, cex=.5, axes=F, xlim=c(-1.5,1.5), ylim=c(-1.5,1.5),
     col=c("darkgoldenrod1", "coral", "darkred")[eval_seurat[[1]]$`5`$fit.seurat$assignments])
axis(1); axis(2)
abline(v=0,a=0,b=1,h=0, lty=2, lwd=.5)

par(mfrow=c(1,1), mar=c(5,4,3,1))
seurat.pca <- prcomp(cbind(eval_seurat[[1]]$`5`$fit.seurat$G2M,
                           eval_seurat[[1]]$`5`$fit.seurat$S), scale=TRUE)
seurat.cell_times_est <- as.numeric(coord2rad(cbind(seurat.pca$x[,1],seurat.pca$x[,2])))
(seurat.pca$sdev^2)/sum(seurat.pca$sdev^2)
plot(seurat.pca$x[,1],seurat.pca$x[,2], pch=16, cex=.5,
     xlim=c(-4,4),ylim=c(-4,4),
     col=c("darkgoldenrod1", "coral", "darkred")[eval_seurat[[1]]$`5`$fit.seurat$assignments],
     axes=F, xlab="PC1 (83%)", ylab="PC2 (17%)")
axis(1); axis(2)
abline(v=0,h=0, lty=2, lwd=.5)
par(new=TRUE)
library(circular)
plot(as.circular(eval_seurat[[1]]$`5`$fit.seurat$cell_times_est[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G1"]),
     stack=T, shrink=1, cex=.5, col="darkgoldenrod1")
points(as.circular(eval_seurat[[1]]$`5`$fit.seurat$cell_times_est[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="S"]),
  stack=T, shrink=1, cex=.5, col="coral")
points(as.circular(eval_seurat[[1]]$`5`$fit.seurat$cell_times_est[
  eval_seurat[[1]]$`5`$fit.seurat$assignments=="G2M"]),
  stack=T, shrink=1, cex=.5, col="darkred", bins=200)


# seurat phase predicted vs our predicted
source("peco/R/cycle.corr.R")
seurat_rotate <- rotation(eval_seurat[[1]]$`5`$fit.seurat$ref_time,
                          eval_seurat[[1]]$`5`$fit.seurat$pred_time)
peco_rotate <- rotation(eval_seurat[[1]]$`5`$fit.supervised$ref_time,
                        eval_seurat[[1]]$`5`$fit.supervised$pred_time)

par(mfrow=c(1,2), mar=c(5,4,2,1))
plot(x=eval_seurat[[1]]$`5`$fit.seurat$ref_time,
     y=seurat_rotate,axes=F, xlab="Fucci phase", ylab="Seurat phase",
     pch=16,cex=.5, col="gray50")
abline(0,1,lty=1, lwd=.5,col="black")
axis(2,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c("0",expression(pi/2), expression(pi),
              expression(3*pi/2), expression(2*pi)))
axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c("0",expression(pi/2), expression(pi),
              expression(3*pi/2), expression(2*pi)))
plot(x=eval_seurat[[1]]$`5`$fit.supervised$ref_time,
     y=peco_rotate, xlab="Fucci phase", ylab="Predicted phase", axes=F,
     pch=16,cex=.5, col="gray50")
abline(0,1,lty=1, lwd=.5,col="black")
axis(2,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c("0",expression(pi/2), expression(pi),
              expression(3*pi/2), expression(2*pi)))
axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
     labels=c("0",expression(pi/2), expression(pi),
              expression(3*pi/2), expression(2*pi)))

# par(mfrow=c(2,1))
# hist(peco_rotate, axes=F, nclass=30)
# axis(2);
# axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#      labels=c("0",expression(pi/2), expression(pi),
#               expression(3*pi/2), expression(2*pi)))
# hist(seurat_rotate, axes=F, nclass=30)
# axis(2);
# axis(1,at=c(0,pi/2, pi, 3*pi/2, 2*pi),
#      labels=c("0",expression(pi/2), expression(pi),
#               expression(3*pi/2), expression(2*pi)))
#

# compare seurat and peco fit
fits_seurat <- lapply(1:5, function(i) {
  fit_g <- fit.trendfilter.generic(eval_seurat[[1]]$`5`$fit.supervised$Y[i,
                                 order(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift)])
  # fun_g <- approxfun(x=as.numeric(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
  #   order(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift)]),
  #                    y=as.numeric(fit_g$trend.yy), rule=2)
  # pred_g <- fun_g(as.numeric(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
  #   order(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift)]))
  # return(pred_g)
  return(fit_g$trend.yy)
})
names(fits_seurat) <- rownames(eval_seurat[[1]]$`5`$fit.supervised$Y)





# fucci intensity and expression by phase
par(mfrow=c(2,5), mar=c(3,4,3,1))
#labs <- rownames(eval_seurat[[1]]$`5`$fit.supervised$Y_reordered)
genes <- c("ENSG00000170312","ENSG00000175063",
          "ENSG00000131747", "ENSG00000198518", "ENSG00000197061")
labs <- c("CDK1", "UBE2C", "TOP2A", "HIST1H4E", "HIST1H4C")
for (i in 1:5) {
  if (i>1) { ylab <- ""} else { ylab <- "Normalized log2CPM"}
  ii <- which(rownames(eval_seurat[[1]]$`5`$fit.supervised$Y_reordered) == genes[i])
  plot(x=eval_seurat[[1]]$`5`$fit.supervised$cell_times_reordered,
       y=eval_seurat[[1]]$`5`$fit.supervised$Y_reordered[ii,], col="gray50",
       xlab="Fucci phase",
       ylab=ylab, axes=F, cex=.5, pch=16,
       main = labs[i])
  points(x=eval_seurat[[1]]$`5`$fit.supervised$cell_times_reordered,
         y=eval_seurat[[1]]$`5`$fit.supervised$mu_reordered[ii,], col="blue", pch=16, cex=.7)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
}
for (i in 1:5) {
  if (i>1) { ylab <- ""} else { ylab <- "Normalized log2CPM"}
  tim <- eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift[
    order(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift)]
  ii <- which(rownames(eval_seurat[[1]]$`5`$fit.supervised$Y_reordered) == genes[i])
  plot(x=tim,
       y=eval_seurat[[1]]$`5`$fit.supervised$Y[ii,
         order(eval_seurat[[1]]$`5`$fit.seurat$pred_time_shift)], col="gray50",
       xlab="Fucci phase",
       ylab=ylab, axes=F, cex=.5, pch=16,
       main = labs[i])
  points(x=tim,
         y=fits_seurat[[ii]], col="blueviolet", pch=16, cex=.7)
  axis(2);
  axis(1,at=c(0,pi/2,pi,3*pi/2,2*pi), labels=c(0,expression(pi/2), expression(pi),
                                               expression(3*pi/2), expression(2*pi)))
  abline(h=0, col="black", lty=1, lwd=.7)
}






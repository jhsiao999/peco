[![Travis-CI Build
Status](https://travis-ci.com/jhsiao999/peco.svg?branch=master)](https://travis-ci.com/jhsiao999/peco)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/jhsiao999/peco?branch=master&svg=true)](https://ci.appveyor.com/project/jhsiao999/peco)
[![CircleCI build
status](https://circleci.com/gh/jhsiao999/peco.svg?style=svg)](https://circleci.com/gh/jhsiao999/peco)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# peco

**peco** is an R package for **P**r**E**dicting **C**ell cycle
pr**O**gression in a continuum using scRNA-seq data.

## Installation

To install and load the package, run:

``` r
install.packages("devtools")
library(devtools)
install_github("jhsiao999/peco")
```

`peco` uses `SingleCellExperiment` class objects.

``` r
library(peco)
library(SingleCellExperiment)
library(doParallel)
library(foreach)
```

## Overview

`peco` is a supervised approach for PrEdicting cell cycle phase in a
COntinuum using single-cell RNA sequencing data. The R package provides
functions to build training dataset and also functions to use existing
training data to predict cell cycle on a continuum.

Our work demonstrated that peco is able to predict continuous cell cylce
phase using a small set of cylcic genes: *CDK1*, *UBE2C*, *TOP2A*,
*HISTH1E*, and *HISTH1C* (identified as cell cycle marker genes in
studies of yeast (\[Spellman et al., 1998\]\[spellman\]) and HeLa cells
(\[Whitfield et al., 2002\]\[whitfield\])).

Below we provide two use cases. Vignette 1 shows how to use the
built-training dataset to predict continuous cell cycle. Vignette 2
shows how to make a training datast and build a predictor using training
data.

Users can also view the vigenettes via `browseVignettes("peco")`.

## Vignette

### Training data

`training_human` stores built-in training data of 101 significant cyclic
genes. Below are the slots contained in `training_human`:

  - `predict.yy`: a gene by sample matrix (101 by 888) that stores
    predict cyclic expression values.
  - `cellcycle_peco_reordered`: cell cycle phase in a unit circle
    (angle), ordered from 0 to 2\(pi\)
  - `cellcycle_function`: lists of 101 function corresponding to the top
    101 cyclic genes identified in our dataset
  - `sigma`: standard error associated with cyclic trends of gene
    expression
  - `pve`: proportion of variance explained by the cyclic trend

<!-- end list -->

``` r
data("training_human")
```

`peco` is integrated with `SingleCellExperiment` object in Bioconductor.
Below shows an example of inputting `SingleCellExperiment` object to
perform cell cycle phase prediction.

`sce_top101genes` includes 101 genes and 888 single-cell samples and one
assay slot of `counts`.

``` r
data("sce_top101genes")
assays(sce_top101genes)
```

    ## List of length 1
    ## names(1): counts

Transform the expression values to quantile-normalizesd
counts-per-million values. `peco` uses the `cpm_quantNormed` slot as
input data for predictions.

``` r
sce_top101genes <- data_transform_quantile(sce_top101genes)
```

    ## computing on 2 cores

``` r
assays(sce_top101genes)
```

    ## List of length 3
    ## names(3): counts cpm cpm_quantNormed

### Predict cell cycle phase

Apply the prediction model using function `cycle_npreg_outsample`.

``` r
model_101genes_predict <- cycle_npreg_outsample(
    Y_test=sce_top101genes,
    sigma_est=training_human$sigma[rownames(sce_top101genes),],
    funs_est=training_human$cellcycle_function[rownames(sce_top101genes)],
    method.trend="trendfilter",
    ncores=1,
    get_trend_estimates=FALSE)
```

`peco` adds the predict cell cycle phase to the `colData` slot.
sce\_top101genes

``` r
head(colData(model_101genes_predict)$cellcycle_peco)
```

    ## 20170905-A01 20170905-A02 20170905-A03 20170905-A06 20170905-A07 
    ##     1.099557     4.680973     2.670354     4.303982     4.052655 
    ## 20170905-A08 
    ##     1.413717

Visualize results of prediction for one gene. Below we choose CDK1
(“ENSG00000170312”). Because CDK1 is a known cell cycle gene, this
visualization serves as a sanity check for the results of fitting. The
fitted function `training_human$cellcycle_function[[1]]` was obtained
from our training data.

``` r
plot(y=assay(sce_top101genes,"cpm_quantNormed")["ENSG00000170312",],
     x=colData(sce_top101genes)$theta_shifted, main = "CDK1",
     ylab = "quantile normalized expression")
points(y=training_human$cellcycle_function[["ENSG00000170312"]](seq(0,2*pi, length.out=100)),
       x=seq(0,2*pi, length.out=100), col = "blue", pch =16)
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Visualize cyclic expression trend

Visualize results of prediction for the top 10 genesone genes. Use
`fit_cyclical_many` to estimate cyclic function based on the input data.

``` r
# predicted cell time in the input data
theta_predict = colData(sce_top101genes)$theta_shifted
names(theta_predict) = rownames(colData(sce_top101genes))

# expression values of 10 genes in the input data
yy_input = assay(sce_top101genes,"cpm_quantNormed")[1:6,]

# apply trendfilter to estimate cyclic gene expression trend
fit_cyclic <- fit_cyclical_many(Y=yy_input, 
                                theta=theta_predict)
```

    ## computing on 2 cores

``` r
par(mfrow=c(2,3))
for (i in 1:6) {
plot(y=yy_input[i,],
     x=fit_cyclic$cellcycle_peco_ordered, main = names(fit_cyclic)[i],
     ylab = "quantile normalized expression")
points(y=fit_cyclic$cellcycle_function[[i]](seq(0,2*pi, length.out=100)),
       x=seq(0,2*pi, length.out=100), col = "blue", pch =16)
}
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

-----

#### Contact

Please contact me at [joyce.hsiao1@gmail.com](joyce.hsiao1@gmail.com)
for questions on the package or the methods.

#### How to cite

> Hsiao, C. J., Tung, P., Blischak, J. D., Burnett, J., Dey, K. K.,
> Barr, A. K., Stephens, M., and Gilad, Y. (2018). [Characterizing and
> inferring quantitative cell-cycle phase in single-cell RNA-seq data
> analysis.](https://doi.org/10.1101/526848) bioRxiv
> <doi:10.1101/526848>

#### License

Copyright (c) 2019-2020, Joyce Hsiao.

All source code and software in this repository are made available under
the terms of the [GNU General Public
License](https://www.gnu.org/licenses/gpl-3.0.en.html). See file
[LICENSE](LICENSE) for the full text of the license.

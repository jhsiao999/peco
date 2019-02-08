# peco

**peco** is a R package for predicting cell cycle progression in a
continuum using scRNA-seq data.

## Installation 

To install the package, follow the commands below:

```
install.packages("devtools")
library(devtools)
install_github('jhsiao999/peco')
```

To load the package

```
library(peco)
```

## Vignettes

I've included an example of using `peco` to predict cell cycle phase
using single-cell RNA-seq data. Run the following command to view the
vignette.

```
browseVignettes("peco")
```

## Contact

Please contact me at [joyce.hsiao1@gmail.com](joyce.hsiao1@gmail.com)
for questions on the package or the methods.

## How to cite

> Hsiao, C. J., Tung, P., Blischak, J. D., Burnett, J., Dey, K. K. ,
> Barr, K., Stephens, M., and Gilad, Y. (2018). Characterizing and
> inferring quantitative cell-cycle phase in single-cell RNA-seq data
> analysis. bioRxiv

## Developer notes

Alternatively, to install and test the peco package, run the following
commands in the command-line shell:

```bash
R CMD build --resave-data peco
R CMD check --as-cran peco_0.1.3.tar.gz
R CMD INSTALL peco_0.1.3.tar.gz
```

## Licenses

This package is distributed under [GPL - General Public License (>= 2)]

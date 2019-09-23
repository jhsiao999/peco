# peco

[![Travis-CI Build Status](https://travis-ci.com/jhsiao999/peco.svg?branch=master)](https://travis-ci.com/jhsiao999/peco)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/jhsiao999/peco?branch=master&svg=true)](https://ci.appveyor.com/project/jhsiao999/peco)
[![CircleCI build status](https://circleci.com/gh/jhsiao999/peco.svg?style=svg)](https://circleci.com/gh/jhsiao999/peco)

**peco** is an R package for **P**r**E**dicting **C**ell cycle pr**O**gression in a
continuum using scRNA-seq data.

## Installation 

To install and load the package, run:

```R
install.packages("devtools")
library(devtools)
install_github("jhsiao999/peco")
library(peco)
```

## Vignette

I've included an example of using `peco` to predict cell cycle phase
using single-cell RNA-seq data. Run the following command to view the
vignette.

```R
browseVignettes("peco")
```

## Contact

Please contact me at [joyce.hsiao1@gmail.com](joyce.hsiao1@gmail.com)
for questions on the package or the methods.

## How to cite

> Hsiao, C. J., Tung, P., Blischak, J. D., Burnett, J., Dey, K. K.,
> Barr, A. K., Stephens, M., and Gilad, Y. (2018). [Characterizing and
> inferring quantitative cell-cycle phase in single-cell RNA-seq data
> analysis.](https://doi.org/10.1101/526848) bioRxiv doi:10.1101/526848

## License

Copyright (c) 2018-2019, Joyce Hsiao.

All source code and software in this repository are made available
under the terms of the [GNU General Public
License](https://www.gnu.org/licenses/gpl-3.0.en.html). See
file [LICENSE](LICENSE) for the full text of the license.

## Developer notes

Alternatively, to install and test the peco package, run the following
commands in the command-line shell:

```bash
R CMD build --resave-data peco
R CMD INSTALL peco_0.99.4.tar.gz
```

<!--- R CMD check --as-cran peco_0.99.0.tar.gz---!>


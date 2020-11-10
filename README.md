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
pr**O**gression in a continuum using scRNA-seq data. **peco**
provides functions to build a training data set and predict cell cycle
on a continuum.

## Installation

To install and load the package, run:

``` r
install.packages("devtools")
devtools::install_github("jhsiao999/peco")
library(peco)
```

For for a detailed illustration of `peco`, see the
[vignette][vignette].

Contact
-------

Please contact Joyce Hsiao at
[joyce.hsiao1@gmail.com](joyce.hsiao1@gmail.com) for questions on the
package or the methods.

How to cite
-----------

> Hsiao, C. J., Tung, P., Blischak, J. D., Burnett, J., Dey, K. K.,
> Barr, A. K., Stephens, M., and Gilad, Y. (2020). Characterizing and
> inferring quantitative cell cycle phase in single-cell RNA-seq data
> analysis. *Genome Biology*, 30(4): 611-621,
> [doi:10.1101/gr.247759.11][paper]

License
-------

Copyright (c) 2019-2020, Joyce Hsiao.

All source code and software in this repository are made available
under the terms of the [GNU General Public License][gpl]. See file
[LICENSE](LICENSE) for the full text of the license.

[gpl]: https://www.gnu.org/licenses/gpl-3.0.en.html
[paper]: doi.org/10.1101/gr.247759.118
[vignette]: https://jhsiao999.github.io/peco/articles/peco.html

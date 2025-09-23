## Dapar Toolshed

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check.yaml](https://github.com/edyp-lab/DaparToolshed/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/edyp-lab/DaparToolshed/actions/workflows/check-standard.yaml)
[![R-CMD-check-bioc](https://github.com/edyp-lab/DaparToolshed/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/edyp-lab/DaparToolshed/actions/workflows/check-bioc.yml)
[![codecov.io](https://codecov.io/github/edyp-lab/DaparToolshed/coverage.svg?branch=master)](https://codecov.io/github/edyp-lab/DaparToolshed?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![CRAN status](https://www.r-pkg.org/badges/version/DaparToolshed)](https://CRAN.R-project.org/package=DaparToolshed)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/edyp-lab/DaparToolshed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/edyp-lab/DaparToolshed/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
> Evolving `DAPAR` package towards Shiny modules and data structures from the 
packages MultiAssayExperiment and SummarizedExperiment




`DaparToolshed` is a package which provides all the necessary functions to 
analyze quantitative data from label-free proteomics experiments bases on the 
MultiAssayExperiment data structure.
Contrarily to most other similar R packages, it is endowed with rich and 
user-friendly graphical interfaces, so that no programming skill is 
required (see `Prostar` package).




### What is DaparToolshed?

`DaparToolshed` is a [Bioconductor
package](http://bioconductor.org/packages/omXplore) that provides all the 
necessary functions to analyze quantitative data from label-free proteomics 
experiments bases on the `MultiAssayExperiment` data structure.
Contrarily to most other similar R packages, it is endowed with rich and 
user-friendly graphical interfaces, so that no programming skill is 
required (see `Prostar` package).



### Getting started

See the
[DaparToolshed introduction](https://edyp-lab.github.io/DaparToolshed/articles/DaparToolshed.html)
to get started.



### License

The `DaparToolshed` code is provided under a permissive [Artistic 2.0
license](https://opensource.org/licenses/Artistic-2.0). The
documentation, including the manual pages and the vignettes, are
distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/4.0/).


# Installation

To install this package, start R (version "4.3") and enter:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DaparToolshed")
```

This will also install dependencies.

It is also possible to install `DaparToolshed` from Github:

```
library(devtools)
install_github('edyp-lab/DaparToolshed')

```

For older versions of R, please refer to the appropriate Bioconductor release.



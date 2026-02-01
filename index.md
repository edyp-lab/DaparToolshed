# Installation

## Dapar Toolshed

`DaparToolshed` is a package which provides all the necessary functions
to analyze quantitative data from label-free proteomics experiments
bases on the MultiAssayExperiment data structure. Contrarily to most
other similar R packages, it is endowed with rich and user-friendly
graphical interfaces, so that no programming skill is required (see
`Prostar` package).

### What is DaparToolshed?

`DaparToolshed` is a [Bioconductor
package](http://bioconductor.org/packages/omXplore) that provides all
the necessary functions to analyze quantitative data from label-free
proteomics experiments bases on the `MultiAssayExperiment` data
structure. Contrarily to most other similar R packages, it is endowed
with rich and user-friendly graphical interfaces, so that no programming
skill is required (see `Prostar` package).

### Getting started

See the [DaparToolshed
introduction](https://edyp-lab.github.io/DaparToolshed/articles/DaparToolshed.html)
to get started.

### License

The `DaparToolshed` code is provided under a permissive [Artistic 2.0
license](https://opensource.org/licenses/Artistic-2.0). The
documentation, including the manual pages and the vignettes, are
distributed under a [CC BY-SA
license](https://creativecommons.org/licenses/by-sa/4.0/).

To install this package, start R (version “4.3”) and enter:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("DaparToolshed")

This will also install dependencies.

It is also possible to install `DaparToolshed` from Github:

    library(devtools)
    install_github('edyp-lab/DaparToolshed')

For older versions of R, please refer to the appropriate Bioconductor
release.

# Filter a peptide assay on the basis of its adjacency matrix.

These functions filters (delete) peptides of an assay, applying a
function on peptides and proteins. They can be used alone but the usual
usage is to create an instance of a class FunctionFilter and to pass it
to the function filterFeaturesOneSE in order to create a new assay,
embedded into the QFeatures object.

## Usage

``` r
AdjMatFilters()

allPeptides(object, ...)

specPeptides(object, ...)

subAdjMat_specificPeptides(X)

sharedPeptides(object, ...)

subAdjMat_sharedPeptides(X)

topnFunctions()

topnPeptides(object, fun, top)

subAdjMat_topnPeptides(X, qData, fun, top)
```

## Arguments

- object:

  An object of class `SummarizedExperiment`

- ...:

  Additional arguments

- X:

  xxx

- fun:

  A [`list()`](https://rdrr.io/r/base/list.html) of additional
  parameters

- top:

  A `integer(1)` which is the number of xxx

- qData:

  xxx

## Value

NA

## Details

This function builds an intermediate matrix with scores for each peptide
based on 'fun' parameter. Once this matrix is built, one select the 'n'
peptides which have the higher score

The list of filter functions is given by `adjMatFilters()`:

- `specPeptides()`: returns a new assay of class `SummazizedExperiment`
  with only specific peptides;

- `sharedpeptides()`: returns a new assay of class
  `SummazizedExperiment` with only shared peptides;

- `opnPeptides()`: returns a new assay of class `SummazizedExperiment`
  with only the 'n' peptides which best satisfies the condition. The
  condition is represented by functions which calculates a score for
  each peptide among all samples. The list of these functions is given
  by `topnFunctions()`:

- `rowMedians()`: xxx;

- [`rowMeans()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html):
  xxx;

- [`rowSums()`](https://rdrr.io/pkg/Matrix/man/colSums-methods.html):
  xxx;

## See also

The QFeatures-filtering-oneSE man page for the class `FunctionFilter`.

## Author

Samuel Wieczorek

## Examples

``` r
library(Matrix)
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
library(QFeatures)
#------------------------------------------------
# This function will keep only specific peptides
#------------------------------------------------

f1 <- FunctionFilter("specPeptides", list())

#------------------------------------------------
# This function will keep only shared peptides
#------------------------------------------------

f2 <- FunctionFilter("sharedPeptides", list())

#------------------------------------------------
# This function will keep only the 'n' best peptides
# w.r.t the quantitative sum of each peptides among
# all samples
#------------------------------------------------

f3 <- FunctionFilter("topnPeptides", fun = "rowSums", top = 2)

#------------------------------------------------------
# To run the filter(s) on the dataset, use [xxx()]
# IF several filters must be used, store them in a list
#------------------------------------------------------

data(subR25pept)
lst.filters <- list()
lst.filters <- append(lst.filters, f1)
lst.filters <- append(lst.filters, f3)

subR25prot <- filterFeaturesOneSE(
    object = subR25pept,
    i = 1,
    name = "filtered",
    filters = lst.filters
)
```

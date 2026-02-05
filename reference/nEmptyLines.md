# Number of empty lines in the data

This function counts the number of empty lines (all elements are equal
to NA).

## Usage

``` r
nEmptyLines(df)
```

## Arguments

- df:

  A `data.frame`.

## Value

A `integer(1)`

## Author

Samuel Wieczorek

## Examples

``` r
library(QFeatures)
data(subR25prot)
nEmptyLines(SummarizedExperiment::assay(subR25prot, 1))
#> [1] 0
```

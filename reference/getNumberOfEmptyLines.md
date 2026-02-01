# Returns the number of empty lines in the data

Returns the number of empty lines in a matrix.

## Usage

``` r
getNumberOfEmptyLines(qData)
```

## Arguments

- qData:

  A matrix corresponding to the quantitative data.

## Value

An integer

## Author

Samuel Wieczorek

## Examples

``` r
library(QFeatures)
data(subR25prot)
qData <- assay(subR25prot[[1]])
getNumberOfEmptyLines(qData)
#> [1] 0
```

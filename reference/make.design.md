# Builds the design matrix

Builds the design matrix

## Usage

``` r
make.design(sTab)
```

## Arguments

- sTab:

  The data.frame which correspond to the
  [`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html) function
  of package `MSnbase`.

## Value

A design matrix

## Author

Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
make.design(SummarizedExperiment::colData(subR25pept))
#>      ConditionA ConditionB
#> [1,]          1          0
#> [2,]          1          0
#> [3,]          1          0
#> [4,]          0          1
#> [5,]          0          1
#> [6,]          0          1
```

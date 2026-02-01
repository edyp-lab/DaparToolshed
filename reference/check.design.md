# Check if the design is valid

Check if the design is valid

## Usage

``` r
check.design(sTab)
```

## Arguments

- sTab:

  The data.frame which correspond to the
  [`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html) function
  of package `MSnbase`.

## Value

A boolean

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
check.design(SummarizedExperiment::colData(subR25pept)[, seq(3)])
#> $valid
#> [1] TRUE
#> 
#> $warn
#> NULL
#> 
```

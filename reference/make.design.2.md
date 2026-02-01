# Builds the design matrix for designs of level 2

Builds the design matrix for designs of level 2

## Usage

``` r
make.design.2(sTab)
```

## Arguments

- sTab:

  The data.frame which correspond to the
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  function of package `MSnbase`.

## Value

A design matrix

## Author

Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
make.design.2(SummarizedExperiment::colData(subR25pept))
#>                Condition1.RepBio1 Condition1.RepBio2 Condition1.RepBio3
#> Intensity_C_R1                  1                  0                  0
#> Intensity_C_R2                  0                  1                  0
#> Intensity_C_R3                  0                  0                  1
#> Intensity_D_R1                  0                  0                  0
#> Intensity_D_R2                  0                  0                  0
#> Intensity_D_R3                  0                  0                  0
#>                Condition2.RepBio4 Condition2.RepBio5 Condition2.RepBio6
#> Intensity_C_R1                  0                  0                  0
#> Intensity_C_R2                  0                  0                  0
#> Intensity_C_R3                  0                  0                  0
#> Intensity_D_R1                  1                  0                  0
#> Intensity_D_R2                  0                  1                  0
#> Intensity_D_R3                  0                  0                  1

```

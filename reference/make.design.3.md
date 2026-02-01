# Builds the design matrix for designs of level 3

Builds the design matrix for designs of level 3

## Usage

``` r
make.design.3(sTab)
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
sTab <- cbind(SummarizedExperiment::colData(subR25pept), Tech.Rep = 1:6)
make.design.3(sTab)
#>                Condition1.RepBio1.RepTech1 Condition1.RepBio2.RepTech2
#> Intensity_C_R1                           1                           0
#> Intensity_C_R2                           0                           1
#> Intensity_C_R3                           0                           0
#> Intensity_D_R1                           0                           0
#> Intensity_D_R2                           0                           0
#> Intensity_D_R3                           0                           0
#>                Condition1.RepBio3.RepTech3 Condition2.RepBio4.RepTech4
#> Intensity_C_R1                           0                           0
#> Intensity_C_R2                           0                           0
#> Intensity_C_R3                           1                           0
#> Intensity_D_R1                           0                           1
#> Intensity_D_R2                           0                           0
#> Intensity_D_R3                           0                           0
#>                Condition2.RepBio5.RepTech5 Condition2.RepBio6.RepTech6
#> Intensity_C_R1                           0                           0
#> Intensity_C_R2                           0                           0
#> Intensity_C_R3                           0                           0
#> Intensity_D_R1                           0                           0
#> Intensity_D_R2                           1                           0
#> Intensity_D_R3                           0                           1

```

# Builds a complete color palette for the conditions given in argument

xxxx

## Usage

``` r
GetColorsForConditions(conds, pal = NULL)
```

## Arguments

- conds:

  The extended vector of samples conditions

- pal:

  A vector of HEX color code that form the basis palette from which to
  build the complete color vector for the conditions.

## Value

NA

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
GetColorsForConditions(design.qf(subR25pept)$Condition)
#> [1] "#E41A1C" "#E41A1C" "#E41A1C" "#377EB8" "#377EB8" "#377EB8"
GetColorsForConditions(design.qf(subR25pept)$Condition, ExtendPalette(2))
#> [1] "#E41A1C" "#E41A1C" "#E41A1C" "#377EB8" "#377EB8" "#377EB8"
```

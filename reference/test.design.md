# Check if xxxxxx

Check if xxxxxx

## Usage

``` r
test.design(tab)
```

## Arguments

- tab:

  A data.frame which correspond to xxxxxx

## Value

A list of two items

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
test.design(SummarizedExperiment::colData(subR25pept)[, seq(3)])
#> $valid
#> [1] FALSE
#> 
#> $warn
#> [1] "The value 25fmol in column 'Condition' is not correctly set.\n"
#> [2] "The value 25fmol in column 'Condition' is not correctly set.\n"
#> [3] "The value 25fmol in column 'Condition' is not correctly set.\n"
#> [4] "The value 10fmol in column 'Condition' is not correctly set.\n"
#> [5] "The value 10fmol in column 'Condition' is not correctly set.\n"
#> [6] "The value 10fmol in column 'Condition' is not correctly set.\n"
#> 
```

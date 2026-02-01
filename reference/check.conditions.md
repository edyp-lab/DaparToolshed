# Check if the design is valid

Check if the design is valid

## Usage

``` r
check.conditions(conds)
```

## Arguments

- conds:

  A `vector` containing the conditions.

## Value

A `list` including : "valid" : Wether the conditions are valid or not.
"warn" : A message describing the issue if the conditions ar not valid.

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
check.conditions(design.qf(subR25pept)$Condition)
#> $valid
#> [1] TRUE
#> 
#> $warn
#> NULL
#> 
```

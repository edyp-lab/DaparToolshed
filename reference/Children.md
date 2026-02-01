# Names of all chidren of a node

xxx

## Usage

``` r
Children(level, parent = NULL)
```

## Arguments

- level:

  A string designing the type of entity/pipeline. Available values are:
  `peptide`, `protein`

- parent:

  xxx

## Value

A vector

## Examples

``` r
Children('protein', 'Missing')
#> [1] "Missing POV" "Missing MEC"
Children('protein', 'Missing POV')
#> character(0)
Children('protein', c('Missing POV', 'Missing MEC'))
#> character(0)
Children('protein', c('Missing', 'Missing POV', 'Missing MEC'))
#> [1] "Missing POV" "Missing MEC"
```

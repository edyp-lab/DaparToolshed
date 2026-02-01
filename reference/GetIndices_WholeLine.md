# Search lines which respects query on all their elements.

This function looks for the lines where each element respect the query.

## Usage

``` r
GetIndices_WholeLine(metacell.mask)
```

## Arguments

- metacell.mask:

  xxx

## Value

xxx

## Examples

``` r
data(subR25pept)
level <- 'peptide'
pattern <- "Missing POV"
metacell.mask <- match.metacell(metadata = qMetacell(subR25pept[[1]]), 
pattern = pattern, level = level)
ind <- GetIndices_WholeLine(metacell.mask)
```

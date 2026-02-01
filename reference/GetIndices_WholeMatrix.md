# Search lines which respects request on one or more conditions.

This function looks for the lines that respect the request in either all
conditions or at least one condition.

## Usage

``` r
GetIndices_WholeMatrix(metacell.mask, op = "==", percent = FALSE, th = 0)
```

## Arguments

- metacell.mask:

  xxx

- op:

  String for operator to use. List of operators is available with
  'SymFilteringOperators()'.

- percent:

  A boolean to indicate whether the threshold represent an absolute
  value (percent = FALSE) or a percentage (percent=TRUE).

- th:

  A floating number which is in the interval `[0, 1]`

## Value

xxx

## Examples

``` r
data(subR25pept)
level <- 'peptide'
pattern <- "Missing"
metacell.mask <- match.metacell(
metadata = qMetacell(subR25pept[[1]]), pattern = pattern, level = level)
percent <- FALSE
th <- 3
op <- ">="
ind <- GetIndices_WholeMatrix(metacell.mask, op, percent, th)
```

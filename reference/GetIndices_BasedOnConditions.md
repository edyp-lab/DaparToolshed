# Search lines which respects request on one or more conditions.

This function looks for the lines that respect the request in either all
conditions or at least one condition.

## Usage

``` r
GetIndices_BasedOnConditions(metacell.mask, type, conds, percent, op, th)
```

## Arguments

- metacell.mask:

  xxx

- type:

  Available values are:

  - 'AllCond' (the query is valid in all the conditions),

  - 'AtLeatOneCond' (the query is valid in at leat one condition.

- conds:

  xxx

- percent:

  xxx

- op:

  String for operator to use. List of operators is available with the
  function 'SymFilteringOperators()'.

- th:

  The theshold to apply

## Value

xxx

## Examples

``` r
data(subR25pept)
level <- typeDataset(subR25pept[[1]])
pattern <- 'Missing'
metacell.mask <- match.metacell(
metadata=qMetacell(subR25pept[[1]]), pattern=pattern, level=level)
type <- 'AllCond'
conds <- design.qf(subR25pept)$Condition
op <- '>='
th <- 0.5
percent <- TRUE
ind <- GetIndices_BasedOnConditions(metacell.mask, type, conds, percent, op, th)
```

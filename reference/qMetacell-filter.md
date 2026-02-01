# Search lines which respects request on one or more conditions.

This function looks for the lines that respect the request in either all
conditions or at least one condition.

## Usage

``` r
SymFilteringOperators()

qMetacellFilteringScope()

qMetacellWholeMatrix(
  object,
  cmd,
  pattern,
  percent = "Percentage",
  th,
  operator
)

qMetacellWholeLine(object, cmd, pattern)

qMetacellOnConditions(
  object,
  cmd,
  mode,
  pattern,
  conds,
  percent = "Percentage",
  operator,
  th
)
```

## Arguments

- object:

  An instance of the class `SummarizedExperiment`

- cmd:

  A `character(1)` indicating the action to perform. Either "keep" or
  "delete".

- pattern:

  A [`character()`](https://rdrr.io/r/base/character.html) indicating
  the tag pattern of interest.

- percent:

  A `character(1)` indicating whether the threshold represent an
  absolute value ("Count") or a percentage ("Percentage").

- th:

  The threshold to apply

- operator:

  String for operator to use. List of operators is available with
  'SymFilteringOperators()'.

- mode:

  A `character(1)` indicating how the task is performed. Either
  "AllCond" or "AtLeastOneCond".

- conds:

  A vector of conditions in the dataset.

## Value

A vector of operators

NA

A vector of filtering scopes

NA

NA

NA

## Examples

``` r
SymFilteringOperators()
#>    lessthan        less greaterthan     greater       equal   different 
#>        "<="         "<"        ">="         ">"        "=="        "!=" 

data(subR25prot)
obj <- subR25prot[[1]]
level <- typeDataset(obj)
pattern <- "Missing"
mask <- match.metacell(
    metadata = qMetacell(obj),
    pattern = pattern,
    level = level
)
percent <- FALSE
th <- 3
op <- ">="
cmd <- 'delete' 
ind <- qMetacellWholeMatrix(obj, cmd, pattern, percent, th, op)

data(subR25prot)
ind <- qMetacellWholeLine(obj, cmd, pattern)

conds <- design.qf(subR25prot)$Condition
op <- ">="
th <- 0.5
percent <- "Percentage"
mode <- "AllCond"
ind <- qMetacellOnConditions(obj, cmd, mode, pattern, conds, percent, op, th)

qMetacellFilteringScope()
#> [1] "None"           "WholeLine"      "WholeMatrix"    "AllCond"       
#> [5] "AtLeastOneCond"

data(subR25prot)
obj <- subR25prot[[1]]
```

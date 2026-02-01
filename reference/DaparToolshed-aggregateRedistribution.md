# Aggregate an assay's quantitative features with shared peptide redistribution

This function aggregates the quantitative features of an assay, applying
a summarization function (`fun`) to sets of features. The `fcol`
variable name points to a rowData column that defines how to group the
features during aggregate. This variable has to be an adjacency matrix.
This function uses
[`DaparToolshed::inner.aggregate.iter()`](https://edyp-lab.github.io/DaparToolshed/reference/DaparToolshed-aggregate.md)
to aggregate quantitative data.

The list of agregation methods can be obtained with the function
[`aggregateMethods()`](https://edyp-lab.github.io/DaparToolshed/reference/DaparToolshed-aggregate.md).
This function compiles both methods from the packages `DaparToolshed`
and `QFeatures`.

## Usage

``` r
aggregateRedistribution(object, ...)

# S4 method for class 'QFeatures'
aggregateRedistribution(
  object,
  i,
  name = "newAssay",
  fcol,
  init.method = "Mean",
  method = "Mean",
  ponderation = "Global",
  n = NULL,
  uniqueiter = FALSE,
  max_iter = 500
)

# S4 method for class 'SummarizedExperiment'
aggregateRedistribution(
  object,
  fcol,
  init.method = "Mean",
  method = "Mean",
  ponderation = "Global",
  n = NULL,
  uniqueiter = FALSE,
  conds,
  max_iter = 500
)
```

## Arguments

- object:

  An instance of class `QFeatures` or `SummarizedExperiment`

- ...:

  Additional parameters.

- i:

  The index or name of the assay which features will be aggregated the
  create the new assay.

- name:

  A `character(1)` naming the new assay. Default is `newAssay`. Note
  that the function will fail if there's already an assay with `name`.

- fcol:

  A `character(1)` naming a rowdata variable (of assay `i` in case of a
  `QFeatures`) defining how to aggregate the features of the assay. This
  variable is a (possibly sparse) matrix. See below for details.

- init.method:

  A function used for initializing the aggregation. Available functions
  are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. See
  [`DaparToolshed::inner.aggregate.iter()`](https://edyp-lab.github.io/DaparToolshed/reference/DaparToolshed-aggregate.md)
  for details.

- method:

  A function used for the aggregation. Available functions are `Sum`,
  `Mean`, `Median` or `medianPolish`. See
  [`DaparToolshed::inner.aggregate.iter()`](https://edyp-lab.github.io/DaparToolshed/reference/DaparToolshed-aggregate.md)
  for details.

- ponderation:

  A `character(1)` defining what to consider to create the coefficient
  for redistribution of shared peptides. Available values are `Global`
  (default), `Condition` or `Sample`.

- n:

  A `numeric(1)` specifying the number of peptides to use for each
  protein. If `NULL`, all peptides are considered.

- uniqueiter:

  A `boolean` indication if there should be only 1 iteration or not.

- max_iter:

  A `numeric(1)` setting the maximum number of iteration.

- conds:

  A [`character()`](https://rdrr.io/r/base/character.html) vector which
  is the names of conditions.

## Value

A `QFeatures` object with an additional assay or a
`SummarizedExperiment` object (or subclass thereof).

## Details

This function uses
[`DaparToolshed::inner.aggregate.iter()`](https://edyp-lab.github.io/DaparToolshed/reference/DaparToolshed-aggregate.md)
to aggregate quantitative data.

## Iterative aggregation function

xxxxxx xxxxx

## Quantitative metadata aggregation

The function to aggregate the quantitative metadata is `aggQmetadat()`

## See also

The *QFeatures* vignette provides an extended example and the
*Aggregation* vignette, for a complete quantitative proteomics data
processing pipeline.

## Examples

``` r
## ---------------------------------------
## An example QFeatures with PSM-level data
## ---------------------------------------
# \donttest{
data(subR25prot)
library(SummarizedExperiment)
subR25prot
#> An instance of class QFeatures (type: bulk) with 2 sets:
#> 
#>  [1] original: SummarizedExperiment with 100 rows and 6 columns 
#>  [2] logAssay: SummarizedExperiment with 100 rows and 6 columns 

## Aggregate peptides into proteins
## using the adjacency matrix
feat1 <- aggregateRedistribution(object = subR25prot,
i = 1,
name = 'aggregated',
fcol = 'adjacencyMatrix',
init.method = 'Mean',
method = 'Mean')
#> Error in .aggregateRedistribution(object, fcol, init.method, method, ponderation,     n, uniqueiter, conds, max_iter): 'fcol' must refer to a sparse matrix.
feat1
#> Error: object 'feat1' not found

assay(feat1[[1]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'feat1' not found
assay(feat1[[2]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'feat1' not found
aggcounts(feat1[[2]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'aggcounts': object 'feat1' not found
assay(feat1[[3]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assay': object 'feat1' not found
aggcounts(feat1[[3]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'aggcounts': object 'feat1' not found
rowData(feat1[[2]])
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'rowData': object 'feat1' not found
# }
```

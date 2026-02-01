# Sets the metacell dataframe for datasets which are from Proline software

In the quantitative columns, a missing value is identified by no value
rather than a value equal to 0.

In these datasets, the metacell info is computed from the 'PSM count'
columns.

Conversion rules Initial conversion rules for proline
\|————–\|—————–\|—–\| \| Quanti \| PSM count \| Tag \|
\|————–\|—————–\|—–\| \| == 0 \| N.A. \| whatever \| 2.0 \| \| \> 0 \|
\> 0 \| 1.1 \| \| \> 0 \| == 0 \| 1.2 \| \| \> 0 \| unknown col \| 1.0
\| \|————–\|—————–\|—–\|

## Usage

``` r
Metacell_proline(qdata, conds, df = NULL, level = NULL)
```

## Arguments

- qdata:

  An object of class `MsnSet`

- conds:

  A 1-col dataframe with the condition associated to each sample.

- df:

  A dataframe with the same dimension as qdata containing the metacell.

- level:

  A string designing the type of entity/pipeline. Available values are:
  `peptide`, `protein`

## Value

NA

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
conds <- design.qf(subR25pept)$Condition
qdata <- SummarizedExperiment::assay(subR25pept[[2]])
df <- Metacell_proline(qdata, conds, level = "peptide")

```

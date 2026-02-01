# Sets the metacell dataframe

Initial conversion rules for maxquant \|————\|———————–\|——–\| \| Quanti
\| Identification \| Tag \| \|————\|———————–\|——–\| \| == 0 \| whatever
\| 2.0 \| \| \> 0 \| 'By MS/MS' \| 1.1 \| \| \> 0 \| 'By matching' \|
1.2 \| \| \> 0 \| unknown col \| 1.0 \| \|————\|———————–\|——–\|

## Usage

``` r
Metacell_maxquant(qdata, conds, df = NULL, level = NULL)
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
df2 <- Metacell_maxquant(qdata, conds, level = "peptide")
```

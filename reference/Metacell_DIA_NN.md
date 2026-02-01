# Sets the metacell dataframe for datasets which are from Dia-NN software

Actually, this function uses the generic function to generate metacell
info

## Usage

``` r
Metacell_DIA_NN(qdata, conds, df, level = NULL)
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
df <- Metacell_DIA_NN(qdata, conds, df, level = "peptide")

```

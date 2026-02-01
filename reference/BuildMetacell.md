# Metacell function

Create metacell

## Usage

``` r
BuildMetacell(from = NULL, level, qdata = NULL, conds = NULL, df = NULL)
```

## Arguments

- from:

  A string designing the software used, either "maxquant", "proline" or
  "DIA-NN"

- level:

  A string designing the type of entity/pipeline. Available values are:
  `peptide`, `protein`

- qdata:

  A matrix of quantitative data

- conds:

  A 1-col dataframe with the condition associated to each sample.

- df:

  A data.frame which contains the type of identification of the
  entities. It must have the same dimensions as `qData`.

## Value

NA

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
conds <- design.qf(subR25pept)$Condition
qdata <- SummarizedExperiment::assay(subR25pept[[2]])
metacell1 <- BuildMetacell('maxquant', 'peptide', qdata, conds)
metacell2 <- BuildMetacell('proline', 'peptide', qdata, conds)
```

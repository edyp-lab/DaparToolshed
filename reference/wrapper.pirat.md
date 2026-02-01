# Missing values imputation using Pirat

This method is a wrapper to the function `pipeline_llkimpute()` of the
package `Pirat` adapted to an object of class `QFeatures` of
`SummarizedExperiment`.

## Usage

``` r
wrapper.pirat(data, adjmat, rnas_ab = NULL, adj_rna_pg = NULL, ...)
```

## Arguments

- data:

  An object of class `QFeatures` or `SummarizedExperiment`. If data is
  of class `QFeatures`, the last assay will be imputed.

- adjmat:

  Adjacency matrix corresponding to the `SummarizedExperiment` or the
  last assay of `QFeatures`.

- rnas_ab:

  Transcriptomic data with sample as row, used only if extension = 'T'.

- adj_rna_pg:

  Adjacency matrix of rna (rows) and peptides or precursors (columns),
  used only if extension = 'T'.

- ...:

  Additional arguments to pass to `my_pipeline_llkimpute()`

## Value

`QFeatures` including a new assay with imputed data or
`SummarizedExperiment` with imputed data.

## Author

Manon Gaudin

## Examples

``` r
data(subR25pept)

# Delete whole empty lines
filter_emptyline <- FunctionFilter("qMetacellWholeLine", cmd = 'delete', pattern = 'Missing MEC')
subR25pept <- filterFeaturesOneSE(object = subR25pept, i = length(subR25pept), name = "Filtered", 
              filters = list(filter_emptyline))

subR25pept <- wrapper.pirat(data = subR25pept, i = length(subR25pept), extension = "base")
#> Error in wrapper.pirat(data = subR25pept, i = length(subR25pept), extension = "base"): 'adjmat' is required.
```

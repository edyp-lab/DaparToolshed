# Normalisation for QFeatures

This method is a wrapper that provides several methods to normalize
quantitative data from objects of class `QFeatures` or
`SummarizedExperiment`.

They are organized in six main families : GlobalQuantileAlignement,
sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn For the first
family, there is no type. For the five other families, two type
categories are available : "Overall" which means that the value for each
protein (ie line in the expression data tab) is computed over all the
samples ; "within conditions" which means that the value for each
protein (ie line in the
[`SummarizedExperiment::assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
data tab) is computed condition by condition. The available methods are
described in
[`normalizeMethods()`](https://edyp-lab.github.io/DaparToolshed/reference/normalization_methods.md).

## Usage

``` r
normalizeFunction(
  obj,
  method,
  conditions = NULL,
  type = "overall",
  subset.norm = NULL,
  quantile = 0.15,
  scaling = FALSE,
  span = 0.7
)
```

## Arguments

- obj:

  An object of class `QFeatures` or `SummarizedExperiment`. If data is
  of class `QFeatures`, the last assay will be normalized.

- method:

  Define the normalization method used :
  ``` "GlobalQuantileAlignment"``,  ```"QuantileCentering"`, `"MeanCentering"`, `"SumByColumns"`, `"LOESS"`or`"vsn"\`.

- conditions:

  A vector of conditions in the dataset. If not provided, the vector
  `"Condition"` from the column metadata will be used.

- type:

  "overall" (shift all the sample distributions at once) or "within
  conditions" (shift the sample distributions within each condition at a
  time).

- subset.norm:

  A vector of index indicating rows to be used for normalization

- quantile:

  A float that corresponds to the quantile used to align the data.

- scaling:

  A boolean that indicates if the variance of the data have to be forced
  to unit (variance reduction) or not.

- span:

  xxx

## Value

`QFeatures` including a new assay with normalized data or
`SummarizedExperiment` with normalized data.

## Author

Manon Gaudin

## Examples

``` r
data(subR25pept)
normalized <- normalizeFunction(subR25pept, method = 'GlobalQuantileAlignment')
```

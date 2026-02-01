# Push p-values based on metacell tags

This function allows to push p-values to 1 based on metacell tags.

## Usage

``` r
pushpvalue(
  obj,
  pvalue,
  scope = "WholeMatrix",
  pattern = "Imputed MEC",
  percent = TRUE,
  threshold = 1,
  conditions = NULL,
  operator = ">=",
  level = NULL,
  value = 1.00000000001
)
```

## Arguments

- obj:

  An object of class `QFeatures` or `SummarizedExperiment`. If data is
  of class `QFeatures`, the last assay will be used.

- pvalue:

  A vector of p-values.

- scope:

  A string for scope to use. Available values are "WholeLine",
  "WholeMatrix", "AllCond" and "AtLeastOneCond".

- pattern:

  A vector of tag to use.

- percent:

  A boolean to indicate whether the threshold represent an absolute
  value (percent = FALSE) or a percentage (percent = TRUE).

- threshold:

  A value that corresponds to the threshold value. Either an integer if
  percent = FALSE, or a float between 0 and 1 of percent = TRUE.

- conditions:

  A vector of conditions in the dataset. If not provided, the vector
  `"Condition"` from the column metadata will be used.

- operator:

  A string for operator to use. Available operators are "\<=", "\<",
  "\>=", "\>", "==" and "!=".

- level:

  A string for dataset type. Either "peptide" or "protein" If not
  provided, the string obtained from `typeDataset(obj)` will be used.

- value:

  A float, value to assign to the pushed p-value. By default, the value
  is set slightly above 1 to be able to differentiate the pushed value.

## Value

A vector with pushed p-values.

## Author

Manon Gaudin

## Examples

``` r
data(subR25prot)
obj <- subR25prot
# Simulate imputation
obj <- NAIsZero(obj, 1)
obj <- NAIsZero(obj, 2)
allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
pushpvalue(obj, allComp$P_Value[, 1], scope = "WholeMatrix", pattern = c("Missing MEC", "Missing POV"), percent = TRUE, threshold = 0.5, operator = ">=",)
#>   [1] 3.943475e-01 9.212880e-01 1.592125e-06 3.994247e-02 8.736050e-01
#>   [6] 6.701855e-02 2.268774e-01 1.418284e-02 9.136290e-01 7.068714e-06
#>  [11] 1.463797e-06 4.786901e-04 2.715090e-05 6.943010e-05 6.776124e-07
#>  [16] 6.266918e-04 2.337633e-04 2.929802e-05 8.377510e-04 2.738104e-06
#>  [21] 6.067336e-05 7.481319e-06 4.204920e-06 2.139632e-04 8.078146e-05
#>  [26] 9.891898e-04 8.078612e-03 6.453347e-06 3.660800e-07 5.923487e-07
#>  [31] 3.710237e-06 4.894113e-04 7.085489e-07 2.167121e-07 1.809192e-05
#>  [36] 1.409477e-04 1.028317e-05 1.372522e-05 6.543750e-07 1.504810e-05
#>  [41] 2.905776e-05 1.059255e-06 3.489724e-06 7.841412e-06 7.564365e-06
#>  [46] 1.730482e-04 1.413031e-05 7.778103e-04 4.037015e-06 6.621161e-05
#>  [51] 4.979948e-05 7.158380e-03 4.779760e-07 7.170332e-06 1.556121e-04
#>  [56] 7.757997e-06 5.127698e-05 5.013309e-01 7.751411e-01 3.339172e-01
#>  [61] 2.225911e-01 9.862644e-01 2.039747e-01 6.941527e-01 4.613080e-01
#>  [66] 8.786207e-02 8.092450e-01 2.602074e-01 9.002982e-01 4.689636e-01
#>  [71] 1.318567e-01 9.607034e-01 7.834765e-01 1.369905e-01 3.188140e-01
#>  [76] 2.714059e-01 8.786232e-02 8.453993e-01 4.930609e-09 9.735355e-01
#>  [81] 4.792999e-01 7.621764e-01 4.784515e-01 3.804386e-01 9.491756e-01
#>  [86] 4.926189e-01 2.673247e-01 2.271292e-01 3.105061e-02 9.803894e-01
#>  [91] 7.677740e-02 4.646400e-01 2.308136e-01 5.871518e-01 3.932452e-01
#>  [96] 2.213727e-01 2.123834e-02 2.639901e-01 7.720213e-01 5.249533e-01
```

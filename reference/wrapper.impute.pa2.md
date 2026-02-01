# Missing values imputation from a `SummarizedExperiment` object

This method is a wrapper to the function
[`impute.pa2()`](https://edyp-lab.github.io/DaparToolshed/reference/impute.pa2.md)
adapted to objects of class `SummarizedExperiment`.

## Usage

``` r
wrapper.impute.pa2(
  obj,
  design,
  q.min = 0,
  q.norm = 3,
  eps = 0,
  distribution = "unif"
)
```

## Arguments

- obj:

  An object of class `SummarizedExperiment`.

- design:

  A data.frame containing the columns "quantCols" corresponding to the
  samples name and "Condition" to the condition of each sample.

- q.min:

  A quantile value of the observed values allowing defining the maximal
  value which can be generated. This maximal value is defined by the
  quantile q.min of the observed values distribution minus eps. Default
  is 0 (the maximal value is the minimum of observed values minus eps).

- q.norm:

  A quantile value of a normal distribution allowing defining the
  minimal value which can be generated. Default is 3 (the minimal value
  is the maximal value minus qn\*median(sd(observed values)) where sd is
  the standard deviation of a row in a condition).

- eps:

  A value allowing defining the maximal value which can be generated.
  This maximal value is defined by the quantile q.min of the observed
  values distribution minus eps. Default is 0.

- distribution:

  The type of distribution used. Values are `unif` (default) or `beta`.

## Value

The object `obj` which has been imputed

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
# \donttest{
utils::data(subR25pept)
design <- design.qf(subR25pept)
subR25pept <- wrapper.impute.pa2(subR25pept[[1]], design)
#> Warning: data length is not a multiple of split variable
#> Error: subscript contains invalid names
# }
```

# Missing values imputation from a `MSnSet` object

This method is a variation to the function `impute.pa()` from the
package `imp4p`.

## Usage

``` r
impute.pa2(
  tab,
  conditions,
  q.min = 0,
  q.norm = 3,
  eps = 0,
  distribution = "unif"
)
```

## Arguments

- tab:

  An object of class `MSnSet`.

- conditions:

  A vector of conditions in the dataset

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

  The type of distribution used. Values are unif or beta.

## Value

The object `obj` which has been imputed

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
library(QFeatures)
utils::data(subR25pept)
qdata <- assay(subR25pept[[1]])
conds <- design.qf(subR25pept)$Condition
obj.imp <- impute.pa2(qdata, conds, distribution = "beta")
```

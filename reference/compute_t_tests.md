# xxxxxx

xxxxxx

## Usage

``` r
compute_t_tests(
  obj,
  i = 1,
  contrast = c("OnevsOne", "OnevsAll"),
  type = c("Student", "Welch")
)
```

## Arguments

- obj:

  A matrix of quantitative data, without any missing values.

- i:

  xxx

- contrast:

  Indicates if the test consists of the comparison of each biological
  condition versus each of the other ones (contrast=1; for example
  H0:"C1=C2" vs H1:"C1!=C2", etc.) or each condition versus all others
  (contrast=2; e.g. H0:"C1=(C2+C3)/2" vs H1:"C1!=(C2+C3)/2", etc. if
  there are three conditions).

- type:

  xxxxx

## Value

A list of two items : logFC and P_Value; both are dataframe. The first
one contains the logFC values of all the comparisons (one column for one
comparison), the second one contains the pvalue of all the comparisons
(one column for one comparison). The names of the columns for those two
dataframes are identical and correspond to the description of the
comparison.

## Author

Florence Combes, Samuel Wieczorek

## Examples

``` r
library(SummarizedExperiment)
data(subR25prot)
obj <- subR25prot
filter <- FunctionFilter('qMetacellOnConditions',
  cmd = 'delete',
  mode = 'AtLeastOneCond',
  pattern = c("Missing POV", "Missing MEC"),
  conds = design.qf(obj)$Condition,
  percent = TRUE,
  th = 0.8,
  operator = '>')
obj <- filterFeaturesOneSE(obj, name = "Filtered", filters = list(filter))
ttest <- compute_t_tests(obj, 3)

ttest <- compute_t_tests(obj, 3)
```

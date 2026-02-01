# Wrapper for One-way Anova statistical test

Wrapper for One-way Anova statistical test

## Usage

``` r
wrapperClassic1wayAnova(obj, i, with_post_hoc = "No", post_hoc_test = "No")
```

## Arguments

- obj:

  An object of class `QFeatures`.

- i:

  xxx

- with_post_hoc:

  a character string with 2 possible values: "Yes" and "No" (default)
  saying if function must perform a Post-Hoc test or not.

- post_hoc_test:

  character string, possible values are "No" (for no test; default
  value) or TukeyHSD" or "Dunnett". See details of
  [`postHocTest()`](https://edyp-lab.github.io/DaparToolshed/reference/postHocTest.md)
  function to choose the appropriate one.

## Value

A list of two dataframes. First one called "logFC" contains all pairwise
comparisons logFC values (one column for one comparison) for each
analysed feature (Except in the case without post-hoc testing, for which
NAs are returned.); The second one named "P_Value" contains the
corresponding p-values.

## Details

This function allows to perform a 1-way Analysis of Variance. Also
computes the post-hoc tests if the `with_post_hoc` parameter is set to
yes. There are two possible post-hoc tests: the Tukey Honest Significant
Differences (specified as "TukeyHSD") and the Dunnett test (specified as
"Dunnett").

## See also

postHocTest()

## Author

Hélène Borges

## Examples

``` r
# \donttest{
data(subR25prot)
obj <- subR25prot
filter <- FunctionFilter('qMetacellOnConditions',
  cmd = 'delete',
  mode = 'AtLeastOneCond',
  pattern = c("Missing POV", "Missing MEC"),
  conds = design.qf(obj)$Condition,
  percent = TRUE,
  th = 0.8,
  operator = '<')

obj <- filterFeaturesOneSE(obj, filter)

anovatest <- wrapperClassic1wayAnova(obj, 2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 't': contrasts can be applied only to factors with 2 or more levels
# }
```

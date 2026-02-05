# Post-hoc tests for classic 1-way ANOVA

This function allows to compute a post-hoc test after a 1-way ANOVA
analysis. It expects as input an object obtained with the function
`classic1wayAnova`. The second parameter allows to choose between 2
different post-hoc tests: the Tukey Honest Significant Differences
(specified as "TukeyHSD") and the Dunnett test (specified as "Dunnett").

## Usage

``` r
postHocTest(aov_fits, post_hoc_test = "TukeyHSD")
```

## Arguments

- aov_fits:

  a list containing aov fitted model objects

- post_hoc_test:

  a character string indicating which post-hoc test to use. Possible
  values are "TukeyHSD" or "Dunnett". See details for what to choose
  according to your experimental design.

## Value

a list of 2 dataframes: first one called "LogFC" contains all pairwise
comparisons logFC values (one column for one comparison) for each
analysed feature; The second one named "P_Value" contains the
corresponding pvalues.

## Details

This is a function allowing to realise post-hoc tests for a set of
proteins/peptides for which a classic 1-way anova has been performed
with the function `classic1wayAnova`. Two types of tests are currently
available: The Tukey HSD's test and the Dunnett's test. Default is
Tukey's test. The Tukey HSD's test compares all possible pairs of means,
and is based on a studentized range distribution. Here is used the
[`TukeyHSD()`](https://rdrr.io/r/stats/TukeyHSD.html) function, which
can be applied to balanced designs (same number of samples in each
group), but also to midly unbalanced designs. The Dunnett's test
compares a single control group to all other groups. Make sure the
factor levels are properly ordered.

## Author

Hélène Borges

## Examples

``` r
# \donttest{
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
qdata <- SummarizedExperiment::assay(obj[[3]])
conds <- design.qf(obj)$Condition
anova_tests <- apply(qdata, 1, classic1wayAnova, conditions = as.factor(conds))
anova_tests <- t(anova_tests)

names(anova_tests) <- rownames(qdata)
pht <- postHocTest(aov_fits = anova_tests)
# }
```

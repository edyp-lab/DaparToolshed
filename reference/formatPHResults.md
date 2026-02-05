# Extract logFC and raw pvalues from multiple post-hoc models summaries

Extract logFC and raw pvalues from multiple post-hoc models summaries

## Usage

``` r
formatPHResults(post_hoc_models_summaries)
```

## Arguments

- post_hoc_models_summaries:

  a list of summaries of post-hoc models.

## Value

a list of 2 dataframes containing the logFC values and pvalues for each
comparison.

## Author

Helene Borges

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
tms <- lapply(
  anova_tests,
  function(x) {
    summary(multcomp::glht(x,
      linfct = multcomp::mcp(conditions = "Tukey")
    ),
      test = multcomp::adjusted("none")
    )
   }
)
res <- formatPHResults(tms)
# }
```

# Function to perform a One-way Anova statistical test on a MsnBase dataset

Function to perform a One-way Anova statistical test on a MsnBase
dataset

## Usage

``` r
classic1wayAnova(current_line, conditions)
```

## Arguments

- current_line:

  The line currently treated from the quantitative data to perform the
  ANOVA

- conditions:

  The conditions represent the different classes of the studied factor

## Value

A named vector containing all the different values of the aov model

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

anova_tests <- t(apply(assay(obj[[2]]), 1, classic1wayAnova,
conditions = as.factor(design.qf(obj)$Condition)
))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': contrasts can be applied only to factors with 2 or more levels
# }
```

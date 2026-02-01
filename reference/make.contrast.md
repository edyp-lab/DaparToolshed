# Builds the contrast matrix

Builds the contrast matrix

## Usage

``` r
make.contrast(design, condition, contrast = 1, design.level = 1)
```

## Arguments

- design:

  The data.frame which correspond to the
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  function of package `MultiAssayExperiment`.

- condition:

  xxxxx

- contrast:

  An integer that Indicates if the test consists of the comparison of
  each biological condition versus each of the other ones (Contrast=1;
  for example H0:"C1=C2" vs H1:"C1!=C2", etc.) or each condition versus
  all others (Contrast=2; e.g. H0:"C1=(C2+C3)/2" vs H1:"C1!=(C2+C3)/2",
  etc. if there are three conditions).

- design.level:

  xxx

## Value

A contrast matrix

## Author

Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
design <- make.design(SummarizedExperiment::colData(subR25pept))
conds <- design.qf(subR25pept)$Condition
make.contrast(design, conds)
#> [1] "( ConditionA )/ 1 -( ConditionB )/ 1"
```

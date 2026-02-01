# Computes a hierarchical differential analysis

Computes a hierarchical differential analysis

## Usage

``` r
limmaCompleteTest(qData, sTab, comp.type = "OnevsOne")
```

## Arguments

- qData:

  A matrix of quantitative data, without any missing values.

- sTab:

  A dataframe of experimental design (design.qf()).

- comp.type:

  A string that corresponds to the type of comparison. Values are:
  'anova1way', 'OnevsOne' and 'OnevsAll'; default is 'OnevsOne'.

## Value

A list of two dataframes : logFC and P_Value. The first one contains the
logFC values of all the comparisons (one column for one comparison), the
second one contains the pvalue of all the comparisons (one column for
one comparison). The names of the columns for those two dataframes are
identical and correspond to the description of the comparison.

## Author

Hélène Borges, Thomas Burger, Quentin Giai-Gianetto, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
qData <- as.matrix(SummarizedExperiment::assay(subR25pept[[2]]))
sTab <- SummarizedExperiment::colData(subR25pept)
limma <- limmaCompleteTest(qData, sTab, comp.type = "anova1way")
#> Warning: Partial NA coefficients for 7 probe(s)
```

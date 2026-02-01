# Identification of differentially abundant peptide/protein

This function allows to identify differentially abundant peptide/protein

## Usage

``` r
is.differential(pvalue, logFC, thpvalue, thlogFC)
```

## Arguments

- pvalue:

  A vector of p-values.

- logFC:

  A vector of logFC.

- thpvalue:

  A float indicating the p-value threshold.

- thlogFC:

  A float indicating the logFC threshold.

## Value

A vector indicating which peptide/protein is differentially abundant (1)
or not (0).

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
is.differential(allComp$P_Value[, 1], allComp$logFC[, 1], 0.05, 0.5)
#> Error in is.differential(allComp$P_Value[, 1], allComp$logFC[, 1], 0.05,     0.5): object 'pval' not found
```

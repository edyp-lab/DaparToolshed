# Imputation of peptides having no values in a biological condition.

This method is a wrapper to the function `impute.mle()` of the package
`imp4p` adapted to an object of class `SummarizedExperiment`. It does
not impute MEC missing values.

## Usage

``` r
wrapper.impute.mle(obj, grp)
```

## Arguments

- obj:

  An object of class `SummarizedExperiment`.

- grp:

  A vector of conditions in the dataset.

## Value

The `SummarizedExperiment::assay(obj)` matrix with imputed values
instead of missing values.

## Author

Samuel Wieczorek

## Examples

``` r
utils::data(subR25pept)
level <- 'peptide'
# Delete whole empty lines
metacell.mask <- DaparToolshed::match.metacell(
qMetacell(subR25pept[[1]]), 
c("Missing POV", "Missing MEC"), level)
indices <- GetIndices_WholeMatrix(metacell.mask, op = ">=", th = 1)
grp <- design.qf(subR25pept)$Condition
subR25pept <- wrapper.impute.mle(subR25pept[[2]], grp)
```

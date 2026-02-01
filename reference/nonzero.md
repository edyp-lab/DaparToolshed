# Retrieve the indices of non-zero elements in sparse matrices

This function retrieves the indices of non-zero elements in sparse
matrices of class dgCMatrix from package Matrix. This function is
largely inspired from the package `RINGO`.

## Usage

``` r
nonzero(x)
```

## Arguments

- x:

  A sparse matrix of class dgCMatrix

## Value

A two-column matrix

## Author

Samuel Wieczorek

## Examples

``` r
library(Matrix)
mat <- Matrix(c(0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1),
    nrow = 5, byrow = TRUE, sparse = TRUE)
res <- nonzero(mat)
```

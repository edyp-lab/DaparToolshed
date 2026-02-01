# splits an adjacency matrix into specific and shared

Method to split an adjacency matrix into specific and shared

## Usage

``` r
splitAdjacencyMat(X)
```

## Arguments

- X:

  An adjacency matrix

## Value

A list of two adjacency matrices

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
ll <- splitAdjacencyMat(X)
```

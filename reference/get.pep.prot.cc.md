# Build the list of connex composant of the adjacency matrix

Build the list of connex composant of the adjacency matrix

## Usage

``` r
get.pep.prot.cc(X)
```

## Arguments

- X:

  An adjacency matrix

## Value

A list of CC

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
ll <- get.pep.prot.cc(X)
```

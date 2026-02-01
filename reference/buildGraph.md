# Display a CC

Display a CC

## Usage

``` r
buildGraph(The.CC, X)
```

## Arguments

- The.CC:

  A cc (a list)

- X:

  xxxxx

## Value

A plot

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
ll <- get.pep.prot.cc(X)
g <- buildGraph(ll[[1]], X)
```

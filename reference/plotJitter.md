# Jitter plot of CC

Jitter plot of CC

## Usage

``` r
plotJitter(list.of.cc = NULL)
```

## Arguments

- list.of.cc:

  List of cc such as returned by the function get.pep.prot.cc

## Value

A plot

## Author

Thomas Burger

## Examples

``` r
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[1]])
ll <- get.pep.prot.cc(X)
plotJitter(ll)

```

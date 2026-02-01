# Applies aov() on a vector of protein abundances using the design derived from the sample names (simple aov wrapper)

Applies aov() on a vector of protein abundances using the design derived
from the sample names (simple aov wrapper)

## Usage

``` r
OWAnova(current_protein, conditions)
```

## Arguments

- current_protein:

  a real vector

- conditions:

  the list of groups the protein belongs to

## Value

See aov()

## Author

Thomas Burger

## Examples

``` r
protein_abundance <- rep(rnorm(3, mean= 18, sd=2), each=3) + rnorm(9)
groups <- c(rep("group1",3),rep("group2",3),rep("group3",3))
OWAnova(protein_abundance,groups)
#> Call:
#>    aov(formula = intensities ~ conditions, data = NULL)
#> 
#> Terms:
#>                 conditions Residuals
#> Sum of Squares    4.045624  1.090825
#> Deg. of Freedom          2         6
#> 
#> Residual standard error: 0.426385
#> Estimated effects may be unbalanced
```

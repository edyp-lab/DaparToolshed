# iteratively applies OWAnova() on the features of an MSnSet object

iteratively applies OWAnova() on the features of an MSnSet object

## Usage

``` r
applyAnovasOnProteins(obj, i)
```

## Arguments

- obj:

  a QFreatures object

- i:

  xxx '

## Value

a list of linear models

## Author

Thomas Burger

## Examples

``` r
data(subR25prot)
applyAnovasOnProteins(subR25prot[1:5,], 1)
#>      CON__A2I7N1;CON__A2I7N0 CON__P00761
#> [1,] aov,13                  aov,13     
#>      P02768upsedyp|ALBU_HUMAN_upsedyp;CON__P02768-1 CON__P04264 CON__P07477
#> [1,] aov,13                                         aov,13      aov,13     
#> attr(,"names")
#> [1] "CON__A2I7N1;CON__A2I7N0"                       
#> [2] "CON__P00761"                                   
#> [3] "P02768upsedyp|ALBU_HUMAN_upsedyp;CON__P02768-1"
#> [4] "CON__P04264"                                   
#> [5] "CON__P07477"                                   
```

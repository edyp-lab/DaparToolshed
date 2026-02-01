# Number of each metacell tags

Number of each metacell tags

## Usage

``` r
GetNbTags(obj)
```

## Arguments

- obj:

  A instance of the class `SummarizedExperiment`

## Value

An integer

## Examples

``` r
data(subR25prot)
GetNbTags(subR25prot[[1]])
#>                 Any          Quantified Quant. by direct id  Quant. by recovery 
#>                   0                   0                 550                  14 
#>             Missing         Missing POV         Missing MEC             Imputed 
#>                   0                  24                  12                   0 
#>         Imputed POV         Imputed MEC       Combined tags 
#>                   0                   0                   0 
```

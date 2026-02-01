# Returns the possible number of values in lines in the data

Returns the possible number of values in lines in the data

## Usage

``` r
getListNbValuesInLines(object, conds, type = "WholeMatrix")
```

## Arguments

- object:

  An object of class `QFeatures`

- conds:

  xxxx

- type:

  WholeMatrix, AllCond or AtLeastOneCond

## Value

An integer

## Author

Samuel Wieczorek, Enora Fremy

## Examples

``` r
data(subR25prot)
res <- getListNbValuesInLines(subR25prot[[1]])
```

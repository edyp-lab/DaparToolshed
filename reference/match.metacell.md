# Similar to the function `is.na()` but focused on the equality with the paramter 'type'.

Similar to the function [`is.na()`](https://rdrr.io/r/base/NA.html) but
focused on the equality with the paramter 'type'.

## Usage

``` r
match.metacell(metadata, pattern = NULL, level)
```

## Arguments

- metadata:

  A data.frame

- pattern:

  The value to search in the dataframe

- level:

  A string designing the type of entity/pipeline. Available values are:
  `peptide`, `protein`

## Value

A boolean dataframe

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
metadata <- qMetacell(subR25pept[[1]])
m <- match.metacell(metadata, pattern = "Missing", level = "peptide")
m <- match.metacell(metadata, pattern = 'Missing POV', level = "peptide")
m <- match.metacell(metadata, pattern = c('Missing', 'Missing POV'), level = "peptide")
```

# Version number of Prostar suite

This function gives the version number of the packages of the Prostar
suite and which propose data processing. This information can be useful
if the user wants to publish its works or to rerun a data processing
pipeline tin a given set of conditions. The packages which are concerned
are `Prostar`, `DaparToolshed` and `DaparToolshedData`

## Usage

``` r
ProstarVersions()
```

## Value

A `list(3)`

## Author

Samuel Wieczorek

## Examples

``` r
ProstarVersions()
#> $Prostar
#> [1] NA
#> 
#> $DaparToolshed
#> [1] "0.99.24"
#> 
#> $DaparToolshedData
#> [1] NA
#> 
```

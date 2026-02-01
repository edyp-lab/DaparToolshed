# Standardize names

Replace ".", ' ', '-' in
[`character()`](https://rdrr.io/r/base/character.html) by '\_' to be
compliant with functions of `Shinyjs`, `Shiny`

## Usage

``` r
ReplaceSpecialChars(x)
```

## Arguments

- x:

  A [`character()`](https://rdrr.io/r/base/character.html) to be
  processed

## Value

A [`character()`](https://rdrr.io/r/base/character.html) of the same
length as 'x' with modified names.

## Author

Samuel Wieczorek

## Examples

``` r
ReplaceSpecialChars(c("foo.1", "foo-2", "foo 3"))
#> [1] "foo_1" "foo-2" "foo_3"
```

# Convert a list to unnumbered HTML list

xxx

## Usage

``` r
ConvertListToHtml(ll)
```

## Arguments

- ll:

  A [`list()`](https://rdrr.io/r/base/list.html) of
  [`character()`](https://rdrr.io/r/base/character.html)

## Value

HTML

## Examples

``` r
ConvertListToHtml(list('foo1', 'foo2', 'foo3'))
#> [1] "<ul><li>foo1</li> <li>foo2</li> <li>foo3</li></ul>"
```

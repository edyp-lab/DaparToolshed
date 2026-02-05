# Customised resetZoom Button of highcharts plots

Customised resetZoom Button of highcharts plots

## Usage

``` r
customChart(
  hc,
  chartType = "scatter",
  zoomType = "None",
  width = 0,
  height = 0
)
```

## Arguments

- hc:

  A highcharter object

- chartType:

  The type of the plot

- zoomType:

  The type of the zoom (one of "x", "y", "xy", "None")

- width:

  The width of the plot

- height:

  The height of the plot

## Value

A highchart plot

## Author

Samuel Wieczorek

## Examples

``` r
if (interactive()) {
    library(highcharter)
    hc <- highchart()
    hc_chart(hc, type = "line")
    hc_add_series(hc, data = c(29, 71, 40))
    customChart(hc)
}

NULL
#> NULL
```

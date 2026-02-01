# Exports a `QFeatures` object to a Excel file.

This function exports an instance of the class `QFeatures` to a Excel
file. The resulting file is composed of four sheets:

- `quantitative data` which contains the content of
  [`assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object with a color code for each cell w.r.t. to cell quantitative
  metadata.

- `metadata` which is the content of
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with only one-dimensionnal data (i.e. the adjacencyMatrix and the
  qMetacell slots are not part of the sheet),

- `exp. design` which is the content of
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
  Each condition in the table is colored with a different color,

- `quantitative metadata` which is the content of
  [`qMetacell()`](https://edyp-lab.github.io/DaparToolshed/reference/QFeatures-accessors.md).
  There is a color code for the different tags.

xxx

## Usage

``` r
write2excel(object, ...)

# S4 method for class 'QFeatures'
write2excel(object, i = NULL, filename = "newFile", writeColdData = TRUE, ...)

# S4 method for class 'SummarizedExperiment'
write2excel(object, filename, exp.design, writeColData = TRUE, ...)

.write2excel(object, filename, exp.design, writeColData = TRUE)

addColors(wb, n, tags, colors)
```

## Arguments

- object:

  xxx

- ...:

  xxx

- i:

  xxx

- filename:

  xxx

- writeColdData:

  xxx

- exp.design:

  xxx

- writeColData:

  A boolean

- wb:

  A workbook

- n:

  A `integer(1)` which is the number of sheet in the workbook.

- tags:

  xxx

- colors:

  A [`character()`](https://rdrr.io/r/base/character.html) which
  contains the HEX code for colors. The size of this vector must be the
  same as the number of tags.

## Value

A Excel file.

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25prot)

#---------------------------------------
# Export the whole dataset
#---------------------------------------

write2excel(subR25prot, filename = "foo")
#> NULL
unlink('foo.xls')
write2excel(subR25prot, 1, "foo")
#> [1] "foo.xlsx"
unlink('foo.xls')
```

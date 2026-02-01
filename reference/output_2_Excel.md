# This function exports a data.frame to a Excel file.

This function exports a `MSnSet` data object to a Excel file. Each of
the three data.frames in the `MSnSet` object (ie experimental data,
phenoData and metaData are respectively integrated into separate sheets
in the Excel file).

The colored cells in the experimental data correspond to the original
missing values which have been imputed.

## Usage

``` r
readExcel(file, sheet = NULL)

listSheets(file)

write_Assay_To_Excel(wb, obj, i, n)

WriteHistory(wb, obj, n)

Write_SamplesData_to_Excel(wb, obj, n)

Write_RowData(wb, obj, i, n)

write.excel(obj, filename)
```

## Arguments

- file:

  The name of the Excel file.

- sheet:

  The name of the sheet

- wb:

  xxxx

- obj:

  xxx

- i:

  xxx

- n:

  xxx

- filename:

  A character string for the name of the Excel file.

## Value

A Excel file (.xlsx)

## Author

Samuel Wieczorek

## Examples

``` r
# \donttest{
library(QFeatures)
data(subR25prot)
df <- assay(subR25prot[[1]])
tags <- qMetacell(subR25prot[[1]])
colors <- list(
    "Missing POV" = "lightblue",
    "Missing MEC" = "orange",
    "Quant. by recovery" = "lightgrey",
    "Quant. by direct id" = "white",
    "Combined tags" = "red"
)
write.excel(subR25prot, filename = "toto.xlsx")
#> [1] "toto.xlsx"


data(subR25pept)
write.excel(subR25pept, "foo.xlsx")
#> [1] "foo.xlsx"
# }

```

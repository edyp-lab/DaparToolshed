# Delete the lines in the matrix of intensities and the metadata table given their indices.

Delete the lines in the matrix of intensities and the metadata table
given their indices.

## Usage

``` r
GetIndices_FunFiltering(
  obj,
  conds,
  level,
  pattern = NULL,
  type = NULL,
  percent,
  op,
  th
)
```

## Arguments

- obj:

  An object of class `SummarizedExperiment` containing quantitative
  data.

- conds:

  xxx

- level:

  A vector of integers which are the indices of lines to delete.

- pattern:

  A string to be included in the `SummarizedExperiment` object for log.

- type:

  xxx

- percent:

  xxx

- op:

  xxx

- th:

  xxx

## Value

An instance of class `SummarizedExperiment` that have been filtered.

## Author

Samuel Wieczorek

## Examples

``` r
if (FALSE) { # interactive()
NA
}
```

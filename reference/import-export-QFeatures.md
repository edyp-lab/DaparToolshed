# Creates an object of class `QFeatures` from text file.

Creates an object of class `QFeatures` from a single tabulated-like file
for quantitative and meta-data and a dataframe for the samples
description.

## Usage

``` r
createQFeatures(
  data = NULL,
  file = "myDataset",
  sample,
  indQData,
  keyId = "AutoID",
  indexForMetacell = NULL,
  logData = FALSE,
  force.na = TRUE,
  typeDataset,
  parentProtId = NULL,
  analysis = "foo",
  description = NULL,
  processes = NULL,
  typePipeline = NULL,
  name.pipeline = NULL,
  software = NULL,
  name = "original"
)
```

## Arguments

- data:

  The name of a tab-separated file that contains the data.

- file:

  A `character(1)`. The name of a file xxx

- sample:

  A dataframe describing the samples (in lines).

- indQData:

  A vector of string where each element is the name of a column in
  designTable that have to be integrated in the
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  table of the `QFeatures` object.

- keyId:

  A `character(1)` or `numeric(1)` which is the indice of the column
  containing the ID of entities (peptides or proteins)

- indexForMetacell:

  They must be in the same order as the samples in the experimental
  design

- logData:

  xxx

- force.na:

  A `boolean` that indicates if the '0' and 'NaN' values of quantitative
  values must be replaced by 'NA' (Default is FALSE)

- typeDataset:

  A string that indicates whether the dataset is about

- parentProtId:

  A `character(1)` For peptide entities, a string which is the name of a
  column in rowData. It contains the id of parent proteins and is used
  to generate adjacency matrix and process to aggregation.

- analysis:

  A `character(1)` which is the name of the MS study.

- description:

  A text which describes the study.

- processes:

  A vector of A [`character()`](https://rdrr.io/r/base/character.html)
  which contains the name of processes which has already been run on the
  data. Default is 'original'.

- typePipeline:

  A `character(1)` The type of pipeline used with this dataset. The list
  of predefined pipelines in DaparToolshed can be obtained with the
  function `pipelines()`. Default value is NULL

- name.pipeline:

  A string

- software:

  A `character(1)`

- name:

  A `character(1)` which is the name of the assay in the QFeatures
  object. Default is 'original'

## Value

An instance of class `QFeatures`.

## Author

Samuel Wieczorek, Manon Gaudin

## Examples

``` r
NULL
#> NULL
```

# List of metacell tags

This function gives the list of metacell tags available.

- onlyPresent: In this case, the function gives the tags found in a
  dataset. In addition, and w.r.t to the hierarchy of tags, if all
  leaves of a node are present, then the tag corresponding to this node
  is added.

These names are common to all assays contained in the object. This is
why they are stored in the global metadata. This function is used
whenever it i s necessary to (re)detect MEC and POV (new dataset or when
post processing protein qMetacell after aggregation)

## Usage

``` r
paramshistory(object, ...)

# S4 method for class 'QFeatures'
paramshistory(object, i, slotName = "paramshistory")

# S4 method for class 'SummarizedExperiment'
paramshistory(object, slotName = "paramshistory")

paramshistory(object, i, slotName = "paramshistory") <- value

GetMetacellTags(object, ...)

# S4 method for class 'QFeatures'
GetMetacellTags(object, i, ...)

# S4 method for class 'SummarizedExperiment'
GetMetacellTags(object, ...)

# S4 method for class 'data.frame'
GetMetacellTags(object, ...)

qMetacell(object, ...)

# S4 method for class 'QFeatures'
qMetacell(object, i)

# S4 method for class 'SummarizedExperiment'
qMetacell(object)

qMetacell(object, i, slotName = "qMetacell") <- value

GetUniqueTags(object, ...)

# S4 method for class 'QFeatures'
GetUniqueTags(object, i)

# S4 method for class 'SummarizedExperiment'
GetUniqueTags(object)

.GetMetadataSlot(object, slotName = NULL)

.GetRowdataSlot(object, slotName = NULL)

ConnectedComp(object, ...)

# S4 method for class 'QFeatures'
ConnectedComp(object, i, slotName = "ConnectedComp")

# S4 method for class 'SummarizedExperiment'
ConnectedComp(object, slotName = "ConnectedComp")

ConnectedComp(object, i, slotName = "ConnectedComp") <- value

typeDataset(object, ...)

# S4 method for class 'QFeatures'
typeDataset(object, i, slotName = "typeDataset")

# S4 method for class 'SummarizedExperiment'
typeDataset(object, slotName = "typeDataset")

typeDataset(object, i, slotName = "typeDataset") <- value

idcol(object, ...)

# S4 method for class 'QFeatures'
idcol(object, i, slotName = "idcol")

# S4 method for class 'SummarizedExperiment'
idcol(object, slotName = "idcol")

idcol(object, i, slotName = "idcol") <- value

parentProtId(object, ...)

# S4 method for class 'QFeatures'
parentProtId(object, i, slotName = "parentProtId")

# S4 method for class 'SummarizedExperiment'
parentProtId(object, slotName = "parentProtId")

parentProtId(object, i, slotName = "parentProtId") <- value

filename(object, ...)

# S4 method for class 'QFeatures'
filename(object, slotName = "filename")

filename(object, slotName = "filename") <- value

analysis(object, ...)

# S4 method for class 'QFeatures'
analysis(object, i, slotName = "analysis")

# S4 method for class 'SummarizedExperiment'
analysis(object, slotName = "analysis")

analysis(object, i, slotName = "analysis") <- value

version(object, ...)

# S4 method for class 'QFeatures'
version(object, slotName = "version")

version(object, slotName = "version") <- value

design.qf(object, ...)

# S4 method for class 'QFeatures'
design.qf(object, slotName = "design")

design.qf(object, slotName = "design") <- value

mainAssay(object)

HypothesisTest(object, ...)

# S4 method for class 'QFeatures'
HypothesisTest(object, i, slotName = "HypothesisTest")

# S4 method for class 'SummarizedExperiment'
HypothesisTest(object, slotName = "HypothesisTest")

HypothesisTest(object, i) <- value

DifferentialAnalysis(object, ...)

# S4 method for class 'QFeatures'
DifferentialAnalysis(object, i, slotName = "DifferentialAnalysis")

# S4 method for class 'SummarizedExperiment'
DifferentialAnalysis(object, slotName = "DifferentialAnalysis")

DifferentialAnalysis(object, i) <- value

names_metacell(object, ...)

# S4 method for class 'QFeatures'
names_metacell(object, i, slotName = "names_metacell")

# S4 method for class 'SummarizedExperiment'
names_metacell(object, slotName = "names_metacell")

names_metacell(object, i, slotName = "names_metacell") <- value
```

## Arguments

- object:

  xxx

- ...:

  Additional parameters

- i:

  The index or name of the assays to extract the quantitative metadata
  from. All must have a rowdata variable named as `slotName`

- slotName:

  xxx

- value:

  xxx

## Value

A vector of tags.

NA

NA

NA

NA

NA

NA

NA

NA

## Details

Additional slots for rowdata of a `SummarizedExperiment` object:

- qMetacell: xxx

Additional slots for Metadata for a `QFeatures` object:

- xxx: xxxx

Additional slots for Metadata for a `SummarizedExperiment` object:

- qMetacell: xxxx

- parentProtId: xxx

- idcol: xxxx

- typeDataset: xxx

## Quantitative metadata

Default slotName is "qMetacell". The value is an adjacency matrix with
row and column names. The matrix will be coerced to compressed,
column-oriented sparse matrix (class `dgCMatrix`) as defined in the
`Matrix` package, as generaled by the
[`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
constructor.

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25pept)
GetMetacellTags(subR25pept, 1, level="peptide")
#> [1] "Quantified"          "Quant. by direct id" "Quant. by recovery" 
#> [4] "Missing"             "Missing POV"         "Missing MEC"        
#> [7] "Imputed"             "Imputed POV"         "Imputed MEC"        
GetMetacellTags(subR25pept, 1, level="peptide", onlyPresent=TRUE)
#> [1] "Quant. by direct id" "Quant. by recovery"  "Missing POV"        
#> [4] "Missing MEC"         "Quantified"          "Missing"            

data(subR25pept)
design.qf(subR25pept)
```

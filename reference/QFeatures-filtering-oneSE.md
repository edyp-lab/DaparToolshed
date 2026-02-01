# Filter features of one SE based on their rowData

The `filterFeaturesOneSE` methods enables users to filter features based
on a variable in their `rowData`. It is directly inspired of the
function `filterFeature` of the package `QFeatures`. The first
difference is that the filter only applies to one `SummarizedExperiment`
contained in the object rather than applying on all the SE. This method
generates a new `SummarizedExperiment` object which is added to the
`QFeatures` object. If the SE on which the filter applies is the last
one of the object, then a new xxxx. If it is not the last one, the new
SE is added and all the further SE are deleted. The features matching
the. The filters can be provided as instances of class
`AnnotationFilter` (see the package `QFeatures`) or of class
`FunctionFilter` (see below).

## Usage

``` r
FunctionFilter(name, ...)

filterFeaturesOneSE(object, ...)

# S4 method for class 'QFeatures'
filterFeaturesOneSE(object, i, name = "newAssay", filters)
```

## Arguments

- name:

  A `character(1)` naming the new assay. Default is `newAssay`. Note
  that the function will fail if there's already an assay with `name`.

- ...:

  Additional arguments

- object:

  An instance of class `QFeatures` or `SummarizedExperiment`.

- i:

  The index or name of the assay which features will be filtered the
  create the new assay.

- filters:

  A [`list()`](https://rdrr.io/r/base/list.html) containing instances of
  class `AnnotationFilter` or `FunctionFilter`

## Value

A filtered `QFeature` object

## Function filters

The function filters are filters as defined in the `DaparToolshed`
package. Each filter is defined by a name (which is the name of a
function) and a list which contains the parameters passed to the
function. Those filters can be created with the `FunctionFilter`
constructor.

Those functions are divided into two main categories:

- the one that filter on one rowData feature,

- the one based on a two-dimensional information such as the adjacency
  matrix

for the first category, all filters of class
[AnnotationFilter::AnnotationFilter](https://rdrr.io/pkg/AnnotationFilter/man/AnnotationFilter.html)
can be used as they are used in `QFeatures`

For the second category, the package `DaparToolshed` provides filter
functions based either on the adjacency matrix:

- 
- 
- 

Or based on the quantitative metadata (identification):

- 
- 
- 

## Author

Samuel Wieczorek

## Examples

``` r
data(feat1, package = 'QFeatures')
## ----------------------------------------
## Creating function filters
## ----------------------------------------

#FunctionFilter('FUN',
#               param1 = 'value_of_param1',
#               param2 = 'value_of_param2')

FunctionFilter('qMetacellWholeLine',
               cmd = 'delete',
               pattern = 'imputed POV')
#> An object of class "FunctionFilter"
#> Slot "name":
#> [1] "qMetacellWholeLine"
#> 
#> Slot "params":
#> $cmd
#> [1] "delete"
#> 
#> $pattern
#> [1] "imputed POV"
#> 
#> 

## ----------------------------------------------------------------
## Filter the last assay to keep only specific peptides. This filter
## only applies on peptide dataset.
## ----------------------------------------------------------------

spec.filter <- FunctionFilter('specPeptides', list())
## using a user-defined character filter
filterFeaturesOneSE(feat1, list(FunctionFilter('specPeptides', list())))
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 


## ----------------------------------------------------------------
## Filter the last assay to keep only specific peptides and topn 
## peptides. The two filters are run sequentially.
## ----------------------------------------------------------------

lst.filters <- list(FunctionFilter('specPeptides', list()))
lst.filters <- append(lst.filters,
FunctionFilter('topnPeptides',
fun = 'rowSums',
top = 2))
filterFeaturesOneSE(feat1, lst.filters)
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 

## ----------------------------------------------------------------
## Filter the last assay to delete peptides where, in at least one 
## condition, there is less than 80% of samples marked as 'imputed POV'
## ----------------------------------------------------------------

filter <- FunctionFilter('qMetacellOnConditions',
cmd = 'delete',
mode = 'AtLeastOneCond',
pattern = 'imputed POV',
conds = SummarizedExperiment::colData(feat1)$Condition,
percent = TRUE,
th = 0.8,
operator = '<')

 filterFeaturesOneSE(feat1, filter)
#> An instance of class QFeatures (type: bulk) with 1 set:
#> 
#>  [1] psms: SummarizedExperiment with 10 rows and 2 columns 

```

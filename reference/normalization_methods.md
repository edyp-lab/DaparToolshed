# Normalisation

Provides several methods to normalize quantitative data from a
`SummarizedExperiment` object. They are organized in six main families :
GlobalQuantileAlignement, sumByColumns, QuantileCentering,
MeanCentering, LOESS, vsn For the first family, there is no type. For
the five other families, two type categories are available : "Overall"
which means that the value for each protein (ie line in the expression
data tab) is computed over all the samples ; "within conditions" which
means that the value for each protein (ie line in the
[`SummarizedExperiment::assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
data tab) is computed condition by condition.

## Usage

``` r
normalizeMethods(target = "all")

GlobalQuantileAlignment(qData)

SumByColumns(qData, conds = NULL, type = NULL, subset.norm = NULL)

QuantileCentering(
  qData,
  conds = NULL,
  type = "overall",
  subset.norm = NULL,
  quantile = 0.15
)

MeanCentering(
  qData,
  conds,
  type = "overall",
  subset.norm = NULL,
  scaling = FALSE
)

vsn(qData, conds, type = NULL)

LOESS(qData, conds, type = "overall", span = 0.7)
```

## Arguments

- target:

  Category of normalization method to show. Either "all", "withTracking"
  or "withoutTracking".

- qData:

  A data.frame with quantitative data to normalize.

- conds:

  A [`character()`](https://rdrr.io/r/base/character.html) vector which
  is the names of conditions for each sample in the dataset.

- type:

  "overall" (shift all the sample distributions at once) or "within
  conditions" (shift the sample distributions within each condition at a
  time).

- subset.norm:

  A vector of index indicating rows to be used for normalization.

- quantile:

  A float that corresponds to the quantile used to align the data.

- scaling:

  A boolean that indicates if the variance of the data have to be forced
  to unit (variance reduction) or not.

- span:

  A float between 0 and 1 indicating the span of loess smoothing window.

## Value

xxx

A normalized numeric matrix

A normalized numeric matrix

A normalized numeric matrix

A normalized numeric matrix

A normalized numeric matrix

A normalized numeric matrix

## Author

Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier, Enora
Fremy

## Examples

``` r
## Get the list of methods
normalizeMethods()
#> [1] "GlobalQuantileAlignment" "SumByColumns"           
#> [3] "QuantileCentering"       "MeanCentering"          
#> [5] "LOESS"                   "vsn"                    


data(subR25pept)
qData <- SummarizedExperiment::assay(subR25pept[[1]])
conds <- design.qf(subR25pept)$Condition



#normalized <- GlobalQuantileAlignment(qData)

normalized <- SumByColumns(qData, conds,
    type = "within conditions",
    subset.norm = seq_len(10)
)

normalized <- QuantileCentering(
SummarizedExperiment::assay(subR25pept), conds,
type = "within conditions", subset.norm = seq_len(10)
)

normalized <- MeanCentering(qData, conds, type = "overall")

# normalized <- vsn(qData, conds, type = "overall")

normalized <- LOESS(qData, conds, type = "overall")
```

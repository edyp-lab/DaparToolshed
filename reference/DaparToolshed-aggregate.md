# Aggregate an assay's quantitative features

This function aggregates the quantitative features of an assay, applying
a summarization function (`fun`) to sets of features. The `fcol`
variable name points to a rowData column that defines how to group the
features during aggregate. This variable has to be an adjacency matrix.
This function uses
[`QFeatures::aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
to aggregate quantitative data.

The list of agregation methods can be obtained with the function
`aggregateMethods()`. This function compiles both methods from the
packages `DaparToolshed` and `QFeatures`.

Aggregate the quantitative metadata tag.

This function aggregate both quantitative and rowdata from the last
assay contained in a `QFeatures`. Note that the function assumes that
the intensities in the QFeatures are already log-transformed.

This function creates a column for the protein dataset after aggregation
by using the previous peptide dataset.

Aggregation of rowData of a `QFeatures` assay.

Aggregate the metadata

xxx

This function computes the number of proteins that are only defined by
specific peptides, shared peptides or a mixture of two.

This function computes the number of peptides used to aggregate
proteins.

Method to compute the number of quantified peptides used for aggregating
each protein

Method to compute the detailed number of quantified peptides used for
aggregating each protein

Method to compute the detailed number of quantified peptides for each
protein

Method to create a plot with proteins and peptides on a MSnSet object
(peptides)

This function aggregate quantitative data using a method of
redistribution of shared peptides. Intensity of shared peptides are
redistributed proportionally to each protein. Note that the function
assumes that the intensities are not log-transformed.

Aggregation using sum method.

Aggregation using mean method.

Aggregation using median method.

Aggregation using medianPolish method. Note that this method is
parallelized to be more efficient.

Aggregation using robustSummary method.

## Usage

``` r
aggregateFeatures4Prostar(object, ...)

# S4 method for class 'QFeatures'
aggregateFeatures4Prostar(
  object,
  i,
  fcol,
  name = "newAssay",
  fun = MsCoreUtils::robustSummary,
  shared = TRUE,
  n = NULL,
  ...
)

# S4 method for class 'SummarizedExperiment'
aggregateFeatures4Prostar(
  object,
  fcol,
  fun = MsCoreUtils::robustSummary,
  conds,
  shared = TRUE,
  n = NULL,
  ...
)

aggQmetacell(qMeta, X, level, conds)

aggregateMethods()

RunAggregation(
  qf,
  includeSharedPeptides = "Yes_As_Specific",
  operator = "Mean",
  considerPeptides = "allPeptides",
  adjMatrix = "adjacencyMatrix",
  ponderation = "Global",
  n = NULL,
  aggregated_col = NULL,
  max_iter = 500
)

BuildColumnToProteinDataset(peptideData, matAdj, columnName, proteinNames)

Add_Aggregated_rowData(obj, col, i.agg)

metacell_agg(aggregatedSE, originalSE, adj_mat, conds, protname_order)

select_topn(pepData, X, n = 10, funpept = "Mean")

getProteinsStats(X)

CountPep(X)

GetNbPeptidesUsed(pepData, X)

GetDetailedNbPeptidesUsed(pepData, X)

GetDetailedNbPeptides(X)

GraphPepProt(mat)

ExtractUniquePeptides(X)

inner.aggregate.iter(
  pepData,
  X,
  init.method = "Mean",
  method = "Mean",
  n = NULL,
  uniqueiter = FALSE,
  topn_fun = "Mean",
  max_iter = 500
)

inner.sum(pepData, X)

inner.mean(pepData, X)

inner.median(pepData, X)

inner.medianpolish(pepData, X)

inner.robustsummary(pepData, X)
```

## Arguments

- object:

  An instance of class `QFeatures` or `SummarizedExperiment`

- ...:

  Additional parameters passed the `fun`.

- i:

  The index or name of the assay which features will be aggregated the
  create the new assay.

- fcol:

  A `character(1)` naming a rowdata variable (of assay `i` in case of a
  `QFeatures`) defining how to aggregate the features of the assay. This
  variable is a (possibly sparse) matrix. See below for details.

- name:

  A `character(1)` naming the new assay. Default is `newAssay`. Note
  that the function will fail if there's already an assay with `name`.

- fun:

  A function used for quantitative feature aggregation. See details for
  examples.

- shared:

  A `boolean` indication if shared peptides should be considered. If
  `TRUE`, shared peptides

- n:

  A `numeric(1)` specifying the number of peptides to use for each
  protein. If `NULL`, all peptides are considered.

- conds:

  A [`character()`](https://rdrr.io/r/base/character.html) vector which
  is the names of conditions for each sample in the dataset.

- qMeta:

  A `matrix` with quantitative metadata tag.

- X:

  A `matrix` acting as an adjacency matrix.

- level:

  A `character(1)` which is the type of dataset

- qf:

  An instance of class
  [QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html).
  The last assay contained in `qf` will be aggregated. Intensities are
  assumed to already be log-transformed.

- includeSharedPeptides:

  How shared peptides are handled. Either `Yes_As_Specific` (default),
  `Yes_Iterative_Redistribution`, `Yes_Simple_Redistribution` or `No`.
  See below for details.

- operator:

  A function used for quantitative feature aggregation. Available
  functions are `Sum`, `Mean`, `Median`, `medianPolish` or
  `robustSummary`. See below for details.

- considerPeptides:

  A `character(1)` defining what peptide to consider. Available values
  are `allPeptides` (default) and `topN`.

- adjMatrix:

  A `character(1)` naming a rowdata variable from the last assay of `qf`
  containing an adjacency matrix.

- ponderation:

  A `character(1)` defining what to consider to create the coefficient
  for redistribution of shared peptides. Available values are `Global`
  (default), `Condition` or `Sample`.

- aggregated_col:

  A [`character()`](https://rdrr.io/r/base/character.html) of column
  names from rowdata to be aggregated.

- max_iter:

  A `numeric(1)` setting the maximum number of iteration.

- peptideData:

  A data.frame of meta data of peptides. It is the rowData of the
  SummarizedExperiment object.

- matAdj:

  The adjacency matrix used to agregate the peptides data.

- columnName:

  The name(s) of the column in Biobase::rowData(peptides_MSnset) that
  the user wants to keep in the new protein data.frame.

- proteinNames:

  The names of the protein in the new dataset (i.e. rownames)

- obj:

  An instance of class
  [QFeatures::QFeatures](https://rformassspectrometry.github.io/QFeatures/reference/QFeatures-class.html).

- col:

  A [`character()`](https://rdrr.io/r/base/character.html) of column
  names from rowdata to be aggregated.

- i.agg:

  A `numeric(1)` indicating the index of the assay to which add the
  aggregated rowData, using the previous assay's rowData.

- aggregatedSE:

  An instance of class
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  containing the aggregated data.

- originalSE:

  An instance of class
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  containing the non-aggregated data.

- adj_mat:

  An adjacency matrix.

- protname_order:

  A [`character()`](https://rdrr.io/r/base/character.html) vector with
  the protein name in order.

- pepData:

  A `matrix` containing the peptide intensities. Note that the function
  assume that data is already log-transformed.

- funpept:

  A function used for determining a peptide's value. Available functions
  are `Sum`, `Mean` or `Median`.

- mat:

  An adjacency matrix.

- init.method:

  A function used for initializing the aggregation. Available functions
  are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. See
  below for details.

- method:

  A function used for the aggregation. Available functions are `Sum`,
  `Mean`, `Median`, `medianPolish` or `robustSummary`. See below for
  details.

- uniqueiter:

  A bole

- topn_fun:

  A function used to determine how to choose the top n peptides.
  Available functions are `Sum`, `Mean` or `Median`. See below for
  details.

## Value

A `QFeatures` object with an additional assay or a
`SummarizedExperiment` object (or subclass thereof).

NA

A QFeatures with an aggregated assay added.

A vector

An instance of `QFeatures` class with aggregated rowData in specified
assay.

A `SummarizedExperiment` containing the aggregated data.

An adjacency matrix with only the top n peptides selected.

A list

A vector of boolean which is the adjacency matrix but with NA values if
they exist in the intensity matrix.

A data.frame

A list of two items

A data.frame

A histogram

A `matrix` containing the aggregated values.

A `matrix` containing the aggregated values.

A `matrix` containing the aggregated values.

A `matrix` containing the aggregated values.

A `matrix` containing the aggregated values.

A `matrix` containing the aggregated values.

## Details

This function uses
[`QFeatures::aggregateFeatures()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
to aggregate quantitative data.

Aggregation of quantitative data is performed using aggregateFeatures,
or inner. aggregate.iter if `Yes_Iterative_Redistribution` or
`Yes_Simple_Redistribution` is selected.

The handling of shared peptide is as follow :

- `Yes_As_Specific` : Shared peptides are used multiple times. Each
  peptide is duplicated as many times as the number of proteins in which
  they are present, and thus are considered as if they are specific to
  each protein.

- `Yes_Simple_Redistribution` : Intensity of shared peptides are
  redistributed proportionally to each protein. See
  `inner.aggregate.iter` for more information.

- `Yes_Iterative_Redistribution` : Intensity of shared peptides are
  redistributed proportionally to each protein. See
  `inner.aggregate.iter` for more information.

- `No` : No shared peptides are used. If a peptide contained only shared
  peptides, its intensity is set as 0 for every sample.

Available functions are :

- `Sum` : base::colSums()\] or base::rowSums() if
  `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.

- `Mean` : base::colMeans()\] or base::rowMeans() if
  `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.

- `Median` : matrixStats::mcolMedians()\] or matrixStats::rowMedians()
  if `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.

- `medianPolish` : MsCoreUtils::medianPolish().

- `robustSummary` : MsCoreUtils::robustSummary().

Available functions are :

- `Sum` : base::rowSums()

- `Mean` : base::rowMeans()

- `Median` : matrixStats::rowMedians()

- `medianPolish` : MsCoreUtils::medianPolish(), not available for
  `topn_fun`. Note that this method takes significantly more time than
  the others, and is parallelized to be more efficient.

- `robustSummary` : MsCoreUtils::robustSummary(), not available for
  `topn_fun`. Note that this method takes significantly more time than
  the others, and is parallelized to be more efficient.

## Iterative aggregation function

xxxxxx xxxxx

## Quantitative metadata aggregation

The function to aggregate the quantitative metadata is `aggQmetadat()`.

## See also

The *QFeatures* vignette provides an extended example and the
*Aggregation* vignette, for a complete quantitative proteomics data
processing pipeline.

## Author

Samuel Wieczorek, Manon Gaudin

Samuel Wieczorek

Manon Gaudin

Alexia Dorffer

Alexia Dorffer, Samuel Wieczorek

## Examples

``` r
NULL
#> NULL
data(subR25pept)
qMeta <- qMetacell(subR25pept, 1)
X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
level <- typeDataset(subR25pept[[1]])
conds <- SummarizedExperiment::colData(subR25pept)$Condition
aggQmeta <- aggQmetacell(qMeta, X, level, conds)

data(subR25pept)

# Remove empty lines
filter_emptyline <- FunctionFilter("qMetacellWholeLine", cmd = 'delete', pattern = 'Missing MEC')
subR25pept <- filterFeaturesOneSE(object = subR25pept, i = length(subR25pept), name = "Filtered",
              filters = list(filter_emptyline))
# Remove proteins with no peptide associated in adjacency matrix
indx <- which(Matrix::colSums(SummarizedExperiment::rowData(subR25pept[[length(subR25pept)]])$adjacencyMatrix) != 0)
SummarizedExperiment::rowData(subR25pept[[length(subR25pept)]])$adjacencyMatrix <- SummarizedExperiment::rowData(subR25pept[[length(subR25pept)]])$adjacencyMatrix[, indx]
  
obj.agg <- RunAggregation(subR25pept, "Yes_As_Specific", "Sum", "allPeptides",
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Adding aggregated metadata
obj.agg <- RunAggregation(subR25pept, "Yes_As_Specific", "Mean", "allPeptides",
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Adding aggregated metadata
obj.agg <- RunAggregation(subR25pept, "Yes_As_Specific", "Sum", "topN", n = 4,
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Adding aggregated metadata
obj.agg <- RunAggregation(subR25pept, "Yes_As_Specific", "Mean", "topN", n = 4,
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
#> Adding aggregated metadata

obj.agg <- RunAggregation(subR25pept, "No", "Sum", "allPeptides")
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.
obj.agg <- RunAggregation(subR25pept, "No", "Sum", "topN", n = 4)
#> Aggregating data
#> Your quantitative data contain missing values. Please read the relevant
#> section(s) in the aggregateFeatures manual page regarding the effects
#> of missing values on data aggregation.

obj.agg <- RunAggregation(subR25pept, "Yes_Simple_Redistribution", "Sum", "allPeptides",
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Adding aggregated metadata
obj.agg <- RunAggregation(subR25pept, "Yes_Iterative_Redistribution", "Sum", "topN", n = 4,
aggregated_col = c("Sequence", "Mass"))
#> Aggregating data
#> Adding aggregated metadata

library(QFeatures)
#> Loading required package: MultiAssayExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: ‘QFeatures’
#> The following object is masked from ‘package:base’:
#> 
#>     sweep

data(subR25pept)
protID <- parentProtId(subR25pept[[2]])
X <- QFeatures::adjacencyMatrix(subR25pept[[2]])

X.split <- DaparToolshed::splitAdjacencyMat(X)
X.shared <- X.split$Xshared
X.unique <- X.split$Xspec


#adjacencyMatrix(subR25pept[[2]]) <- X.unique
#rowdata.pep <- rowData(subR25pept[[2]])


# subR25pept <- aggregateFeatures4Prostar(
#   object = subR25pept,
#   i = length(subR25pept),
#   name = 'aggregated',
#   fcol = 'adjacencyMatrix',
#   fun = 'colSumsMat')
# 
# 
# .names <- "Sequence"
# 
# proteinNames <- rownames(subR25pept[[2]])
# data <- rowData(subR25pept[[1]])
# 
# new.col <- BuildColumnToProteinDataset(
#   peptideData = rowData(subR25pept[[1]]), 
#   matAdj = adjacencyMatrix(subR25pept[[2]]), 
#   columnName = "Sequence",
#   proteinNames = rownames(subR25pept[[2]]))
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
X.topn <- select_topn(assay(subR25pept[[2]]), X, n = 3)

data(subR25pept)
obj.last <- subR25pept[[2]]
X <- BuildAdjacencyMatrix(subR25pept[[2]])
getProteinsStats(X)
#> $nbPeptides
#> [1] 100
#> 
#> $nbSpecificPeptides
#> [1] 97
#> 
#> $nbSharedPeptides
#> [1] 3
#> 
#> $nbProt
#> [1] 96
#> 
#> $protOnlyUniquePep
#>  [1] "1005" "1017" "103"  "106"  "1060" "1073" "1079" "1101" "111"  "1172"
#> [11] "1210" "1238" "1246" "1261" "1268" "129"  "1421" "1433" "1437" "1469"
#> [21] "150"  "152"  "1652" "1685" "1731" "1733" "1745" "1781" "1791" "182" 
#> [31] "1845" "186"  "1941" "198"  "1987" "199"  "2"    "2008" "204"  "2081"
#> [41] "210"  "2112" "212"  "2165" "2169" "220"  "2232" "225"  "2258" "2260"
#> [51] "2287" "2292" "231"  "2312" "2347" "237"  "239"  "240"  "250"  "254" 
#> [61] "272"  "292"  "311"  "359"  "364"  "378"  "379"  "397"  "40"   "400" 
#> [71] "414"  "423"  "426"  "439"  "440"  "471"  "484"  "534"  "558"  "573" 
#> [81] "599"  "660"  "662"  "683"  "704"  "706"  "800"  "828"  "858"  "973" 
#> 
#> $protOnlySharedPep
#> [1] "1121" "115"  "116"  "1315" "2342" "944" 
#> 
#> $protMixPep
#> character(0)
#> 


data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
CountPep(X)
#> 100 x 96 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 96 column names ‘1005’, ‘1017’, ‘103’ ... ]]
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . 1 . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . 1 1 . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . 1 . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . 1 . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . 1 . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . 1 . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . 1
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . 1 . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . 1 . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . 1 . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . 1 . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . 1 . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . 1 . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          1 . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . 1 . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . 1 . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . 1 .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . 1 . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . 1 . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . 1 . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . 1 . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . 1 . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . 1 . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               1 . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . 1 . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . 1 . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . 1 . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . 1 . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . 1
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . 1 .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . 1 . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . 1 . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . 1 . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . 1 . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . 1 . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . 1 . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . 1 . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . 1 . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . 1 . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    1 . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . 1 . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . 1 .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . 1 . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . 1 . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . 1 . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . 1 . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . 1 . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . 1 . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . 1 . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . 1 . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . 1 . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . 1 . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . 1 . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . 1
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . 1 . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . 1 . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . 1 . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . 1 . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . 1 . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . 1 . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . 1 . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . 1 . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . 1 . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . 1 . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . 1
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . 1 . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . 1 . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . 1 . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . 1 . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . 1 . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . 1 . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . 1 .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . 1 . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . 1 . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . 1 . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . 1 . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               1 . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . 1 . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . 1 . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . 1 . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . 1 . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . 1 . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . 1 .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . 1 . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . 1 . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . 1 . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . 1 . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . 1 . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . 1 . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . 1 . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . 1 . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . 1
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . 1 . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . 1 . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . 1 . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . 1 . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . 1 . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          1 . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . 1 . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . 1 . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . 1 . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                           
#> AAAAQDEITGDGTTTVVCIVGEIIR                .
#> AAADAISDIEIK                             .
#> AAADAISDIEIKDSK                          .
#> AAAEEQAKR                                .
#> AAAEGVANIHIDEATGEMVSK                    .
#> AAAEYEKGEYETAISTINDAVEQGR                .
#> AAAHSSIKEYDQAVK                          .
#> AAAPGIQIVAGEGFQSPIEDR                    .
#> AAAPTVVFIDEIDSIAK                        .
#> AACIVQNGIATWFPIAVTK                      .
#> AADAIIIK                                 .
#> AADETAAAFYPSK                            .
#> AADIINIAK                                .
#> AADIPVVGNAAGHSNDWFDIK                    .
#> AADTPETSDAVHTEQKPEEEKETIQEE              .
#> AAEAATTDITYR                             .
#> AAEAGETGAATSATEGDNNNNTAAGDK              .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             .
#> AAEEADADAEIADEEDAIHDEI                   .
#> AAEIDVINDPK                              .
#> AAEIIIENR                                .
#> AAEIIISDQDNVIPK                          .
#> AAENASNAIAETR                            .
#> AAESIISIANVPDGDSR                        .
#> AAFDEDGNISNVK                            .
#> AAFISAIVGK                               .
#> AAFNGVTFK                                .
#> AAFNYQFDSIIEHSEK                         .
#> AAFTECCQAADK                             .
#> AAFTYIINDPEIAK                           .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK .
#> AAGANVDNVWADVYAK                         .
#> AAGFIIEK                                 .
#> AAGGVVEIIA                               .
#> AAGIFVSTSSYGGGQESTVK                     .
#> AAGITAAYAR                               .
#> AAGIVDIATVISTSAYIIESK                    .
#> AAGKEIGDFEDISTENEK                       .
#> AAIANVYEYR                               .
#> AAIDFYTK                                 .
#> AAIEAGAFEAVTSNHWAEGGK                    .
#> AAIEDGPVTAENISSETAR                      .
#> AAIEDGWVPGK                              .
#> AAIEEIVK                                 .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            .
#> AAIGSSPINFPSSSQR                         .
#> AAIIACAAEYIQK                            .
#> AAIIGSIGSIFK                             .
#> AAIINQYFAQAYK                            .
#> AAIISSGNVK                               .
#> AAINDPAKAPIIINNIIDSGIR                   .
#> AAIQTYIPK                                .
#> AAISFGAKPEEQK                            .
#> AAITDFER                                 .
#> AAITIIQFDGTGTR                           .
#> AAIVQIDATPFR                             .
#> AAIYAIHSIGCK                             .
#> AANAIKDIYGWTQTSIDDYPIK                   .
#> AANAPVYVTR                               .
#> AANIAHDNQTTVEAYK                         .
#> AANIGGVAVSGIEMAQNSQR                     .
#> AANQGAIPPDISIIVK                         .
#> AANQTASSIVDFYNAIGDDEEEK                  .
#> AANSHRIIDIQESQANCSHFFIEPIK               .
#> AAPIVDDEETEFDIYNSK                       .
#> AAPSPISHVAEIR                            .
#> AAQDVWNR                                 .
#> AAQIGFNTACVEK                            .
#> AAQIGSSFIAQIK                            .
#> AAQQQWGNDFYK                             .
#> AASDAIPPASPK                             .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              .
#> AASKPFIETFICGR                           .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          .
#> AASSINRVDTIR                             .
#> AASYADKINTPEIWSQIGTAQIDGIR               .
#> AATIIPQFVGIK                             .
#> AATIISNII                                .
#> AATVVATSDCIIWAIDR                        .
#> AAVDCECEFQNIEHNEK                        .
#> AAVEEGIIPGGGTAIVK                        .
#> AAVPFNREQIESVIR                          .
#> AAVSGKPYFFFGSDSAPHPVQNK                  .
#> AAWWSPTGDYIAFIK                          .
#> AAYAIGGIGSGFANNEK                        .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    .
#> AAYGDISDEEEK                             .
#> AAYSYMFDSIR                              .
#> ACAAQTNATFIK                             .
#> ACDTSNDNFPIQYDGSK                        .
#> ACGIFSGYPDTFK                            1
#> ACGIIISEER                               .
#> ACGVSRPVIAASITTNDASAIK                   .
#> ACNFQFPEIAYPGK                           .
#> ACPVGNEAGVTTSIR                          .
#> ACVIVVSDIK                               .
#> ACVVYGGSPIGNQIR                          .
#> ADAEWVQSTASK                             .
#> ADASGEGVEDEASGVHK                        .

# \donttest{
library(QFeatures)
data(subR25pept)
X <- BuildAdjacencyMatrix(data(subR25pept)[[2]])
#> Error in data(subR25pept)[[2]]: subscript out of bounds
GetNbPeptidesUsed(assay(subR25pept), X)
#> 96 x 6 Matrix of class "dgeMatrix"
#>      Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2
#> 1005              1              1              1              1              1
#> 1017              1              1              1              1              1
#> 103               1              1              0              1              1
#> 106               1              1              1              1              1
#> 1060              1              1              1              1              1
#> 1073              1              1              0              0              1
#> 1079              1              1              1              1              1
#> 1101              1              1              1              1              1
#> 111               1              1              1              1              1
#> 1121              1              1              1              1              1
#> 115               1              1              1              1              1
#> 116               1              1              1              1              1
#> 1172              1              1              0              1              1
#> 1210              1              1              0              1              1
#> 1238              0              0              1              0              0
#> 1246              1              1              1              1              1
#> 1261              1              1              1              1              1
#> 1268              1              1              1              1              1
#> 129               1              1              1              1              1
#> 1315              1              1              1              1              1
#> 1421              0              0              0              0              0
#> 1433              2              2              1              1              2
#> 1437              1              0              0              0              0
#> 1469              1              1              1              1              1
#> 150               1              1              1              1              1
#> 152               1              1              1              1              1
#> 1652              0              0              0              0              0
#> 1685              1              1              1              1              1
#> 1731              1              1              1              1              1
#> 1733              1              1              1              1              1
#> 1745              1              1              1              0              1
#> 1781              1              0              1              1              0
#> 1791              1              1              1              1              1
#> 182               1              1              1              1              1
#> 1845              0              1              0              1              1
#> 186               1              1              1              1              1
#> 1941              1              1              1              1              1
#> 198               0              0              0              0              1
#> 1987              1              1              1              1              1
#> 199               1              1              1              1              1
#> 2                 1              1              1              1              1
#> 2008              1              1              1              1              1
#> 204               1              1              1              0              1
#> 2081              0              1              1              1              1
#> 210               2              2              2              2              1
#> 2112              2              2              2              2              2
#> 212               1              1              1              1              1
#> 2165              1              1              1              1              1
#> 2169              1              0              1              0              1
#> 220               1              1              1              1              1
#> 2232              1              1              1              1              1
#> 225               1              1              1              1              1
#> 2258              2              2              2              2              2
#> 2260              1              1              1              1              1
#> 2287              0              0              0              0              0
#> 2292              1              1              1              1              1
#> 231               1              1              1              1              1
#> 2312              1              0              1              1              0
#> 2342              1              1              1              1              1
#> 2347              1              0              0              0              0
#> 237               1              1              1              1              1
#> 239               0              1              0              0              1
#> 240               1              1              1              1              1
#> 250               1              1              1              1              1
#> 254               2              2              2              2              2
#> 272               1              1              1              1              1
#> 292               1              1              1              1              1
#> 311               0              1              0              0              0
#> 359               1              1              1              1              1
#> 364               1              1              1              1              1
#> 378               1              1              1              1              1
#> 379               1              1              1              1              1
#> 397               1              1              1              1              1
#> 40                0              0              0              1              1
#> 400               1              1              1              0              1
#> 414               1              1              1              1              1
#> 423               1              1              1              1              1
#> 426               2              2              2              2              2
#> 439               1              1              1              1              1
#> 440               0              1              1              0              0
#> 471               1              1              1              1              1
#> 484               1              1              1              1              1
#> 534               0              1              0              0              0
#> 558               1              1              1              1              1
#> 573               2              2              2              2              2
#> 599               1              1              1              1              1
#> 660               1              0              0              1              0
#> 662               0              0              0              0              0
#> 683               1              1              1              1              1
#> 704               1              0              1              1              0
#> 706               1              0              1              1              0
#> 800               1              1              1              0              0
#> 828               1              1              1              1              1
#> 858               1              1              0              1              1
#> 944               1              1              1              1              1
#> 973               1              1              1              1              1
#>      Intensity_D_R3
#> 1005              0
#> 1017              1
#> 103               0
#> 106               1
#> 1060              1
#> 1073              1
#> 1079              1
#> 1101              1
#> 111               1
#> 1121              1
#> 115               1
#> 116               1
#> 1172              0
#> 1210              1
#> 1238              0
#> 1246              1
#> 1261              1
#> 1268              1
#> 129               1
#> 1315              1
#> 1421              0
#> 1433              2
#> 1437              0
#> 1469              1
#> 150               1
#> 152               1
#> 1652              0
#> 1685              1
#> 1731              1
#> 1733              1
#> 1745              1
#> 1781              1
#> 1791              1
#> 182               1
#> 1845              1
#> 186               1
#> 1941              1
#> 198               0
#> 1987              1
#> 199               1
#> 2                 1
#> 2008              1
#> 204               1
#> 2081              1
#> 210               2
#> 2112              2
#> 212               1
#> 2165              1
#> 2169              1
#> 220               1
#> 2232              1
#> 225               1
#> 2258              2
#> 2260              1
#> 2287              0
#> 2292              1
#> 231               1
#> 2312              0
#> 2342              1
#> 2347              0
#> 237               1
#> 239               0
#> 240               1
#> 250               1
#> 254               2
#> 272               1
#> 292               1
#> 311               1
#> 359               1
#> 364               1
#> 378               1
#> 379               1
#> 397               1
#> 40                1
#> 400               0
#> 414               1
#> 423               1
#> 426               2
#> 439               1
#> 440               1
#> 471               1
#> 484               1
#> 534               0
#> 558               1
#> 573               2
#> 599               1
#> 660               0
#> 662               0
#> 683               1
#> 704               1
#> 706               0
#> 800               0
#> 828               1
#> 858               1
#> 944               1
#> 973               1
# }

library(SummarizedExperiment)
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
ll.n <- GetDetailedNbPeptidesUsed(assay(subR25pept), X)

data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
n <- GetDetailedNbPeptides(X)

data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
GraphPepProt(X)


data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
ExtractUniquePeptides(X)
#> 100 x 96 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 96 column names ‘1005’, ‘1017’, ‘103’ ... ]]
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . 1 . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . 1 . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . 1 . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . 1 . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . 1 . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . 1
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . 1 . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . 1 . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . 1 . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . 1 . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . 1 . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          1 . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . 1 . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . 1 . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . 1 .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . 1 . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . 1 . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . 1 . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . 1 . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . 1 . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . 1 . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . 1 . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . 1 . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . 1 . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . 1 . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . 1
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . 1 .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . 1 . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . 1 . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . 1 . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . 1 . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . 1 . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . 1 . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . 1 . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . 1 . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . 1 . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    1 . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . 1 . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . 1 .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . 1 . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . 1 . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . 1 . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . 1 . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . 1 . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . 1 . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . 1 . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . 1 . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . 1 . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . 1 . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . 1 . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . 1
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . 1 . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . 1 . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . 1 . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . 1 . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . 1 . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . 1 . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . 1 . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . 1 . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . . . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . 1 . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . 1 . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . . . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . . .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . . . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . 1
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . 1 . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . . . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . . . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . . . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . 1 . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . . . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . . . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . 1 . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . 1 . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . 1 . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . . . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . 1 .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . 1 . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . 1 . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . 1 . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . 1 . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               1 . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . . . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . 1 . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . 1 . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . . . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . 1 . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . . . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . . . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . . . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . . . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . . . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . . . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                                                               
#> AAAAQDEITGDGTTTVVCIVGEIIR                . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIK                             . . . . . . . . . . . . . . . . . . .
#> AAADAISDIEIKDSK                          . . . . . . . . . . . . . . . . . . .
#> AAAEEQAKR                                . . . . . . 1 . . . . . . . . . . . .
#> AAAEGVANIHIDEATGEMVSK                    . . . . . . . . . . . . . . . . . . .
#> AAAEYEKGEYETAISTINDAVEQGR                . . . . . . . . . . . . . . . . . . .
#> AAAHSSIKEYDQAVK                          . . . . . . . . . . . . . . . . . . .
#> AAAPGIQIVAGEGFQSPIEDR                    . . . . . . . . . . . . . . . . . . .
#> AAAPTVVFIDEIDSIAK                        . . . . . . . . . 1 . . . . . . . . .
#> AACIVQNGIATWFPIAVTK                      . . . . . . . . . . . . . . . . . 1 .
#> AADAIIIK                                 . . . . . . . . . . . . . . . . . . .
#> AADETAAAFYPSK                            . . . . . . . . . . . . 1 . . . . . .
#> AADIINIAK                                . . . . . . . . . . . . . . . . . . .
#> AADIPVVGNAAGHSNDWFDIK                    . . . . . . . . . . . . . . . . . . .
#> AADTPETSDAVHTEQKPEEEKETIQEE              . . . . . . . . . . . . . . . . . . .
#> AAEAATTDITYR                             . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDK              . . . . . . . . . . . . . . . . . . .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             . . . . . . . . . . . . . . . . . . .
#> AAEEADADAEIADEEDAIHDEI                   . . . . . . . . . . . . . . . . . . .
#> AAEIDVINDPK                              . . . . . . . . . . . . . . . . . . .
#> AAEIIIENR                                . . . . . . . . . . . . . . . . . . .
#> AAEIIISDQDNVIPK                          . . . . . . . . . . . . . . . . . . .
#> AAENASNAIAETR                            . . . . . . . . . . . . . . . . . . .
#> AAESIISIANVPDGDSR                        . . . . . . . . . . . . . . . . . . .
#> AAFDEDGNISNVK                            . . . . . . . . . . . . . . . . . . .
#> AAFISAIVGK                               . . . . . . . . . . . . . . . . . . .
#> AAFNGVTFK                                . 1 . . . . . . . . . . . . . . . . .
#> AAFNYQFDSIIEHSEK                         . . . . . . . . . . . . . . . . . . .
#> AAFTECCQAADK                             . . . . . . . . . . . . . . . . . . .
#> AAFTYIINDPEIAK                           . . . . . . . . . . . . . . . 1 . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  . . . . . . . . 1 . . . . . . . . . .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK . . . . . . . . 1 . . . . . . . . . .
#> AAGANVDNVWADVYAK                         . . . . . . . . . . . . . . . . . . .
#> AAGFIIEK                                 . . . . . . . . . . . . . . . . . . .
#> AAGGVVEIIA                               . . . . . . . . . . . . . . . . . . .
#> AAGIFVSTSSYGGGQESTVK                     . . . . . . . 1 . . . . . . . . . . .
#> AAGITAAYAR                               . . . . . . . . . . . . . . . . . . .
#> AAGIVDIATVISTSAYIIESK                    . . . . . . . . . . . . . . . . . . .
#> AAGKEIGDFEDISTENEK                       . . . . . . . . . . . . . . . . . . .
#> AAIANVYEYR                               . . . . . . . . . . . . . 1 . . . . .
#> AAIDFYTK                                 . . . . . . . . . . . . . . . . . . .
#> AAIEAGAFEAVTSNHWAEGGK                    . . . . . . . . . . . . . . . . . . .
#> AAIEDGPVTAENISSETAR                      . . . . . . . . . . . . . . . . . . .
#> AAIEDGWVPGK                              . . . . . . . . . . . . . . . . . . .
#> AAIEEIVK                                 . . . . . . . . . . . . . . . . . . .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            . . . . . . . . . . . . . . . . . . .
#> AAIGSSPINFPSSSQR                         . . . . . . . . . . 1 . . . . . . . .
#> AAIIACAAEYIQK                            . . . . . . . . . . . . . . . . . . .
#> AAIIGSIGSIFK                             . . . . . . . . . . . . . . . . . . .
#> AAIINQYFAQAYK                            . . . . . . . . . . . . . . . . . . .
#> AAIISSGNVK                               . . . . . . . . . . . . . . . . . . .
#> AAINDPAKAPIIINNIIDSGIR                   . . . . . . . . . . . . . . . . . . .
#> AAIQTYIPK                                . . . . . . . . . . . . . . . . . . .
#> AAISFGAKPEEQK                            . . . . . . . . . . . . . . . . . . .
#> AAITDFER                                 . . . . . . . . . . . . . . . . . . .
#> AAITIIQFDGTGTR                           . . . . . . . . . . . . . . . . . . .
#> AAIVQIDATPFR                             . . . . . . . . . . . . . . . . . . .
#> AAIYAIHSIGCK                             . . . . . . . . . . . . . . . . . . .
#> AANAIKDIYGWTQTSIDDYPIK                   . . . . . . . . . . . . . . . . . . .
#> AANAPVYVTR                               . . . . . . . . . . . . . . . . . . .
#> AANIAHDNQTTVEAYK                         . . . . . . . . . . . . . . . . . . .
#> AANIGGVAVSGIEMAQNSQR                     . . . . . . . . . . . . . . . . . . .
#> AANQGAIPPDISIIVK                         . . . . . . . . . . . . . . . . . . .
#> AANQTASSIVDFYNAIGDDEEEK                  . . . . . . . . . . . . . . . . . . .
#> AANSHRIIDIQESQANCSHFFIEPIK               . . . . . . . . . . . . . . . . . . .
#> AAPIVDDEETEFDIYNSK                       . . . . 1 . . . . . . . . . . . . . .
#> AAPSPISHVAEIR                            . . . . . . . . . . . . . . . . . . .
#> AAQDVWNR                                 . . . . . . . . . . . . . . . . . . .
#> AAQIGFNTACVEK                            . . . . . . . . . . . . . . . . . . .
#> AAQIGSSFIAQIK                            . . . . . . . . . . . . . . . . . . .
#> AAQQQWGNDFYK                             . 1 . . . . . . . . . . . . . . . . .
#> AASDAIPPASPK                             . . . . . . . . . . . . . . . . . . .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              . . . . . . . . . . . . . . . . . . .
#> AASKPFIETFICGR                           . . . . . . . . . . . . . . . . . . .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          . . . . . . . . . . . . . . . . . . .
#> AASSINRVDTIR                             . . . . . . . . . . . . . . . . . . .
#> AASYADKINTPEIWSQIGTAQIDGIR               . . . . . 1 . . . . . . . . . . . . .
#> AATIIPQFVGIK                             . . . . . . . . . . . . . . . . . . .
#> AATIISNII                                . . . . . . . . . . . . . . . . . . .
#> AATVVATSDCIIWAIDR                        . . . . . . . . . . . . . . . . . . .
#> AAVDCECEFQNIEHNEK                        . . . . . . . . . . . . . . . . . . .
#> AAVEEGIIPGGGTAIVK                        . . 1 . . . . . . . . . . . . . . . .
#> AAVPFNREQIESVIR                          . . . . . . . . . . . . . . . . . . .
#> AAVSGKPYFFFGSDSAPHPVQNK                  . . . 1 . . . . . . . . . . . . . . .
#> AAWWSPTGDYIAFIK                          1 . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEK                        . . . . . . . . . . . . . . . . . . .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    . . . . . . . . . . . . . . . . . . .
#> AAYGDISDEEEK                             . . . . . . . . . . . . . . . . . . .
#> AAYSYMFDSIR                              . . . . . . . . . . . . . . . . . . .
#> ACAAQTNATFIK                             . . . . . . . . . . . . . . . . 1 . .
#> ACDTSNDNFPIQYDGSK                        . . . . . . . . . . . . . . . . . . .
#> ACGIFSGYPDTFK                            . . . . . . . . . . . . . . . . . . .
#> ACGIIISEER                               . . . . . . . . . . . 1 . . . . . . .
#> ACGVSRPVIAASITTNDASAIK                   . . . . . . . . . . . . . . . . . . .
#> ACNFQFPEIAYPGK                           . . . . . . . . . . . . . . . . . . .
#> ACPVGNEAGVTTSIR                          . . . . . . . . . . . . . . . . . . .
#> ACVIVVSDIK                               . . . . . . . . . . . . . . 1 . . . .
#> ACVVYGGSPIGNQIR                          . . . . . . . . . . . . . . . . . . .
#> ADAEWVQSTASK                             . . . . . . . . . . . . . . . . . . .
#> ADASGEGVEDEASGVHK                        . . . . . . . . . . . . . . . . . . .
#>                                           
#> AAAAQDEITGDGTTTVVCIVGEIIR                .
#> AAADAISDIEIK                             .
#> AAADAISDIEIKDSK                          .
#> AAAEEQAKR                                .
#> AAAEGVANIHIDEATGEMVSK                    .
#> AAAEYEKGEYETAISTINDAVEQGR                .
#> AAAHSSIKEYDQAVK                          .
#> AAAPGIQIVAGEGFQSPIEDR                    .
#> AAAPTVVFIDEIDSIAK                        .
#> AACIVQNGIATWFPIAVTK                      .
#> AADAIIIK                                 .
#> AADETAAAFYPSK                            .
#> AADIINIAK                                .
#> AADIPVVGNAAGHSNDWFDIK                    .
#> AADTPETSDAVHTEQKPEEEKETIQEE              .
#> AAEAATTDITYR                             .
#> AAEAGETGAATSATEGDNNNNTAAGDK              .
#> AAEAGETGAATSATEGDNNNNTAAGDKK             .
#> AAEEADADAEIADEEDAIHDEI                   .
#> AAEIDVINDPK                              .
#> AAEIIIENR                                .
#> AAEIIISDQDNVIPK                          .
#> AAENASNAIAETR                            .
#> AAESIISIANVPDGDSR                        .
#> AAFDEDGNISNVK                            .
#> AAFISAIVGK                               .
#> AAFNGVTFK                                .
#> AAFNYQFDSIIEHSEK                         .
#> AAFTECCQAADK                             .
#> AAFTYIINDPEIAK                           .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQK  .
#> AAFVAIGNTYGFITPNIWAEQPIPVSPIDIYSDEASAQKK .
#> AAGANVDNVWADVYAK                         .
#> AAGFIIEK                                 .
#> AAGGVVEIIA                               .
#> AAGIFVSTSSYGGGQESTVK                     .
#> AAGITAAYAR                               .
#> AAGIVDIATVISTSAYIIESK                    .
#> AAGKEIGDFEDISTENEK                       .
#> AAIANVYEYR                               .
#> AAIDFYTK                                 .
#> AAIEAGAFEAVTSNHWAEGGK                    .
#> AAIEDGPVTAENISSETAR                      .
#> AAIEDGWVPGK                              .
#> AAIEEIVK                                 .
#> AAIGCIESIIIAQDAQAWNNTYDINVTPK            .
#> AAIGSSPINFPSSSQR                         .
#> AAIIACAAEYIQK                            .
#> AAIIGSIGSIFK                             .
#> AAIINQYFAQAYK                            .
#> AAIISSGNVK                               .
#> AAINDPAKAPIIINNIIDSGIR                   .
#> AAIQTYIPK                                .
#> AAISFGAKPEEQK                            .
#> AAITDFER                                 .
#> AAITIIQFDGTGTR                           .
#> AAIVQIDATPFR                             .
#> AAIYAIHSIGCK                             .
#> AANAIKDIYGWTQTSIDDYPIK                   .
#> AANAPVYVTR                               .
#> AANIAHDNQTTVEAYK                         .
#> AANIGGVAVSGIEMAQNSQR                     .
#> AANQGAIPPDISIIVK                         .
#> AANQTASSIVDFYNAIGDDEEEK                  .
#> AANSHRIIDIQESQANCSHFFIEPIK               .
#> AAPIVDDEETEFDIYNSK                       .
#> AAPSPISHVAEIR                            .
#> AAQDVWNR                                 .
#> AAQIGFNTACVEK                            .
#> AAQIGSSFIAQIK                            .
#> AAQQQWGNDFYK                             .
#> AASDAIPPASPK                             .
#> AASIIYGIGFSTEAQQQPTNSFSGGWR              .
#> AASKPFIETFICGR                           .
#> AASSIDPIVTDYAVGYFNHISGITFDAVQSK          .
#> AASSINRVDTIR                             .
#> AASYADKINTPEIWSQIGTAQIDGIR               .
#> AATIIPQFVGIK                             .
#> AATIISNII                                .
#> AATVVATSDCIIWAIDR                        .
#> AAVDCECEFQNIEHNEK                        .
#> AAVEEGIIPGGGTAIVK                        .
#> AAVPFNREQIESVIR                          .
#> AAVSGKPYFFFGSDSAPHPVQNK                  .
#> AAWWSPTGDYIAFIK                          .
#> AAYAIGGIGSGFANNEK                        .
#> AAYAIGGIGSGFANNEKEIVDICNVAFSSSPQVIVEK    .
#> AAYGDISDEEEK                             .
#> AAYSYMFDSIR                              .
#> ACAAQTNATFIK                             .
#> ACDTSNDNFPIQYDGSK                        .
#> ACGIFSGYPDTFK                            1
#> ACGIIISEER                               .
#> ACGVSRPVIAASITTNDASAIK                   .
#> ACNFQFPEIAYPGK                           .
#> ACPVGNEAGVTTSIR                          .
#> ACVIVVSDIK                               .
#> ACVVYGGSPIGNQIR                          .
#> ADAEWVQSTASK                             .
#> ADASGEGVEDEASGVHK                        .

data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
qdata.agg <- inner.aggregate.iter(assay(subR25pept[[2]]), X)

library(QFeatures)
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[2]])
i.sum <- inner.sum(assay(subR25pept[[2]]), X)

# \donttest{
library(QFeatures)
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept)
#> Error in BuildAdjacencyMatrix(subR25pept): inherits(obj.pep, "SummarizedExperiment") is not TRUE
i.mean <- inner.mean(assay(subR25pept), X)
# }

# \donttest{
library(QFeatures)
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept)
#> Error in BuildAdjacencyMatrix(subR25pept): inherits(obj.pep, "SummarizedExperiment") is not TRUE
i.mean <- inner.median(assay(subR25pept[[2]]), X)
# }

# \donttest{
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept)
#> Error in BuildAdjacencyMatrix(subR25pept): inherits(obj.pep, "SummarizedExperiment") is not TRUE
i.mean <- inner.medianpolish(assay(subR25pept[[2]]), X)
# }

# \donttest{
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept)
#> Error in BuildAdjacencyMatrix(subR25pept): inherits(obj.pep, "SummarizedExperiment") is not TRUE
i.mean <- inner.robustSummary(assay(subR25pept[[2]]), X)
#> Error in inner.robustSummary(assay(subR25pept[[2]]), X): could not find function "inner.robustSummary"
# }
```

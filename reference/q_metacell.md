# Quantitative metadata vocabulary for entities

This function gives the vocabulary used for the quantitative metadata of
each entity in each condition.

This function is based on the qMetacell dataframe to look for either
missing values (used to update an initial dataset) or imputed values
(used when post processing protein qMetacell after aggregation)

In the quantitative columns, a missing value is identified by no value
rather than a value equal to 0. Conversion rules Quanti Tag NA or 0 NA

Update the quantitative metadata information of missing values that were
imputed

Gives all the tags of the metadata vocabulary containing the pattern
(parent and all its children).

Agregation rules for the cells quantitative metadata of peptides. Please
refer to the qMetacell.def vocabulary in `qMetacell.def()`

## Usage

``` r
metacell.def(level)

custom_metacell_colors()

Set_POV_MEC_tags(obj, conds)

Set_POV_MEC_tags2(conds, df, level)

Metacell_generic(qdata, conds, level)

UpdateMetacellAfterImputation(object, ...)

# S4 method for class 'SummarizedExperiment'
UpdateMetacellAfterImputation(object, from, to, ...)

search.metacell.tags(pattern, level, depth = "1")

metacombine(met, level)
```

## Arguments

- level:

  A string designing the type of entity/pipeline. Available values are:
  `peptide`, `protein`

- obj:

  An object of class `SummarizedExperiment`

- conds:

  A 1-col dataframe with the condition associated to each sample.

- df:

  An object of class `SummarizedExperiment`

- qdata:

  A matrix of quantitative data

- object:

  An object of class `SummarizedExperiment`

- ...:

  Additional parameters

- from:

  xxx

- to:

  xxx

- pattern:

  The string to search.

- depth:

  Either "0", "1" or "\*".

- met:

  Metacells

## Value

A data.frame containing the different tags and corresponding colors for
the level given in parameter

A list

An instance of class `QFeatures`.

NA

NA

NA

## Glossary

Peptide-level vocabulary

\|– 'Any' \| \| \| \|– 1.0 'Quantified' \| \| \| \| \| \|– 1.1 "Quant.
by direct id" (color 4, white) \| \| \| \| \| \|– 1.2 "Quant. by
recovery" (color 3, lightgrey) \| \| \| \|– 2.0 "Missing" (no color) \|
\| \| \| \| \|– 2.1 "Missing POV" (color 1) \| \| \| \| \| \|– 2.2
'Missing MEC' (color 2) \| \| \| \|– 3.0 'Imputed' \| \| \| \| \| \|–
3.1 'Imputed POV' (color 1) \| \| \| \| \| \|– 3.2 'Imputed MEC' (color
2)

Protein-level vocabulary: \|– 'Any' \| \| \| \|– 1.0 'Quantified' \| \|
\| \| \| \|– 1.1 "Quant. by direct id" (color 4, white) \| \| \| \| \|
\|– 1.2 "Quant. by recovery" (color 3, lightgrey) \| \| \| \|– 2.0
"Missing" \| \| \| \| \| \|– 2.1 "Missing POV" (color 1) \| \| \| \| \|
\|– 2.2 'Missing MEC' (color 2) \| \| \| \|– 3.0 'Imputed' \| \| \| \|
\| \|– 3.1 'Imputed POV' (color 1) \| \| \| \| \| \|– 3.2 'Imputed MEC'
(color 2) \| \| \| \|– 4.0 'Combined tags' (color 3bis, lightgrey)

## Conversion to the glossary

A generic conversion

Conversion for Proline datasets

Conversion from Maxquant datasets

## Basic agreagtion

Agregation of non imputed values (2.X) with quantitative values

|                      |
|----------------------|
| (1.0, 1.X, 3.0, 3.X) |
| Not possible         |
| —————————-           |
|                      |

|                                                                  |
|------------------------------------------------------------------|
| Agregation of different types of missing values (among 2.1, 2.2) |
|                                                                  |

- Agregation of 2.1 peptides between each other gives a missing value
  non imputed (2.0)

- Agreagtion of 2.2 peptides between each other givesa missing value non
  imputed (2.0)

- Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed
  (2.0) \|—————————-

|                                                                            |
|----------------------------------------------------------------------------|
| Agregation of a mix of quantitative values (among 1.0, 1.1, 1.2, 3.0, 3.X) |
|                                                                            |

- if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, then
  the final metadata is set the this tag

- if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, then
  the final metadata is set to 1.0

- if the set of metacell to agregate is a mix of 3.X and 3.0, then the
  final metadata is set to 3.0

- if the set of metacell to agregate is a mix of 3.X and 3.0 and other
  (X.X), then the final metadata is set to 4.0 \|—————————-

## Post processing

Update metacell with POV/MEC status for the categories 2.0 and 3.0 TODO

## Author

Thomas Burger, Samuel Wieczorek

Samuel Wieczorek

## Examples

``` r
metacell.def('protein')
#>                   node     parent   color
#> 1                  Any            #FFFFFF
#> 2           Quantified        Any #0A31D0
#> 3  Quant. by direct id Quantified #6178D9
#> 4   Quant. by recovery Quantified #B9C4F2
#> 5              Missing        Any #CF8205
#> 6          Missing POV    Missing #E5A947
#> 7          Missing MEC    Missing #F1CA8A
#> 8              Imputed        Any #A40C0C
#> 9          Imputed POV    Imputed #E34343
#> 10         Imputed MEC    Imputed #F59898
#> 11       Combined tags        Any #1E8E05
metacell.def('peptide')
#>                   node     parent   color
#> 1                  Any            #FFFFFF
#> 2           Quantified        Any #0A31D0
#> 3  Quant. by direct id Quantified #6178D9
#> 4   Quant. by recovery Quantified #B9C4F2
#> 5              Missing        Any #CF8205
#> 6          Missing POV    Missing #E5A947
#> 7          Missing MEC    Missing #F1CA8A
#> 8              Imputed        Any #A40C0C
#> 9          Imputed POV    Imputed #E34343
#> 10         Imputed MEC    Imputed #F59898

library(QFeatures)
data(subR25prot)
conds <- design.qf(subR25prot)$Condition
df <- Set_POV_MEC_tags(subR25prot[[2]], conds)

library(SummarizedExperiment)
data(subR25pept)
conds <- design.qf(subR25pept)$Condition
qdata <- assay(subR25pept[[2]])
df <- Metacell_generic(qdata, conds, 'peptide')

data(subR25prot)
subR25prot[[2]] <- UpdateMetacellAfterImputation(subR25prot[[2]], 'Missing', 'Imputed')

search.metacell.tags('Missing POV', 'peptide')
#> [1] "Missing POV"
search.metacell.tags('Quantified', 'peptide')
#> [1] "Quantified"          "Quant. by direct id" "Quant. by recovery" 

# \donttest{
ll <- omXplore::metacell.def('peptide')$node
for (i in 1:length(ll))
test <- lapply(combn(ll, i, simplify = FALSE), 
function(x) tag <- metacombine(x, 'peptide'))
# }
```

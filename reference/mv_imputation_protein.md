# Finds the LAPALA

Methods available are:

- wrapper.impute.detQuant(): This method is a wrapper of the function
  `impute.detQuant()` for objects of class `MSnSet`

- wrapper.impute.KNN(): Can impute only POV missing values. This method
  is a wrapper for objects of class `QFeatures` and imputes missing
  values with a fixed value. This function imputes the missing values
  condition by condition.

- wrapper.impute.slsa(): Imputation of peptides having no values in a
  biological condition. This method is a wrapper to the function
  `impute.slsa()` of the package `imp4p` adapted to an object of class
  `MSnSet`.

- wrapper.impute.fixedValue(): This method is a wrapper to objects of
  class `MSnSet` and imputes missing values with a fixed value.

- wrapper.impute.pa(): Imputation of peptides having no values in a
  biological condition. This method is a wrapper to the function
  `impute.pa` of the package `imp4p` adapted to an object of class
  `MSnSet`.

## Usage

``` r
findMECBlock(obj, grp)

reIntroduceMEC(obj, grp, MECIndex)

wrapper.impute.KNN(obj = NULL, grp, K)

wrapper.impute.fixedValue(obj, grp, fixVal = 0, na.type)

wrapper.impute.pa(obj = NULL, grp, q.min = 0.025)

wrapper.impute.detQuant(obj, qval = 0.025, factor = 1, na.type)

getQuantile4Imp(qdata, qval = 0.025, factor = 1)

wrapper.impute.slsa(obj = NULL, design = NULL)
```

## Arguments

- obj:

  An object of class `QFeatures`.

- grp:

  xxx

- MECIndex:

  A data.frame that contains index of MEC (see findMECBlock)

- K:

  the number of neighbors.

- fixVal:

  A float.

- na.type:

  A string which indicates the type of missing values to impute.
  Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.

- q.min:

  Same as the function `impute.pa()` in the package `imp4p`

- qval:

  An expression set containing quantitative values of various replicates

- factor:

  A scaling factor to multiply the imputation value with

- qdata:

  xxx

- design:

  xxx

## Value

A data.frame containing the indexes of LAPALA

A list of two vectors, respectively containing the imputation values and
the rescaled imputation values

## Utilities functions

- findMECBlock(): xxx

- reIntroduceMEC(): xxx

- getQuantile4Imp(): Quantile imputation value definition. This method
  returns the q-th quantile of each column of an expression set, up to a
  scaling factor

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25prot)
obj <- subR25prot[[2]]
grp <- design.qf(subR25prot)$Condition
lapala <- findMECBlock(obj, grp)
na.type = c("Missing POV", "Missing MEC")
obj.imp.pov <- wrapper.impute.detQuant(obj, na.type = na.type)
obj.imp.pov <- reIntroduceMEC(obj, grp, lapala)

obj.imp.pov <- wrapper.impute.KNN(obj, grp, 3)
#> Warning: 2 rows with more than 99 % entries missing;
#>  mean imputation used for these rows
#> Warning: 2 rows with more than 99 % entries missing;
#>  mean imputation used for these rows

obj.imp.pov <- wrapper.impute.fixedValue(obj, grp, 0.001, na.type = "Missing POV")
obj.imp.mec <- wrapper.impute.fixedValue(obj, grp, 0.001, na.type = "Missing MEC")
obj.imp.na <- wrapper.impute.fixedValue(
obj, grp, 0.001, 
na.type = c("Missing MEC", "Missing POV"))

obj.imp.pov <- wrapper.impute.pa(obj, grp)

qdata <- SummarizedExperiment::assay(obj)
quant <- getQuantile4Imp(qdata)


```

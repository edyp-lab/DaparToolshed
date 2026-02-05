# Missing values imputation using the LSimpute algorithm.

This method is a wrapper to the function `impute.mi()` of the package
`imp4p` adapted to an object of class `SummarizedExperiment`.

## Usage

``` r
wrapper.dapar.impute.mi(
  obj,
  design,
  nb.iter = 3,
  nknn = 15,
  selec = 600,
  siz = 500,
  weight = 1,
  ind.comp = 1,
  progress.bar = FALSE,
  x.step.mod = 300,
  x.step.pi = 300,
  nb.rei = 100,
  method = 4,
  gridsize = 300,
  q = 0.95,
  q.min = 0,
  q.norm = 3,
  eps = 0,
  methodi = "slsa",
  lapala = TRUE,
  distribution = "unif"
)
```

## Arguments

- obj:

  An object of class `SummarizedExperiment`.

- design:

  xxx

- nb.iter:

  Same as the function `mi.mix` in the package `imp4p`

- nknn:

  Same as the function `mi.mix` in the package `imp4p`

- selec:

  Same as the function `mi.mix` in the package `imp4p`

- siz:

  Same as the function `mi.mix` in the package `imp4p`

- weight:

  Same as the function `mi.mix` in the package `imp4p`

- ind.comp:

  Same as the function `mi.mix` in the package `imp4p`

- progress.bar:

  Same as the function `mi.mix` in the package `imp4p`

- x.step.mod:

  Same as the function `estim.mix` in the package `imp4p`

- x.step.pi:

  Same as the function `estim.mix` in the package `imp4p`

- nb.rei:

  Same as the function `estim.mix` in the package `imp4p`

- method:

  Same as the function `estim.mix` in the package `imp4p`

- gridsize:

  Same as the function `estim.mix` in the package `imp4p`

- q:

  Same as the function `mi.mix` in the package `imp4p`

- q.min:

  Same as the function `impute.pa` in the package `imp4p`

- q.norm:

  Same as the function `impute.pa` in the package `imp4p`

- eps:

  Same as the function `impute.pa` in the package `imp4p`

- methodi:

  Same as the function `mi.mix` in the package `imp4p`

- lapala:

  xxxxxxxxxxx

- distribution:

  The type of distribution used. Values are `unif` (default) or `beta`.

## Value

The `Biobase::exprs(obj)` matrix with imputed values instead of missing
values.

## Author

Samuel Wieczorek

## Examples

``` r
# \donttest{
utils::data(subR25pept)
design <- design.qf(subR25pept)
obj.imp.na <- wrapper.dapar.impute.mi(subR25pept[[2]], design, nb.iter = 1, lapala = TRUE)
#> Warning: tab.imp contains missing values
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R1
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R2
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R3
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R1
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R2
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R3
#> 
#>  Iterations: 
#> 
#>  1 / 1  - Imputation MNAR OK - 
#> Imputation MCAR in progress - 
#> 3 % - 6 % - 9 % - 12 % - 15 % - 18 % - 21 % - 24 % - 27 % - 30 % - 33 % - 36 % - 39 % - 42 % - 45 % - 48 % - 51 % - 54 % - 57 % - 60 % - 63 % - 66 % - 69 % - 72 % - 75 % - 78 % - 81 % - 84 % - 87 % - 90 % - 93 % - 96 % - 99 % - 
obj.imp.pov <- wrapper.dapar.impute.mi(subR25pept[[2]], design, nb.iter = 1, lapala = FALSE)
#> Warning: tab.imp contains missing values
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R1
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R2
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_C_R3
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R1
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R2
#> Warning: Few missing values are found as imputed in the input imputed matrix (<20) in sample Intensity_D_R3
#> 
#>  Iterations: 
#> 
#>  1 / 1  - Imputation MNAR OK - 
#> Imputation MCAR in progress - 
#> 3 % - 6 % - 9 % - 12 % - 15 % - 18 % - 21 % - 24 % - 27 % - 30 % - 33 % - 36 % - 39 % - 42 % - 45 % - 48 % - 51 % - 54 % - 57 % - 60 % - 63 % - 66 % - 69 % - 72 % - 75 % - 78 % - 81 % - 84 % - 87 % - 90 % - 93 % - 96 % - 99 % - 
# }
```

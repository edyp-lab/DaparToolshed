# Computes the adjusted p-values separately on contrast using CP4P

Computes the adjusted p-values separately on contrast using CP4P

## Usage

``` r
separateAdjPval(x, pval.threshold = 1.05, method = 1)
```

## Arguments

- x:

  a proteins x contrasts dataframe of (raw) p-values

- pval.threshold:

  all the p-values above the threshold are not considered. Default is
  1.05 (which is equivalent to have no threshold). Applying a threshold
  nearby 1 can be instrumental to improve the uniformity under the null,
  notably in case of upstream mutliple contrat correction (for
  experienced users only)

- method:

  a method to estimate pi_0, see CP4P

## Value

a proteins x contrasts table of adjusted p-values

## Author

Thomas Burger

## Examples

``` r
data(subR25prot)
obj <- subR25prot[1:5,]
separateAdjPval(
testAnovaModels(
applyAnovasOnProteins(obj, 1), "TukeyHSD")$P_Value)
#> Procedure of Benjamini-Hochberg is used. pi0 is fixed to 1.
#>           V1_pval
#> [1,] 0.7209545661
#> [2,] 0.8916825703
#> [3,] 0.0001321402
#> [4,] 0.1253349546
#> [5,] 0.8916825703
```

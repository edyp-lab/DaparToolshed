# Computes the adjusted p-values on all the stacked contrasts using CP4P

Computes the adjusted p-values on all the stacked contrasts using CP4P

## Usage

``` r
globalAdjPval(x, pval.threshold = 1.05, method = 1, display = TRUE)
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

  method a method to estimate pi_0, see CP4P

- display:

  if T, a calibration plot is diplayed using CP4P

## Value

a proteins x contrasts table of adjusted p-values

## Author

Thomas Burger

## Examples

``` r
data(subR25prot)
obj <- subR25prot[1:5,]
globalAdjPval(testAnovaModels(
applyAnovasOnProteins(obj, 1), "TukeyHSD")$P_Value)

#> Procedure of Benjamini-Hochberg is used. pi0 is fixed to 1.
#>                                                     V1_pval
#> CON__A2I7N1.CON__A2I7N0                        0.7209545661
#> CON__P00761                                    0.8916825703
#> P02768upsedyp.ALBU_HUMAN_upsedyp.CON__P02768.1 0.0001321402
#> CON__P04264                                    0.1253349546
#> CON__P07477                                    0.8916825703
```

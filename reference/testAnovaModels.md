# Applies a statistical test on each element of a list of linear models

Applies a statistical test on each element of a list of linear models

## Usage

``` r
testAnovaModels(aov_fits, test = "Omnibus")
```

## Arguments

- aov_fits:

  a list of linear models, such as those outputted by
  applyAnovasOnProteins

- test:

  a character string among "Omnibus", "TukeyHSD", "TukeySinglestep",
  "TukeyStepwise", "TukeyNoMTC", "DunnettSinglestep", "DunnettStepwise"
  and "DunnettNoMTC". "Omnibus" tests the all-mean equality, the Tukey
  tests compares all pairs of means and the Dunnet tests compare all the
  means to the first one. For multiple tests (Dunnet's or Tukey's) it is
  possible to correct for multiplicity (either with single-step or
  step-wise FWER) or not. All the Tukey's and Dunnet's tests use the
  multcomp package expect for "TukeyHSD" which relies on the stats
  package. "TukeyHSD" and "TukeyStepwise" gives similar results.

## Value

a list of 2 tables (p-values and fold-changes, respecively)

## Author

Thomas Burger

## Examples

``` r
data(subR25prot)
obj <- subR25prot[1:5,]
testAnovaModels(applyAnovasOnProteins(obj, 1))
#> $logFC
#>                                                anova_1way_logFC
#> CON__A2I7N1;CON__A2I7N0                                      NA
#> CON__P00761                                                  NA
#> P02768upsedyp|ALBU_HUMAN_upsedyp;CON__P02768-1               NA
#> CON__P04264                                                  NA
#> CON__P07477                                                  NA
#> 
#> $P_Value
#>                                                anova_1way_pval
#> CON__A2I7N1;CON__A2I7N0                           4.325725e-01
#> CON__P00761                                       8.916825e-01
#> P02768upsedyp|ALBU_HUMAN_upsedyp;CON__P02768-1    2.052406e-05
#> CON__P04264                                       5.013337e-02
#> CON__P07477                                       8.704512e-01
#> 
```

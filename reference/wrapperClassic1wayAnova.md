# Wrapper for One-way Anova statistical test

Wrapper for One-way Anova statistical test

## Usage

``` r
wrapperClassic1wayAnova(obj, i, with_post_hoc = "No", post_hoc_test = "No")
```

## Arguments

- obj:

  An object of class `QFeatures`.

- i:

  xxx

- with_post_hoc:

  a character string with 2 possible values: "Yes" and "No" (default)
  saying if function must perform a Post-Hoc test or not.

- post_hoc_test:

  character string, possible values are "No" (for no test; default
  value) or TukeyHSD" or "Dunnett". See details of
  [`postHocTest()`](https://edyp-lab.github.io/DaparToolshed/reference/postHocTest.md)
  function to choose the appropriate one.

## Value

A list of two dataframes. First one called "logFC" contains all pairwise
comparisons logFC values (one column for one comparison) for each
analysed feature (Except in the case without post-hoc testing, for which
NAs are returned.); The second one named "P_Value" contains the
corresponding p-values.

## Details

This function allows to perform a 1-way Analysis of Variance. Also
computes the post-hoc tests if the `with_post_hoc` parameter is set to
yes. There are two possible post-hoc tests: the Tukey Honest Significant
Differences (specified as "TukeyHSD") and the Dunnett test (specified as
"Dunnett").

## See also

postHocTest()

## Author

Hélène Borges

## Examples

``` r
# \donttest{
data(subR25prot)
obj <- subR25prot
filter <- FunctionFilter('qMetacellOnConditions',
  cmd = 'delete',
  mode = 'AtLeastOneCond',
  pattern = c("Missing POV", "Missing MEC"),
  conds = design.qf(obj)$Condition,
  percent = TRUE,
  th = 0.8,
  operator = '<')

obj <- filterFeaturesOneSE(obj, filter)

anovatest <- wrapperClassic1wayAnova(obj, 2)
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       21.32981       21.29457       21.36592       21.31562       20.89947 
#> Intensity_D_R3 
#>       21.37284 
#> [1] 21.32981 21.29457 21.36592 21.31562 20.89947 21.37284
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       32.49615       32.44441       32.50197       32.51763       32.42795 
#> Intensity_D_R3 
#>       32.51065 
#> [1] 32.49615 32.44441 32.50197 32.51763 32.42795 32.51065
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       31.33698       31.27985       31.41931       29.80146       29.74882 
#> Intensity_D_R3 
#>       29.80845 
#> [1] 31.33698 31.27985 31.41931 29.80146 29.74882 29.80845
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       25.64031       26.04866       26.20568       24.84398       25.47695 
#> Intensity_D_R3 
#>       25.45457 
#> [1] 25.64031 26.04866 26.20568 24.84398 25.47695 25.45457
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.29776       29.36039       29.21422       29.33841       29.32371 
#> Intensity_D_R3 
#>       29.23913 
#> [1] 29.29776 29.36039 29.21422 29.33841 29.32371 29.23913
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       23.88920       22.45831       23.30215       22.58146       21.66536 
#> Intensity_D_R3 
#>       22.13327 
#> [1] 23.88920 22.45831 23.30215 22.58146 21.66536 22.13327
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       25.50761       25.09225       25.75728       23.14074             NA 
#> Intensity_D_R3 
#>       23.24802 
#> [1] 25.50761 25.09225 25.75728 23.14074       NA 23.24802
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       22.91019       23.40860       23.23999       22.60914       22.36805 
#> Intensity_D_R3 
#>       22.68710 
#> [1] 22.91019 23.40860 23.23999 22.60914 22.36805 22.68710
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.48653       26.22159       26.37310       26.32929       26.27489 
#> Intensity_D_R3 
#>       26.44597 
#> [1] 26.48653 26.22159 26.37310 26.32929 26.27489 26.44597
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.80292       26.68525       26.84326       24.83697       25.01137 
#> Intensity_D_R3 
#>       25.05821 
#> [1] 26.80292 26.68525 26.84326 24.83697 25.01137 25.05821
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.53202       29.46227       29.49766       28.01600       27.88763 
#> Intensity_D_R3 
#>       27.92360 
#> [1] 29.53202 29.46227 29.49766 28.01600 27.88763 27.92360
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.94224       29.82996       29.90139       28.36328       27.80378 
#> Intensity_D_R3 
#>       28.45505 
#> [1] 29.94224 29.82996 29.90139 28.36328 27.80378 28.45505
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.25609       29.29103       28.93267       27.28577       27.30643 
#> Intensity_D_R3 
#>       27.36294 
#> [1] 29.25609 29.29103 28.93267 27.28577 27.30643 27.36294
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.28859       26.90873       27.03999       25.40581       25.34614 
#> Intensity_D_R3 
#>       25.12522 
#> [1] 27.28859 26.90873 27.03999 25.40581 25.34614 25.12522
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.18028       29.11304       29.18449       27.66221       27.71967 
#> Intensity_D_R3 
#>       27.68953 
#> [1] 29.18028 29.11304 29.18449 27.66221 27.71967 27.68953
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.10099       30.13907       29.90024       28.03407       27.33501 
#> Intensity_D_R3 
#>       28.24500 
#> [1] 30.10099 30.13907 29.90024 28.03407 27.33501 28.24500
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.40410       29.53918       29.58327       27.36553       27.99911 
#> Intensity_D_R3 
#>       27.65581 
#> [1] 29.40410 29.53918 29.58327 27.36553 27.99911 27.65581
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.83886       29.77338       30.00156       28.17188       28.32710 
#> Intensity_D_R3 
#>       28.42009 
#> [1] 29.83886 29.77338 30.00156 28.17188 28.32710 28.42009
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       24.93116       25.20165       24.86377       23.98416       23.54769 
#> Intensity_D_R3 
#>       23.52085 
#> [1] 24.93116 25.20165 24.86377 23.98416 23.54769 23.52085
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.19836       28.15291       28.11600       26.62937       26.46952 
#> Intensity_D_R3 
#>       26.48577 
#> [1] 28.19836 28.15291 28.11600 26.62937 26.46952 26.48577
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.82122       27.72488       28.17902       25.78786       25.60108 
#> Intensity_D_R3 
#>       25.90814 
#> [1] 27.82122 27.72488 28.17902 25.78786 25.60108 25.90814
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.69856       27.62066       27.58025       26.14579       25.93391 
#> Intensity_D_R3 
#>       26.07502 
#> [1] 27.69856 27.62066 27.58025 26.14579 25.93391 26.07502
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.53866       29.43313       29.61507       27.96780       27.85757 
#> Intensity_D_R3 
#>       27.86638 
#> [1] 29.53866 29.43313 29.61507 27.96780 27.85757 27.86638
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.77621       28.59623       28.39693       27.26419       26.99550 
#> Intensity_D_R3 
#>       26.89227 
#> [1] 28.77621 28.59623 28.39693 27.26419 26.99550 26.89227
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.77831       27.68686       27.91456       25.91734       25.98581 
#> Intensity_D_R3 
#>       25.51008 
#> [1] 27.77831 27.68686 27.91456 25.91734 25.98581 25.51008
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.54252       28.52955       28.77850       27.28215       27.52761 
#> Intensity_D_R3 
#>       26.90185 
#> [1] 28.54252 28.52955 28.77850 27.28215 27.52761 26.90185
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.06950       26.05167       25.90626       22.19096       23.97745 
#> Intensity_D_R3 
#>       24.13265 
#> [1] 26.06950 26.05167 25.90626 22.19096 23.97745 24.13265
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.77408       28.61491       28.77533       27.21771       27.20955 
#> Intensity_D_R3 
#>       27.31337 
#> [1] 28.77408 28.61491 28.77533 27.21771 27.20955 27.31337
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       31.44801       31.45865       31.50697       29.90268       29.84244 
#> Intensity_D_R3 
#>       29.87807 
#> [1] 31.44801 31.45865 31.50697 29.90268 29.84244 29.87807
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       31.79750       31.75031       31.79254       30.29603       30.23112 
#> Intensity_D_R3 
#>       30.30715 
#> [1] 31.79750 31.75031 31.79254 30.29603 30.23112 30.30715
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       31.28361       31.22297       31.34154       29.60459       29.54609 
#> Intensity_D_R3 
#>       29.40787 
#> [1] 31.28361 31.22297 31.34154 29.60459 29.54609 29.40787
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.80090       29.29748       29.38412       27.16116       27.44380 
#> Intensity_D_R3 
#>       27.83201 
#> [1] 29.80090 29.29748 29.38412 27.16116 27.44380 27.83201
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.44391       30.34477       30.42642       28.47804       28.57315 
#> Intensity_D_R3 
#>       28.57924 
#> [1] 30.44391 30.34477 30.42642 28.47804 28.57315 28.57924
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.92267       29.88297       29.91755       28.37323       28.36378 
#> Intensity_D_R3 
#>       28.38319 
#> [1] 29.92267 29.88297 29.91755 28.37323 28.36378 28.38319
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.29474       28.13467       28.39681       26.64031       26.62742 
#> Intensity_D_R3 
#>       26.77681 
#> [1] 28.29474 28.13467 28.39681 26.64031 26.62742 26.77681
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.75764       30.26941       30.80139       28.77712       28.68853 
#> Intensity_D_R3 
#>       28.76721 
#> [1] 30.75764 30.26941 30.80139 28.77712 28.68853 28.76721
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.81617       27.74951       27.77228       25.98325       26.25034 
#> Intensity_D_R3 
#>       26.16503 
#> [1] 27.81617 27.74951 27.77228 25.98325 26.25034 26.16503
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.44892       28.38952       28.49720       26.50416       26.55385 
#> Intensity_D_R3 
#>       26.80415 
#> [1] 28.44892 28.38952 28.49720 26.50416 26.55385 26.80415
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.05279       29.02794       29.08922       27.54760       27.62296 
#> Intensity_D_R3 
#>       27.57824 
#> [1] 29.05279 29.02794 29.08922 27.54760 27.62296 27.57824
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.05918       30.02366       30.09068       28.46435       28.25228 
#> Intensity_D_R3 
#>       28.54208 
#> [1] 30.05918 30.02366 30.09068 28.46435 28.25228 28.54208
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.91734       30.78986       31.05006       29.35178       29.36913 
#> Intensity_D_R3 
#>       29.53206 
#> [1] 30.91734 30.78986 31.05006 29.35178 29.36913 29.53206
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.11540       28.08294       28.10164       26.63062       26.51084 
#> Intensity_D_R3 
#>       26.53301 
#> [1] 28.11540 28.08294 28.10164 26.63062 26.51084 26.53301
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.14114       30.10899       30.15037       28.61038       28.42258 
#> Intensity_D_R3 
#>       28.57264 
#> [1] 30.14114 30.10899 30.15037 28.61038 28.42258 28.57264
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.00911       29.03007       29.07523       27.75939       27.72540 
#> Intensity_D_R3 
#>       27.57708 
#> [1] 29.00911 29.03007 29.07523 27.75939 27.72540 27.57708
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.08982       27.16442       27.28710       25.71916       25.71466 
#> Intensity_D_R3 
#>       25.79764 
#> [1] 27.08982 27.16442 27.28710 25.71916 25.71466 25.79764
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.33070       28.23370       28.18386       26.96235       26.49462 
#> Intensity_D_R3 
#>       26.56344 
#> [1] 28.33070 28.23370 28.18386 26.96235 26.49462 26.56344
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       28.23945       28.18281       28.33398       26.76661       26.51134 
#> Intensity_D_R3 
#>       26.60950 
#> [1] 28.23945 28.18281 28.33398 26.76661 26.51134 26.60950
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.04729       26.28995       25.95790       25.16250       24.81787 
#> Intensity_D_R3 
#>       25.09648 
#> [1] 26.04729 26.28995 25.95790 25.16250 24.81787 25.09648
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.26843       29.31361       29.39991       27.71332       27.65294 
#> Intensity_D_R3 
#>       27.81592 
#> [1] 29.26843 29.31361 29.39991 27.71332 27.65294 27.81592
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.96851       27.10410       27.19150       25.62408       25.88697 
#> Intensity_D_R3 
#>       25.78584 
#> [1] 26.96851 27.10410 27.19150 25.62408 25.88697 25.78584
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.71017       29.67528       29.83415       28.25696       27.87496 
#> Intensity_D_R3 
#>       28.08811 
#> [1] 29.71017 29.67528 29.83415 28.25696 27.87496 28.08811
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.82723       26.80243       26.57433       25.58695       24.13108 
#> Intensity_D_R3 
#>       24.06357 
#> [1] 26.82723 26.80243 26.57433 25.58695 24.13108 24.06357
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       29.58048       29.56727       29.60419       28.08542       28.06974 
#> Intensity_D_R3 
#>       28.14107 
#> [1] 29.58048 29.56727 29.60419 28.08542 28.06974 28.14107
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       26.67317       26.59178       26.58104       24.88869       25.11128 
#> Intensity_D_R3 
#>       25.03788 
#> [1] 26.67317 26.59178 26.58104 24.88869 25.11128 25.03788
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.96939       27.79216       28.04020       25.83810       25.77211 
#> Intensity_D_R3 
#>       26.31403 
#> [1] 27.96939 27.79216 28.04020 25.83810 25.77211 26.31403
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       30.23352       30.41195       30.46985       28.62972       28.72923 
#> Intensity_D_R3 
#>       28.71148 
#> [1] 30.23352 30.41195 30.46985 28.62972 28.72923 28.71148
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       27.54156       27.58126       27.65772       26.06069       26.34484 
#> Intensity_D_R3 
#>       26.04053 
#> [1] 27.54156 27.58126 27.65772 26.06069 26.34484 26.04053
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       25.30271       25.14981             NA             NA             NA 
#> Intensity_D_R3 
#>       25.97520 
#> [1] 25.30271 25.14981       NA       NA       NA 25.97520
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>       21.97376       20.95156       23.60399       21.84518       22.41160 
#> Intensity_D_R3 
#>       22.98021 
#> [1] 21.97376 20.95156 23.60399 21.84518 22.41160 22.98021
#> Intensity_C_R1 Intensity_C_R2 Intensity_C_R3 Intensity_D_R1 Intensity_D_R2 
#>             NA             NA             NA       21.37597             NA 
#> Intensity_D_R3 
#>             NA 
#> [1]       NA       NA       NA 21.37597       NA       NA
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': error in evaluating the argument 'x' in selecting a method for function 't': contrasts can be applied only to factors with 2 or more levels
# }
```

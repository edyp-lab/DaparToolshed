# Heuristic to choose the value of the hyperparameter (fudge factor) used to regularize the variance estimator in the likelihood ratio statistic

\#' fudge2LRT: heuristic to choose the value of the hyperparameter
(fudge factor) used to regularize the variance estimator in the
likelihood ratio statistic (as implemented in samLRT). We follow the
heuristic described in \[1\] and adapt the code of the fudge2 function
in the siggene R package. \[1\] Tusher, Tibshirani and Chu, Significance
analysis of microarrays applied to the ionizing radiation response, PNAS
2001 98: 5116-5121, (Apr 24).

## Usage

``` r
fudge2LRT(
  lmm.res.h0,
  lmm.res.h1,
  cc,
  n,
  p,
  s,
  alpha = seq(0, 1, 0.05),
  include.zero = TRUE
)
```

## Arguments

- lmm.res.h0:

  a vector of object containing the estimates (used to compute the
  statistic) under H0 for each connected component. If the fast version
  of the estimator was used (as implemented in this package), lmm.res.h0
  is a vector containing averages of squared residuals. If a fixed
  effect model was used, it is a vector of lm objects and if a mixed
  effect model was used it is a vector or lmer object.

- lmm.res.h1:

  similar to lmm.res.h0, a vector of object containing the estimates
  (used to compute the statistic) under H1 for each protein.

- cc:

  a list containing the indices of peptides and proteins belonging to
  each connected component.

- n:

  the number of samples used in the test

- p:

  the number of proteins in the experiment

- s:

  a vector containing the maximum likelihood estimate of the variance
  for the chosen model. When using the fast version of the estimator
  implemented in this package, this is the same thing as the input
  lmm.res.h1. For other models (e.g. mixed models) it can be obtained
  from samLRT.

- alpha:

  A vector of proportions used to build candidate values for the
  regularizer. We use quantiles of s with these proportions. Default to
  seq(0, 1, 0.05)

- include.zero:

  logical value indicating if 0 should be included in the list of
  candidates. Default to TRUE.

## Value

(same as the fudge2 function of siggene): s.zero: the value of the fudge
factor s0. alpha.hat: the optimal quantile of the 's' values. If s0=0,
'alpha.hat' will not be returned. vec.cv: the vector of the coefficients
of variations. Following Tusher et al. (2001), the optimal 'alpha'
quantile is given by the quantile that leads to the smallest CV of the
modified test statistics. msg: a character string summarizing the most
important information about the fudge factor.

## Author

Thomas Burger, Laurent Jacob

## Examples

``` r
NULL
#> NULL
```

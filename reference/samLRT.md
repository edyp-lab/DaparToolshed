# xxxxxx

This function computes a regularized version of the likelihood ratio
statistic. The regularization adds a user-input fudge factor s1 to the
variance estimator. This is straightforward when using a fixed effect
model (cases 'numeric' and 'lm') but requires some more care when using
a mixed model.

## Usage

``` r
samLRT(lmm.res.h0, lmm.res.h1, cc, n, p, s1)
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

- s1:

  the fudge factor to be added to the variance estimate

## Value

llr.sam: a vector of numeric containing the regularized log likelihood
ratio statistic for each protein. s: a vector containing the maximum
likelihood estimate of the variance for the chosen model. When using the
fast version of the estimator implemented in this package, this is the
same thing as the input lmm.res.h1. lh1.sam: a vector of numeric
containing the regularized log likelihood under H1 for each protein.
lh0.sam: a vector of numeric containing the regularized log likelihood
under H0 for each connected component. sample.sizes: a vector of numeric
containing the sample size (number of biological samples times number of
peptides) for each protein. This number is the same for all proteins
within each connected component.

## Author

Thomas Burger, Laurent Jacob

## Examples

``` r
NULL
#> NULL
```

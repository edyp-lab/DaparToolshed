# Computes the adjusted p-values

This function is a wrapper to the function adjust.p from the `cp4p`
package. It returns the adjusted p-values corresponding to the p-values
of the differential analysis. The adjusted p-values is computed with the
function `p.adjust`{stats}.

## Usage

``` r
diffAnaComputeAdjustedPValues(pval, pi0Method = 1)
```

## Arguments

- pval:

  The result (p-values) of the differential analysis

- pi0Method:

  The parameter pi0.method of the method adjust.p in the package `cp4p`

## Value

The computed adjusted p-values

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25prot)
obj <- subR25prot
# Simulate imputation
obj <- NAIsZero(obj, 1)
obj <- NAIsZero(obj, 2)
allComp <- limmaCompleteTest(
SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), 
comp.type="OnevsOne")
diffAnaComputeAdjustedPValues(pval = allComp$P_Value[, 1])
#> Procedure of Benjamini-Hochberg is used. pi0 is fixed to 1.
#>   [1] 5.257966e-01 9.697768e-01 1.447386e-05 7.396754e-02 9.495707e-01
#>   [6] 1.218519e-01 3.494295e-01 2.780948e-02 9.697768e-01 3.409309e-05
#>  [11] 1.447386e-05 1.112298e-03 9.362379e-05 1.928614e-04 8.856861e-06
#>  [16] 1.392649e-03 5.565793e-04 9.450975e-05 1.782449e-03 2.281753e-05
#>  [21] 1.784511e-04 3.409309e-05 2.628075e-05 5.218615e-04 2.183283e-04
#>  [26] 2.060812e-03 1.615722e-02 3.409309e-05 8.856861e-06 8.856861e-06
#>  [31] 2.628075e-05 1.112298e-03 8.856861e-06 8.856861e-06 6.461401e-05
#>  [36] 3.709151e-04 4.284653e-05 5.434734e-05 8.856861e-06 5.573369e-05
#>  [41] 9.450975e-05 1.176950e-05 2.628075e-05 3.409309e-05 3.409309e-05
#>  [46] 4.326204e-04 5.434734e-05 1.690892e-03 2.628075e-05 1.891760e-04
#>  [51] 1.553848e-04 1.460894e-02 8.856861e-06 3.409309e-05 3.990054e-04
#>  [56] 3.409309e-05 1.553848e-04 6.113792e-01 8.803107e-01 4.637739e-01
#>  [61] 3.494295e-01 9.862644e-01 3.343848e-01 8.166503e-01 5.991249e-01
#>  [66] 1.514868e-01 8.991611e-01 3.874271e-01 9.680626e-01 5.991249e-01
#>  [71] 2.234859e-01 9.862644e-01 8.803107e-01 2.283175e-01 4.490338e-01
#>  [76] 3.877227e-01 1.514868e-01 9.290102e-01 4.930609e-07 9.862644e-01
#>  [81] 5.991249e-01 8.803107e-01 5.991249e-01 5.211488e-01 9.862644e-01
#>  [86] 6.081715e-01 3.874271e-01 3.494295e-01 5.858606e-02 9.862644e-01
#>  [91] 1.371025e-01 5.991249e-01 3.497176e-01 6.989902e-01 5.257966e-01
#>  [96] 3.494295e-01 4.084296e-02 3.874271e-01 8.803107e-01 6.324739e-01
```

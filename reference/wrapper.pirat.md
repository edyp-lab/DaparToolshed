# Missing values imputation using Pirat

This method is a wrapper to the function `pipeline_llkimpute()` of the
package `Pirat` adapted to an object of class `QFeatures` of
`SummarizedExperiment`.

## Usage

``` r
wrapper.pirat(data, adjmat, rnas_ab = NULL, adj_rna_pg = NULL, ...)
```

## Arguments

- data:

  An object of class `QFeatures` or `SummarizedExperiment`. If data is
  of class `QFeatures`, the last assay will be imputed.

- adjmat:

  Adjacency matrix corresponding to the `SummarizedExperiment` or the
  last assay of `QFeatures`.

- rnas_ab:

  Transcriptomic data with sample as row, used only if extension = 'T'.

- adj_rna_pg:

  Adjacency matrix of rna (rows) and peptides or precursors (columns),
  used only if extension = 'T'.

- ...:

  Additional arguments to pass to `my_pipeline_llkimpute()`

## Value

`QFeatures` including a new assay with imputed data or
`SummarizedExperiment` with imputed data.

## Author

Manon Gaudin

## Examples

``` r
data(subR25pept)

# Delete whole empty lines
filter_emptyline <- FunctionFilter("qMetacellWholeLine", cmd = 'delete', pattern = 'Missing MEC')
subR25pept <- filterFeaturesOneSE(object = subR25pept, i = length(subR25pept), name = "Filtered",
              filters = list(filter_emptyline))

subR25pept <- wrapper.pirat(data = subR25pept,
adjmat = SummarizedExperiment::rowData(subR25pept[[length(subR25pept)]])$adjacencyMatrix,
extension = "base")
#> Starting Python environment...
#> Installing pyenv ...
#> Done! pyenv has been installed to '/home/runner/.local/share/r-reticulate/pyenv/bin/pyenv'.
#> Using Python: /home/runner/.pyenv/versions/3.10.19/bin/python3.10
#> Creating virtual environment '/home/runner/.cache/R/basilisk/1.22.0/Pirat/1.4.4/envPirat' ... 
#> + /home/runner/.pyenv/versions/3.10.19/bin/python3.10 -m venv /home/runner/.cache/R/basilisk/1.22.0/Pirat/1.4.4/envPirat
#> Done!
#> Installing packages: pip, wheel, setuptools
#> + /home/runner/.cache/R/basilisk/1.22.0/Pirat/1.4.4/envPirat/bin/python -m pip install --upgrade pip wheel setuptools
#> Installing packages: 'torch==2.0.0', 'numpy==1.24'
#> + /home/runner/.cache/R/basilisk/1.22.0/Pirat/1.4.4/envPirat/bin/python -m pip install --upgrade --no-user 'torch==2.0.0' 'numpy==1.24'
#> Virtual environment '/home/runner/.cache/R/basilisk/1.22.0/Pirat/1.4.4/envPirat' successfully created.
#> Warning: NaNs produced
#> Warning: NaNs produced
#> 
#> Call:
#> stats::lm(formula = log(probs) ~ m_ab_sorted[seq(length(probs))])
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3960 -0.3014  0.2577  0.4306  0.9575 
#> 
#> Coefficients:
#>                                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)                     11.21070    1.42188   7.884 2.74e-11 ***
#> m_ab_sorted[seq(length(probs))] -0.56918    0.06006  -9.477 3.10e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual standard error: 0.6962 on 71 degrees of freedom
#> Multiple R-squared:  0.5585, Adjusted R-squared:  0.5523 
#> F-statistic: 89.81 on 1 and 71 DF,  p-value: 3.104e-14
#> 
```

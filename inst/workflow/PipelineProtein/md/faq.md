### FAQ table of contents


* [Why does the table in experimental design blink during edition?](#why-does-the-table-in-experimental-design-blink-during-edition)
* [How to build a valid experimental design?](#how-to-build-a-valid-experimental-design)
* [Why do the items of the contextual menus for plots remain 'undefined'?](#why-do-the-items-of-the-contextual-menus-for-plots-remain-undefined)
* [Why does my volcano plot look so aligned?](#why-does-my-volcano-plot-look-so-aligned)
* [How to recover differential analysis results?](#how-to-recover-differential-analysis-results)
* [Why isn't it possible to adjust the logFC threshold on each differential analysis comparison?](#why-isnt-it-possible-to-adjust-the-logfc-threshold-on-each-differential-analysis-comparison)
* [How does "Push p-values" works?](#how-does-push-p-values-works)
* [How should I tune the calibration plot options?](#how-should-i-tune-the-calibration-plot-options)
* [Why are the adjusted p-values hidden in the tabular outputs of the differential analysis?](#why-are-the-adjusted-p-values-hidden-in-the-tabular-outputs-of-the-differential-analysis)
* [How the cells metadata tags are structured?](#how-the-cells-metadata-tags-are-structured)
* [How the cells metadata tags are aggregated?](#How-the-cells-metadata-tags-are-aggregated)


### Why does the table in experimental design blink during edition?

When you edit the experimental design (during converting a text file to MSnset or during the update of the design),
it may happen that the cells begin to blink in a random order. Then, no more operation is possible in the table. 
This happens if you edit the cells too fast with respect to the table update speed. We apologize for this caveat : 
this is a known bug of the package used to provide the table. No patch is available yet. The only workaround is to close then reopen Prostar.


### How to build a valid experimental design? 
In Prostar, the differential analysis is devoted to the processing of hierarchical unpaired experimental designs. 
However, in former versions, this was not explicit enough, so that users with paired samples could used Prostar with wrong assumptions. 
To clear this out, we have changed the experimental design construction step so that its explicitly appears unpaired.

As a result, the samples must now be numbered as in the following example:

* Condition 1: 1 - 2 - 3 - 4,
* Condition 2: 5 - 6 - 7 - 8

As opposed to:

* Condition 1: 1 - 2 - 3 - 4,
* Condition 2: 1 - 2 - 3 - 4

Which, depending on the context, could suggest that the 8 samples comes only from 4 different biological subjects, 
and thus leading to paired tests - for instance, patients that are compared between Before (Condition 1) and After (Condition 2) some treatment.

However, one should note that even if the experimental design now looks different, 
this is just due to a numbering convention, and the statistical test is not impacted.


### Why do the items of the contextual menus for plots remain 'undefined'?
This happens if the version of the package 'highcharter' is less or equal to 0.5.0. 
To fix this issue, you should install the devel version of the package by typing the following command in a R console: 
```R
devtools::install_github('jbkunst/highcharter')
```

### Why does my volcano plot look so aligned?

![](https://raw.githubusercontent.com/prostarproteomics/Prostar/master/inst/ProstarApp/www/images/dfPriorIssue.png "Volcano plot resulting from a Limma issue")

In very uncommun situations, one may obtain a bowl shape volcano plot such as depicted above. This is due to using Limma on a dataset for which it is not adapted: Briefly, the numerical values in the quantitative matrix appears to have a repetitive pattern that prevent Limma routines to compute the number of degrees of freedom of the Chi2 distribution on which the protein variances should be fitted. As a result, Limma returns a result directly proportional to the fold-change, and the p-values are none-informative. In such cases, which are fortunately extremely odd, we advise to replace Limma test by a classical t-test.

### How to recover differential analysis results?

From Prostar 1.14, the differential analysis results are not exported anymore when using the "Export fo file" functionality, regardless the format (MSnset, excel or zipped CSV). This is due to the separate management of the "data mining" and "data processing" outputs. As a result, after performing the differential analysis, the results must be downloaded thanks to the devoted buttons (otherwise, they will be lost when closing the Prostar session). However, all the p-value computatations (from the "hypothesis testing" menu) can be exported and recovered from one session to another one.

### Why isn't it possible to adjust the logFC threshold on each differential analysis comparison?
Shortly, because from a statistical viewpoint, doing so roughly amount to FDR cheating. We have observed that numerous practitioners use the logFC threshold as a way to discard some proteins on the volcano plot, so that other proteins of interests appear more strikingly. In addition to be an uncontrolled and subjective way of sorting the proteins regardless of p-values, it has an important side effect on FDR computation: FDR computation requires a sufficiently large amount of proteins "below" the horizontal threshold on the volcano plot. However, all the proteins filtered out because of a too low logFC are not considered in FDR computation, so that tuning the logFC threshold to a too high value (so as to fine tune the protein selection) may lead to a spurious FDR. Finally, in case of more than two conditions (say A, B and C) , it would not make sense to defined differently (i.e. with a different logFC threshold) differentialy abundant proteins across various comparisons (e.g. when comparing AvsB, BvsC and CvsA). For all these reasons, we advise Prostar user to define once and for all the logFC threshold of each proteomics experiment (to a minimal value, such that below the threshold, a protein cannot be interesting froma biological viewpoint, because de FC cannot be properly exploited). More detailed explanations can be found in the following articles:   

* <a href="http://prabig-prostar.univ-lyon1.fr/Articles/OnTheMissuseOfFudgeFactorInProteomics_FV.pdf" target="_blank">_Q. Giai Gianetto, Y. Coute, C. Bruley and T. Burger_. Uses and misuses of the fudge factor in quantitative discovery proteomics. **Proteomics**, 16(14):1955-60, 2016.</a>  
* <a href="http://prabig-prostar.univ-lyon1.fr/Articles/fdrtuto.pdf" target="_blank">_T. Burger_. Gentle introduction to the statistical foundations of false discovery rate in quantitative proteomics. **Journal of Proteome Research**, 17(1):12-22, 2017.</a>
* <a href="https://sites.google.com/site/thomasburgerswebpage/download/5tips-VF-accepted.pdf" target="_blank">_S. Wieczorek, Q. Giai Gianetto, T. Burger_. Five simple yet essential steps to correctly estimate the rate of false differentially abundant proteins in mass spectrometry analyses.  **Journal of Proteomics**, vol 207, p. 103441, In Press, 2019</a> (Section 2).


### How does "Push p-values" works?
This is now clearly explained and illustrated in <a href="https://www.bioconductor.org/packages/release/bioc/vignettes/Prostar/inst/doc/Prostar_UserManual.pdf" target="_blank">Prostar user manual</a> as well as in the following article (Section 3):

<a href="https://sites.google.com/site/thomasburgerswebpage/download/5tips-VF-accepted.pdf" target="_blank">_S. Wieczorek, Q. Giai Gianetto, T. Burger_. Five simple yet essential steps to correctly estimate the rate of false differentially abundant proteins in mass spectrometry analyses.  **Journal of Proteomics**, vol 207, p. 103441, In Press, 2019.</a>

### How should I tune the calibration plot options?
This is now clearly explained and illustrated in <a href="https://www.bioconductor.org/packages/release/bioc/vignettes/Prostar/inst/doc/Prostar_UserManual.pdf" target="_blank">Prostar user manual</a> as well as in the following article (Sections 4, 5 and 6):

<a href="https://sites.google.com/site/thomasburgerswebpage/download/5tips-VF-accepted.pdf" target="_blank">_S. Wieczorek, Q. Giai Gianetto, T. Burger_. Five simple yet essential steps to correctly estimate the rate of false differentially abundant proteins in mass spectrometry analyses.  **Journal of Proteomics**, vol 207, p. 103441, In Press, 2019.</a>


### Why are the adjusted p-values hidden in the tabular outputs of the differential analysis?
Adjusted p-values (a. k. a. q-values) are often misunderstood. Because of their names, Prostar users often assume that an adjusted p-value is a protein-level piece of information, which quantitative value differs from the raw p-value, because of a mere "mathematical correction factor", yet with essentially the same interpretation (that of a probabilistic quantification related to the DA status of the protein). However, this does not hold: an adjusted p-value is a list-related piece of information, which is exactly equivalent to an FDR: The adjusted p-value of the N-th protein with the smallest p-value equates the FDR when the list of putative DA proteins is cut to a length of N proteins. Consequently, a same protein with a same p-value, may have different adjusted p-values depending on the other p-values of the datasets (those which are smaller). This makes the adjusted p-value interpretation a bit touchy. In addition, displaying the adjusted p-value information in a tabular sheet were rows could be ordered in different ways (such as Excel spreadsheets for instance) is confusing: If the DA proteins are not sorted by increasing or decreasing p-values, their adjusted p-values cannot be related to the set of DA proteins they apply to. Thus, displaying them is both non-informative and error-prone. After having observed (or having received questions from) Prostar users who were misled by the adjusted p-value information, we have decided to remove it from the tabular outputs. An alternative would be to display the adjusted p-values, but with a different names, easier to interpret, such as "minimum FDR threshold at which the protein would be selected". However, displaying this information promotes p-value hacking, as users may be tempted to spot their proteins of interest and tune the FDR accordingly. As we search for a better solution, any suggestion or comment is welcomed!

### How the cells metadata tags are structured?

#### Prostar prior to 1.30.0

Peptide-level vocabulary

```
──Any
  │── 1.0 Quantitative Value
  │    │── 1.1 Identified
  │    │── 1.2 Recovered
  │
  │── 2.0 Missing value
  │    │── 2.1 Missing POV
  │    │── 2.2 Missing MEC
  │
  │── 3.0 Imputed value
  │    │── 3.1 Imputed POV
  │    │── 3.2 Imputed MEC

```

Protein-level vocabulary:
```
──Any
  │── 1.0 Quantitative Value
  │    │── 1.1 Identified
  │    │── 1.2 Recovered
  │
  │── 2.0 Missing value
  │    │── 2.1 Missing POV
  │    │── 2.2 Missing MEC
  │
  │── 3.0 Imputed value
  │    │── 3.1 Imputed POV
  │    │── 3.2 Imputed MEC
  │
  │── 4.0 Combined value
```

#### Prostar from 1.31.2


<img src="https://github.com/prostarproteomics/Prostar_website/raw/master/docs/img/metacellTags.png"   width="600"/>

### How the cells metadata tags are aggregated?

A set of rules define how the cells metadata tags are aggregated:
Each example below represent a set of three peptides which are aggregated into one protein.

Note: Before an aggregation, a dataset does not contain:

* x.0 tags: 'Quantified' (1.0), 'Missing' (2.0) nor 'Imputed' (3.0)
* 'Combined tags' tag

These tags are specific to the aggregation process and appear only in aggregated proteins.

**Rule 1: Aggregation of a unique type of quantified values (among 1.1, 1.2)**. 
If the type of all the peptides to agregate is either 1.1 or 1.2, then the resulting metadata tag is set to the corresponding tag

**Rule 2: Aggregation of a unique type of missing values (among 2.1, 2.2)**.
The aggregation of 2.1 peptides between each other gives a generic missing value (2.0)

**Rule 3: Aggregation of a unique type of imputed values (among 3.3, 3.2)**.
If the type of all the peptides to agregate is either 3.0, 3.1 or 3.2, then the resulting metadata tag is set to an  'Imputed' value (3.0).


**Rule 4: Aggregation of a mix of quantified values (1.x)**.
If the set of cell matadata tags to agregate is a mix of 1.x, then the final metadata is set to 1.0.


**Rule 5: Aggregation of a mix of missing values (2.x)**.
If the set of metacell to agregate is a mix of 2.x, then the final metadata is set to 'Missing' (2.0)


**Rule 6: Aggregation of a mix of imputed values (3.x)**.
If the set of cells metadata tags to agregate is a mix of 3.x, then the final metadata is set to 'Imputed' (3.0)


**Rule 7: Aggregation of a mix of missing values (2.x) with any of quantitative and/or imputed values (1.x, 3.x)**.
This case is not possible.


**Rule 8: Agregation of a mix of quantitative values (1.x) and imputed values (3.x)**.
If the set of metacell to agregate is a mix of 3.X and 3.0 and other (1.X), then the final metadata is set to 
'Combined tags' (4.0).

<img src="https://github.com/prostarproteomics/Prostar_website/raw/master/docs/img/metacellTags2.png"   width="800"/>


Once the cell metadata tags have been aggregated, the final step of the process consist in analyzing missing values in order to identify 'missing POV' and 'missing MEC'. Thus, at the end of the agregation process, a dataset may contains

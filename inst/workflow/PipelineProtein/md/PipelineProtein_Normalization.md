---
title: "Workflow Protein: Normalization"
abstract: >
  This page describes the computations done with Imputation.

output:
    BiocStyle::html_document:
        highlight: tango
        toc: false
        number_sections: FALSE
        css: style.css
date: 30 mai 2024
link-citations: true
---

To correct batch or group effects.

### Methods available

**GlobalQuantileAlignment**

This method proposes a normalization of important magnitude that should be cautiously used. It proposes to align the quantiles of all the replicates as described in [Other ref. 1]; practically it amounts to replace abundances by order statistics.

**QuantileCentering**

These methods propose to shift the sample distributions (either all of them at once, or within each condition at a time) to align a specific quantile: the median (under the assumption that up-regulations and down-regulations are equally frequent), the 15% quantile (under the assumption that the signal/noise ratio is roughly the same in all the samples), or any other user's choice.

**MeanCentering**

These methods propose to shift the sample distributions (either all of them at once, or within each condition at a time) to align their means. It is also possible to force unit variance (or not).


**SumByColumns**

These methods propose normalizations of important magnitude that should be cautiously used. It operates on the original scale (not the log2 one) and propose to normalize each abundance by the total abundance of the sample (so as to focus on the analyte proportions among each sample).


**LOESS**

This method proposes to apply a cyclic LOESS [Other ref. 4, 5] normalization to the data (either all of them at once, or on each condition independently). It relates to a combination of multiple regression models. The user can tune the regression span (an higher span smooths the fit more, while a lower span captures more trends).


**vsn**

This method proposes to apply the Variance Stabilization Normalization method [Other ref. 6] to the data (either all of them at once, or on each condition independently). No specific parameters required.

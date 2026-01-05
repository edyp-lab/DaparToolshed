#' @title Computes the FDR corresponding to the p-values of the
#' differential analysis
#' 
#' @description 
#' This function returns the FDR corresponding to the p-values of the differential 
#' analysis.
#'
#' @param adj.pvals The adjusted p-values of the differential analysis 
#'
#' @return The computed FDR value (floating number)
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' NULL
#'
#' @export
#'
diffAnaComputeFDR <- function(adj.pvals) {
  BH.fdr <- max(adj.pvals, na.rm = TRUE)
  return(BH.fdr)
}




#' @title Computes the adjusted p-values
#' 
#' @description 
#' This function is a wrapper to the function adjust.p from the `cp4p` package.
#'  It returns the adjusted p-values corresponding to the p-values of the differential 
#' analysis. The adjusted p-values is computed with the function \code{p.adjust}\{stats\}.
#'
#' @param pval The result (p-values) of the differential analysis 
#'
#' @param pi0Method The parameter pi0.method of the method adjust.p in the 
#' package \code{cp4p}
#'
#' @return The computed adjusted p-values
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
#' diffAnaComputeAdjustedPValues(pval = allComp$P_Value[, 1])
#'
#' @export
#'
diffAnaComputeAdjustedPValues <- function(pval, 
                                          pi0Method = 1) {
  pkgs.require('cp4p')
  
  padj <- cp4p::adjust.p(pval, pi0Method)
  return(padj$adjp[, 2])
}


#' @title Performs a calibration plot on an \code{SummarizedExperiment} object,
#' calling the \code{cp4p} package functions.
#' 
#' @description 
#' This function is a wrapper to the calibration.plot method of the
#' \code{cp4p} package for use with \code{SummarizedExperiment} objects.
#'
#' @param vPVal A dataframe that contains quantitative data.
#'
#' @param pi0Method A vector of the conditions (one condition per sample).
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
#' wrapperCalibrationPlot(allComp$P_Value[, 1])
#'
#' @export
#'
wrapperCalibrationPlot <- function(vPVal, pi0Method = "pounds") {
  if (is.null(vPVal)) {
    return(NULL)
  }
  pkgs.require('cp4p')
  
  p <- cp4p::calibration.plot(vPVal, pi0.method = pi0Method)
  
  return(p)
}




#' @title Plots a histogram of p-values
#'
#' @param pval_ll A vector of the p-values.
#'
#' @param bins A integer indicating the number of cells for the histogram
#'
#' @param pi0 A `float` between 0 and 1 corresponding to the proportion of true null hypotheses.
#'
#' @return A histogram of the p-values with pi0 
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
#' histPValue_HC(allComp$P_Value[1])
#'
#' @export
#' @import highcharter
#'
histPValue_HC <- function(pval_ll, bins = 80, pi0 = 1) {
  pkgs.require('graphics')
  
  h <- graphics::hist(sort(unlist(pval_ll)), breaks = bins, plot = FALSE)
  
  # serieInf <- sapply(h$density, function(x) min(pi0, x))
  # serieSup <- sapply(h$density, function(x) max(0, x - pi0))
  
  serieInf <- vapply(h$density, function(x) min(pi0, x), numeric(1))
  serieSup <- vapply(h$density, function(x) max(0, x - pi0), numeric(1))
  
  hc <- highcharter::highchart() %>%
    hc_chart(type = "column") %>%
    hc_add_series(data = serieSup, name = "p-value density") %>%
    hc_add_series(data = serieInf, name = "p-value density") %>%
    hc_title(text = "P-value histogram") %>%
    hc_legend(enabled = FALSE) %>%
    hc_colors(c("green", "red")) %>%
    hc_xAxis(title = list(text = "P-value"), categories = h$breaks) %>%
    hc_yAxis(
      title = list(text = "Density"),
      plotLines = list(
        list(
          color = "blue", 
          width = 2, 
          value = pi0, 
          zIndex = 5)
      )
    ) %>%
    hc_tooltip(
      headerFormat = "",
      pointFormat = "<b> {series.name} </b>: {point.y} ",
      valueDecimals = 2
    ) %>%
    my_hc_ExportMenu(filename = "histPVal") %>%
    hc_plotOptions(
      column = list(
        groupPadding = 0,
        pointPadding = 0,
        borderWidth = 0
      ),
      series = list(
        stacking = "normal",
        animation = list(duration = 100),
        connectNulls = TRUE,
        marker = list(enabled = FALSE)
      )
    ) %>%
    hc_add_annotation(
      labelOptions = list(
        backgroundColor = "transparent",
        verticalAlign = "top",
        y = -30,
        borderWidth = 0,
        x = 20,
        style = list(
          fontSize = "1.5em",
          color = "blue"
        )
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = 80,
            y = pi0
          ),
          text = paste0("pi0=", pi0)
        )
      )
    )
  return(hc)
}


#' @title Push p-values based on metacell tags
#' 
#' @description 
#' This function allows to push p-values to 1 based on metacell tags.
#'
#' @param obj An object of class \code{QFeatures} or \code{SummarizedExperiment}.
#'            If data is of class \code{QFeatures}, the last assay will be used.
#' @param pvalue A vector of p-values.
#' @param scope A string for scope to use. 
#'              Available values are "WholeLine", "WholeMatrix", "AllCond" and "AtLeastOneCond".
#' @param pattern A vector of tag to use.
#' @param percent A boolean to indicate whether the threshold represent an 
#'                absolute value (percent = FALSE) or a percentage (percent = TRUE).
#' @param threshold A value that corresponds to the threshold value. 
#'                  Either an integer if percent = FALSE, or a float between 0 and 1 of percent = TRUE. 
#' @param conditions A vector of conditions in the dataset. 
#'                   If not provided, the vector `"Condition"` from the column metadata will be used.
#' @param operator A string for operator to use. 
#'                 Available operators are "<=", "<", ">=", ">", "==" and "!=".
#' @param level A string for dataset type. Either "peptide" or "protein"
#'              If not provided, the string obtained from `typeDataset(obj)` will be used.
#' 
#' @return A vector with pushed p-values.
#'
#' @author Manon Gaudin
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
#' pushpvalue(obj, allComp$P_Value[, 1], scope = "WholeMatrix", pattern = c("Missing MEC", "Missing POV"), percent = TRUE, threshold = 0.5, operator = ">=",)
#'
#' @export
#'
pushpvalue <- function(obj,
                       pvalue,
                       scope = "WholeMatrix",
                       pattern = "Imputed MEC",
                       percent = TRUE,
                       threshold = 1,
                       conditions = NULL,
                       operator = ">=",
                       level = NULL){
  if (missing(obj))
    stop("'obj' is required.")
  if (missing(pvalue))
    stop("'pvalue' is required.")
  
  if (inherits(obj, 'QFeatures')){
    i <- length(obj)
    obj.se <- obj[[i]]
    SummarizedExperiment::colData(obj.se) <- SummarizedExperiment::colData(obj)
  } else if (inherits(obj, 'SummarizedExperiment')){
    obj.se <- obj
  } else { stop("'obj' must be of class 'SummarizedExperiment' or 'QFeatures'.") }
  
  if (percent) {
    thresholdtype <- "Percentage"
  } else {
    thresholdtype <- "Count"
  }
  
  if (is.null(conditions)){
    conditions <- SummarizedExperiment::colData(obj.se)$Condition
  }
  
  if (is.null(level)){
    level <- typeDataset(obj.se)
  }
  
  # Get indices
  indices <- GetIndices_FunFiltering(obj.se, 
                                     conds = conditions,
                                     level = level,
                                     pattern = pattern,
                                     type = scope,
                                     percent = thresholdtype, 
                                     op = operator, 
                                     th = threshold)
  
  # Push p-values
  pvalue[indices] <- 1
  return(pvalue)
}

#' @title Identification of differentially abundant peptide/protein
#' 
#' @description 
#' This function allows to identify differentially abundant peptide/protein
#'
#' @param pvalue A vector of p-values.
#' @param logFC A vector of logFC.
#' @param thpvalue A float indicating the p-value threshold.
#' @param thlogFC A float indicating the logFC threshold.
#' 
#' @return A vector indicating which peptide/protein is differentially abundant (1) or not (0). 
#'
#' @author Manon Gaudin
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
#' is.differential(allComp$P_Value[, 1], allComp$logFC[, 1], 0.05, 0.5)
#'
#' @export
#'

is.differential <- function(pvalue,
                           logFC,
                           thpvalue,
                           thlogFC){
  differentialList <- rep(0, length(pval))
  logpval <- -log10(pval)
  
  signifItems <- intersect(which(logpval >= thpvalue),
                           which(abs(logFC) >= thlogFC))
  differentialList[signifItems] <- 1
  
  return(differentialList)
}
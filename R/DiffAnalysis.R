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
#' allComp <- limmaCompleteTest(
#' SummarizedExperiment::assay(obj[[length(obj)]]), design_qf(obj), 
#' comp.type="OnevsOne")
#' diffAnaComputeAdjustedPValues(pval = allComp$P_Value[, 1])
#'
#' @export
#'
diffAnaComputeAdjustedPValues <- function(pval, 
                                          pi0Method = 1) {
  pkgsRequire('cp4p')
  
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
#' allComp <- limmaCompleteTest(
#' SummarizedExperiment::assay(obj[[length(obj)]]), 
#' design_qf(obj), 
#' comp.type="OnevsOne")
#' wrapperCalibrationPlot(allComp$P_Value[, 1])
#'
#' @export
#'
wrapperCalibrationPlot <- function(vPVal, pi0Method = "pounds") {
  if (is.null(vPVal)) {
    return(NULL)
  }
  pkgsRequire('cp4p')
  
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
#' allComp <- limmaCompleteTest(
#' SummarizedExperiment::assay(obj[[length(obj)]]), 
#' design_qf(obj), 
#' comp.type="OnevsOne")
#' histPValue_HC(allComp$P_Value[1])
#'
#' @export
#' @import plotly
#'
histPValue_HC <- function(pval_ll, bins = 80, pi0 = 1) {
  pkgsRequire('graphics')
  
  h <- graphics::hist(sort(unlist(pval_ll)), breaks = bins, plot = FALSE)
  
  serieInf <- vapply(h$density, function(x) min(pi0, x), numeric(1))
  serieSup <- vapply(h$density, function(x) max(0, x - pi0), numeric(1))
  
  df_plot <- data.frame(
    mids = h$mids,
    serieInf = serieInf,
    serieSup = serieSup
  )
  
  p <- plot_ly(df_plot, x = ~mids) |>
    add_trace(y = ~serieInf, type = 'bar', name = 'Below pi0', marker = list(color = 'red')) |>
    add_trace(y = ~serieSup, type = 'bar', name = 'Above pi0', marker = list(color = 'green')) |>
    plotly::layout(
      barmode = 'stack',
      showlegend = FALSE,
      bargap = 0,
      title = list(text = "P-value histogram", x = 0.5),
      xaxis = list(title = "P-value"),
      yaxis = list(title = "Density"),
      shapes = list(
        list(
          type = "line",
          x0 = min(df_plot$mids),
          x1 = max(df_plot$mids),
          y0 = pi0,
          y1 = pi0,
          line = list(color = "blue", width = 2)
        )
      ),
      annotations = list(
        list(
          x = max(df_plot$mids) * 0.8,
          y = pi0+5,
          text = paste0("pi0 = ", pi0),
          xref = "x",
          yref = "y",
          showarrow = FALSE,
          font = list(size = 16, color = "blue"),
          xanchor = "left",
          yanchor = "bottom"
        )
      )
    )
  
  return(p)
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
#' @param value A float, value to assign to the pushed p-value. 
#'              By default, the value is set slightly above 1 to be able to differentiate the pushed value.
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
#' allComp <- limmaCompleteTest(SummarizedExperiment::assay(
#' obj[[length(obj)]]), 
#' design_qf(obj), 
#' comp.type="OnevsOne")
#' pushpvalue(obj, 
#' allComp$P_Value[, 1], 
#' scope = "WholeMatrix", 
#' pattern = c("Missing MEC", "Missing POV"), 
#' percent = TRUE, 
#' threshold = 0.5, 
#' operator = ">=",)
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
                       level = NULL,
                       value = 1.00000000001){
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
  pvalue[indices] <- value
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
#' @return A vector indicating which peptide/protein is differentially 
#' abundant (1) or not (0). 
#'
#' @author Manon Gaudin
#'
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' allComp <- limmaCompleteTest(
#' SummarizedExperiment::assay(obj[[length(obj)]]), 
#' design_qf(obj), 
#' comp.type="OnevsOne")
#' isDifferential(allComp$P_Value[, 1], allComp$logFC[, 1], 0.05, 0.5)
#'
#' @export
#'

isDifferential <- function(pvalue,
                           logFC,
                           thpvalue,
                           thlogFC){
  differentialList <- rep(0, length(pvalue))
  logpval <- -log10(pvalue)
  
  signifItems <- intersect(which(logpval >= thpvalue),
                           which(abs(logFC) >= thlogFC))
  differentialList[signifItems] <- 1
  
  return(differentialList)
}
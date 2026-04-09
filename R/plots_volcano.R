#' @title Volcanoplot of the differential analysis
#'
#' @description
#' #' Plots an interactive volcanoplot after the differential analysis.
#' Typically, the log of Fold Change is represented on the X-axis and the
#' log10 of the p-value is drawn on the Y-axis. When the \code{th_pval}
#' and the \code{th_logfc} are set, two lines are drawn respectively on
#' the y-axis and the X-axis to visually distinguish between differential and
#' non differential data. With the use of the package plotly, a
#' customizable tooltip appears when the user put the mouse's pointer over
#' a point of the scatter plot.
#'
#'
#' @param df A dataframe which contains the following slots :
#' x : a vector of the log(fold change) values of the differential analysis,
#' y : a vector of the p-value values returned by the differential analysis.
#' index : a vector of the rowanmes of the data.
#' This dataframe must has been built with the option stringsAsFactors set
#' to FALSE. There may be additional slots which will be used to show
#' informations in the tooltip. The name of these slots must begin with the
#' prefix "tooltip_". It will be automatically removed in the plot.
#' @param th_pval A floating number which represents the p-value that
#' separates differential and non-differential data.
#' @param th_logfc A floating number which represents the log of the
#' Fold Change that separates differential and non-differential data.
#' @param conditions A list of the names of condition 1 and 2 used for the
#' differential analysis.
#' @param pal A `list` containing 2 color to use for the plot 
#' 
#' 
#' @return An interactive volcanoplot
#' @author Samuel Wieczorek
#' @examples
#' data(subR25prot)
#' obj <- subR25prot
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' allComp <- limmaCompleteTest(
#' SummarizedExperiment::assay(obj[[length(obj)]]), 
#' design_qf(obj), 
#' comp.type="OnevsOne")
#' df <- data.frame(
#' x = allComp$logFC[[1]],
#' y = -log10(allComp$P_Value[[1]]),
#' index = as.character(rownames(obj[[1]]))
#' )
#' tooltipSlot <- c("Fasta_headers", "Sequence_length")
#' df <- cbind(df, SummarizedExperiment::rowData(obj[[1]])[, tooltipSlot])
#' colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
#' if (ncol(df) > 3) {
#'   colnames(df)[4:ncol(df)] <- paste0("tooltip_", 
#'                                      colnames(df)[4:ncol(df)])
#' }
#' cond <- c("25fmol", "10fmol")
#' diffAnaVolcanoplot_rCharts(
#'   df,
#'   th_pval = 2.5,
#'   th_logfc = 1,
#'   conditions = cond
#' )
#'
#' @export
#' 
#' @import plotly
#'
#'
diffAnaVolcanoplot_rCharts <- function(
    df,
    th_pval = 1e-60,
    th_logfc = 0,
    conditions = NULL,
    pal = NULL) {
  stopifnot(inherits(df, "data.frame"))
  pkgsRequire('magrittr')
  
  xtitle <- paste("log2 ( mean(", conditions[2], ") / mean(",
                  conditions[1], ") )",
                  sep = ""
  )
  
  if (is.null(pal)) {
    pal <- list(In = "orange", Out = "gray")
  } else if (length(pal) != 2) {
    warning("The palette must be a list of two items: In and Out.
              Set to default.")
    pal <- list(In = "orange", Out = "gray")
  }
  
  g <- x <- y <- NULL
  df <- cbind(df,
              g = ifelse(df$y >= th_pval & abs(df$x) >= th_logfc, "g1", "g2")
  )
  
  
  i_tooltip <- which(startsWith(colnames(df), "tooltip"))
  txt_tooltip <- apply(df[, i_tooltip, drop = FALSE], 1, function(row) {
    paste(
      paste0("<b>", gsub("tooltip_", "", colnames(df)[i_tooltip]), "</b>: ", row),
      collapse = "<br>"
    )
  })
  
  leftBorder <- data.frame(
    x = c(min(df$x), -th_logfc, -th_logfc),
    y = c(th_pval, th_pval, max(df$y))
  )
  rightBorder <- data.frame(
    x = c(max(df$x), th_logfc, th_logfc),
    y = c(th_pval, th_pval, max(df$y))
  )
  
  title <- NULL
  title <- paste0(conditions[1], "_vs_", conditions[2])
  
  p <- plot_ly()
  
  # Scatter principal
  p <- p |> add_trace(
    data = df,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "markers",
    color = ~g,
    colors = c(pal$In, pal$Out),
    text = ~txt_tooltip,  # pour le tooltip
    hoverinfo = "text",
    showlegend = FALSE
  )
  
  # Lignes verticales
  p <- p |> add_trace(
    data = leftBorder,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "lines",
    line = list(color = "grey", dash = "dash"),
    showlegend = FALSE
  )
  
  p <- p |> add_trace(
    data = rightBorder,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "lines",
    line = list(color = "grey", dash = "dash"),
    showlegend = FALSE
  )
  
  p <- p |> plotly::layout(
    margin = list(t = 60, b = 60),
    title = list(
      text = title,
      x = 0.5,
      xanchor = "center",
      font = list(size = 20, color = "black")
    ),
    xaxis = list(
      title = "logFC",
      zeroline = TRUE,
      zerolinecolor = "grey",
      zerolinewidth = 1
    ),
    yaxis = list(
      title = "-log10(pValue)",
      rangemode = "tozero"
    ),
    hovermode = "closest"
  )
  
  return(p)
}

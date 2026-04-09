
#'
#' @title Density plots of logFC values
#' 
#' @description 
#' This function show the density plots of Fold Change (the same as calculated
#' by limma) for a list of the comparisons of conditions in a differential
#' analysis.
#'
#' @param df_logFC A dataframe that contains the logFC values
#'
#' @param th_logFC The threshold on log(Fold Change) to
#' distinguish between differential and non-differential data
#'
#' @param pal A `vector` of HEX codes for colors.
#'
#' @return A plotly density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' library(SummarizedExperiment)
#' data(subR25pept)
#' # Simulate missing value imputation
#' SummarizedExperiment::assay(subR25pept[[1]])[which(is.na(assay(subR25pept[[1]])))] <- 0
#' 
#' qData <- as.matrix(SummarizedExperiment::assay(subR25pept[[1]]))
#' sTab <- SummarizedExperiment::colData(subR25pept)
#' limma <- limmaCompleteTest(qData, sTab)
#' pal <- ExtendPalette(2, "Dark2")
#' hc_logFC_DensityPlot(limma$logFC, th_logFC = 1, pal = pal)
#'
#' @export
#' 
#' @import plotly
#' @importFrom stats density
#' 
#'
hc_logFC_DensityPlot <- function(
    df_logFC,
    th_logFC = 0,
    pal = NULL
) {
  
  if (th_logFC < 0) {
    warning("The parameter 'th_logFC' must be positive or equal to zero.")
    return(NULL)
  }
  
  if (is.null(df_logFC) || ncol(df_logFC) == 0) {
    return(NULL)
  }
  
  if (is.null(pal)) {
    warning("Color palette set to default.")
    pal <- ExtendPalette(ncol(df_logFC), "Paired")
  } else if (length(pal) != ncol(df_logFC)) {
    warning("The color palette has not the same dimension as the 
              number of samples")
    pal <- ExtendPalette(ncol(df_logFC), "Paired")
  }
  
  nValues <- nrow(df_logFC) * ncol(df_logFC)
  nInf <- sum(df_logFC <= -th_logFC)
  nSup <- sum(df_logFC >=  th_logFC)
  nInside <- sum(abs(df_logFC) < th_logFC)
  
  
  p <- plot_ly()
  
  maxY.inf <- 0
  maxY.inside <- 0
  maxY.sup <- 0
  minX <- Inf
  maxX <- -Inf
  
  
  for (i in seq_len(ncol(df_logFC))) {
    
    tmp <- density(df_logFC[, i], na.rm = TRUE)
    
    minX <- min(minX, tmp$x)
    maxX <- max(maxX, tmp$x)
    
    maxY.inf <- max(maxY.inf, max(tmp$y[tmp$x <= -th_logFC], 0))
    maxY.inside <- max(maxY.inside, max(tmp$y[tmp$x > -th_logFC & tmp$x < th_logFC], 0))
    maxY.sup <- max(maxY.sup, max(tmp$y[tmp$x >=  th_logFC], 0))
    
    p <- p |> add_lines(
      x = tmp$x,
      y = tmp$y,
      name = colnames(df_logFC)[i],
      line = list(color = pal[i]),
      hovertemplate = paste0("<b>", colnames(df_logFC)[i], "</b><br>",
                             "y: %{y:.2f}<extra></extra>"),
      showlegend = TRUE
    )
  }
  
  p <- p |> plotly::layout(
    title = "log(FC) repartition",
    margin = list(t = 60, b = 60),
    xaxis = list(title = "log(FC)"),
    yaxis = list(title = "Density"),
    legend = list(
      orientation = "h", 
      x = 0, 
      y = -0.15, 
      xanchor = "left",
      yanchor = "top"
    ),
    shapes = list(
      list(
        type = "rect",
        x0 = -th_logFC,
        x1 = th_logFC,
        y0 = 0,
        y1 = 1,
        xref = "x",
        yref = "paper",
        fillcolor = "lightgrey",
        opacity = 0.5,
        line = list(width = 0)
      )
    )
  )
  
  if (th_logFC > 0) {
    p <- p |> add_annotations(
      x = 10,
      y = maxY.inside-0.1,
      text = sprintf("n Filtered out = %d<br>(%.2f%%)", nInside, 100*nInside/nValues),
      showarrow = FALSE,
      arrowhead = 2,
      ax = 40,
      ay = -40,
      font = list(size = 18)
    )
  }
  
  if (th_logFC >= minX) {
    p <- p |> add_annotations(
      x = mean(c(minX, -th_logFC)),
      y = maxY.inf+0.1,
      text = sprintf("nInf = %d<br>(%.2f%%)", nInf, 100*nInf/nValues),
      showarrow = FALSE,
      font = list(color = "blue",
                  size = 18)
    )
  }
  
  if (th_logFC <= maxX) {
    p <- p |> add_annotations(
      x = mean(c(maxX, th_logFC)),
      y = maxY.sup+0.1,
      text = sprintf("nSup = %d<br>(%.2f%%)", nSup, 100*nSup/nValues),
      showarrow = FALSE,
      font = list(color = "blue",
                  size = 18)
    )
  }
  
  return(p)
}

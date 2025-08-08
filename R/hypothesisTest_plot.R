
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
#' @param pal xxx
#'
#' @return A highcharts density plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept
#' # Simulate missing value imputation
#' assay(obj[[1]])[which(is.na(assay(obj[[1]])))] <- 0
#' assay(obj[[2]])[which(is.na(assay(obj[[2]])))] <- 0
#' 
#' qData <- as.matrix(assay(obj[[2]]))
#' sTab <- SummarizedExperiment::colData(obj)
#' limma <- limmaCompleteTest(qData, sTab)
#' pal <- ExtendPalette(2, "Dark2")
#' hc_logFC_DensityPlot(limma$logFC, th_logFC = 1, pal = pal)
#'
#' @export
#' 
#' @import highcharter
#' 
#'
hc_logFC_DensityPlot <- function(
    df_logFC,
    th_logFC = 0,
    pal = NULL) {
  
  pkgs.require(c("stats", "RColorBrewer", "grDevices"))
  
  if (th_logFC < 0) {
    warning("The parameter 'th_logFC' must be positive or equal to zero.")
    return(NULL)
  }
  
  
  hc <- highcharter::highchart() %>%
    hc_title(text = "log(FC) repartition") %>%
    my_hc_chart(chartType = "spline", zoomType = "x") %>%
    hc_legend(enabled = TRUE) %>%
    hc_xAxis(
      title = list(text = "log(FC)"),
      plotBands = list(
        list(
          from = -th_logFC, 
          to = th_logFC, 
          color = "lightgrey")
      ),
      plotLines = list(
        list(
          color = "grey", 
          width = 2, 
          value = 0, 
          zIndex = 5
        )
      )
    ) %>%
    hc_yAxis(title = list(text = "Density")) %>%
    hc_tooltip(
      headerFormat = "",
      pointFormat = "<b> {series.name} </b>: {point.y} ",
      valueDecimals = 2
    ) %>%
    my_hc_ExportMenu(filename = "densityplot") %>%
    hc_plotOptions(
      series = list(
        animation = list(duration = 100),
        connectNulls = TRUE,
        marker = list(enabled = FALSE)
      )
    )
  
  if (is.null(df_logFC) || ncol(df_logFC) == 0) {
    return(hc)
  }
  
  
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- ExtendPalette(ncol(df_logFC), "Paired")
  } else {
    if (length(pal) != ncol(df_logFC)) {
      warning("The color palette has not the same dimension as the 
                number of samples")
      myColors <- ExtendPalette(ncol(df_logFC), "Paired")
    }
  }
  
  nValues <- nrow(df_logFC) * ncol(df_logFC)
  nInf <- length(which(df_logFC <= -th_logFC))
  nSup <- length(which(df_logFC >= th_logFC))
  nInside <- length(which(abs(df_logFC) < th_logFC))
  hc <- hc %>%
    hc_colors(myColors)
  
  maxY.inf <- NULL
  maxY.inside <- NULL
  maxY.sup <- NULL
  minX <- NULL
  maxX <- NULL

  
  for (i in seq(ncol(df_logFC))) {
    tmp <- stats::density(df_logFC[, i], na.rm = FALSE)
    ind <- tmp$y[which(tmp$x <= -th_logFC)]
    maxY.inf <- max(maxY.inf, ifelse(length(ind) == 0, 0, ind))
    .ind1 <- which(tmp$x > -th_logFC)
    .ind2 <- which(tmp$x < th_logFC)
    maxY.inside <- max(maxY.inf, tmp$y[intersect(.ind1, .ind2)])
    ind <- tmp$y[which(tmp$x > th_logFC)]
    maxY.sup <- max(
      maxY.sup, 
      ifelse(length(ind) == 0, tmp$y[length(tmp$y)], ind)
    )
    minX <- min(minX, tmp$x)
    maxX <- max(maxX, tmp$x)
    
    
    hc <- hc_add_series(hc,
      data.frame(x = tmp$x, y = tmp$y),
      name = colnames(df_logFC)[i]
    )
  }
  
  ## add annotations
  if (th_logFC > 0) {
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape = "connector",
        backgroundColor = "lightgrey",
        # verticalAlign = 'bottom',
        align = "left",
        # distance=0,
        style = list(
          fontSize = "1.5em",
          textOutline = "1px white"
        ),
        borderWidth = 0,
        x = 20
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = 0,
            y = maxY.inside
          ),
          text = paste0("n Filtered out = ", 
            nInside, "<br>(", 
            round(100 * nInside / nValues, digits = 2), "%)")
        )
      )
    )
  }
  if (th_logFC >= minX) {
    hc <- hc %>%
      hc_add_annotation(
        labelOptions = list(
          shape = "connector",
          backgroundColor = "rgba(255,255,255,0.5)",
          verticalAlign = "top",
          borderWidth = 0,
          crop = TRUE,
          style = list(
            color = "blue",
            fontSize = "1.5em",
            textOutline = "1px white"
          ),
          y = -10
        ),
        labels = list(
          list(
            point = list(
              xAxis = 0,
              yAxis = 0,
              x = mean(c(minX, -th_logFC)),
              y = maxY.inf
            ),
            text = paste0("nInf = ", nInf, "<br>(", 
              round(100 * nInf / nValues, digits = 2), ")%")
          )
        )
      )
  }
  
  if (th_logFC <= maxX) {
    hc <- hc %>% hc_add_annotation(
      labelOptions = list(
        shape = "connector",
        backgroundColor = "blue",
        verticalAlign = "top",
        borderWidth = 0,
        style = list(
          color = "blue",
          fontSize = "1.5em",
          textOutline = "1px white"
        ),
        y = -5
      ),
      labels = list(
        list(
          point = list(
            xAxis = 0,
            yAxis = 0,
            x = mean(c(maxX, th_logFC)),
            y = maxY.sup
          ),
          text = paste0("nSup = ", nSup, "<br>(", 
            round(100 * nSup / nValues, digits = 2), ")%")
        )
      )
    )
  }
  
  
  
  return(hc)
}

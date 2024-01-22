#' 
#' #'
#' #' @title Volcanoplot of the differential analysis
#' #'
#' #' @description Plots an interactive volcanoplot after the differential 
#' #' analysis. Typically, the log of Fold Change is represented on the X-axis 
#' #' and the log10 of the p-value is drawn on the Y-axis. When 
#' #' the \code{threshold_pVal} and the \code{threshold_logFC} are set, two 
#' #' lines are drawn respectively on the y-axis and the X-axis to visually 
#' #' distinguish between differential and non differential data. With the use 
#' #' of the package Highcharter, a customizable tooltip appears when the user 
#' #' put the mouse's pointer over a point of the scatter plot.
#' #'
#' #' @param df A dataframe which contains the following slots :
#' #' x : a vector of the log(fold change) values of the differential analysis,
#' #' y : a vector of the p-value values returned by the differential analysis.
#' #' index : a vector of the rownames of the data.
#' #' This dataframe must has been built with the option stringsAsFactors set
#' #' to FALSE. There may be additional slots which will be used to show
#' #' informations in the tooltip. The name of these slots must begin with the
#' #' prefix "tooltip_". It will be automatically removed in the plot.
#' #'
#' #' @param threshold_pVal A floating number which represents the p-value that
#' #' separates differential and non-differential data.
#' #'
#' #' @param threshold_logFC A floating number which represents the log of the
#' #' Fold Change that separates differential and non-differential data.
#' #'
#' #' @param conditions A list of the names of condition 1 and 2 used for the
#' #' differential analysis.
#' #'
#' #' @param clickFunction A string that contains a JavaScript function used to
#' #' show info from slots in df. The variable this.index refers to the slot
#' #' named index and allows to retrieve the right row to show in the tooltip.
#' #'
#' #' @param palette xxx
#' #'
#' #' @param swap A boolean that indicates if the conditions have been swaped
#' #'
#' #' @return An interactive volcanoplot
#' #'
#' #' @author Samuel Wieczorek, Enora Fremy
#' #'
#' #' @examples
#' #' \donttest{
#' #' library(QFeatures)
#' #' utils::data(Exp1_R25_pept, package = "DaparViz")
#' #' obj <- QFeatures::filterNA(Exp1_R25_pept, i = "original") 
#' #' qData <- SummarizedExperiment::assay(obj[["original"]])
#' #' pData <- as.data.frame(colData(obj))
#' #' data <- limma.complete.test(qData, pData)
#' #' df <- data.frame(x = data[, 1], 
#' #' y = -log10(data[, 2]), index = as.character(rownames(qData)))
#' #' colnames(df) <- c("x", "y", "index")
#' #' tooltipSlot <- c("Sequence", "Score")
#' #' df <- cbind(df, rowData(obj[["original"]])[tooltipSlot])
#' #' colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
#' #' if (ncol(df) > 3) {
#' #'     colnames(df)[4:ncol(df)] <-
#' #'         paste("tooltip_", colnames(df)[4:ncol(df)], sep = "")
#' #' }
#' #' hc_clickFunction <- JS("function(event) {
#' #' Shiny.onInputChange('eventPointClicked', 
#' #' [this.index]+'_'+ [this.series.name]);}")
#' #' cond <- unique(colData(obj)[["Condition"]])
#' #' diffAnaVolcanoplot_rCharts(df, 2.5, 1, cond, hc_clickFunction)
#' #' }
#' #'
#' #' @export
#' #'
#' #' @import highcharter
#' #'
#' #' @rdname differential-analysis
#' #'
#' diffAnaVolcanoplot_rCharts <- function(df,
#'     threshold_pVal = 1e-60,
#'     threshold_logFC = 0,
#'     conditions = NULL,
#'     clickFunction = NULL,
#'     palette = NULL,
#'     swap = FALSE) {
#'     xtitle <- paste("log2 ( mean(", 
#'         conditions[2], ") / mean(", 
#'         conditions[1], ") )", sep = "")
#' 
#'     if (is.null(clickFunction)) {
#'         clickFunction <-
#'             JS("function(event) {Shiny.onInputChange('eventPointClicked', 
#'                 [this.index]+'_'+ [this.series.name]);}")
#'     }
#' 
#'     if (is.null(palette)) {
#'         palette <- list(
#'             In = "orange",
#'             Out = "gray"
#'         )
#'     }
#'     .newCol <- NULL
#'     if(df$y >= threshold_pVal & abs(df$x) >= threshold_logFC) 
#'         .newCol <- "g1"
#'     else
#'         .newCol <- "g2"
#'     df <- cbind(df, g = .newCol)
#' 
#'     i_tooltip <- which(startsWith(colnames(df), "tooltip"))
#'     txt_tooltip <- NULL
#'     for (i in i_tooltip) {
#'         t <- txt_tooltip <- paste(txt_tooltip, "<b>", gsub("tooltip_", "",
#'             colnames(df)[i],
#'             fixed = TRUE
#'         ),
#'         " </b>: {point.", colnames(df)[i], "} <br> ",
#'         sep = ""
#'         )
#'     }
#' 
#'     leftBorder <- data.frame(
#'         x = c(min(df$x), -threshold_logFC, -threshold_logFC),
#'         y = c(threshold_pVal, threshold_pVal, max(df$y))
#'     )
#'     rightBorder <- data.frame(
#'         x = c(max(df$x), threshold_logFC, threshold_logFC),
#'         y = c(threshold_pVal, threshold_pVal, max(df$y))
#'     )
#' 
#' 
#'     title <- NULL
#'     if (isTRUE(swap)) {
#'         title <- paste0(conditions[2], "_vs_", conditions[1])
#'     } else {
#'         title <- paste0(conditions[1], "_vs_", conditions[2])
#'     }
#' 
#' 
#'     h1 <- highchart() %>%
#'         hc_add_series(data = df, 
#'             type = "scatter", 
#'             hcaes(df$x, df$y, group = g)) %>%
#'         hc_colors(c(palette$In, palette$Out)) %>%
#'         dapar_hc_chart(
#'             zoomType = "xy", 
#'             chartType = "scatter") %>%
#'         hc_title(
#'             text = title,
#'             margin = 20, align = "center",
#'             style = list(
#'                 size = 20, 
#'                 color = "black", 
#'                 useHTML = TRUE)
#'         ) %>%
#'         hc_legend(enabled = FALSE) %>%
#'         hc_yAxis(title = list(text = "-log10(pValue)")) %>%
#'         hc_xAxis(
#'             title = list(text = "logFC"),
#'             plotLines = list(
#'                 list(
#'                     color = "grey", 
#'                     width = 1, 
#'                     value = 0, 
#'                     zIndex = 5)
#'                 )
#'         ) %>%
#'         hc_tooltip(headerFormat = "", pointFormat = txt_tooltip) %>%
#'         hc_plotOptions(
#'             line = list(
#'                 marker = list(enabled = FALSE),
#'                 dashStyle = "Dash"
#'             ),
#'             series = list(
#'                 animation = list(duration = 100),
#'                 cursor = "pointer",
#'                 point = list(events = list(
#'                     click = clickFunction
#'                 ))
#'             )
#'         ) %>%
#'         dapar_hc_ExportMenu(filename = "volcanoplot") %>%
#'         hc_add_series(data = leftBorder, type = "line", color = "grey") %>%
#'         hc_add_series(data = rightBorder, type = "line", color = "grey")
#' 
#'     return(h1)
#' }

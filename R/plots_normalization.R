

#' @title Builds a plot from a dataframe. Same as compareNormalizationD but
#' uses the library \code{highcharter}
#' 
#' @description 
#' Plot to compare the quantitative proteomics data before and after
#' normalization using the package \code{highcharter}
#'
#'
#' @param qDataBefore A dataframe that contains quantitative data before
#' normalization.
#'
#' @param qDataAfter A dataframe that contains quantitative data after
#' normalization.
#'
#' @param keyId xxx
#'
#' @param conds A vector of the conditions (one condition
#' per sample).
#'
#' @param pal xxx
#'
#' @param subset.view xxx
#'
#' @param n An integer that is equal to the maximum number of displayed points.
#' This number must be less or equal to the size
#' of the dataset. If it is less than it, it is a random selection
#'
#' @param type scatter or line
#'
#' @return A plot
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package="DaparToolshedData")
#' obj <- Exp1_R25_prot
#' qDataBefore <- SummarizedExperiment::assay(obj[[length(obj)]])
#' conds <- omXplore::get_group(obj)
#' id <- rowData(obj[[length(obj)]])[, omXplore::get_colID(obj[[length(obj)]])]
#' # pal <- ExtendPalette(2)
#' qDataAfter <- GlobalQuantileAlignment(qDataBefore)
#'
#' n <- 1
#' compareNormalizationD_HC(
#' qDataBefore = qDataBefore,
#' qDataAfter = qDataAfter, 
#' keyId = id, 
#' pal = pal, 
#' n = n,
#' subset.view = seq_len(n),
#' conds = conds)
#'
#' @import highcharter
#' @importFrom utils str
#'
#' @export
#'
compareNormalizationD_HC <- function(
    qDataBefore,
    qDataAfter,
    keyId = NULL,
    conds = NULL,
    pal = NULL,
    subset.view = NULL,
    n = 100,
    type = "scatter") {
  
  pkgs.require('RColorBrewer')
  
  if (is.null(conds)) {
    warning("'conds' is null.")
    return(NULL)
  }
  if (n < 0 || n > nrow(qDataBefore)){
    warning("'n' must be a positive integer not null and less than the total number
      of entities. Set to default value: 0.2")
    n <- 100
  }
  
  if (is.null(keyId)) {
    keyId <- seq_len(length(qDataBefore))
  }
  
  if (!is.null(subset.view) && length(subset.view) > 0) {
    keyId <- keyId[subset.view]
    if (nrow(qDataBefore) > 1) {
      if (length(subset.view) == 1) {
        qDataBefore <- t(qDataBefore[subset.view, ])
        qDataAfter <- t(qDataAfter[subset.view, ])
      } else {
        qDataBefore <- qDataBefore[subset.view, ]
        qDataAfter <- qDataAfter[subset.view, ]
      }
      n <- 100
    }
  }
  
  if (!match(type, c("scatter", "line"))) {
    warning("'type' must be equal to 'scatter' or 'line'.")
    return(NULL)
  }
  
  # if (is.null(n)) {
  #     n <- seq_len(nrow(qDataBefore))
  # } else {
  # if (n > nrow(qDataBefore)) {
  #     warning("'n' is higher than the number of rows of datasets. 
  #     Set to number of rows.")
  #     n <- nrow(qDataBefore)
  # }
  # Truncate dataset
  ind <- sample(seq_len(nrow(qDataBefore)), min(n, length(subset.view)))
  keyId <- keyId[ind]
  if (nrow(qDataBefore) > 1) {
    if (length(ind) == 1) {
      qDataBefore <- t(qDataBefore[ind, ])
      qDataAfter <- t(qDataAfter[ind, ])
    } else {
      qDataBefore <- qDataBefore[ind, ]
      qDataAfter <- qDataAfter[ind, ]
    }
  }
  #}
  
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- GetColorsForConditions(conds, 
      ExtendPalette(length(unique(conds))))
  } else {
    if (length(pal) != length(unique(conds))) {
      warning("The color palette has not the same dimension as 
                the number of samples")
      myColors <- GetColorsForConditions(conds, 
        ExtendPalette(length(unique(conds))))
    } else {
      myColors <- GetColorsForConditions(conds, pal)
    }
  }
  
  x <- qDataBefore
  y <- qDataAfter / qDataBefore
  
  ## Colors definition
  legendColor <- unique(myColors)
  txtLegend <- unique(conds)
  
  
  series <- list()
  for (i in seq_len(length(conds))) {
    series[[i]] <- list(
      name = colnames(x)[i],
      data = highcharter::list_parse(data.frame(
        x = x[, i],
        y = y[, i],
        name = keyId
      ))
    )
  }
  
  h1 <- highchart() %>%
    dapar_hc_chart(chartType = type) %>%
    hc_add_series_list(series) %>%
    hc_colors(myColors) %>%
    hc_tooltip(headerFormat = "", pointFormat = "Id: {point.name}") %>%
    dapar_hc_ExportMenu(filename = "compareNormalization")
  h1
}

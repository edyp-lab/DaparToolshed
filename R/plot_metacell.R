

#' @title Bar plot of missing values per lines using highcharter
#' 
#' @description 
#' This method plots a bar plot which represents the distribution of the
#' number of missing values (NA) per lines (ie proteins).
#'
#' @param obj An instance of the class `obj`
#' @param pattern xxx
#' @param detailed 'value' or 'percent'
#' @param indLegend xxx
#' @param showValues A logical that indicates whether numeric values should be
#' drawn above the bars.
#' @return A bar plot
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examplesIf interactive()
#' data(ft_na)
#' grp <- omXplore::get_group(ft_na)
#' metacellPerLinesHisto_HC(ft_na[[1]], group = grp, pattern = "Missing POV")
#' metacellPerLinesHisto_HC(ft_na[[1]])
#' metacellPerLinesHisto_HC(ft_na[[1]], pattern = "Quantified")
#' metacellPerLinesHisto_HC(ft_na[[1]], pattern = "Quant. by direct id")
#' metacellPerLinesHisto_HC(ft_na[[1]], pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellPerLinesHisto_HC(ft_na[[1]], pattern = pattern)
#' 
#' 
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]], pattern = "Missing POV")
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]])
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]], pattern = "Quantified")
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]], pattern = "Quant. by direct id")
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]], pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellPerLinesHistoPerCondition_HC(ft_na[[1]], pattern = pattern)
#' 
#' 
#' metacellHisto_HC(ft_na[[1]], pattern = "Missing POV")
#' metacellHisto_HC(ft_na[[1]])
#' metacellHisto_HC(ft_na[[1]], pattern = "Quantified")
#' metacellHisto_HC(ft_na[[1]], pattern = "Quant. by direct id")
#' metacellHisto_HC(ft_na[[1]], pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellHisto_HC(ft_na[[1]], pattern = pattern)
#' 
#' 
#' 
#' @name metacell-plots
#'
NULL


#' @rdname metacell-plots
#' @export
#' @import highcharter
#'
metacellPerLinesHisto_HC <- function(obj,
  group,
  pattern = NULL,
  detailed = FALSE,
  indLegend = "auto",
  showValues = FALSE) {
  stopifnot(inherits(obj, "SummarizedExperiment"))
  
  if(missing(pattern) || is.null(pattern))
    return(NULL)
  

  if (identical(indLegend, "auto")) {
    indLegend <- seq.int(from = 2, to = length(group))
  }
  
  mask <- match.metacell(get_metacell(obj), 
                         pattern = pattern, 
                         level = omXplore::get_type(obj))
  
  
  if (length(mask) == 0)
    return(NULL)
  
  NbNAPerRow <- rowSums(mask)
  
  nb.col <- dim(assay(obj))[2]
  nb.na <- NbNAPerRow
  temp <- table(NbNAPerRow)
  nb.na2barplot <- rep(0, ncol(assay(obj)))
  
  for (i in seq_len(length(temp))) {
    nb.na2barplot[as.integer(names(temp)[i])] <- temp[i]
  }
  
  
  df <- data.frame(
    y = nb.na2barplot,
    y_percent = round(100 * nb.na2barplot / dim(assay(obj))[1], digits = 2)
  )
  
  myColors <- rep("lightgrey", nrow(df))
  
  h1 <- highchart() %>%
    hc_title(text = paste0("Nb of lines with (", paste0(pattern, collapse=', '), ") tags")) %>%
    hc_add_series(data = df, type = "column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions(
      column = list(stacking = "normal"),
      animation = list(duration = 100)
    ) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = row.names(df), 
             title = list(
               text = paste0("Nb of (", paste0(pattern, collapse=', '), ") tags in a line")
             )
    ) %>%
    my_hc_ExportMenu(filename = "missingValuesPlot1") %>%
    hc_tooltip(
      enabled = TRUE,
      headerFormat = "",
      pointFormat = paste0("{point.y} lines<br>
                ({point.y_percent}% of all lines)")
    )
  
  return(h1)
}







#' @rdname metacell-plots
#' @export
#'
metacellPerLinesHistoPerCondition_HC <- function(obj,
  group,
  pattern = NULL,
  indLegend = "auto",
  showValues = FALSE,
  pal = NULL) {
  stopifnot(inherits(obj, "SummarizedExperiment"))
  
  if(missing(pattern) || is.null(pattern))
    return(NULL)
  
  u_conds <- unique(group)
  nbConds <- length(u_conds)
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- GetColorsForConditions(group, ExtendPalette(nbConds))
  } else {
    if (length(pal) != nbConds) {
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- GetColorsForConditions(group, ExtendPalette(nbConds))
    } else {
      myColors <- pal
    }
  }
  
  if (identical(indLegend, "auto")) {
    indLegend <- seq.int(from = 2, to = length(group))
  }
  
  
  ncolMatrix <- max(unlist(lapply(
    u_conds,  function(x)
      length(which(group == x))
  )))
  
  
  mask <- match.metacell(omXplore::get_metacell(obj), 
                         pattern = pattern, 
                         level = omXplore::get_type(obj))
  
  if (length(mask) == 0)
    return(NULL)
  
  ll.df <- list()
  for (i in u_conds)
    {
    df <- as.data.frame(matrix(rep(0, 2 * (1 + nbConds)),
                               nrow = 1 + nbConds,
                               dimnames = list(
                                 seq(seq.int(from=0, to=(nbConds))),
                                 c("y", "y_percent")
                               )))
    rownames(df) <- seq.int(from = 0, to = (nrow(df) - 1))
    ll.df[[i]] <- df
    nSample <- length(which(group == i))
    t <- NULL
    if (nSample == 1) {
      t <- table(as.integer(mask[, which(group == i)]))
    } else {
      t <- table(rowSums(mask[, which(group == i)]))
    }
    #browser()
    df[as.integer(names(t)) + 1, "y"] <- t
    df[as.integer(names(t)) + 1, "y_percent"] <- round(100 * t / nrow(assay(obj)), digits = 2)
    
    ll.df[[i]] <- df
    }
  
  h1 <- highchart() %>%
    hc_title(text = paste0("Nb of lines containing (", 
                           paste0(pattern, collapse=', '), ") tags (condition-wise)")) %>%
    my_hc_chart(chartType = "column") %>%
    hc_plotOptions(
      column = list(stacking = ""),
      dataLabels = list(enabled = FALSE),
      animation = list(duration = 100)
    ) %>%
    hc_colors(unique(myColors)) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = seq.int(from = 0, to = ncolMatrix), 
             title = list(text = paste0("Nb of (", paste0(pattern, collapse=', '), 
                                        ") tags in each line (condition-wise)"))) %>%
    my_hc_ExportMenu(filename = "missingValuesPlot_2") %>%
    hc_tooltip(
      headerFormat = "",
      pointFormat = "{point.y} lines<br>({point.y_percent}% of all lines)"
    )
  
  for (i in seq_len(nbConds)) {
    h1 <- h1 %>% hc_add_series(data = ll.df[[u_conds[i]]])
  }
  
  return(h1)
}






#' @rdname metacell-plots
#' @import highcharter
#'
#' @examples
#' data(Exp1_R25_pept, package = 'DaparToolshedData')
#' vData <- convert2viz(Exp1_R25_pept)
#' pattern <- "Missing POV"
#' pal <- ExtendPalette(2, "Dark2")
#' metacellHisto_HC(vData@ll.obj[[1]], pattern, showValues = TRUE, pal = pal)
#'
#' @export
#'
metacellHisto_HC <- function(obj,
  group = NULL,
  pattern = NULL,
  indLegend = "auto",
  showValues = FALSE,
  pal = NULL) {
  stopifnot(inherits(obj, "SummarizedExperiment"))
  
  if(missing(pattern) || is.null(pattern) || missing(group) || is.null(group))
    return(NULL)
  
  
  u_conds <- unique(group)
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- GetColorsForConditions(group, ExtendPalette(length(u_conds)))
  } else {
    if (length(pal) != length(u_conds)) {
      warning("The color palette has not the same dimension as the number of samples")
      myColors <- GetColorsForConditions(group, ExtendPalette(length(u_conds)))
    } else {
      myColors <- GetColorsForConditions(group, pal)
    }
  }
  
  if (identical(indLegend, "auto")) {
    indLegend <- seq.int(from=2, to = length(group))
  }
  
  
  
  mask <- match.metacell(omXplore::get_metacell(obj), 
                         pattern = pattern, 
                         level = omXplore::get_type(obj))
  if (length(mask) == 0)
    return(NULL)
  
  NbNAPerCol <- colSums(mask)
  
  df <- data.frame(
    y = NbNAPerCol,
    y_percent = round(100 * NbNAPerCol / nrow(mask), digits = 2)
  )

  
  h1 <- highchart() %>%
    my_hc_chart(chartType = "column") %>%
    hc_title(text = paste0("Nb of (", paste0(pattern, collapse=', '), ") tags by replicate")) %>%
    hc_add_series(df, type = "column", colorByPoint = TRUE) %>%
    hc_colors(myColors) %>%
    hc_plotOptions(
      column = list(stacking = "normal"),
      animation = list(duration = 100)
    ) %>%
    hc_legend(enabled = FALSE) %>%
    hc_xAxis(categories = group, title = list(text = "Replicates")) %>%
    my_hc_ExportMenu(filename = "missingValuesPlot_3") %>%
    hc_tooltip(
      headerFormat = "",
      pointFormat = "{point.y} lines<br>({point.y_percent}% of all lines)"
    )
  
  return(h1)
}





#' @title Heatmap of missing values from an object of class `obj`
#' @description 
#' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class `obj` and
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to vizualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' @param obj An object of class `SummarizedExperiment`.
#' @param pattern xxx
#'
#' @return A heatmap
#' @author Samuel Wieczorek, Alexia Dorffer
#'
#' @rdname metacell-plots
#' @export
#'
#'
wrapper.mvImage <- function(obj, 
  group,
  pattern = "Missing MEC") {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  
  indices <- which(apply(match.metacell(
    omXplore::get_metacell(obj), 
    pattern, 
    level = omXplore::get_type(obj)),
    1, sum) > 0)
  
  if (length(indices) == 0) {
    warning("The dataset contains no Missing value on Entire Condition. 
            So this plot is not available.")
    return(NULL)
  } else if (length(indices) == 1) {
    warning("The dataset contains only one Missing value on Entire 
        Condition. Currently, Prostar does not handle such dataset to build 
        the plot. As it has no side-effects on the results, you can continue 
        your imputation.")
    return(NULL)
  }
  
  mvImage(assay(obj)[indices, ], group)
}





#' @title Heatmap of missing values
#' @description 
#' #' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class `MsnSet` and
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to vizualize easily the different
#' number of missing values. A white square is plotted for missing values.
#' 
#' @return A heatmap
#' @author Samuel Wieczorek, Thomas Burger
#' @examples
#' data(Exp1_R25_pept, package = 'DaparToolshedData')
#' vData <- convert2viz(Exp1_R25_pept)
#' mvImage(vData@ll.obj[[1]])
#'
#' @export
#'
#' @rdname metacell-plots
#'
mvImage <- function(obj, group) {
  
  pkgs.require(c('grDevices', 'stats'))
  
  ### build indices of conditions
  indCond <- list()
  ConditionNames <- unique(group)
  for (i in ConditionNames) {
    indCond <- append(indCond, list(which(i == group)))
  }
  indCond <- stats::setNames(indCond, as.list(c("cond1", "cond2")))
  
  nNA1 <- apply(as.matrix(assay(obj)[, indCond$cond1]), 1, 
                function(x) sum(is.na(x)))
  nNA2 <- apply(as.matrix(assay(obj)[, indCond$cond2]), 1, 
                function(x) sum(is.na(x)))
  o <- order(((nNA1 + 1)^2) / (nNA2 + 1))
  exprso <- assay(obj)[o, ]
  
  for (i in seq_len(nrow(exprso))) {
    k <- order(exprso[i, indCond$cond1])
    exprso[i, rev(indCond$cond1)] <- exprso[i, k]
    .temp <- mean(exprso[i, rev(indCond$cond1)], na.rm = TRUE)
    exprso[i, which(!is.na(exprso[i, indCond$cond1]))] <- .temp
    
    k <- order(exprso[i, indCond$cond2])
    exprso[i, indCond$cond2] <- exprso[i, k + length(indCond$cond1)]
    .temp <- mean(exprso[i, indCond$cond2], na.rm = TRUE)
    exprso[i, length(indCond$cond1) +
             which(!is.na(exprso[i, indCond$cond2]))] <- .temp
  }
  
  
  heatmapForMissingValues(exprso,
                          col = grDevices::colorRampPalette(c("yellow", "red"))(100),
                          key = TRUE,
                          srtCol = 0,
                          labCol = group,
                          ylab = "Peptides / proteins",
                          main = "MEC heatmap"
  )
  
  # heatmap_HC(exprso,col = colfunc(100),labCol=conds)
}



#' @title Distribution of Observed values with respect to intensity values
#' 
#' @description 
#' This method shows density plots which represents the repartition of
#' Partial Observed Values for each replicate in the dataset.
#' The colors correspond to the different conditions (slot Condition in in the
#' dataset of class `MsnSet`).
#' The x-axis represent the mean of intensity for one condition and one
#' entity in the dataset (i. e. a protein)
#' whereas the y-axis count the number of observed values for this entity
#' and the considered condition.
#'
#' @param obj xxx
#' @param pal The different colors for conditions
#' @param pattern xxx
#' @param typeofMV xxx
#' @param title The title of the plot
#'
#' @import highcharter
#'
#' @return Density plots
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package = 'DaparToolshedData')
#' pal <- ExtendPalette(length(unique(get_group(Exp1_R25_pept))), "Dark2")
#' pattern <- "Missing MEC"
#' hc_mvTypePlot2(Exp1_R25_pept[[1]], group = get_group(Exp1_R25_pept), pattern = pattern, pal = pal)
#'
#' @import highcharter
#' @rdname metacell-plots
#'
#' @export
#'
hc_mvTypePlot2 <- function(obj,
  group,
  pal = NULL,
  pattern,
  typeofMV = NULL,
  title = NULL) {
  
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  pkgs.require('stats')
  
  
  qdata <- SummarizedExperiment::assay(obj)
  myColors <- pal
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- ExtendPalette(length(unique(group)))
  } else {
    if (length(pal) != length(unique(group))) {
      warning("The color palette has not the same dimension as the 
                number of samples")
      myColors <- ExtendPalette(length(unique(group)))
    }
  }
  
  mTemp <- nbNA <- nbValues <- matrix(
    rep(0, nrow(assay(obj)) * length(unique(group))),
    nrow = nrow(assay(obj)),
    dimnames = list(NULL, unique(group))
  )
  dataCond <- data.frame()
  ymax <- 0
  series <- list()
  myColors <- NULL
  j <- 1
  
  
  for (iCond in unique(group)) {
    if (length(which(group == iCond)) == 1) {
      mTemp[, iCond] <- qdata[, which(group == iCond)]
      nbNA[, iCond] <- as.integer(
        match.metacell(omXplore::get_metacell(obj)[, which(group == iCond)],
                       pattern = pattern,
                       level = omXplore::get_type(obj))
      )
      
      .op1 <- length(which(group == iCond))
      .op2 <- nbNA[, iCond]
      nbValues[, iCond] <- .op1 - .op2
    } else {
      .qcond <- which(group == iCond)
      mTemp[, iCond] <- apply(assay(obj)[, .qcond], 1, mean, na.rm = TRUE)
      
      nbNA[, iCond] <- rowSums(
        match.metacell(omXplore::get_metacell(obj)[, .qcond],
                       pattern = pattern,
                       level = omXplore::get_type(obj))
      )
      
      nbValues[, iCond] <- length(.qcond) - nbNA[, iCond]
    }
    
    
    for (i in seq_len(length(which(group == iCond)))) {
      data <- mTemp[which(nbValues[, iCond] == i), iCond]
      tmp <- NULL
      if (length(data) >= 2) {
        tmp <- stats::density(mTemp[which(nbValues[, iCond] == i), iCond])
        tmp$y <- tmp$y + i
        if (max(tmp$y) > ymax) {
          ymax <- max(tmp$y)
        }
      }
      series[[j]] <- tmp
      myColors <- c(myColors, pal[which(unique(group) == iCond)])
      j <- j + 1
    }
  }
  
  
  hc <- highchart(type = "chart") %>%
    hc_title(text = title) %>%
    my_hc_chart(chartType = "spline", zoomType = "xy") %>%
    hc_legend(align = "left", verticalAlign = "top", layout = "vertical"
    ) %>%
    hc_xAxis(title = list(text = "Mean of intensities")) %>%
    hc_yAxis(title = list(text = "Number of quantity values per condition"),
             tickInterval = 0.5) %>%
    hc_tooltip(
      headerFormat = "",
      pointFormat = "<b> {series.name} </b>: {point.y} ",
      valueDecimals = 2
    ) %>%
    my_hc_ExportMenu(filename = paste0(pattern, "_distribution")) %>%
    hc_plotOptions(
      series = list(
        showInLegend = TRUE,
        animation = list(duration = 100),
        connectNulls = TRUE,
        marker = list(enabled = FALSE)
      )
    )
  
  for (i in seq_len(length(series))) {
    hc <- hc_add_series(hc,
                        data = list_parse(data.frame(cbind(
                          x = series[[i]]$x,
                          y = series[[i]]$y
                        ))),
                        showInLegend = FALSE,
                        color = myColors[i],
                        name = group[i]
    )
  }
  
  # add three empty series for the legend entries. Change color and marker 
  # symbol
  for (c in seq_len(length(unique(group)))) {
    hc <- hc_add_series(hc,
                        data = data.frame(),
                        name = unique(group)[c],
                        color = pal[c],
                        marker = list(symbol = "circle"),
                        type = "line"
    )
  }
  
  hc
  return(hc)
}

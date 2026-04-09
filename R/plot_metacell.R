

#' @title Bar plot of missing values per lines using plotly
#' 
#' @description 
#' This method plots a bar plot which represents the distribution of the
#' number of missing values (NA) per lines (ie proteins).
#'
#' @param obj An instance of the class `QFeatures`
#' @param group A vector 
#' @param pattern A `character()` indicating the tag pattern of interest. 
#' @param detailed 'value' or 'percent'
#' @param indLegend A `vector()` of `integers`
#' @param showValues A logical that indicates whether numeric values should be
#' drawn above the bars.
#' @return A bar plot
#' 
#' @author Florence Combes, Samuel Wieczorek
#' 
#' @examplesIf interactive()
#' data(subR25prot)
#' grp <- design_qf(subR25prot)$Condition
#' metacellPerLinesHisto_HC(subR25prot[[1]], group = grp, pattern = "Missing POV")
#' metacellPerLinesHisto_HC(subR25prot[[1]])
#' metacellPerLinesHisto_HC(subR25prot[[1]], group = grp, pattern = "Quantified")
#' metacellPerLinesHisto_HC(subR25prot[[1]], group = grp, pattern = "Quant. by direct id")
#' metacellPerLinesHisto_HC(subR25prot[[1]], group = grp, pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellPerLinesHisto_HC(subR25prot[[1]], group = grp, pattern = pattern)
#' 
#' 
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]], group = grp, pattern = "Missing POV")
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]])
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]], group = grp, pattern = "Quantified")
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]], group = grp, pattern = "Quant. by direct id")
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]], group = grp, pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellPerLinesHistoPerCondition_HC(subR25prot[[1]], group = grp, pattern = pattern)
#' 
#' 
#' metacellHisto_HC(subR25prot[[1]], group = grp, pattern = "Missing POV")
#' metacellHisto_HC(subR25prot[[1]])
#' metacellHisto_HC(subR25prot[[1]], group = grp, pattern = "Quantified")
#' metacellHisto_HC(subR25prot[[1]], group = grp, pattern = "Quant. by direct id")
#' metacellHisto_HC(subR25prot[[1]], group = grp, pattern = "Quant. by recovery")
#' pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
#' metacellHisto_HC(subR25prot[[1]], group = grp, pattern = pattern)
#' 
#' 
#' 
#' @name metacell-plots
#'
NULL


#' @rdname metacell-plots
#' @export
#' @import plotly
#' @import QFeatures
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
  
  mask <- matchMetacell(qMetacell(obj), 
                         pattern = pattern, 
                         level = typeDataset(obj))
  
  
  if (length(mask) == 0)
    return(NULL)
  
  NbNAPerRow <- rowSums(mask)
  
  nb.col <- dim(SummarizedExperiment::assay(obj))[2]
  nb.na <- NbNAPerRow
  temp <- table(NbNAPerRow)
  nb.na2barplot <- rep(0, ncol(SummarizedExperiment::assay(obj)))
  
  for (i in seq_len(length(temp))) {
    nb.na2barplot[as.integer(names(temp)[i])] <- temp[i]
  }
  
  
  df <- data.frame(
    y = nb.na2barplot,
    y_percent = round(100 * nb.na2barplot / dim(SummarizedExperiment::assay(obj))[1], digits = 2)
  )
  
  myColors <- rep("lightgrey", nrow(df))
  
  p <- plot_ly(
    x = row.names(df),
    y = df[[1]],
    type = 'bar',
    marker = list(color = myColors),
    text = paste0(df[[1]], " lines<br>(",
                  round(df[[1]] / sum(df[[1]]) * 100, 1), "% of all lines)"),
    textposition = 'none',
    hoverinfo = 'text'
  ) %>%
    plotly::layout(
      margin = list(t = 60, b = 60),
      title = paste0("Nb of lines with (", paste0(pattern, collapse=', '), ") tags"),
      xaxis = list(title = paste0("Nb of (", paste0(pattern, collapse=', '), ") tags in a line")),
      yaxis = list(title = "Count"),
      showlegend = FALSE
    )
  
  return(p)
}







#' @rdname metacell-plots
#' @export
#' @import QFeatures
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
  } else if (length(pal) != nbConds) {
    warning("The color palette has not the same dimension as the number of samples")
    myColors <- GetColorsForConditions(group, ExtendPalette(nbConds))
  } else {
    myColors <- pal
  }
  
  ncolMatrix <- max(sapply(u_conds, function(x) sum(group == x)))
  
  mask <- matchMetacell(qMetacell(obj), pattern = pattern, level = typeDataset(obj))
  if (length(mask) == 0) return(NULL)
  
  df_all <- data.frame()
  for (i in u_conds) {
    nSample <- sum(group == i)
    t <- if (nSample == 1) table(as.integer(mask[, which(group == i)])) else table(rowSums(mask[, which(group == i)]))
    
    x_vals <- 0:ncolMatrix
    y_vals <- rep(0, length(x_vals))
    y_vals[as.integer(names(t)) + 1] <- t
    y_percent <- round(100 * y_vals / nrow(SummarizedExperiment::assay(obj)), 2)
    
    df_cond <- data.frame(
      x = x_vals,
      y = y_vals,
      y_percent = y_percent,
      condition = i
    )
    
    df_all <- rbind(df_all, df_cond)
  }
  df_all$condition <- factor(df_all$condition, levels = u_conds)
  
  p <- plot_ly(df_all, 
               x = ~x, 
               y = ~y, 
               color = ~condition, 
               colors = myColors,
               type = 'bar', 
               text = if (showValues) ~y else ~paste0(condition, " : ", y, " lines (", y_percent, "%)"),
               textposition = 'none', 
               hoverinfo = 'text') |>
    plotly::layout(
      barmode = 'group', 
      title = paste0("Nb of lines containing (", paste0(pattern, collapse=', '), ") tags (condition-wise)"),
      xaxis = list(title = paste0("Nb of (", paste0(pattern, collapse=', '), ") tags in each line")),
      yaxis = list(title = ""),
      margin = list(t = 60, b = 60),
      showlegend = FALSE
    )
  
  return(p)
}






#' @rdname metacell-plots
#' @import plotly
#'
#' @examples
#' data(subR25pept)
#' pattern <- "Missing POV"
#' pal <- ExtendPalette(2, "Dark2")
#' metacellHisto_HC(subR25pept[[1]], pattern, showValues = TRUE, pal = pal)
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
  } else if (length(pal) != length(u_conds)) {
    warning("The color palette has not the same dimension as the number of samples")
    myColors <- GetColorsForConditions(group, ExtendPalette(length(u_conds)))
  } else {
    myColors <- GetColorsForConditions(group, pal)
  }
  
  if (identical(indLegend, "auto")) {
    indLegend <- seq.int(from=2, to = length(group))
  }
  
  
  
  mask <- matchMetacell(qMetacell(obj), 
                        pattern = pattern, 
                        level = typeDataset(obj))
  if (length(mask) == 0)
    return(NULL)
  
  NbNAPerCol <- colSums(mask)
  
  df <- data.frame(
    y = NbNAPerCol,
    y_percent = round(100 * NbNAPerCol / nrow(mask), digits = 2),
    group = group
  )
  
  df$sample <- rownames(df)
  df$group <- factor(df$group, levels = u_conds)
  
  p <- plot_ly(
    data = df,
    x = ~sample,
    y = ~y,
    type = "bar",
    color = ~group,
    colors = myColors,
    text = if (showValues) ~y else ~paste0(group, " : ", y, " lines (", y_percent, "%)"),
    textposition = 'none', 
    hoverinfo = 'text'
  ) %>%
    layout(
      title = paste0("Nb of (", paste0(pattern, collapse=', '), ") tags by replicate"),
      xaxis = list(
        title = "Replicates",
        tickmode = "array",
        tickvals = df$sample, 
        ticktext = df$group 
      ),
      yaxis = list(title = ""),
      bargap = 0.2,
      margin = list(t = 60, b = 60),
      showlegend = FALSE
    )
  
  return(p)
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
#' @param pattern A `character()` indicating the tag pattern of interest. 
#'
#' @return A heatmap
#' @author Samuel Wieczorek, Alexia Dorffer
#'
#' @rdname metacell-plots
#' @export
#' @import QFeatures
#'
wrapperMVImage <- function(obj, 
  group = NULL,
  pattern = "Missing MEC") {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  
  if (is.null(group))
    return(NULL)
  
  indices <- which(apply(matchMetacell(
    qMetacell(obj), 
    pattern, 
    level = typeDataset(obj)),
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
  
  mvImage(SummarizedExperiment::assay(obj)[indices, ], group)
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
#' data(subR25pept)
#' mvImage(subR25pept[[1]], design_qf(subR25pept)$Condition)
#'
#' @export
#' 
#' @import QFeatures
#' @importFrom stats setNames
#' @rdname metacell-plots
#' 
#'
mvImage <- function(obj, group) {
  
  pkgsRequire(c('grDevices'))
  
  ### build indices of conditions
  indCond <- list()
  ConditionNames <- unique(group)
  for (i in ConditionNames) {
    indCond <- append(indCond, list(which(i == group)))
  }
  indCond <- stats::setNames(indCond, as.list(c("cond1", "cond2")))
  
  nNA1 <- apply(as.matrix(SummarizedExperiment::assay(obj)[, indCond$cond1]), 1, 
                function(x) sum(is.na(x)))
  nNA2 <- apply(as.matrix(SummarizedExperiment::assay(obj)[, indCond$cond2]), 1, 
                function(x) sum(is.na(x)))
  o <- order(((nNA1 + 1)^2) / (nNA2 + 1))
  exprso <- SummarizedExperiment::assay(obj)[o, ]
  
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
#' @param obj An instance of the class `QFeatures`
#' @param pal The different colors for conditions
#' @param pattern A `character()` indicating the tag pattern of interest. 
#' @param title The title of the plot
#'
#' @import plotly
#'
#' @return Density plots
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' pal <- ExtendPalette(length(unique(design_qf(subR25pept)$Condition)), "Dark2")
#' pattern <- "Missing MEC"
#' hc_mvTypePlot2(subR25pept[[1]], 
#' group = design_qf(subR25pept)$Condition, 
#' pattern = pattern, pal = pal)
#'
#' @import plotly
#' @rdname metacell-plots
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom stats density
#'
hc_mvTypePlot2 <- function(obj,
                           group,
                           pal = NULL,
                           pattern,
                           title = NULL) {
  
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  
  
  qdata <- SummarizedExperiment::assay(obj)
  myColors <- NULL
  if (is.null(pal)) {
    warning("Color palette set to default.")
    myColors <- ExtendPalette(length(unique(group)))
  } else if (length(pal) != length(unique(group))) {
    warning("The color palette has not the same dimension as the 
              number of samples")
    myColors <- ExtendPalette(length(unique(group)))
  } else {
    myColors <- pal
  }
  
  mTemp <- nbNA <- nbValues <- matrix(
    rep(0, nrow(SummarizedExperiment::assay(obj)) * length(unique(group))),
    nrow = nrow(SummarizedExperiment::assay(obj)),
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
        matchMetacell(qMetacell(obj)[, which(group == iCond)],
                      pattern = pattern,
                      level = typeDataset(obj))
      )
      
      .op1 <- length(which(group == iCond))
      .op2 <- nbNA[, iCond]
      nbValues[, iCond] <- .op1 - .op2
    } else {
      .qcond <- which(group == iCond)
      mTemp[, iCond] <- apply(SummarizedExperiment::assay(obj)[, .qcond], 1, mean, na.rm = TRUE)
      
      nbNA[, iCond] <- rowSums(
        matchMetacell(qMetacell(obj)[, .qcond],
                      pattern = pattern,
                      level = typeDataset(obj))
      )
      
      nbValues[, iCond] <- length(.qcond) - nbNA[, iCond]
    }
    
    
    for (i in seq_len(length(which(group == iCond)))) {
      data <- mTemp[which(nbValues[, iCond] == i), iCond]
      tmp <- NULL
      if (length(data) >= 2) {
        tmp <- stats::density(mTemp[which(nbValues[, iCond] == i), iCond], na.rm = TRUE)
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
  
  
  p <- plot_ly() |> 
    plotly::layout(
      title = list(text = title),
      xaxis = list(title = "Mean of intensities"),
      yaxis = list(title = "Number of quantity values per condition",
                   tickmode = "linear",
                   dtick = 0.5),
      legend = list(orientation = "v", x = 0, y = 1),
      hoverlabel = list(namelength = -1)
    )
  
  for (i in seq_len(length(series))) {
    if (!is.null(series[[i]])){
      df <- data.frame(
        x = series[[i]]$x,
        y = series[[i]]$y
      )
      
      p <- p |> plotly::add_trace(
        data = df,
        x = ~x, y = ~y,
        type = "scatter",
        mode = "lines",
        line = list(shape = "spline", color = myColors[i]),
        showlegend = FALSE,
        name = group[i],
        hovertemplate = paste("<b>", group[i], "</b>: %{y:.2f}<extra></extra>")
      )
    }
  }
  
  unique_groups <- unique(group)
  
  for (c in seq_along(unique_groups)) {
    p <- p |> add_trace(
      x = 0, y = 0,
      type = "scatter",
      mode = "lines+markers",
      marker = list(symbol = "circle",
                    size = 8, 
                    color = pal[c]),
      line = list(color = pal[c]),
      name = unique_groups[c],
      showlegend = TRUE,
      visible = "legendonly",
      hoverinfo = "none"
    )
  }
  
  p <- p |> plotly::layout(legend = list(itemclick = FALSE, 
                                         itemdoubleclick = FALSE,
                                         orientation = "h", 
                                         x = 0, 
                                         y = -0.15, 
                                         xanchor = "left",
                                         yanchor = "top"))
  
  return(p)
}

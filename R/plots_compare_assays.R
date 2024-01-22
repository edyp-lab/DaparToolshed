
#' @title Compare two assays
#'
#' @description
#' This plot compares the quantitative proteomics data between two assays.
#' It can be used for example to compare the effect of the normalization 
#' process.
#'
#' The comparison is made with the division operator.
#'
#' @param obj An instance of the class `DaparViz`.
#' @param i A numeric matrix containing quantitative data after 
#' normalization.
#' @param j A numeric matrix containing quantitative data after 
#' normalization
#' @param info xxx
#' @param pal.name xxx
#' @param subset.view xxx
#' @param n xxx
#' @param type The type of plot. Available values are 'scatter' (default) 
#' or 'line'
#' @param FUN xxx
#'
#' @return A plot
#'
#' @author Samuel Wieczorek, Enora Fremy
#'
#' @examples
#' data(vData_ft)
#' plotCompareAssays(vData_ft, 1, 2, n = 5)
#'
#' @import highcharter
#' @importFrom tibble as_tibble
#' @importFrom utils str
#'
#' @export
#'
#' @rdname compare-assays
#'
plotCompareAssays <- function(obj,
                              i,
                              j,
                              info = NULL,
                              pal.name = NULL,
                              subset.view = NULL,
                              n = 100,
                              type = "scatter",
                              FUN = NULL) {
  
  pkgs.require('highcharter')
    if (missing(obj)) {
        stop("'vList' is missing")
    }
    stopifnot(inherits(obj, "list"))
    
    qdata1 <- obj[[i]]@qdata
    qdata2 <- obj[[j]]@qdata
    conds <- obj[[i]]@conds
    
   
    if (nrow(qdata1) != nrow(qdata2) || 
        ncol(qdata1) != ncol(qdata2) ) {
        stop("Assays must have the same dimensions.")
    }

    if (ncol(qdata1) != length(conds)) {
        stop("The length of 'conds' must be equal to the number of columns of
    'qdata1' and 'qdata2'")
    }


    if (is.null(info)) {
        info <- rep(NA, nrow(qdata1))
    } else if (length(info) != nrow(qdata1)) {
        stop("'info' must have the same length as the number of rows of
         'qdata1' and 'qdata2'")
    }


    stopifnot(type %in% c("scatter", "line"))

    # Reduce the assays to the entities given by 'subset.view'
    if (!is.null(subset.view) && length(subset.view) > 0) {
        if (nrow(qdata1) > 1) {
            info <- info[subset.view]

            if (length(subset.view) == 1) {
                qdata1 <- t(qdata1[subset.view, ])
                qdata2 <- t(qdata2[subset.view, ])
            } else {
                qdata1 <- qdata1[subset.view, ]
                qdata2 <- qdata2[subset.view, ]
            }
        }
    }


    if (is.null(n) || n > nrow(qdata1)) {
        warning("'n' is NULL or higher than the number of rows of datasets.
      Set to number of rows.")
        n <- nrow(qdata1)
    }



    # Reduce the assays to 'n' entities so as to make a lighter plot
    if (nrow(qdata1) > 1) {
        ind <- sample(seq_len(nrow(qdata1)), n)
        info <- info[ind]
        if (length(ind) == 1) {
            qdata1 <- t(qdata1[ind, ])
            qdata2 <- t(qdata2[ind, ])
        } else {
            qdata1 <- qdata1[ind, ]
            qdata2 <- qdata2[ind, ]
        }
    }

    myColors <- ExtendPalette(length(unique(conds)), pal.name)

    # Compare by using the division between assays
    x <- qdata1
    y <- qdata2 / qdata1

    series <- list()
    for (i in seq_len(length(conds))) {
        series[[i]] <- list(
            name = colnames(x)[i],
            data = highcharter::list_parse(
                data.frame(x = x[, i],
                           y = y[, i],
                           name = info)
                )
        )
    }

    h1 <- highcharter::highchart() %>%
        customChart(chartType = type) %>%
        highcharter::hc_add_series_list(series) %>%
      highcharter::hc_colors(myColors) %>%
        customExportMenu(fname = "compareAssays")

    if (!all.equal(info, rep(NA, length(info)))) {
        h1 <- h1 %>%
          highcharter::hc_tooltip(headerFormat = "", 
                                  pointFormat = "Id: {point.name}")
    } else {
        h1 <- h1 %>% highcharter::hc_tooltip(enabled = FALSE)
    }


    h1
}

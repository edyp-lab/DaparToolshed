

#' @title Standardize names
#'
#' @description
#'
#' Replace ".", ' ', '-' in `character()` by '_' to be compliant
#' with functions of `Shinyjs`, `Shiny`
#'
#' @param x A `character()` to be processed
#'
#' @return A `character()` of the same length as 'x' with modified
#' names.
#'
#' @author Samuel Wieczorek
#'
#' @export
#'
#' @examples
#'
#' ReplaceSpecialChars(c("foo.1", "foo-2", "foo 3"))
#'
ReplaceSpecialChars <- function(x) {
    if (is.null(x)) {
        return(x)
    }

    for (char in c(".", " ")) {
        x <- gsub(char, "_", x, fixed = TRUE)
    }

    x
}


#' @title Version number of Prostar suite
#'
#' @description
#'
#' This function gives the version number of the packages of the
#' Prostar suite and which propose data processing. This information
#' can be useful if the user wants to publish its works or to
#' rerun a data processing pipeline tin a given set of conditions.
#' The packages which are concerned are `Prostar`, `DaparToolshed`
#' and `DaparToolshedData`
#'
#' @return A `list(3)`
#'
#' @rdname ProstarVersions
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' ProstarVersions()
#'
#' @export
#'
ProstarVersions <- function() {
    v.Prostar <- v.DaparToolshed <- v.DaparToolshedData <- NA

    # tryCatch(
    #     {
    #         find.package("Prostar.2.0")
    #         v.Prostar <- Biobase::package.version("Prostar.2.0")
    #     },
    #     error = function(e) v.Prostar <- NA
    # )

    tryCatch(
        {
            find.package("DaparToolshed")
          v.DaparToolshed <- Biobase::package.version("DaparToolshed")
        },
        error = function(e) v.DaparToolshed <- NA
    )

    tryCatch(
        {
            find.package("DaparToolshedData")
          v.DaparToolshed <- Biobase::package.version("DaparToolshedData")
        },
        error = function(e) v.DaparToolshedData <- NA
    )


    list(
        Prostar = v.Prostar,
        DaparToolshed = v.DaparToolshed,
        DaparToolshedData = v.DaparToolshedData
    )
}










#' @title Number of empty lines in the data
#'
#' @description
#'
#' This function counts the number of empty lines (all elements
#' are equal to NA).
#'
#' @param df A `data.frame`.
#'
#' @return A `integer(1)`
#'
#' @export
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' library(QFeatures)
#' data(subR25prot)
#' nEmptyLines(SummarizedExperiment::assay(subR25prot, 1))
#'
nEmptyLines <- function(df) {
    sum(apply(is.na(as.matrix(df)), 1, all))
}






#' @title Similar to the function `is.na()` but focused on the equality 
#' with the paramter 'type'.
#'
#' @param data A data.frame
#'
#' @param type The value to search in the dataframe
#'
#' @return A boolean dataframe
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' library(QFeatures)
#' data(subR25prot)
#' obj <- subR25prot[[1]]
#' data <- qMetacell(obj)
#' is.OfType(as.data.frame(data), "MEC")
#'
#' @export
#'
is.OfType <- function(data, type) {
    stopifnot(inherits(data, "data.frame"))
    return(type == data)
}



#' @title xxx
#' @description xxx
#' @param set1 xxx
#' @param set2 xxx
#' 
#' @return xxx
#' 
#' @examples
#' is.subset('a', letters)
#' is.subset(c('a', 'c', 't'), letters)
#' is.subset(c('a', 3, 't'), letters)
#' is.subset(3, letters)
#' 
#' @export
#' 
is.subset <- function(set1, set2)
  length(intersect(set1, set2)) > 0 && length(set1) == length(intersect(set1, set2))






#' @title Returns the possible number of values in lines in the data
#'
#' @param object An object of class `QFeatures`
#'
#' @param conds xxxx
#'
#' @param type WholeMatrix, AllCond or AtLeastOneCond
#'
#' @return An integer
#'
#' @author Samuel Wieczorek, Enora Fremy
#'
#' @examples
#' data(subR25prot)
#' res <- getListNbValuesInLines(subR25prot[[1]])
#'
#' @export
#'
#' @importFrom S4Vectors sort
#' @import QFeatures
#'
getListNbValuesInLines <- function(object, conds, type = "WholeMatrix") {
    if (is.null(object)) {
        return(NULL)
    }
    stopifnot(inherits(object, "SummarizedExperiment"))
    if (!(type %in% c("None", "WholeLine", "WholeMatrix", 
        "AtLeastOneCond", "AllCond"))) {
        stop("'type' is not one of: 'None', 
            'WholeMatrix', 'AtLeastOneCond', 'AllCond'")
    }

    data <- as.data.frame(SummarizedExperiment::assay(object))
    ll <- switch(type,
      WholeLine = NULL,
      WholeMatrix = seq(0, ncol(data)),
      AllCond = seq(0, min(unlist(lapply(unique(conds), 
        function(x) length(which(conds == x)))))),
      AtLeastOneCond = seq(0, min(unlist(lapply(unique(conds), 
        function(x) length(which(conds == x))))))
    )

    return(sort(ll))
}







#' @title Retrieve the indices of non-zero elements in sparse matrices
#'
#' @description
#'
#' This function retrieves the indices of non-zero elements in sparse matrices
#' of class dgCMatrix from package Matrix. This function is largely inspired 
#' from the package `RINGO`.
#'
#' @param x A sparse matrix of class dgCMatrix
#'
#' @return A two-column matrix
#'
#' @author Samuel Wieczorek
#'
#' @rdname nonzero
#'
#' @examples
#' library(Matrix)
#' mat <- Matrix(c(0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1),
#'     nrow = 5, byrow = TRUE, sparse = TRUE)
#' res <- nonzero(mat)
#'
#' @export
#'
nonzero <- function(x) {
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Please install Matrix: BiocManager::install('Matrix')")
    }

    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0)) {
        return(matrix(0,
            nrow = 0, ncol = 2,
            dimnames = list(character(0), c("row", "col"))
        ))
    }
    res <- cbind(x@i + 1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}



#' @title Convert a list to unnumbered HTML list
#' @description xxx
#' @param ll A `list()` of `character()`
#' @rdname ConvertListToHtml
#' @export
#' 
#' @examples 
#' ConvertListToHtml(list('foo1', 'foo2', 'foo3'))
#' @return HTML
#' 
ConvertListToHtml <- function(ll) {
    if (length(ll) == 0) {
        return("-")
    }

    paste0(
        "<ul>",
        paste0(
            lapply(
                ll,
                function(x) {
                    paste0(
                        "<li>",
                        paste0(x, collapse = " "),
                        "</li>"
                    )
                }
            ),
            collapse = " "
        ),
        "</ul>"
    )
}




#' @title Customised contextual menu of highcharts plots
#'
#' @param hc A highcharter object
#'
#' @param filename The filename under which the plot has to be saved
#'
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' library("highcharter")
#' hc <- highchart()
#' hc_chart(hc, type = "line")
#' hc_add_series(hc, data = c(29, 71, 40))
#' my_hc_ExportMenu(hc, filename = "foo")
#'
#' @export
#' @return A contextual menu for highcharts plots.
#'
my_hc_ExportMenu <- function(hc, filename) {
  hc_exporting(hc,
               enabled = TRUE,
               filename = filename,
               buttons = list(
                 contextButton = list(
                   menuItems = list("downloadPNG", "downloadSVG", "downloadPDF")
                 )
               )
  )
}






#' @title Customised resetZoomButton of highcharts plots
#'
#' @param hc A highcharter object
#'
#' @param chartType The type of the plot
#'
#' @param zoomType The type of the zoom (one of "x", "y", "xy", "None")
#'
#' @return A highchart plot
#'
#' @author Samuel Wieczorek
#' 
#' @import highcharter
#'
#' @examples
#' library("highcharter")
#' hc <- highchart()
#' hc_chart(hc, type = "line")
#' hc_add_series(hc, data = c(29, 71, 40))
#' my_hc_ExportMenu(hc, filename = "foo")
#'
#' @export
#'
my_hc_chart <- function(hc, chartType, zoomType = "None") {
  hc |>
    hc_chart(
      type = chartType,
      zoomType = zoomType,
      showAxes = TRUE,
      resetZoomButton = list(
        position = list(
          align = "left",
          verticalAlign = "top"
        )
      )
    )
}






#' @title Customised resetZoomButton of highcharts plots
#'
#' @author Samuel Wieczorek
#' 
#' @param obj xxx
#' @param i xxx
#'
#' @examples
#' library("highcharter")
#' hc <- highchart()
#' hc_chart(hc, type = "line")
#' hc_add_series(hc, data = c(29, 71, 40))
#' my_hc_ExportMenu(hc, filename = "foo")
#'
#' @export
#' @import SummarizedExperiment
#' 
#' @return An instance of the SummarizedExperiment structure
#'
CleanRowData <- function(obj, i){
  stopifnot(inherits(obj, 'QFeatures'))
  
    ind <- match(c('qMetacell','adjacencyMatrix'), names(rowData(obj[[i]])))
    ind <- ind[which(!is.na(ind))]
    if (length(ind) > 0)
    rowData(obj[[i]]) <- rowData(obj[[i]])[, -ind]
  
  return(obj[[i]])
}



#' Returns the number of empty lines in a matrix.
#'
#' @title Returns the number of empty lines in the data
#'
#' @param qData A matrix corresponding to the quantitative data.
#' @author Samuel Wieczorek
#'
#' @examples
#' library(QFeatures)
#' data(subR25prot)
#' qData <- SummarizedExperiment::assay(subR25prot[[1]])
#' getNumberOfEmptyLines(qData)
#'
#' @export
#'
#' @return An integer
#'
#'
#'
getNumberOfEmptyLines <- function(qData) {
  n <- sum(apply(is.na(as.matrix(qData)), 1, all))
  return(n)
}





#' @title Loads packages
#' 
#' @description Checks if a package is available to load it
#' 
#' @param ll.deps A `character()` vector which contains packages names
#' 
#' @examples 
#' NULL
#' 
#' @return NA
#' @rdname pkgs.require
#' @export
#' 
#' @author Samuel Wieczorek
#' 
pkgs.require <- function(ll.deps){
  
  if (!requireNamespace('BiocManager', quietly = TRUE)) {
    stop(paste0("Please run install.packages('BiocManager')"))
  }
  
  lapply(ll.deps, function(x) {
    if (!requireNamespace(x, quietly = TRUE)) {
      stop(paste0("Please install ", x, ": BiocManager::install('", x, "')"))
    }
  })
}
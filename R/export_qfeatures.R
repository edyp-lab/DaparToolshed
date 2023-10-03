#' @title Exports a `QFeatures` object to a Excel file.
#'
#' @description
#'
#' This function exports an instance of the class `QFeatures` to a Excel file.
#' The resulting file is composed of four sheets:
#'
#' - `quantitative data` which contains the content of `assay()` object with a
#' color code for each cell w.r.t. to cell quantitative metadata.
#'
#' - `metadata` which is the content of `rowData()` with only one-dimensionnal
#' data (i.e. the adjacencyMatrix and the qMetacell slots are not part of
#' the sheet),
#'
#' - `exp. design` which is the content of `colData()`. Each condition in the 
#' table is colored with a different color,
#'
#' - `quantitative metadata` which is the content of `qMetacell()`. There is a 
#' color code for the different tags.
#'
#'
#' @param object An object of class `QFeatures`.
#' @param i xxx
#' @param filename A character string for the name of the Excel file.
#' @param exp.design xxx
#' @param ... Additional arguments
#' @param wb A workbook
#' @param n A `integer(1)` which is the number of sheet in the workbook.
#' @param tags xxx
#' @param colors A `character()` which contains the HEX code for colors. 
#' The size of this vector must be the same as the number of tags.
#'
#' @return A Excel file.
#'
#' @name QFeatures-excel
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(ft, package='DaparToolshed')
#'
#' #---------------------------------------
#' # Export the whole dataset
#' #---------------------------------------
#'
#' write2excel(ft, filename = "foo")
#' unlink('foo.xls')
#' write2excel(ft, 1, "foo")
#' unlink('foo.xls')
#'
NULL


#' @exportMethod write2excel
#' @rdname QFeatures-excel
setMethod(
    "write2excel", "QFeatures",
    #' @param object xxx
    #' @param i xxx
    #' @param filename xxx
    #' @param ... xxx
    function(object,
             i = NULL,
             filename = "newFile", ...) {
        if (length(object) == 0) {
            return(NULL)
        }


        if (is.null(i)) {
            # One exports all the QFeatures object
        } else {
            # One exports only one SE
            write2excel(
                object[[i]],
                filename,
                design.qf(object),
                ...
            )
        }
    }
)


#' @exportMethod write2excel
#' @rdname QFeatures-excel
setMethod(
    "write2excel", "SummarizedExperiment",
    
    #' @param object xxx
    #' @param filename xxx
    #' @param exp.design xxx
    #' @param ... xxx
    function(object, filename, exp.design, ...) {
        .write2excel(object, filename, exp.design, ...)
    }
)






#' @title xxx
#' @description xxx
#' 
#' @param object xxx
#' @param filename xxx
#' @param exp.design xxx
#'
#' @rdname QFeatures-excel
#' 
.write2excel <- function(object, filename, exp.design) {
    
  pkgs.require(c("DaparViz", "openxlsx"))
  
    name <- paste0(filename, ".xlsx", sep = "")
    wb <- openxlsx::createWorkbook(name)
    # Write the assay data to the first sheet
    i.sheet <- 1
    openxlsx::addWorksheet(wb, "Quantitative data")
    openxlsx::writeData(wb,
        sheet = i.sheet,
        cbind(
            ID = rowData(object)[, idcol(object)],
            assay(object)
        ),
        rowNames = FALSE
    )


    # Add colors to quantitative table
    mc <- qMetacell.def(typeDataset(object))
    colors <- as.list(stats::setNames(mc$color, mc$node))
    tags <- cbind(
        keyId = rep("identified", nrow(object)),
        qMetacell(object)
    )

    addColors(wb, i.sheet, tags, colors)

    # Write the rowData table to the second sheet
    i.sheet <- 2

    # Write only one-dimensional slots
    # browser()
    openxlsx::addWorksheet(wb, "Exp. design")
    openxlsx::writeData(wb,
        sheet = i.sheet,
        exp.design,
        rowNames = FALSE
    )


    # Add colors for sample data sheet
    u_conds <- unique(exp.design$Condition)
    colors <- stats::setNames(DaparViz::ExtendPalette(length(u_conds)), u_conds)
    colors[["blank"]] <- "white"

    tags <- as.data.frame(exp.design)
    tags[, ] <- "blank"
    tags$Sample.name <- exp.design$Condition
    tags$Condition <- exp.design$Condition

    addColors(wb, i.sheet, tags, colors)

    #
    # ## Add the experimental design to the third sheet

    n <- 3
    oneDim <- which(lapply(rowData(object), is.vector) == 1)
    new.rowData <- rowData(object)[, oneDim]
    openxlsx::addWorksheet(wb, "rowData")
    openxlsx::writeData(wb,
        sheet = n,
        cbind(
            ID = new.rowData[, idcol(object)],
            as.data.frame(new.rowData)
        ),
        rowNames = FALSE
    )


    # Add the qMetacell information
    n <- 4
    new.rowData <- qMetacell(object)
    openxlsx::addWorksheet(wb, "qMetacell")
    openxlsx::writeData(wb,
        sheet = n,
        cbind(
            ID = rowData(object)[, idcol(object)],
            new.rowData
        ),
        rowNames = FALSE
    )


    colors <- as.list(stats::setNames(mc$color, mc$node))
    tags <- cbind(
        keyId = rep("identified", nrow(new.rowData)),
        new.rowData
    )

    tags[, ] <- "identified"
    tags[, 1 + seq_len(ncol(new.rowData))] <- new.rowData

    addColors(wb, n, tags, colors)




    openxlsx::saveWorkbook(wb, name, overwrite = TRUE)
    return(name)
}


#'
#' @rdname QFeatures-excel
#'
addColors <- function(wb, n, tags, colors) {
    unique.tags <- NULL
    if (!is.null(tags) && !is.null(colors)) {
        unique.tags <- unique(as.vector(as.matrix(tags)))
        conds.colors <- sum(unique.tags %in% names(colors))
        conds.colors <- conds.colors == length(unique.tags)

        if (conds.colors) {
            lapply(seq_len(length(colors)), function(x) {
                list.tags <- which(names(colors)[x] == tags, arr.ind = TRUE)
                openxlsx::addStyle(wb,
                    sheet = n,
                    cols = list.tags[, "col"],
                    rows = list.tags[, "row"] + 1,
                    style = openxlsx::createStyle(fgFill = colors[x])
                )
            })
        } else {
            warning("The length of colors vector must be equal to the number 
            of different tags. As is it not the case, colors are ignored")
        }
    }
}

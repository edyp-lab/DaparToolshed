#' @title This function exports a data.frame to a Excel file.
#' 
#' @description 
#' This function exports a \code{MSnSet} data object to a Excel file.
#' Each of the three data.frames in the \code{MSnSet} object (ie experimental
#' data, phenoData and metaData are respectively integrated into separate sheets
#' in the Excel file). 
#' 
#' The colored cells in the experimental data correspond to the original
#' missing values which have been imputed.
#' 
#' @param wb A Workbook object containing a worksheet.
#' @param obj An instance of the class `QFeatures`
#' @param i An integer which is the index of the assay in the QFeatures object
#' @param n The total number of sheets
#' @param filename A character string for the name of the Excel file.
#' @param file The name of the Excel file.
#' @param sheet The worksheet to write to. Can be the worksheet index or name.
#' 
#' @return A Excel file (.xlsx)
#'
#' @author Samuel Wieczorek
#' @name output_2_Excel
#' @examples
#' \donttest{
#' library(QFeatures)
#' data(subR25prot)
#' df <- SummarizedExperiment::assay(subR25prot[[1]])
#' tags <- qMetacell(subR25prot[[1]])
#' colors <- list(
#'     "Missing POV" = "lightblue",
#'     "Missing MEC" = "orange",
#'     "Quant. by recovery" = "lightgrey",
#'     "Quant. by direct id" = "white",
#'     "Combined tags" = "red"
#' )
#' file <- tempfile('toto.xlsx')
#' write.excel(subR25prot, filename = file)
#' unlink(file)
#' 
#' data(subR25pept)
#' file <- tempfile('foo.xlsx')
#' write.excel(subR25pept, file)
#' unlink(file)
#' }
#'
#' 
NULL



#' @rdname output_2_Excel
#'
#' @export
#' 
#'
readExcel <- function(file, sheet = NULL) {
  pkgs.require('readxl')
  
  if(is.null(sheet))
    return(NULL)
  
  data <- NULL
  data <- readxl::read_excel(file, 
    sheet,
    col_types = 'guess')
  
  return(
    as.data.frame(
      data, 
      asIs = TRUE, 
      stringsAsFactors = FALSE
    )
  )
}



#' @rdname output_2_Excel
#' @export
#'
listSheets <- function(file) {
  pkgs.require('openxlsx')
  return(openxlsx::getSheetNames(file))
}





#' @export
#' @rdname output_2_Excel
#' @import QFeatures
#' @importFrom stats setNames
#' 
write_Assay_To_Excel <- function(wb, obj, i, n){
  .name <- paste0(names(obj)[i], '_quantdata')
  
  openxlsx::addWorksheet(wb, .name)
  
  data <- SummarizedExperiment::assay(obj[[i]])
  openxlsx::writeData(wb, 
    sheet = n, 
    cbind(
      ID = rownames(data),
      data), 
    rowNames = FALSE)
  
  
  # Add colors to quantitative table
  mc <- metacell.def(typeDataset(obj[[i]]))
  colors <- as.list(stats::setNames(mc$color, mc$node))
  tags <- cbind(
    keyId = rep("Any", nrow(obj[[i]])),
    qMetacell(obj[[i]])
  )
  

  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors) && !is.null(qMetacell(obj[[i]]))) {
    unique.tags <- unique(as.vector(as.matrix(tags)))
    test <- sum(unique.tags %in% names(colors)) == length(unique.tags)
    if (!isTRUE(test)) {
      warning("The length of colors vector must be equal to the number 
            of different tags. As is it not the case, colors are ignored")
    } else {
       
      lapply(seq_len(length(colors)), function(x) {
        list.tags <- which(names(colors)[x] == tags, arr.ind = TRUE)

        openxlsx::addStyle(wb,
          sheet = n,
          cols = list.tags[, "col"],
          rows = list.tags[, "row"] + 1,
          style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  }
  return(wb)
}



#' @export
#' @rdname output_2_Excel
#' 
WriteHistory <- function(wb, obj, n){
  
  openxlsx::addWorksheet(wb, 'history')
  
  openxlsx::writeData(wb, 
    sheet = n, 
    DaparToolshed::paramshistory(obj[[length(obj)]]), 
    rowNames = FALSE)
  
  return(wb)
}


#' @export
#' @rdname output_2_Excel
#' @import SummarizedExperiment
#' @importFrom stats setNames
#' 
Write_SamplesData_to_Excel <- function(wb, obj, n){
  
  openxlsx::addWorksheet(wb, "Samples Meta Data")
  openxlsx::writeData(wb, 
    sheet = n, 
    SummarizedExperiment::colData(obj), 
    rowNames = FALSE)
  
  
  # Add colors for sample data sheet
  u_conds <- unique(design.qf(obj)$Condition)
  colors <- stats::setNames(ExtendPalette(length(u_conds)), u_conds)
  colors[["blank"]] <- "white"
  
  tags <- as.data.frame(SummarizedExperiment::colData(obj))
  tags[, ] <- "blank"
  tags$quantCols <- design.qf(obj)$Condition
  tags$Condition <- design.qf(obj)$Condition
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)) {
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(
      sum(unique.tags %in% names(colors)) == length(unique.tags))) {
      warning("The length of colors vector must be equal to the number 
            of different tags. As is it not the case, colors are ignored")
    } else {
      lapply(seq_len(length(colors)), function(x) {
        list.tags <- which(names(colors)[x] == tags, arr.ind = TRUE)
        openxlsx::addStyle(wb,
          sheet = n,
          cols = list.tags[, "col"],
          rows = list.tags[, "row"] + 1,
          style = openxlsx::createStyle(fgFill = colors[x])
        )
      })
    }
  }
  
  return(wb)
}



#' @export
#' @rdname output_2_Excel
#' 
Write_RowData <- function(wb, obj, i, n){
  .name <- paste0(names(obj)[i], '_rowdata')
  openxlsx::addWorksheet(wb, .name)

  obj[[i]] <- CleanRowData(obj, i)

  openxlsx::writeData(wb, 
    sheet = n, 
    SummarizedExperiment::rowData(obj[[i]]), 
    rowNames = FALSE)
  
  
  # Add colors to quantitative table
  # mc <- metacell.def(omXplore::get_type(obj[[i]]))
  # colors <- as.list(stats::setNames(mc$color, mc$node))
  # tags <- cbind(
  #   keyId = rep("Quant. by direct id", nrow(obj[[i]])),
  #   omXplore::get_metacell(obj[[i]])
  # )
  # 
  # 
  # unique.tags <- NULL
  # if (!is.null(tags) && !is.null(colors)) {
  #   unique.tags <- unique(as.vector(as.matrix(tags)))
  #   if (!isTRUE(
  #     sum(unique.tags %in% names(colors)) == length(unique.tags))) {
  #     warning("The length of colors vector must be equal to the number 
  #           of different tags. As is it not the case, colors are ignored")
  #   }
  #   if (isTRUE(
  #     sum(unique.tags %in% names(colors)) == length(unique.tags))) {
  #     lapply(seq_len(length(colors)), function(x) {
  #       list.tags <- which(names(colors)[x] == tags, arr.ind = TRUE)
  #       openxlsx::addStyle(wb,
  #         sheet = n,
  #         cols = list.tags[, "col"],
  #         rows = list.tags[, "row"] + 1,
  #         style = openxlsx::createStyle(fgFill = colors[x])
  #       )
  #     })
  #   }
  # }
  
  return(wb)
}



#' @export
#' @rdname output_2_Excel
#' 
#' @import openxlsx
#'
write.excel <- function(obj, filename) {

  if (!inherits(obj, "QFeatures")){
    message('Obj is not a QFeatures')
  return(NULL)
  }

  
  name <- filename
  wb <- openxlsx::createWorkbook(filename)
  
  
  n <- 1

  wb <- WriteHistory(wb, obj, n)
  n <- n + 1
  wb <- Write_SamplesData_to_Excel(wb, obj, n)
 
  for (i in seq(length(obj))){
    n <- n + 1
    wb <- write_Assay_To_Excel(wb, obj, i, n)
    n <- n + 1
    wb <- Write_RowData(wb, obj, i, n)
  }

  ## Add feature Data sheet

  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
  return(filename)
}


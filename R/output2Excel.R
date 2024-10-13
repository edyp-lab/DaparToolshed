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
#' @param wb xxxx
#' @param obj xxx
#' @param i xxx
#' @param n xxx
#' @param df An data.frame
#'
#' @param tags xxx
#'
#' @param colors xxx
#'
#' @param tabname xxx
#'
#' @param filename A character string for the name of the Excel file.
#' @param file The name of the Excel file.
#'
#' @param sheet The name of the sheet
#' 
#' @return A Excel file (.xlsx)
#'
#' @author Samuel Wieczorek
#' @name output_2_Excel
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' df <- Biobase::exprs(Exp1_R25_pept[seq_len(100)])
#' tags <- GetMetacell(Exp1_R25_pept[seq_len(100)])
#' colors <- list(
#'     "Missing POV" = "lightblue",
#'     "Missing MEC" = "orange",
#'     "Quant. by recovery" = "lightgrey",
#'     "Quant. by direct id" = "white",
#'     "Combined tags" = "red"
#' )
#' write.excel(df, tags, colors, filename = "toto")
#' 
#' 
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' write.excel(obj, "foo.xlsx")
#' 
#'
#' 
NULL



#' @rdname output_2_Excel
#'
#' @export
#'
readExcel <- function(file, sheet=NULL) {
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
#'
listSheets <- function(file) {
  pkgs.require('openxlsx')
  return(openxlsx::getSheetNames(file))
}





#' @export
#' @rdname output_2_Excel
#' @import omXplore
#' 
write_Assay_To_Excel <- function(wb, obj, i, n){
  .name <- paste0(names(obj)[i], '_quantdata')
  
  openxlsx::addWorksheet(wb, .name)
  
  data <- assay(obj[[i]])
  openxlsx::writeData(wb, 
    sheet = n, 
    cbind(
      ID = rownames(data),
      data), 
    rowNames = FALSE)
  
  
  # Add colors to quantitative table
  mc <- metacell.def(omXplore::get_type(obj[[i]]))
  colors <- as.list(stats::setNames(mc$color, mc$node))
  tags <- cbind(
    keyId = rep("Quant. by direct id", nrow(obj[[i]])),
    omXplore::get_metacell(obj[[i]])
  )
  
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)) {
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(
      sum(unique.tags %in% names(colors)) == length(unique.tags))) {
      warning("The length of colors vector must be equal to the number 
            of different tags. As is it not the case, colors are ignored")
    }
    if (isTRUE(
      sum(unique.tags %in% names(colors)) == length(unique.tags))) {
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
  
  dfHistory <- data.frame(process = c(), parameter = c(), value = c())

  for (i in seq(length(obj))){

    paramName <- names(paramshistory(obj[[i]]))
    if(is.null(paramName)){
      paramName <- '-'
      paramValue <- '-'
      dfHistory <- rbind(dfHistory, 
        data.frame(process = names(obj)[i], 
          parameter = paramName, 
          value = paramValue)
      )
    } else {
    
      for (x in paramName){
        dfHistory <- rbind(dfHistory, 
          data.frame(process = names(obj)[i], 
            parameter = x, 
            value = unname(unlist(paramshistory(obj[[i]])[[x]]))
            )
        )
      }
    
  }
  
  }

  openxlsx::writeData(wb, 
    sheet = n, 
    dfHistory, 
    rowNames = FALSE)
  
  return(wb)
}


#' @export
#' @rdname output_2_Excel
#' @import omXplore
#' 
Write_SamplesData_to_Excel <- function(wb, obj, n){
  
  openxlsx::addWorksheet(wb, "Samples Meta Data")
  openxlsx::writeData(wb, 
    sheet = n, 
    SummarizedExperiment::colData(obj), 
    rowNames = FALSE)
  
  
  # Add colors for sample data sheet
  u_conds <- unique(omXplore::get_group(obj))
  colors <- stats::setNames(ExtendPalette(length(u_conds)), u_conds)
  colors[["blank"]] <- "white"
  
  tags <- as.data.frame(SummarizedExperiment::colData(obj))
  tags[, ] <- "blank"
  tags$Sample.name <- omXplore::get_group(obj)
  tags$Condition <- omXplore::get_group(obj)
  
  unique.tags <- NULL
  if (!is.null(tags) && !is.null(colors)) {
    unique.tags <- unique(as.vector(as.matrix(tags)))
    if (!isTRUE(
      sum(unique.tags %in% names(colors)) == length(unique.tags))) {
      warning("The length of colors vector must be equal to the number 
            of different tags. As is it not the case, colors are ignored")
    }
    if (isTRUE(
      sum(unique.tags %in% names(colors)) == length(unique.tags))) {
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
  .name <- paste0(names(obj)[i], '_coldata')
  openxlsx::addWorksheet(wb, .name)

  openxlsx::writeData(wb, 
    sheet = n, 
    rowData(obj[[i]]), 
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
write.excel <- function(obj, filename) {
  pkgs.require(c('stats', 'openxlsx'))
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


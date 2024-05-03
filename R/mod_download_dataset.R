#' @title dl
#'
#' @description  A shiny Module.
#' 
#' 
#' @param id internal
#' @param dataIn internal
#' @param extension Available values are `csv` (default), `RData` and `Excel`.
#' @param widget.type Available values are `Button` and `Link` (default).
#' @param name internal
#' @param excel.style xxx
#'
#' @return NA
#'
#' @name download_dataset
#' @examplesIf interactive()
#' data(ft_na)
#' shiny::runApp(download_dataset(ft_na))
#'
NULL


#' @import shiny
#'
#' @rdname download_dataset
#'
#' @export
#'
download_dataset_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3('--- Default donwload dataset tool ---'),
    uiOutput(ns('dl_xl')),
    uiOutput(ns('dl_csv')),
    uiOutput(ns('dl_raw'))
  )
}

#' @rdname download_dataset
#'
#' @export
#'
download_dataset_server <- function(id,
  dataIn = reactive({NULL}),
  extension = c('csv', 'xlsx', 'RData'),
  widget.type = 'Link',
  name = 'foo', 
  excel.style = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      UI_type = NULL,
      export_file_RData = NULL,
      export_file_csv = NULL,
      export_file_xlsx = NULL
    )
    
    observeEvent(dataIn(), ignoreNULL = TRUE,{
      rv$export_file_csv <- tryCatch({
        out.csv <- tempfile(fileext = ".csv")
        write.csv(x = dataIn(), file = out.csv)
        out.csv
      },
        warning = function(w) NULL,
        error = function(e) NULL
      )
      
      rv$export_file_xlsx <- tryCatch({
        out.xlsx <- tempfile(fileext = ".xlsx")
        write.excel(obj = dataIn(), filename = out.xlsx)
        out.xlsx
      },
        warning = function(w) NULL,
        error = function(e) NULL
      )
      
      
      rv$export_file_RData <- tryCatch({
        out.RData <- tempfile(fileext = ".RData")
        saveRDS(dataIn(), file = out.RData)
        out.RData
      },
        warning = function(w) NULL,
        error = function(e) NULL
      )
      
      
      print(rv$export_file_xlsx)
      print(rv$export_file_RData)
      print(rv$export_file_csv)
    })
    
    GetType <- reactive({
      if(length(extension) != length(widget.type)){
        warning("Widget.type is not correctly configured. As one cannot decide, 
            all values are set to default ('Link')")
        rv$UI_type <- rep('Link', length(extension))
      } else {
        rv$UI_type <- widget.type
      }
      
      rv$UI_type
    })
    
    output$dl_csv <- renderUI({
      req('csv' %in% extension)
      req(rv$export_file_csv)
      type <- GetType()[which(extension == 'csv')]
      
      do.call(paste0('download', type),
        list(
          ns("downloadDatacsv"),
          "csv",
          class = if (type=='Button') actionBtnClass else ''
        )
      )
    })
    
    
    output$dl_xl <- renderUI({
      req('xlsx' %in% extension)
      req(rv$export_file_xlsx)
      type <- GetType()[which(extension == 'xlsx')]
      do.call(paste0('download', type),
        list(
          ns("downloadDataExcel"), 
          "xlsx",
          class = if (type=='Button') actionBtnClass else ''
        )
      )
    })
    
    output$dl_raw <- renderUI({
      req('RData' %in% extension)
      req(rv$export_file_RData)
      type <- GetType()[which(extension == 'RData')]
      do.call(paste0('download', type),
        list(
          ns("downloadDataRData"), 
          "RData",
          class = if (type=='Button') actionBtnClass else ''
        )
      )
    })
    
    output$downloadDatacsv <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        file.copy(
          from = rv$export_file_csv,
          to = file
        )
      }
    )
    
    output$downloadDataRData <- downloadHandler(
      filename = function() {
        paste ("data-", Sys.Date(), ".RData", sep = "")
      },
      content = function(file) {
        file.copy(
          from = rv$export_file_RData,
          to = file
        )
      }
    )
    
    output$downloadDataExcel <- downloadHandler(
      filename = function() {
        paste("data-", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(file) {
        file.copy(
          from = rv$export_file_xlsx,
          to = file
        )
      }
    )
  }
  )
}




#' @rdname download_dataset
#'
#' @export
#'
download_dataset <- function(data){
  ui <- download_dataset_ui("dl")
  
  server <- function(input, output, session) {
    
    download_dataset_server("dl",
      dataIn = reactive({data}),
      extension = c('csv', 'xlsx', 'RData')
    )
  }
  
  app <- shiny::shinyApp(ui = ui, server = server)
}

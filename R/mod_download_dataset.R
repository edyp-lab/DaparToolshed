#' @title Download_btns shiny module.
#' 
#' @description A shiny module that shows download buttons in different formats.
#' 
#' @param id internal
#' @param settings xxx
#' @param obj A data.frame
#' @param name xxxx.
#' @param colors xxx
#' @param tags xxx
#'
#'
#' @name download_btns
#' @examplesIf interactive()
#' data(lldata)
#' shiny::runApp(download_btns(lldata))
#'
NULL


#' @rdname download_btns
#'
#' @export
#'
download_btns_ui <- function(id, settings = list()) {
  ns <- NS(id)
  
  tagList(
    downloadButton(ns("download_as_Excel_btn"), "Excel", class = settings$actionBtnClass),
    downloadButton(ns("download_as_csv_btn"), "csv", class = settings$actionBtnClass),
    downloadButton(ns("download_as_RData_btn"), "RData", class = settings$actionBtnClass)
  )
}




#' @rdname download_btns
#' @return NA
#'
#' @export
#'
download_btns_server <- function(id,
  obj = reactive({NULL}), 
  name, 
  colors = reactive({NULL}), 
  tags = reactive({NULL})) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$download_as_csv_btn <- downloadHandler(
      filename = function() {
        paste(name(), "-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.table(obj(), file, sep = ";", row.names = FALSE)
      }
    )
    
    output$download_as_RData_btn <- downloadHandler(
      filename = function() {
        paste ("data-", Sys.Date(), ".RData", sep = "")
      },
      content = function(fname) {
        saveRDS(obj(), file=fname)
      }
    )
    
    output$download_as_Excel_btn <- downloadHandler(
      filename = function() {
        paste(name(), "-", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(file) {
        fname <- paste("temp", Sys.Date(), ".xlsx", sep = "")
        write.excel(
          df = obj(),
          colors = colors(),
          tags = tags(),
          filename = fname
        )
        
        file.rename(fname, file)
      }
    )
  }
  )
}


#' @export
#' @rdname download_btns
#' 
download_btns <- function(obj){
  ui <- fluidPage(
    download_btns_ui(id = 'ex', 
      settings = list(actionBtnClass = "btn-primary"))
  )
  
  server <- function(input, output) {
    
    download_btns_server(id = "ex",
      data = reactive({obj}),
      name = reactive({"myTest"}),
      colors = reactive({NULL}),
      tags = reactive({NULL
      })
    )
  }
  
  shinyApp(ui, server)
}

# Module UI

#' @title   open_dataset_ui and open_dataset_server
#' 
#' @description  A shiny Module.
#' 
#' @param id xxx
#' 
#' @name open_dataset
#'
#' @keywords internal
#' 
#' @examples 
#' if (interactive()){
#' ui <- fluidPage(
#' tagList(
#'   open_dataset_ui("qf_file"),
#'   textOutput('res')
#' )
#' )
#' 
#' server <- function(input, output, session) {
#'   rv <- reactiveValues(
#'     obj = NULL
#'   )
#'   rv$obj <- open_dataset_server("qf_file")
#'   
#'   output$res <- renderText({
#'     rv$obj()
#'     paste0('Names of the datasets: ', names(rv$obj()))
#'   })
#' }
#' 
#' shinyApp(ui, server)
#' }
#' 
NULL




#' @export 
#' @rdname open_dataset
#' @importFrom shiny NS tagList 
#' @import shinyjs
#' 
open_dataset_ui <- function(id){
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "Open file", multiple = FALSE),
    actionButton(ns('load_qf_btn'), 'Load file')
  )
}


#' @rdname open_dataset
#' @export
#' @importFrom shinyjs info
#' @import QFeatures
#' 
open_dataset_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    rv.openqf <- reactiveValues(
      dataRead = NULL,
      dataOut = NULL
    )

    ## -- Open a MSnset File --------------------------------------------
    observeEvent(input$load_qf_btn, ignoreInit = TRUE, {
      input$file
      
      authorizedExtension <- "qf"
      
      warn.wrong.file <- "Warning : this file is not a QFeatures file !
       Please choose another one."
      tryCatch(
        {
          if (length(grep(GetExtension(input$file$datapath),
                          authorizedExtension,ignore.case = TRUE)) == 0) {
            warning(warn.wrong.file)
          } else {
            rv.openqf$dataRead <- readRDS(input$file$datapath)
          }
          
          if (!inherits(rv.openqf$dataRead,'QFeatures')) {
            warning(warn.wrong.file)
          }
          
          rv.openqf$dataOut <- rv.openqf$dataRead
        },
        warning = function(w) {
          shinyjs::info(conditionMessage(w))
          return(NULL)
        },
        error = function(e) {
          shinyjs::info(conditionMessage(e))
          return(NULL)
        },
        finally = {
          # cleanup-code
        }
      )
      
    })

    reactive({rv.openqf$dataOut})
  })
  
}






#----------------------------------------------------

ui <- fluidPage(
  tagList(
    open_dataset_ui("qf_file"),
    textOutput('res')
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    obj = NULL,
    result = NULL
  )
  
 
    rv$result <- open_dataset_server("qf_file")
  
  observeEvent(req(rv$result()), {
    rv$obj <- rv$result()
  })
  
  
  output$res <- renderText({
    rv$obj
    print(rv$obj)
    paste0('Names of the datasets: ', names(rv$obj))
  })
}

shinyApp(ui, server)

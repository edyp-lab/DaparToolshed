#' @title   mod_open_demo_dataset_ui and mod_open_demo_dataset_server
#' 
#' @description  A shiny Module.
#' 
#' @param id xxx
#' 
#' @name open_demo_dataset
#'
#' @keywords internal
#' 
NULL




#' @export 
#' @rdname open_demo_dataset
#' @importFrom shiny NS tagList 
#' @import shinyjs
#' 
open_demoDataset_ui <- function(id){
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    tagList(
      uiOutput(ns("chooseDemoDataset")),
      uiOutput(ns("linktoDemoPdf")),
      shinyjs::disabled(
          actionButton(ns('load_dataset_btn'), 'Load dataset', 
                       class= actionBtnClass))
    )
  )
}


#' @rdname open_demo_dataset
#' 
#' @export
#' @keywords internal
#' 
#' @importFrom BiocGenerics get
#' @importFrom utils data
#' @importFrom shinyjs info
#' @import QFeatures
#' 
open_demoDataset_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    
    
    .package <- 'DaparToolshedData'
    
    rv.openDemo <- reactiveValues(
      dataRead = NULL,
      dataOut = NULL
    )
    
    
    ### function for demo mode
    output$chooseDemoDataset <- renderUI({
      
      pkgs.require(.package)
      
      selectInput(ns("demoDataset"),
                  "Demo dataset",
                  choices = c('None', utils::data(package=.package)$results[,"Item"]),
                  selected = character(0),
                  width='200px')
    })
    
    
    
    observeEvent(req(input$demoDataset != 'None'), {
      nSteps <- 1
      withProgress(message = '',detail = '', value = 0, {
        incProgress(1/nSteps, detail = 'Loading dataset')
        utils::data(list=input$demoDataset, package=.package)
        rv.openDemo$dataRead <- BiocGenerics::get(input$demoDataset)
        if (!inherits(rv.openDemo$dataRead, "QFeatures")) {
          shinyjs::info("Warning : this file is not a QFeatures file ! 
                      Please choose another one.")
          return(NULL)
        }
        shinyjs::toggleState('load_dataset_btn', condition = !is.null(rv.openDemo$dataRead))
      }) # End withProgress
    }) # End observeEvent
    
    
    observeEvent(input$load_dataset_btn, {
      rv.openDemo$dataOut <- rv.openDemo$dataRead
     })
    
    output$linktoDemoPdf <- renderUI({
      req(input$demoDataset)
      
    })
    
    reactive({rv.openDemo$dataOut })
  })
  
}




###################################################################
##                                                               ##
##                                                               ##
###################################################################

library(shiny)
library(DaparToolshedData)
library(shinyjs)

ui <- fluidPage(
  tagList(
    open_demoDataset_ui("demo"),
    htmlOutput('res')
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    obj = NULL
  )
  
  rv$obj <- open_demoDataset_server("demo")
  
  output$res <- renderText({
    rv$obj()
    HTML(paste0(tags$br(tags$strong('Names of the datasets: ')), names(rv$obj())))
  })
}

shinyApp(ui = ui, server = server)

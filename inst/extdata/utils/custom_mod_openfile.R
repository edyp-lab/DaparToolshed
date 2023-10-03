#' @title Read a file
#' 
#' @description This function is the implementation of the empty readFiles
#' function (used in MagellantNTK)
#' 
#' 
#' @param name xxx
#' @param path xxx
#' 
#' @export
#' 
readFile <- function(name, path){
  
  ext <- unlist(strsplit(name, '.', fixed=TRUE))[2]
  object <- readRDS(path)
  
  return(object)
  
}





#' @export
custom_openfile_ui <- function(id){
  ns <- NS(id)
  tagList(
    useShinyjs(),
    box(
      title = "Open file", status = "primary", solidHeader = TRUE,
      collapsible = TRUE,
      fileInput(ns('file'), "Select file", multiple = FALSE)
    ),
    
    box(
      title = "Dataset infos", status = "primary", solidHeader = TRUE,
      uiOutput(ns('dataInfosUI'))
    ),
    
    shinyjs::hidden(
      div(id = 'datasetInfosPanel',
      box(
      title = "Dataset infos", status = "primary", solidHeader = TRUE,
      'Work in progress...'
    )
      )
    )
  )
}

#' @export
custom_openfile_server <- function(id, path = reactive({NULL})){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  
  dataOut <- reactiveValues(
    dataset = NULL,
    name = NULL
  )
  
  observeEvent(req(input$file),{
    .extension <- unlist(strsplit(input$file$name, split='.', fixed = TRUE))[2]
    ext_dataset <- c('RData', 'rdata')
    ext_file_to_convert <- c('txt', 'csv', 'xls', 'xlsx')
    dataOut$name <- input$file$name
    
    if(.extension %in% ext_dataset)
     dataOut$dataset <- readFile(input$file$name, input$file$datapath)
    else if (.extension %in% ext_file_to_convert)
      dataOut$dataset <- mod_convert_server(input$file$name, input$file$datapath)
     
    show('datasetInfosPanel')
  })
  
  
  output$dataInfosUI <- renderUI({
    req(dataOut$dataset)
    tagList(
      h3('Dataset infos panel'),
      p(paste0('Nb items: ', length(dataOut$dataset)))
    )
  })
  
  return(
    list(
      data = reactive({dataOut$dataset}),
      name = reactive({dataOut$name})
      )
  )
})
}

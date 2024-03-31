#' @title Wrapper to Convert pipeline
#' 
#' @description
#' These functions are inspired by the functions run_workflow() in the package
#' `MagellanNTK`. They wrap the entire workflow into a single function
#' 
#' @examplesIf interactive()
#' shiny::runApp(Convert())
#' 
#' @name Convert_wrapper
#' 
NULL


#' @rdname Convert_wrapper
#' @export
convert_dataset_ui <- function(id){
  requireNamespace('MagellanNTK')
  ns <- NS(id)
    tagList(
      MagellanNTK::nav_ui(ns('PipelineConvert_Convert'))
    )
}


#' @export
#' @rdname Convert_wrapper
#' 
convert_dataset_server <- function(id){
  
  requireNamespace('MagellanNTK')
  

  path <- system.file('extdata/workflow/PipelineConvert', package = 'DaparToolshed')
  dataIn <- NULL
  tl.layout <- NULL
  mode <- "user"
  
  
  MagellanNTK::source_shinyApp_files()
  
  MagellanNTK::source_wf_files(path)
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    dataOut <- reactiveVal()
    session$userData$workflow.path <- path
    
    #observeEvent(dataIn, {
     # session$userData$workflow.path <- path
      
      dataOut(
        nav_server(id = 'PipelineConvert_Convert',
          dataIn = reactive({data.frame()}),
          tl.layout = tl.layout
        )
      )
    #})
    
    return(reactive({dataOut()}))
    
  })
}





#' @rdname Convert_wrapper
#' @export
convert_dataset <- function() {
  
  ui <- convert_dataset_ui('Convert')
  
  server <- function(input, output, session) {
    
    res <- convert_dataset_server('Convert')
    
    observeEvent(req(res()$dataOut()$trigger), {
      print(res()$dataOut()$value)
    })
  }

  app <- shiny::shinyApp(ui, server)
}

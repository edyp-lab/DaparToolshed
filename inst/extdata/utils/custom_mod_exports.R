
#' @title Customizable shiny module for exporting `QFeatures` objects
#' 
#' @param id The `id` parameter for the Shiny module
#' @param object n instance of class `QFeatures`
#' 
#' @name default_export_plugin


#' @export
#' @rdname default_export_plugin
#' 
mod_export_ui <- function(id){
  ns <- NS(id)
  tagList(
    h3('This is the custom export plugin')
  )
}

#' @export
#' @rdname default_export_plugin
#' 
mod_export_server <- function(id, 
  object,
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})
  ){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    
    
  })
}
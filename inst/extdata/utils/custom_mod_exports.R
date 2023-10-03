
#' @title xxxx
#' 
#' @description xx
#' 
#' @param id xxxx
#' @param object xxx
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
mod_export_server <- function(id, object){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    
    
  })
}
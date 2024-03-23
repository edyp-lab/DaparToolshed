
#' @export
mod_convert_ui <- function(id){
  ns <- NS(id)
  tagList(
    useShinyjs(),
    h3('Convert module')
  )
}

#' @export
mod_convert_server <- function(id, file, path){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    dataOut <- reactiveValues(
      dataset = NULL,
      name = NULL
    )

  })
}

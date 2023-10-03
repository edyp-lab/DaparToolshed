#' @title module foo
#' 
#' @description 
#' xxxxx
#' 
#' @name mod-foo
#'
#' @examples 
#' if (interactive()){
#' data(ft_na)
#' ui <- foo_ui('query')
#' 
#' server <- function(input, output, session) {
#'   
#'   rv <- reactiveValues(
#'     res = NULL
#'   )
#'   ll.tags <- c('None' = 'None', 
#'   qMetadata.def(typeDataset(ft_na[[1]]))$node)
#'   
#'   rv$res <- foo_server('query', 
#'   reset = reactive({NULL}),
#'   is.enabled = reactive({NULL})
#'   )
#'   
#'   observeEvent(rv$res$dataOut()$value, ignoreNULL = TRUE, ignoreInit = TRUE, {
#'     print(rv$res$dataOut()$value)
#'     print(rv$res$dataOut()$widgets)
#'   })
#' }
#' 
#' shinyApp(ui=ui, server=server)
#' }
NULL


#' @param id xxx
#' 
#' @rdname foo
#' 
foo_ui <- function(id){
  
  ns <- NS(id)
  tagList(
    uiOutput(ns('widget1_ui')),
    uiOutput(ns("widget2_ui")),
    uiOutput(ns('widget3_ui')),
    uiOutput(ns('valid_btn_ui'))
  )
}


foo_server <- function(id,
  obj,
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})
) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    widget1 = NULL,
    widget2 = NULL,
    widget3 = NULL
  )
  
  
  rv.custom.default.values <- list()
  
  moduleServer(id,function(input, output, session) {
    ns <- session$ns
    
    # DO NOT MODIFY THIS FUNCTION CALL
    eval(
      str2expression(
        Get_AdditionalModule_Core_Code(
          w.names = names(widgets.default.values),
          rv.custom.names = names(rv.custom.default.values)
        )
      )
    )
    
    
    output$widget1_ui <- renderUI({
      widget <- selectInput(ns('widget1'), 'widget1', choices = 1:3)
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    output$widget2_ui <- renderUI({
      widget <- selectInput(ns('widget2'), 'widget2', choices = 1:3)
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    output$widget3_ui <- renderUI({
      widget <- selectInput(ns('widget3'), 'widget3', choices = 1:3)
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$valid_btn_ui <- renderUI({
      widget <- actionButton(ns("valid_btn"), "Validate")
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    observeEvent(input$valid_btn, {
      dataOut$trigger <- as.numeric(Sys.time())
      dataOut$value <- obj()
      dataOut$widgets <- reactiveValuesToList(rv.widgets)
    })
    
    
    return(reactive({dataOut}))
  })
  
}
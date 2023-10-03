
#' @export
custom_EDA_ui <- function(id){
  ns <- NS(id)
  tagList(
    h3('Custom EDA plugin'),
    uiOutput(ns('chooseArrayUI')),
    plotOutput(ns('test'))
  )
}

#' @export
custom_EDA_server <- function(id, object = reactive({NULL})){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$chooseArrayUI <- renderUI({
      req(object())
      if (length(names(object())) == 0) {
        choices <- list(" " = character(0))
      } else {
        choices <- names(object())
      }
      
      selectInput(ns("chooseArray"), "Array",
                  choices = choices,
                  selected = names(object())[length(object())],
                  width = 200
      )
                  
    })
    
  output$test <- renderPlot({
    req(input$chooseArray)
    plot(object()[[input$chooseArray]][,1])
  })
    
  })
}
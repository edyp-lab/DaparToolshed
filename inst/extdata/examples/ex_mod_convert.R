

ui <- fluidPage(
  nav_ui('Convert')
  
)


server <- function(input, output){
  
  rv <- reactiveValues(
    dataIn = NULL,
    dataOut = NULL
  )
  
  observe({
    rv$dataOut <- nav_server(id = 'Convert',
                             dataIn = reactive({data.frame()}))
  })
  
  observeEvent(rv$dataOut, {
    print('toto')
    browser()
  })
}


shinyApp(ui, server)

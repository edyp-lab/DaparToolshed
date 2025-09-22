
library(shiny)
library(shinyBS)

ui <- fluidPage(
  tagList(
    shinyjs::useShinyjs(),
    actionButton('external_reset', 'Reset'),
    mod_query_metacell_ui('query'),
    uiOutput('res'),
    shinyjs::disabled(actionButton('perform', 'Perform')),
    
  )
)

server <- function(input, output) {
  
  
  
  utils::data(subR25prot)
  
  tmp <- mod_query_metacell_server('query', 
                                   obj = reactive({subR25prot}),
                                   reset = reactive({input$external_reset + input$perform}),
                                   op_names = reactive({c('Push p-value', 'Keep original p-value')})
  )
  
  observeEvent(tmp()$trigger, {
    #print(tmp()$indices)
    shinyjs::toggleState("perform",
                         condition = length(tmp()$indices) > 0
    )
    
  })
  
  output$res <- renderUI({
    p(paste0(tmp()$params$MetacellTag, collapse='\n'))
  })
}

shinyApp(ui = ui, server = server)


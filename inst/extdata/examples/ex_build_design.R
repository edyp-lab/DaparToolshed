library(shiny)
library(highcharter)

################################################################
ui <- tagList(
  mod_buildDesign_ui("designEx"),
  
)

server <- function(input, output, session) {
  design <- mod_buildDesign_server("designEx", c('a', 'a', 'a', 'b', 'b', 'b'))
  
  observeEvent(design()$trigger, {
    print(design()$design)
  })
}
shinyApp(ui = ui, server = server)



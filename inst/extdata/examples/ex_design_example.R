library(shiny)

################################################################
ui <- tagList(
  mod_designExample_ui("designEx")
)

server <- function(input, output, session) {
  mod_designExample_server("designEx", 3)
}
shinyApp(ui = ui, server = server)

library(shiny)
library(shinyjs)

ui <- mod_format_DT_ui("dt")

server <- function(input, output, session) {
  mod_format_DT_server("dt", data = iris )
}

shinyApp(ui = ui, server = server)

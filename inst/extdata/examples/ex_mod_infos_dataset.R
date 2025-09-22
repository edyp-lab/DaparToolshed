
library(QFeatures)
library(omXplore)
library(shinyjs)
library(shiny)

data(subR25prot)

ui <- fluidPage(infos_dataset_ui("mod_info"))

server <- function(input, output, session) {
  infos_dataset_server("mod_info", reactive({subR25prot}))
}

shinyApp(ui, server)


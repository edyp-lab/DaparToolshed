
library(QFeatures)
library(DaparViz)
library(shinyjs)
library(shiny)

data(ft_na, package='DaparViz')

ui <- fluidPage(infos_dataset_ui("mod_info"))

server <- function(input, output, session) {
  infos_dataset_server("mod_info", reactive({ft_na}))
}

shinyApp(ui, server)


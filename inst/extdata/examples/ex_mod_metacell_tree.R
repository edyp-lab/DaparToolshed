library(shiny)


ui <- fluidPage(
  mod_metacell_tree_ui('tree')
)

server <- function(input, output) {
  rv <- reactiveValues(
    tags = NULL
  )
  
  observe({
    rv$tags <- mod_metacell_tree_server('tree', type = reactive({'peptide'}))
  })
  
  
}

if (interactive())
  shinyApp(ui = ui, server = server)

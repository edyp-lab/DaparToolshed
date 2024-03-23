library(highcharter)
library(shinyBS)
library(shiny)
library(DaparViz)

ui <- fluidPage(
  mod_ds_metacell_ui('test')
)

server <- function(input, output) {
  data(Exp1_R25_pept, package='DaparToolshedData')
  vList <- convert2Viz(Exp1_R25_pept)
  vData <- vList[[1]]
  
  rv <- reactiveValues(
    tags = NULL
  )
  #pattern <- c('Missing POV', 'Missing MEC')
   pattern <- NULL
   
  observe({
    rv$tags <- mod_ds_metacell_server('test',
                           obj = reactive({vData}),
                           pal = reactive({NULL}),
                           pattern = reactive({pattern}),
                           showSelect = reactive({is.null(pattern)})
    )
  })
  
}

if (interactive())
  shinyApp(ui = ui, server = server)


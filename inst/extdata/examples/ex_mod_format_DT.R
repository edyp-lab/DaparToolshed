library(omXplore)
library(DT)

ui <- fluidPage(
  DTOutput("dt")
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(selected = NULL)

  data(subR25prot)
  obj <- subR25prot[[1]]
  

  mc <- metacell.def(typeDataset(obj))
  colors <- as.list(setNames(mc$color, mc$node))
  dt_style = list(data = qMetacell(obj), colors = colors)
  
  rv$selected <- format_DT_server("dt", 
    dataIn = reactive({assay(obj)}))
  
   observeEvent(rv$selected(), {
     print(paste0('Selected line(s):', rv$selected()))
   })
}

if (interactive())
  shinyApp(ui = ui, server = server)

library(omXplore)
library(DT)

ui <- format_DT_ui("dt")

server <- function(input, output, session) {
  
  rv <- reactiveValues(selected = NULL)

  data(subR25prot)
  obj <- subR25prot[[1]]
  

  mc <- metacell.def(obj@type)
  colors <- as.list(setNames(mc$color, mc$node))
  
  
  
  dt_style = list(data = obj@metacell, colors = colors)
  
  rv$selected <- format_DT_server("dt", 
                                      data = reactive({obj@qdata}),
                                      data_nostyle = reactive({iris[1:21,1:3]}),
                                      withDLBtns = FALSE,
                                      showRownames = FALSE,
                                      dom = 'Bt',
                                      dt_style = reactive({dt_style}),
                                      filename = "export",
                                      selection = 'single'
                                      )
  
   observeEvent(rv$selected(), {
     print(paste0('Selected line(s):', rv$selected()))
   })
}

if (interactive())
  shinyApp(ui = ui, server = server)

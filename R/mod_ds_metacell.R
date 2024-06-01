#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#'
#' @description
#' xxxx
#'
#' @name metacell-plots
#' 
#' @param id xxx
#' @param obj An instance of the class `SummarizedExperiment`
#' @param pal xxx
#' @param pattern xxx
#' @param showSelect xxx 
#' 
#' @return NA
#'
#' @examplesIf interactive()
#' data(ft_na)
#' grp <- omXplore::get_group(ft_na)
#' shiny::runApp(mod_ds_metacell_Histos(ft_na[[1]], group = grp))
#'
NULL

#' @rdname metacell-plots
#' @export
#' 
mod_ds_metacell_Histos_ui <- function(id) {
    ns <- NS(id)
    tagList(
        shinyjs::useShinyjs(),

        uiOutput(ns('chooseTagUI')),
        fluidRow(
            column(width = 4,
                   highchartOutput(ns("histo")), height = "600px"),
            column(width = 4,
                   highchartOutput(ns("histo_per_lines"))),
            column(width = 4,
                   highchartOutput(ns("histo_per_lines_per_conds")))
            )
        )
}


#' @rdname metacell-plots
#' @export
#' 
mod_ds_metacell_Histos_server <- function(id,
                                   obj = reactive({NULL}),
                                   group = reactive({NULL}),
                                    pal = reactive({NULL}), 
                                   pattern = reactive({NULL}),
                                   showSelect = reactive({TRUE})) {
    moduleServer(id, function(input, output, session) {
            ns <- session$ns
            
            addResourcePath(prefix = "img_ds_metacell", 
                            directoryPath = system.file('images', package='DaparToolshed'))
            
            
            
            
            
            rv <- reactiveValues(
                chooseTag = pattern(),
                showSelect = if(is.null(pattern())) TRUE else showSelect(),
                type = NULL
            )
            
            observeEvent(obj(), {
              #rv$type <- omXplore::get_type(obj())
              rv$type <- metadata(obj())$typeDataset
            })
              
              tmp.tags <- mod_metacell_tree_server('tree', 
                obj = reactive({obj()}))

            
            observeEvent(tmp.tags()$values, ignoreNULL = FALSE, ignoreInit = TRUE,{
                rv$chooseTag <- tmp.tags()$values
            })
            
            
            output$chooseTagUI <- renderUI({
                req(obj())
                tagList(
                  p('Select one or several tag(s) to display statistics about'),
                  mod_metacell_tree_ui(ns('tree'))
                )
             })

            output$histo <- renderHighchart({
               tmp <- NULL
               tmp <- metacellHisto_HC(obj(),
                 group = group(),
                 pattern = rv$chooseTag,
                 pal = pal())
                tmp
            })



            output$histo_per_lines <- renderHighchart({
              tmp <- NULL
               tmp <-
                  metacellPerLinesHisto_HC(obj(),
                    group = group(),
                    pattern = rv$chooseTag,
                    indLegend = seq.int(from = 2, to = length(group))
                    )
                # future(createPNGFromWidget(tmp,pattern))
                # })
                tmp
            })



            output$histo_per_lines_per_conds <- renderHighchart({
               tmp <- NULL
                # isolate({
                # pattern <- paste0(GetCurrentObjName(),".MVplot2")
                tmp <- metacellPerLinesHistoPerCondition_HC(obj(),
                  group = group(),
                  pattern = rv$chooseTag,
                  pal = pal()
                  )
                # future(createPNGFromWidget(tmp,pattern))
                # })
                tmp
            })
        }
    )
}



#' @export
#' @rdname metacell-plots
#' 
mod_ds_metacell_Histos <- function(obj,
  group,
  pal,
  pattern,
  showSelect){
ui <- fluidPage(
  mod_ds_metacell_Histos_ui('test')
)

server <- function(input, output) {
  
  rv <- reactiveValues(
    tags = NULL
  )
  pattern <- NULL
  
  observe({
    rv$tags <- mod_ds_metacell_Histos_server('test',
      obj = reactive({obj}),
      group = reactive({group}),
      pal = reactive({NULL}),
      pattern = reactive({pattern}),
      showSelect = reactive({is.null(pattern)})
    )
  })
  
}

app <- shiny::shinyApp(ui = ui, server = server)
}


#' @title Displays a correlation matrix of the quantitative data of a
#' numeric matrix.
#'
#' @description
#' xxxx
#'
#' @name metacell-plots
#' 
#' @param id xxx
#' @param vizData xxx
#' @param pal xxx
#' @param pattern xxx
#' @param showSelect xxx 
#' 
#' @return NA
#'
#' @example inst/extdata/examples/ex_mod_ds_metacell.R
#'
NULL

#' @rdname metacell-plots
#' @export
#' 
mod_ds_metacell_ui <- function(id) {
    ns <- NS(id)
    tagList(
        shinyjs::useShinyjs(),

        p('Select one or several tag(s) to display statistics about'),
        uiOutput(ns('chooseTagUI')),
        fluidRow(
            column(width = 4,
                   highchartOutput(ns("qMetacell")), height = "600px"),
            column(width = 4,
                   highchartOutput(ns("qMetacell_per_lines"))),
            column(width = 4,
                   highchartOutput(ns("qmetacell_per_lines_per_conds")))
            )
        )
}


#' @rdname metacell-plots
#' @export
#' 
mod_ds_metacell_server <- function(id,
                                   vizData = reactive({NULL}),
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
            
            observe({
              req(vizData())
              
              rv$type <- vizData()@type
            })
              
              tmp.tags <- mod_metacell_tree_server('tree', type = reactive({rv$type}))

            
            observeEvent(tmp.tags()$values, ignoreNULL = FALSE, ignoreInit = TRUE,{
                rv$chooseTag <- tmp.tags()$values
            })
            
            
            output$chooseTagUI <- renderUI({
                req(vizData())
                mod_metacell_tree_ui(ns('tree'))
             })

            output$qMetacell <- renderHighchart({
               tmp <- NULL
               tmp <- metacellHisto_HC(vizData = vizData(),
                                       pattern = rv$chooseTag,
                                       pal = pal()
                                        )
                tmp
            })



            output$qMetacell_per_lines <- renderHighchart({
              tmp <- NULL
               tmp <-
                  metacellPerLinesHisto_HC(vizData = vizData(),
                                           pattern = rv$chooseTag,
                                           indLegend = seq.int(from = 2, to = length(vizData()@conds))
                                           )
                # future(createPNGFromWidget(tmp,pattern))
                # })
                tmp
            })



            output$qmetacell_per_lines_per_conds <- renderHighchart({
               tmp <- NULL
                # isolate({
                # pattern <- paste0(GetCurrentObjName(),".MVplot2")
                tmp <- metacellPerLinesHistoPerCondition_HC(vizData = vizData(),
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




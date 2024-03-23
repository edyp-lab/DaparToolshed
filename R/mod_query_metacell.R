#' @title xxxx
#'
#' @description
#' xxxx
#'
#' @name query-metacell
#' 
#' @param id xxx
#' @param obj xxx
#' @param filters xxx
#' @param val_vs_percent xxx
#' @param operator xxx
#' @param reset xxx
#' @param op_names xxx
#' 
#' @return NA
#'
#' @example inst/extdata/examples/ex_mod_query_metacell.R
#'
NULL

options(shiny.reactlog=TRUE) 

#' @rdname query-metacell
#' @import shiny
#' @export
#' 
mod_query_metacell_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
                column(2, mod_metacell_tree_ui(ns('tree'))),
                column(2, uiOutput(ns("Choose_keepOrRemove_ui"))),
                column(2, uiOutput(ns("choose_metacellFilters_ui"))),
                column(6, uiOutput(ns("MetacellFilters_widgets_set2_ui")))
            ),
        fluidRow(
            column(6, uiOutput(ns("metacellFilter_request_ui"))),
            column(3, uiOutput(ns("show_example_ui"))),
            column(3, uiOutput(ns('showApplyBtn')))
        )
    )
}




#' @rdname query-metacell
#' @import shiny
#' @export
#' 
mod_query_metacell_server <- function(id,
                                      obj = reactive({NULL}),
                                      filters = reactive({NULL}),
                                      val_vs_percent = reactive({NULL}),
                                      operator = reactive({NULL}),
                                      reset = reactive({NULL}), 
                                      op_names = reactive({NULL})
                                      ){
    
  pkgs.require(c("shinyjs", "shinyBS"))    
    
    rv.widgets <- reactiveValues(
        MetacellTag = NULL,
        MetacellFilters = "None",
        KeepRemove = "delete",
        metacell_value_th = 0,
        metacell_percent_th = 0,
        val_vs_percent = "Count",
        metacellFilter_operator = "<="
    )
    
    rv <- reactiveValues(
        tags = NULL, 
        tmp.tags = NULL,
        op_names = setNames(nm = c('delete', 'keep')),
        reset_tree = NULL
    )

    dataOut <- reactiveValues(
        trigger = as.numeric(Sys.time()),
        params = list(),
        query = NULL,
        indices = NULL
    )
    
    GetFiltersScope <- function()
        c("Whole Line" = "WholeLine",
          "Whole matrix" = "WholeMatrix",
          "For every condition" = "AllCond",
          "At least one condition" = "AtLeastOneCond"
        )
    

    moduleServer(id, function(input, output, session) {
            ns <- session$ns

            popover_for_help_server("filterScope_help",
                title = "Scope",
                    content = HTML("To filter the missing values, the choice of
                        the lines to be kept is made by different options:
                               <ul>
                               <li><strong>None</strong>: No filtering, the
                        quantitative data is left unchanged.</li>
                        <li><strong>(Remove) Empty lines</strong>: All the
                        lines with 100% of missing values are filtered
                            out.</li>
                            <li><strong>Whole Matrix</strong>: The lines
                        (across all conditions) which contain less quantitative
                        value than a user-defined threshold are kept;</li>
                        <li><strong>For every condition</strong>: The lines
                        for which each condition contain less quantitative
                        value than a user-defined threshold are deleted;</li>
                        <li><strong>At least one condition</strong>: The
                        lines for which at least one condition contain less
                        quantitative value than a user-defined threshold are
                            deleted.</li>
                               </ul>")
            )


            init_rv_widgets <- function(){
                #print('marqueur 2')
                rv.widgets$MetacellTag <- NULL
                rv.widgets$MetacellFilters <- "None"
                rv.widgets$KeepRemove <- "delete"
                rv.widgets$metacell_value_th <- 0
                rv.widgets$metacell_percent_th <- 0
                rv.widgets$val_vs_percent <- "Count"
                rv.widgets$metacellFilter_operator <- "<="
                
                rv$reset_tree <- as.numeric(Sys.time()) 
                }
            
            
            # observeEvent(op_names, {
            #     rv$op_names <- c("delete", "keep")
            #     if(!is.null(op_names))
            #         rv$op_names <- op_names
            # })
            
            observeEvent(reset(), {
                #print("rororororororo")
                init_rv_widgets()
                rv$reset_tree <- as.numeric(Sys.time()) 
                rv$tags <- NULL
                dataOut$trigger <- as.numeric(Sys.time())
                dataOut$params <- list()
                dataOut$query <- NULL
                dataOut$indices <- NULL
            })
            
            
            
            observeEvent(c(id, obj()), {
                #print('marqueur 1')
                if(is.null(obj())){
                    dataOut$trigger <- as.numeric(Sys.time())
                    dataOut$params<- list(
                        MetacellTag = NULL,
                        KeepRemove = NULL,
                        MetacellFilters = NULL,
                        metacell_percent_th = NULL,
                        metacell_value_th = NULL,
                        val_vs_percent = NULL,
                        metacellFilter_operator = NULL
                        )
                    
                    dataOut$query <- ''
                    dataOut$indices <- NULL
                } else {
                    init_rv_widgets()
                    
                    rv$tags <- NULL
                    dataOut$trigger <- as.numeric(Sys.time())
                    dataOut$params <- list()
                    dataOut$query <- NULL
                    dataOut$indices <- NULL
                }
                
                
            }, priority=1000)

            tmp.tags <- mod_metacell_tree_server('tree',
                                                 obj = reactive({obj()}),
                                                 reset = reactive({rv$reset_tree})
                                                 )
            
            observeEvent(tmp.tags()$values, ignoreNULL = FALSE, ignoreInit = TRUE, {
                #print('marqueur 3')
                rv.widgets$MetacellTag <- tmp.tags()$values
            }, priority = 900)
            

            observeEvent(input$ChooseKeepRemove, {
                rv.widgets$KeepRemove <- input$ChooseKeepRemove
            })
            
            observeEvent(input$ChooseMetacellFilters, {
                rv.widgets$MetacellFilters <- input$ChooseMetacellFilters
            })
            
            observeEvent(input$choose_val_vs_percent, {
                rv.widgets$val_vs_percent <- input$choose_val_vs_percent
            })
            
            observeEvent(input$choose_metacell_value_th, {
                rv.widgets$metacell_value_th <- input$choose_metacell_value_th
            })
            
            observeEvent(input$choose_metacell_percent_th, {
                rv.widgets$metacell_percent_th <- input$choose_metacell_percent_th
            })
            
            observeEvent(input$choose_metacellFilter_operator, {
                rv.widgets$metacellFilter_operator <- input$choose_metacellFilter_operator
            })


            
            output$Choose_keepOrRemove_ui <- renderUI({
                req(rv.widgets$MetacellTag)
                
                if(is.null(op_names()))
                    rv$op <- setNames(nm = c('delete', 'keep'))
                else
                    rv$op <- setNames(c('delete', 'keep'), nm = op_names())
                
                
                radioButtons(ns("ChooseKeepRemove"), "Type of filter operation",
                             choices =  rv$op,
                             selected = rv.widgets$KeepRemove
                             )
            })


            output$choose_metacellFilters_ui <- renderUI({
                req(rv.widgets$MetacellTag)
                selectInput(ns("ChooseMetacellFilters"),
                            popover_for_help_ui(ns("filterScope_help")),
                            choices = c("None" = "None", GetFiltersScope()),
                            selected = rv.widgets$MetacellFilters,
                            width = "200px"
                            )
            })



            output$showApplyBtn <- renderUI({
                req(rv.widgets$MetacellTag)
                req(rv.widgets$MetacellFilters != "None")
                
                actionButton(ns('buildQueryBtn'), 'Apply', class = actionBtnClass)
            })
            
            output$show_example_ui <- renderUI({
                req(rv.widgets$MetacellTag)
                req(rv.widgets$MetacellFilters != "None")

                mod_filtering_example_server(id = "filteringExample",
                                             obj = reactive({obj()}),
                                             indices = reactive({CompileIndices()}),
                                             params = reactive({rv.widgets}),
                                             txt = reactive({WriteQuery()})
                                             )

                mod_filtering_example_ui(ns("filteringExample"))
            })


            output$MetacellFilters_widgets_set2_ui <- renderUI({
                req(rv.widgets$MetacellTag)
                req(!(rv.widgets$MetacellFilters %in% c("None", "WholeLine")))

                popover_for_help_server("choose_val_vs_percent_help",
                    title = paste("#/% of values to ", rv.widgets$KeepRemove),
                        content = "Define xxx"
                )

                tagList(
                    fluidRow(
                        column(4,
                            radioButtons(ns("choose_val_vs_percent"),
                                         popover_for_help_ui(ns("choose_val_vs_percent_help")),
                                         choices = setNames(nm = c("Count", "Percentage")),
                                         selected = rv.widgets$val_vs_percent
                                         )
                        ),
                        column(8,
                            selectInput(ns("choose_metacellFilter_operator"),
                                        "Choose operator()",
                                        choices = setNames(nm = SymFilteringOperators()),
                                        selected = rv.widgets$metacellFilter_operator,
                                        width = "100px"
                                        ),
                            uiOutput(ns("choose_value_ui")),
                            uiOutput(ns("choose_percentage_ui"))
                        )
                    )
                )
            })

            output$choose_value_ui <- renderUI({
                req(rv.widgets$val_vs_percent == "Count")
                req(!(rv.widgets$MetacellFilters %in% c("None", "WholeLine")))

                popover_for_help_server("metacell_value_th_help",
                    title = "Count threshold",
                    content = "Define xxx"
                    )

                tagList(
                    popover_for_help_ui(ns("modulePopover_keepVal")),
                    selectInput(ns("choose_metacell_value_th"),
                                popover_for_help_ui(ns("metacell_value_th_help")),
                        choices = getListNbValuesInLines(obj(), type = rv.widgets$MetacellFilters),
                        selected = rv.widgets$metacell_value_th,
                        width = "150px"
                    )
                )
            })



            output$choose_percentage_ui <- renderUI({
                req(rv.widgets$val_vs_percent == "Percentage")
                req(!(rv.widgets$MetacellFilters %in% c("None", "WholeLine")))

                popover_for_help_server("metacell_percent_th_help",
                    title = "Percentage threshold",
                        content = "Define xxx"
                )
                
                tagList(
                    popover_for_help_ui(ns("modulePopover_keepVal_percent")),
                    sliderInput(ns("choose_metacell_percent_th"),
                                popover_for_help_ui(ns("metacell_percent_th_help")),
                        min = 0,
                        max = 100,
                        step = 1,
                        value = rv.widgets$metacell_percent_th,
                        width = "250px"
                    )
                )
            })


            WriteQuery <- reactive({
                req(rv.widgets$MetacellTag)
                
                query <- ''
                
                # Format the variables to be inserted in text
                operation <- names(rv$op)[ which(rv$op == rv.widgets$KeepRemove)]
                nblines <- length(CompileIndices())
                plural <- ''
                if (nblines > 1)
                    plural <- 's'
                tags <- ''
                if (!is.null(rv.widgets$MetacellTag))
                    tags <- paste0(rv.widgets$MetacellTag, collapse=', ')
                method <- switch(rv.widgets$MetacellFilters,
                                 None = '',
                                 WholeLine = 'entire line',
                                 WholeMatrix = "the whole matrix.",
                                 AllCond = "each condition.",
                                 AtLeastOneCond = "at least one condition.")
                arithmetic_op <- rv.widgets$metacellFilter_operator
                
                threshold <- rv.widgets$metacell_value_th
                if (rv.widgets$val_vs_percent == "Percentage")
                    threshold <- paste(as.character(rv.widgets$metacell_percent_th)," %", sep = "")

                # Builds the query
                query <- switch(rv.widgets$MetacellFilters, 
                            None = {
                                query <- "Query is empty."
                                },
                            WholeLine = {
                                query <- "#operation# #nblines# #datatype##plural# that are only (#tags#)."
                                },
                            WholeMatrix = {
                                query <- "#operation# #nblines# #datatype##plural# where number of (#tags#) data #arithmetic_op# #threshold# in #method#"
                                },
                            AllCond = {
                                query <- "#operation# #nblines# #datatype##plural# where number of (#tags#) data #arithmetic_op# #threshold# in #method#"
                                },
                            AtLeastOneCond = {
                                query <- "#operation# #nblines# #datatype##plural# where number of (#tags#) data #arithmetic_op# #threshold# in #method#"
                                }
                            )
                
                
                query <- gsub('#tags#', tags, query)
                query <- gsub('#operation#', operation, query)
                query <- gsub('#nblines#', nblines, query)
                query <- gsub('#arithmetic_op#', arithmetic_op, query)
                query <- gsub('#threshold#', threshold, query)
                query <- gsub('#method#', method, query)
                query <- gsub('#datatype#', typeDataset(obj()), query)
                query <- gsub('#plural#', plural, query)
                
                
                #browser()
                query
            })

            output$metacellFilter_request_ui <- renderUI({
                txt_summary <- paste("You are about to ", WriteQuery())
                tags$p(style = "font-size: small; text-align : center; color: purple;",
                       txt_summary)
            })

            # Set useless widgets to default values
            observeEvent(rv.widgets$MetacellFilters == "WholeLine",{
                rv.widgets$metacell_percent_th <- 0
                rv.widgets$metacell_value_th <- 0
                rv.widgets$val_vs_percent <- "Percentage"
                },  priority = 1000)


            CompileIndices <- reactive({
                req(obj())
                req(rv.widgets$MetacellTag)
                req(rv.widgets$MetacellFilters != "None")

                
                th <- switch(rv.widgets$val_vs_percent,
                    Percentage = rv.widgets$metacell_percent_th / 100,
                    Count = as.integer(rv.widgets$metacell_value_th)
                )
                
                
                GetIndices_MetacellFiltering(
                    obj = obj(),
                    level = typeDataset(obj()),
                    pattern = rv.widgets$MetacellTag,
                    type = rv.widgets$MetacellFilters,
                    percent = rv.widgets$val_vs_percent == "Percentage",
                    op = rv.widgets$metacellFilter_operator,
                    th = th
                )
            })


            
            observeEvent(input$buildQueryBtn, ignoreInit = TRUE, ignoreNULL = TRUE, {
                
                dataOut$trigger <- as.numeric(Sys.time())
                dataOut$params<- list(
                    MetacellTag = rv.widgets$MetacellTag,
                    KeepRemove = rv.widgets$KeepRemove,
                    MetacellFilters = rv.widgets$MetacellFilters,
                    metacell_percent_th = rv.widgets$metacell_percent_th,
                    metacell_value_th = rv.widgets$metacell_value_th,
                    val_vs_percent = rv.widgets$val_vs_percent,
                    metacellFilter_operator = rv.widgets$metacellFilter_operator
                )
                
                dataOut$query <- WriteQuery()
                dataOut$indices <- CompileIndices()
                # reset
                #init_rv_widgets()
                #rv$tags <- NULL
            })

            reactive({dataOut})
        }
    )
}




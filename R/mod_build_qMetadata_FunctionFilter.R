#' @title Build queries for filtering quantitative metadata
#'
#' @description
#'
#' This function is a shiny module to create a list of queries (instances of the
#' class `FunctionFilter` to filter the quantitative metadata of an instance
#' of the class `SummarizedExperiment`)
#'
#' @name query_qMetacell
#'
#' @return As for all modules used with `MagellanNTK`, the return value is a
#' `list()` of two items:
#' - trigger : xxx
#' - value: In this case, it contains a list() of three slots:
#'   - ll.fun: a list() of instances of the class `FunctionFilter`,
#'   - ll.query: a list of `character()` which describe the queries in natural
#'   language,
#'   - ll.widgets.value: a list of the values of widgets.
#'
#' @examples
#' if (interactive()) {
#'     data(ft_na, package='DaparToolshed')
#'     ui <- mod_build_qMetacell_FunctionFilter_ui("query")
#'
#'     server <- function(input, output, session) {
#'         rv <- reactiveValues(
#'             res = NULL
#'         )
#'         ll.tags <- c(
#'             "None" = "None",
#'             qMetacell.def(typeDataset(ft_na[[1]]))$node
#'         )
#'         rv.custom$res <- mod_build_qMetacell_FunctionFilter_server("query",
#'             obj = reactive({
#'                 ft_na[[1]]
#'             }),
#'             conds = reactive({
#'                 colData(ft_na)$Condition
#'             }),
#'             list_tags = reactive({
#'                 ll.tags
#'             }),
#'             keep_vs_remove = reactive({
#'                 setNames(nm = c("delete", "keep"))
#'             }),
#'             val_vs_percent = reactive({
#'                 setNames(nm = c("Count", "Percentage"))
#'             }),
#'             operator = reactive({
#'                 setNames(nm = SymFilteringOperators())
#'             })
#'         )
#'
#'
#'         observeEvent(rv.custom$res$dataOut()$trigger, ignoreNULL = TRUE, 
#'         ignoreInit = TRUE, {
#'             print(rv.custom$res$dataOut()$fun)
#'         })
#'     }
#'
#'     shinyApp(ui = ui, server = server)
#' }
NULL









#' @param id xxx
#'
#' @export
#'
#' @rdname query_qMetacell
#'
mod_build_qMetacell_FunctionFilter_ui <- function(id) {
    ns <- NS(id)
    tagList(
        div(
            fluidRow(
                column(2, uiOutput(ns("chooseTag_ui"))),
                column(2, uiOutput(ns("chooseKeepRemove_ui"))),
                column(2, uiOutput(ns("chooseScope_ui"))),
                column(6, tagList(
                    uiOutput(ns("qMetacellScope_widgets_set2_ui"))
                ))
            ),
            div(
                style = "display:inline-block; 
                vertical-align: middle; align: center;",
                uiOutput(ns("qMetacellScope_request_ui"))
            ),
            actionButton(ns("BuildFilter_btn"), "Add filter")
        )
    )
}




#' @param id xxx
#' @param obj An instance of the class `SummarizedExperiment`
#' @param conds A `character()` which contains the name of the conditions. The
#' length of this vector must be equal to the number of samples in the assay
#' (i.e. number of columns in àssay(obj))
#' @param list_tags xxx
#' @param keep_vs_remove xxx
#' @param val_vs_percent xxx
#' @param operator xxx
#' @param reset A `ìnteger(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
#' @rdname query_qMetacell
#'
#' @export
#'
mod_build_qMetacell_FunctionFilter_server <- function(id,
    obj,
    conds,
    list_tags = reactive({NULL}),
    keep_vs_remove = reactive({NULL}),
    val_vs_percent = reactive({NULL}),
    operator = reactive({NULL}),
    reset = reactive({NULL}),
    is.enabled = reactive({TRUE})) {

    # Define default selected values for widgets
    # This is only for simple workflows
    widgets.default.values <- list(
        tag = "None",
        scope = "None",
        keepRemove = "delete",
        valueTh = 0,
        percentTh = 0,
        valPercent = "Count",
        operator = "<="
    )

    rv.custom.default.values <- list(
        indices = NULL,
        functionFilter = NULL,
        query = list(),
        fun.list = list(),
        widgets.value = list()
    )

    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        eval(
            str2expression(
                MagellanNTK::Get_AdditionalModule_Core_Code(
                    w.names = names(widgets.default.values),
                    rv.custom.names = names(rv.custom.default.values)
                )
            )
        )



        mod_popover_for_help_server("tag_help",
            title = "Nature of data to filter",
            content = "Define xxx"
        )
        help.txt1 <- "To filter the missing values, the choice of the lines to 
        be kept is made by different options:
    <ul>
    <li><strong>None</strong>: No filtering, the quantitative data is left 
    unchanged.</li>
    <li><strong>(Remove) Empty lines</strong>: All the lines with 100% of 
    missing values are filtered out.</li>
    <li><strong>Whole Matrix</strong>: The lines (across all conditions) 
    which contain less quantitative value than a user-defined threshold are 
    kept;</li>
    <li><strong>For every condition</strong>: The lines for which each 
    condition contain less quantitative value than a user-defined threshold 
    are deleted;</li>
    <li><strong>At least one condition</strong>: The lines for which at least 
    one condition contain less quantitative value than a user-defined 
    threshold are deleted.</li>
    </ul>"

        mod_popover_for_help_server("filterScope_help",
            title = "Scope",
            content = HTML(help.txt1)
        )


        output$chooseTag_ui <- renderUI({
            widget <- selectInput(ns("tag"),
                mod_popover_for_help_ui(ns("tag_help")),
                choices = list_tags(),
                selected = rv.widgets$tag,
                width = "200px"
            )
            MagellanNTK::toggleWidget(widget, is.enabled())
        })

        output$chooseKeepRemove_ui <- renderUI({
            req(rv.widgets$tag != "None")
            widget <- radioButtons(ns("keepRemove"),
                "Type of filter operation",
                choices = keep_vs_remove(),
                selected = rv.widgets$keepRemove
            )
            MagellanNTK::toggleWidget(widget, is.enabled())
        })

        output$chooseScope_ui <- renderUI({
            req(rv.widgets$tag != "None")
            widget <- selectInput(ns("scope"),
                mod_popover_for_help_ui(ns("filterScope_help")),
                choices = c(
                    "None" = "None",
                    "Whole Line" = "WholeLine",
                    "Whole matrix" = "WholeMatrix",
                    "For every condition" = "AllCond",
                    "At least one condition" = "AtLeastOneCond"
                ),
                selected = rv.widgets$scope,
                width = "200px"
            )
            MagellanNTK::toggleWidget(widget, is.enabled())
        })




        output$qMetacellScope_widgets_set2_ui <- renderUI({
            req(!(rv.widgets$scope %in% c("None", "WholeLine")))
            mod_popover_for_help_server("chooseValPercent_help",
                title = paste("#/% of values to ", rv.widgets$keepRemove),
                content = "Define xxx"
            )

            widget1 <- radioButtons(ns("valPercent"),
                mod_popover_for_help_ui(ns("chooseValPercent_help")),
                choices = val_vs_percent(),
                selected = rv.widgets$valPercent
            )

            widget2 <- selectInput(ns("operator"),
                "Choose operator",
                choices = operator(),
                selected = rv.widgets$operator,
                width = "100px"
            )


            tagList(
                fluidRow(
                    column(4, MagellanNTK::toggleWidget(widget1, is.enabled())),
                    column(
                        8, MagellanNTK::toggleWidget(widget2, is.enabled()),
                        uiOutput(ns("value_ui")),
                        uiOutput(ns("percentage_ui"))
                    )
                )
            )
        })


        output$value_ui <- renderUI({
            req(rv.widgets$valPercent == "Count")
            req(!(rv.widgets$scope %in% c("None", "WholeLine")))
            mod_popover_for_help_server("value_th_help",
                title = "Count threshold",
                content = "Define xxx"
            )

            widget <- selectInput(ns("valueTh"),
                mod_popover_for_help_ui(ns("value_th_help")),
                choices = getListNbValuesInLines(
                    object = obj(),
                    conds = conds(),
                    type = rv.widgets$scope
                ),
                selected = rv.widgets$valueTh,
                width = "150px"
            )
            tagList(
                mod_popover_for_help_ui(ns("keepVal_help")),
                MagellanNTK::toggleWidget(widget, is.enabled())
            )
        })

        output$percentage_ui <- renderUI({
            req(rv.widgets$valPercent == "Percentage")
            req(!(rv.widgets$scope %in% c("None", "WholeLine")))

            mod_popover_for_help_server("percentTh_help",
                title = "Percentage threshold",
                content = "Define xxx"
            )
            widget <- sliderInput(ns("percentTh"),
                mod_popover_for_help_ui(ns("percentTh_help")),
                min = 0,
                max = 100,
                step = 1,
                value = rv.widgets$percentTh,
                width = "250px"
            )
            tagList(
                mod_popover_for_help_ui(ns("keepVal_percent_help")),
                MagellanNTK::toggleWidget(widget, is.enabled())
            )
        })


        ## -------------------------------------------------------------------
        ##
        ## This function xxx
        ##
        ## -------------------------------------------------------------------
        WriteQuery <- reactive({
            if (rv.widgets$scope == "None") {
                txt_summary <- "No filtering is processed."
            } else if (rv.widgets$scope == "WholeLine") {
                txt_summary <- paste(
                    rv.widgets$keepRemove,
                    "lines that contain only",
                    rv.widgets$tag
                )
            } else {
                text_method <- switch(rv.widgets$scope,
                    WholeMatrix = "the whole matrix.",
                    AllCond = "each condition.",
                    AtLeastOneCond = "at least one condition."
                )

                if (rv.widgets$valPercent == "Count") {
                    text_threshold <- rv.widgets$valueTh
                } else {
                    text_threshold <- paste(as.character(rv.widgets$percentTh), 
                        " %", sep = "")
                }

                txt_summary <- paste(
                    rv.widgets$keepRemove,
                    " lines where number of ",
                    rv.widgets$tag,
                    " data ",
                    rv.widgets$operator,
                    " ",
                    text_threshold,
                    " in ",
                    text_method
                )
            }
            txt_summary
        })


        output$qMetacellFilter_request_ui <- renderUI({
            tags$p(paste("You are going to ", WriteQuery()),
                style = "font-size: small; text-align : center; color: purple;"
            )
        })


        # Set useless widgets to default values
        observeEvent(rv.widgets$scope == "WholeLine",
            {
                rv.widgets$percentThh <- 0
                rv.widgets$valueTh <- 0
                rv.widgets$valPercent <- "Percentage"
            },
            priority = 1000
        )




        BuildFunctionFilter <- reactive({
            req(obj())
            req(rv.widgets$tag != "None")
            req(rv.widgets$scope != "None")

            th <- switch(rv.widgets$valPercent,
                Percentage = rv.widgets$percentTh / 100,
                Count = as.integer(rv.widgets$valueTh)
            )

            ff <- switch(rv.widgets$scope,
                WholeLine = FunctionFilter("qMetacellWholeLine",
                    cmd = rv.widgets$keepRemove,
                    pattern = rv.widgets$tag
                ),
                WholeMatrix = FunctionFilter("qMetacellWholeMatrix",
                    cmd = rv.widgets$keepRemove,
                    pattern = rv.widgets$tag,
                    percent = rv.widgets$valPercent,
                    th = th,
                    operator = rv.widgets$operator
                ),
                AllCond = FunctionFilter("qMetacellOnConditions",
                    cmd = rv.widgets$keepRemove,
                    mode = rv.widgets$scope,
                    pattern = rv.widgets$tag,
                    conds = conds(),
                    percent = rv.widgets$valPercent,
                    th = th,
                    operator = rv.widgets$operator
                ),
                AtLeastOneCond = FunctionFilter("qMetacellOnConditions",
                    cmd = rv.widgets$keepRemove,
                    mode = rv.widgets$scope,
                    pattern = rv.widgets$tag,
                    conds = conds(),
                    percent = rv.widgets$valPercent,
                    th = th,
                    operator = rv.widgets$operator
                )
            )

            ff
        })


        observeEvent(input$BuildFilter_btn, {
            rv.custom$ll.fun <- append(rv.custom$ll.fun, BuildFunctionFilter())
            rv.custom$ll.query <- append(rv.custom$ll.query, WriteQuery())
            rv.custom$ll.widgets.value <- append(
                rv.custom$ll.widgets.value,
                list(reactiveValuesToList(rv.widgets))
            )


            # Append a new FunctionFilter to the list
            dataOut$trigger <- as.numeric(Sys.time())
            dataOut$ll.fun <- rv.custom$ll.fun
            dataOut$ll.query <- rv.custom$ll.query
            dataOut$ll.widgets.value <- rv.custom$ll.widgets.value
        })



        list(
            trigger = reactive({
                dataOut$trigger
            }),
            ll.fun = reactive({
                dataOut$ll.fun
            }),
            ll.query = reactive({
                dataOut$ll.query
            }),
            ll.widgets.value = reactive({
                dataOut$ll.widgets.value
            })
        )
    })
}

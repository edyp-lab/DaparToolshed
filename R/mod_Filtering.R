
#' @title Module Filtering
#'
#' @description
#'
#' xxxxx
#'
#' @section Step 'Quanti metadata filtering':
#'
#' xxxxxxx
#'
#' @section Step 'Variable filtering':
#'
#' xxxxx
#'
#' @section Step 'Save':
#'
#' xxxxxx
#'
#' @seealso The user manual of the package `MagellanNTK`.
#'
#'
#' @name mod_filtering
#' 
#' @return NA
#'
#' @examples
#' if (interactive()) {
#'     run_workflow("Filtering", verbose = TRUE)
#' }
NULL



#' @param id
#'
#' @rdname mod_filtering
#'
#' @export
#'
mod_Filtering_ui <- function(id) {
    ns <- NS(id)
}


#' @param id A `character(1)` which is the 'id' of the module.
#' @param dataIn An instance of the class `QFeatures`
#' @param steps.enabled A `logical()` which indicates whether each step is
#' enabled or disabled in the UI.
#' @param remoteReset A `logical(1)` which acts asa a remote command to reset
#' the module to its default values. Default is FALSE.
#' @param steps.status A `logical()` which indicates the status of each step
#' which can be either 'validated', 'undone' or 'skipped'.
#' enabled or disabled in the UI.
#' @param current.pos A `interger(1)` which acts as a remote command to make
#'  a step active in the timeline. Default is 1.
#' @param verbose A `logical(1)` to indicates whether run and debug infos must
#' be printed in the console. Default is FALSE.
#'
#' @rdname mod_filtering
#' @importFrom shinyjs toggle hidden
#' @importFrom MagellanNTK Timestamp toggleWidget Get_Workflow_Core_Code
#'
#' @export
#'
mod_Filtering_server <- function(id,
    dataIn = reactive({NULL}),
    steps.enabled = reactive({NULL}),
    remoteReset = reactive({FALSE}),
    steps.status = reactive({NULL}),
    current.pos = reactive({1}),
    verbose = FALSE) {

  pkgs.require("DT")
    
    # This list contains the basic configuration of the process
    config <- MagellanNTK::Config(
        name = "Filtering",
        # Define the type of module
        mode = "process",

        # List of all steps of the process
        steps = c(
            "Description",
            "Quanti metadata filtering",
            "Variable filtering",
            "Save"
        ),

        # A vector of boolean indicating if the steps are mandatory or not.
        mandatory = c(TRUE, FALSE, FALSE, TRUE),
        path_to_md_dir = system.file("md/", package = "DaparToolshed")
    )


    # Define default selected values for widgets
    # This is only for simple workflows
    widgets.default.values <- list(
        Quantimetadatafiltering_tag = "None",
        Quantimetadatafiltering_scope = "None",
        Quantimetadatafiltering_keepRemove = "delete",
        Quantimetadatafiltering_valueTh = 0,
        Quantimetadatafiltering_percentTh = 0,
        Quantimetadatafiltering_valuePercent = 0,
        Quantimetadatafiltering_valPercent = "Value",
        Quantimetadatafiltering_operator = "<=",
        Variablefiltering_cname = "None",
        Variablefiltering_value = NULL,
        Variablefiltering_operator = ""
    )


    rv.custom.default.values <- list(
        funFilter = NULL,
        varFilters = list(),
        varQueries = list(),
        varFilter_DT = data.frame(
            query = "-",
            nbDeleted = "-",
            TotalMainAssay = "-",
            stringsAsFactors = FALSE
        ),
        qMetacell_Filter_SummaryDT = data.frame(
            query = "-",
            nbDeleted = "-",
            TotalMainAssay = "-",
            stringsAsFactors = FALSE
        )
    )

    ### -------------------------------------------------------------###
    ###                                                             ###
    ### ------------------- MODULE SERVER --------------------------###
    ###                                                             ###
    ### -------------------------------------------------------------###
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Insert necessary code which is hosted by MagellanNTK
        # DO NOT MODIFY THIS LINE
        eval(str2expression(Get_Workflow_Core_Code(
            w.names = names(widgets.default.values),
            rv.custom.names = names(rv.custom.default.values)
        )))





        # Add observer to catch the remoteReset so as to
        # manage rv.custom values

        # >>>
        # >>> START ------------- Code for Description UI---------------
        # >>>


        output$Description <- renderUI({
            file <- paste0(config@path_to_md_dir, "/", id, ".md")

            tagList(
                # In this example, the md file is found in the module_examples 
                # directory but with a real app, it should be provided by the 
                # package which contains the UI for the different steps of 
                # the process module. system.file(xxx)

                if (file.exists(file)) {
                    includeMarkdown(file)
                } else {
                    p("No Description available")
                },


                # Used to show some information about the dataset which is 
                # loaded. This function must be provided by the package of the 
                # process module
                uiOutput(ns("datasetDescription_ui")),

                # Insert validation button
                uiOutput(ns("Description_btn_validate_ui"))
            )
        })

        output$datasetDescription_ui <- renderUI({
            # Insert your own code to vizualise some information
            # about your dataset. It will appear once the 'Start' button
            # has been clicked
        })

        output$Description_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Description_btn_validate"),
                "Start",
                class = btn_success_color
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Description"])
        })


        observeEvent(input$Description_btn_validate, {
            rv$dataIn <- dataIn()
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
            rv$steps.status["Description"] <- global$VALIDATED
        })



        #--------------------------------------------------------------
        #
        #             Quantitative metadata filtering UI
        #
        #--------------------------------------------------------------

        output$Quantimetadatafiltering <- renderUI({
            wellPanel(
                # uiOutput for all widgets in this UI
                # This part is mandatory
                # The renderUI() function of each widget is managed by MagellanNTK
                # The dev only have to define a reactive() function for each
                # widget he want to insert
                # Be aware of the naming convention for ids in uiOutput()
                # For more details, please refer to the dev document.
                # mod_build_qMetacell_FunctionFilter_ui(ns('query')),
                DT::dataTableOutput(ns("qMetacell_Filter_DT")),
                uiOutput(ns("Quantimetadatafiltering_buildQuery_ui")),
                uiOutput(ns("example_ui")),
                mod_ds_qMetacell_ui(ns("plots")),
                # Insert validation button
                uiOutput(ns("Quantimetadatafiltering_btn_validate_ui"))
            )
        })



        output$example_ui <- renderUI({
            req(length(rv.custom$funFilter$ll.fun()) > 0)
            req(rv$steps.status["Quantimetadatafiltering"] == 0)

            temp <- filterFeaturesOneSE(
                object = mainAssay(rv$dataIn),
                filters = rv.custom$funFilter$ll.fun()
            )
            mod_filterExample_server(
                id = "filteringExample",
                objBefore = reactive({
                    mainAssay(rv$dataIn)
                }),
                objAfter = reactive({
                    temp
                }),
                query = reactive({
                    rv.custom$funFilter$ll.query()
                })
            )
            widget <- mod_filterExample_ui(ns("filteringExample"))
            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Quantimetadatafiltering"])
        })


        mod_ds_qMetacell_server(
            id = "plots",
            se = reactive({
                mainAssay(rv$dataIn)
            }),
            init.pattern = "missing",
            conds = design.qf(rv$dataIn)$Condition
        )

        output$qMetacell_Filter_DT <- DT::renderDataTable(
            server = TRUE,
            {
                df <- rv.custom$qMetacell_Filter_SummaryDT
                df[, "query"] <- ConvertListToHtml(rv.custom$funFilter$ll.query())
                showDT(df)
            }
        )

        output$Quantimetadatafiltering_buildQuery_ui <- renderUI({
            widget <- mod_build_qMetacell_FunctionFilter_ui(ns("query"))
            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Quantimetadatafiltering"])
        })

        rv.custom$funFilter <- mod_build_qMetacell_FunctionFilter_server(
            id = "query",
            obj = reactive({
                mainAssay(rv$dataIn)
            }),
            conds = reactive({
                design.qf(rv$dataIn)$Condition
            }),
            list_tags = reactive({
                req(rv$dataIn)
                c(
                    "None" = "None",
                    qMetacell.def(typeDataset(mainAssay(rv$dataIn)))$node
                )
            }),
            keep_vs_remove = reactive({
                stats::setNames(nm = c("delete", "keep"))
            }),
            val_vs_percent = reactive({
                stats::setNames(nm = c("Count", "Percentage"))
            }),
            operator = reactive({
                stats::setNames(nm = SymFilteringOperators())
            })
        )


        # Update the list of queries each time a new filter is added
        # observeEvent(rv.custom$funFilter$trigger(), ignoreInit = TRUE, 
        # ignoreNULL = TRUE,{
        #   rv.custom$qMetacell_Filter_SummaryDT$query <- data.frame(
        #     query = rv.custom$funFilter$ll.query(),
        #     nbDeleted = '-',
        #     TotalMainAssay = '-',
        #     stringsAsFactors = FALSE)
        # })

        output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
                "Perform qMetacell filtering",
                class = btn_success_color
            )

            cond <- length(rv.custom$funFilter$ll.fun()) > 0
            cond <- cond && rv$steps.enabled["Quantimetadatafiltering"]
            MagellanNTK::toggleWidget(widget, cond)
        })


        observeEvent(input$Quantimetadatafiltering_btn_validate, {
            rv$dataIn <- filterFeaturesOneSE(
                object = rv$dataIn,
                i = length(rv$dataIn),
                name = "qMetacellFiltered",
                filters = rv.custom$funFilter$ll.fun()
            )

            # Add infos
            nBefore <- nrow(rv$dataIn[[length(rv$dataIn) - 1]])
            nAfter <- nrow(rv$dataIn[[length(rv$dataIn)]])

            rv.custom$qMetacell_Filter_SummaryDT[, "nbDeleted"] <- nBefore - nAfter
            rv.custom$qMetacell_Filter_SummaryDT[, "TotalMainAssay"] <- nrow(mainAssay(rv$dataIn))

            par <- rv.custom$funFilter$ll.widgets.value()
            params(rv$dataIn, length(rv$dataIn)) <- par
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
            rv$steps.status["Quantimetadatafiltering"] <- global$VALIDATED
        })



        # -------------------------------------------------------
        #
        #                  Variable filtering
        #
        # -------------------------------------------------------
        output$Variablefiltering <- renderUI({
            wellPanel(
                DT::dataTableOutput(ns("VarFilter_DT")),
                # Build queries
                uiOutput(ns("Variablefiltering_cname_ui")),
                uiOutput(ns("Variablefiltering_value_ui")),
                uiOutput(ns("Variablefiltering_operator_ui")),
                uiOutput(ns("Variablefiltering_addFilter_btn_ui")),
                # Show example
                uiOutput(ns("Variablefiltering_example_ui")),
                # Process the queries
                uiOutput(ns("Variablefiltering_btn_validate_ui"))
            )
        })


        output$Variablefiltering_example_ui <- renderUI({
            req(length(rv.custom$varFilters) > 0)
            req(rv$steps.status["Variablefiltering"] == 0)

            temp <- filterFeaturesOneSE(
                object = mainAssay(rv$dataIn),
                filters = rv.custom$varFilters
            )

            mod_filterExample_server(
                id = "varFilterExample",
                objBefore = reactive({
                    mainAssay(rv$dataIn)
                }),
                objAfter = reactive({
                    temp
                }),
                query = reactive({
                    rv.custom$varQueries
                })
            )

            widget <- mod_filterExample_ui(ns("varFilterExample"))
            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Variablefiltering"])
        })


        output$VarFilter_DT <- DT::renderDataTable(
            server = TRUE,
            {
                rv.custom$varFilter_DT[, "query"] <- ConvertListToHtml(rv.custom$varQueries)
                showDT(rv.custom$varFilter_DT)
            }
        )

        showDT <- function(df) {
            DT::datatable(df,
                extensions = c("Scroller"),
                escape = FALSE,
                rownames = FALSE,
                options = list(
                    dom = "rt",
                    initComplete = .initComplete(),
                    deferRender = TRUE,
                    bLengthChange = FALSE
                )
            )
        }


        output$Variablefiltering_cname_ui <- renderUI({
            .choices <- c(
                "None",
                colnames(rowData(mainAssay(rv$dataIn)))
            )

            widget <- selectInput(ns("Variablefiltering_cname"),
                "Column name",
                choices = stats::setNames(.choices, nm = .choices),
                width = "300px"
            )

            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Variablefiltering"])
        })


        output$Variablefiltering_operator_ui <- renderUI({
            req(rv.widgets$Variablefiltering_value)
            if (is.na(as.numeric(rv.widgets$Variablefiltering_value))) {
                .operator <- c("==", "!=", "startsWith", "endsWith", "contains")
            } else {
                .operator <- DaparToolshed::SymFilteringOperators()
            }


            widget <- selectInput(ns("Variablefiltering_operator"),
                "operator",
                choices = stats::setNames(nm = .operator),
                width = "100px"
            )
            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Variablefiltering"])
        })

        output$Variablefiltering_value_ui <- renderUI({
            widget <- textInput(ns("Variablefiltering_value"),
                "value",
                width = "100px"
            )
            MagellanNTK::toggleWidget(widget, 
                rv$steps.enabled["Variablefiltering"])
        })


        output$Variablefiltering_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Variablefiltering_btn_validate"),
                "Perform",
                class = btn_success_color
            )
            
            .cond1 <- rv$steps.enabled["Variablefiltering"]
            .cond2 <- length(rv.custom$varFilters) > 0
            MagellanNTK::toggleWidget( widget, .cond1 && .cond2)
        })


        output$Variablefiltering_addFilter_btn_ui <- renderUI({
            widget <- actionButton(
                ns("Variablefiltering_addFilter_btn"),
                "Add filter"
            )
            MagellanNTK::toggleWidget(
                widget,
                rv$steps.enabled["Variablefiltering"]
            )
        })


        observeEvent(input$Variablefiltering_addFilter_btn, {
            type.val <- as.numeric(rv.widgets$Variablefiltering_value)
            if (!is.na(type.val)) {
                value <- type.val
            } else {
                rv.widgets$Variablefiltering_value
            }


            rv.custom$varFilters <- append(
                rv.custom$varFilters,
                VariableFilter(
                    field = rv.widgets$Variablefiltering_cname,
                    value = value,
                    condition = rv.widgets$Variablefiltering_operator
                )
            )
            rv.custom$varQueries <- append(
                rv.custom$varQueries,
                paste0(
                    rv.widgets$Variablefiltering_cname, " ",
                    rv.widgets$Variablefiltering_operator, " ",
                    value
                )
            )
        })



        observeEvent(input$Variablefiltering_btn_validate, {
            rv$dataIn <- filterFeaturesOneSE(
                object = rv$dataIn,
                i = length(rv$dataIn),
                name = "variableFiltered",
                filters = rv.custom$varFilters
            )
            # Add infos
            nBefore <- nrow(rv$dataIn[[length(rv$dataIn) - 1]])
            nAfter <- nrow(mainAssay(rv$dataIn))


            rv.custom$varFilter_DT[, "nbDeleted"] <- nBefore - nAfter
            rv.custom$varFilter_DT[, "TotalMainAssay"] <- nrow(mainAssay(rv$dataIn))

            # Add the parameters values to the new dataset
            params(rv$dataIn[[length(rv$dataIn)]]) <- rv.custom$varQueries

            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
            rv$steps.status["Variablefiltering"] <- global$VALIDATED
        })





        #--------------------------------------------------------
        output$Save <- renderUI({
            wellPanel(
                uiOutput(ns("Save_btn_validate_ui")),
                uiOutput(ns("mod_dl_ui"))
            )
        })


        output$mod_dl_ui <- renderUI({
            req(config@mode == "process")
            req(rv$steps.status["Save"] == global$VALIDATED)
            MagellanNTK::mod_dl_ui(ns("createQuickLink"))
        })


        output$Save_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Save_btn_validate"),
                "Perform",
                class = btn_success_color
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Save"])
        })


        observeEvent(input$Save_btn_validate, {
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
            rv$steps.status["Save"] <- global$VALIDATED
            MagellanNTK::mod_dl_server("createQuickLink", dataIn = reactive({
                rv$dataIn
            }))
        })


        #----------------------------------------------------
        # Insert necessary code which is hosted by MagellanNTK
        # DO NOT MODIFY THIS LINE
        eval(parse(text = MagellanNTK::Module_Return_Func()))
    })
}

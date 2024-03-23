
#' @title Module agregation
#'
#' @description
#' xxxxx
#'
#' @section Step 'Description':
#'
#' xxxxxxx
#'
#' @section Step 'Filter peptides':
#'
#' xxxxxx
#'
#' @section Step 'Agregation':
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
#' @name mod_agregation
#' 
#' @return NA
#'
#' @examples
#' if (interactive()) {
#'     run_workflow("Agregation", verbose = TRUE)
#' }
NULL


#' @param id
#'
#' @rdname mod_agregation
#'
#' @export
#'
mod_Agregation_ui <- function(id) {
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
#'
#'
#'
#' @rdname mod_agregation
#' @importFrom shinyjs toggle hidden
#' @importFrom utils write.table
#' @importFrom MagellanNTK Timestamp toggleWidget Get_Workflow_Core_Code
#'
#' @export
#'
mod_Agregation_server <- function(id,
    dataIn = reactive({NULL}),
    steps.enabled = reactive({NULL}),
    remoteReset = reactive({FALSE}),
    steps.status = reactive({NULL}),
    current.pos = reactive({1}),
    verbose = FALSE) {

  pkgs.require("DT")
    
    
    # This list contains the basic configuration of the process
    config <- list(
        # Define the type of module
        mode = "process",

        # List of all steps of the process
        steps = c("Description", "Filter peptides", "Agregation", "Save"),

        # A vector of boolean indicating if the steps are mandatory or not.
        mandatory = c(TRUE, FALSE, TRUE, TRUE),
        path_to_md_dir = system.file("md/", package = "DaparToolshed")
    )


    # Define default selected values for widgets
    # This is only for simple workflows
    widgets.default.values <- list(
        Filterpeptides_ProteinId = NULL,
        Filterpeptides_useOfShared = NULL,
        Filterpeptides_consider = NULL,
        Filterpeptides_nTopn = NULL,
        Filterpeptides_barplotType = "all",
        Agregation_method = NULL,
        Agregation_filterProtAfterAgregation = FALSE,
        nbPeptides = 0,
        AddMetadata_columnsForProteinDatasetBox = NULL
    )


    rv.custom.default.values <- list(
        temp.agregate = NULL
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
        # Filter peptides step UI
        #--------------------------------------------------------------

        output$Filterpeptides <- renderUI({
            wellPanel(
                # uiOutput for all widgets in this UI
                # This part is mandatory
                # The renderUI() function of each widget is managed by MagellanNTK
                # The dev only have to define a reactive() function for each
                # widget he want to insert
                # Be aware of the naming convention for ids in uiOutput()
                # For more details, please refer to the dev document.
                uiOutput(ns("warningAgregationMethod")),

                # uiOutput(ns("Filterpeptides_ProteinId_ui")),
                uiOutput(ns("Filterpeptides_consider_ui")),
                uiOutput(ns("Filterpeptides_nTopn_ui")),

                # Insert validation button
                uiOutput(ns("Filterpeptides_btn_validate_ui")),

                # Insert additional code

                uiOutput(ns("Filterpeptides_barplotType_ui")),
                highchartOutput(ns("peptideBarplot"), width = "400px")
            )
        })

        output$Filterpeptides_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Filterpeptides_btn_validate"),
                "Run",
                class = btn_success_color
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Filterpeptides"])
        })


        observeEvent(input$Filterpeptides_btn_validate, {
            f <- switch(rv.widgets$Filterpeptides_consider,
                allPeptides = FunctionFilter("allPeptides", list()),
                topnPeptides = FunctionFilter("topnPeptides", 
                    fun = "rowSums", top = 2),
                specPeptides = FunctionFilter("specPeptides", list())
            )

            ll.filters <- list()
            rv$dataIn <- filterFeaturesOneSE(
                object = rv$dataIn,
                i = length(rv$dataIn),
                name = "adjMatfiltered",
                filters = list(f)
            )
            # DO NOT MODIFY THE THREE FOLLOWINF LINES
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
            rv$steps.status["Filterpeptides"] <- global$VALIDATED
        })

        output$Filterpeptides_ProteinId_ui <- renderUI({
            # req (is.null(GetProteinId(rv$dataIn))
            #
            # widget <- selectInput(ns("Agregation_ProteinId"),
            # "Choose the protein ID",
            # choices = c("None", colnames(Biobase::fData(rv$dataIn))),
            # selected = rv.widgets$Agregation_proteinId)
            # toggleWidget(rv$steps.enabled['Agregation'], widget )
        })



        output$Filterpeptides_consider_ui <- renderUI({
            widget <- radioButtons(ns("Filterpeptides_consider"), "Consider",
                choices = DaparToolshed::AdjMatFilters(),
                selected = rv.widgets$Filterpeptides_consider
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Filterpeptides"])
        })




        output$Filterpeptides_nTopn_ui <- renderUI({
            req(rv.widgets$Filterpeptides_consider == "topnPeptides")

            widget <- numericInput(ns("Filterpeptides_nTopn"), "N",
                value = rv.widgets$Filterpeptides_nTopn,
                min = 0,
                step = 1,
                width = "100px"
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Filterpeptides"])
        })






        output$Filterpeptides_barplotType_ui <- renderUI({
            widget <- selectInput(ns("Filterpeptides_barplotType"), 
                "Barplot type",
                choices = c(
                    "All peptides" = "all"
                    # "Only specific peptides" = "onlySpec",
                    # "Only shared peptides" = "onlyShared"
                ),
                selected = rv.widgets$Filterpeptides_barplotType,
                width = "200"
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Filterpeptides"])
        })

        output$peptideBarplot <- renderHighchart({
            req(adjacencyMatrix(last_assay(rv$dataIn)))
            withProgress(message = "Rendering plot, pleast wait...", 
                detail = "", value = 1, {
                X <- adjacencyMatrix(last_assay(rv$dataIn))
                if (is.null(X)) {
                    .last <- last_assay(rv$dataIn)
                    X <- makeAdjacencyMatrix(rowData(.last)[, "Protein_group_IDs"])
                    rownames(X) <- rownames(last_assay(rv$dataIn))
                }

                # Modify the adjacencymatrix to reflect the filter parameters
                # X <- updateAdjacencyMatrix(X, 
                # mode = rv.widgets$Filterpeptides_barplotType)

                GraphPepProt_hc(X)
            })
        })


        #--------------------------------------------------------------
        # Agregation step UI
        #--------------------------------------------------------------

        output$Agregation <- renderUI({
            wellPanel(
                # uiOutput for all widgets in this UI
                # This part is mandatory
                # The renderUI() function of each widget is managed by MagellanNTK
                # The dev only have to define a reactive() function for each
                # widget he want to insert
                # Be aware of the naming convention for ids in uiOutput()
                # For more details, please refer to the dev document.
                uiOutput(ns("Agregation_method_ui")),

                # Insert validation button
                uiOutput(ns("Agregation_btn_validate_ui")),

                # Insert additional code
                uiOutput(ns("ObserverAggregationDone")),
                shinyjs::hidden(downloadButton(ns("downloadAggregationIssues"),
                    "Download issues",
                    class = actionBtnClass
                )),
                hr(),
                DT::DTOutput(ns("aggregationStats"))
            )
        })


        output$Agregation_method_ui <- renderUI({
            widget <- radioButtons(ns("Agregation_method"), "method",
                choices = aggregateMethods(),
                selected = rv.widgets$Filterpeptides_method
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Agregation"])
        })

        # output$Agregation_useOfShared_ui <- renderUI({
        #   mod_popover_for_help_server("modulePopover_useOfShared",
        #   data = list(title="Include shared peptides",
        #   content= HTML(
        #   paste0("<ul><li><strong>No:</strong> 
        #   only protein-specific peptides</li><li><strong>Yes 1:
        #   </strong> shared peptides processed as protein specific
        #   </li><li><strong>Yes 2</strong>: proportional redistribution 
        #   of shared peptides</li></ul>")
        #   )
        #   )
        #   )
        #
        #   widget <- radioButtons(ns("Agregation_useOfShared"),
        #   NULL,
        #   choices = c("No" = "No",
        #   "As protein specific"= "asSpec" ,
        #   "For redistribution" = "forDistrib" ),
        #   selected = rv.widgets$Agregation_useOfShared)
        #
        #   tagList(
        #     mod_popover_for_help_ui(ns("modulePopover_useOfShared")),
        #     toggleWidget(rv$steps.enabled['Agregation'], widget )
        #   )
        #
        # })


        output$warningAgregationMethod <- renderUI({
            req(rv$dataIn)
            m <- match.qMetacell(qMetacell(rv$dataIn, length(rv$dataIn)),
                pattern = "missing",
                level = "peptide"
            )

            if (length(which(m)) > 0) {
                tags$p(
                    style = "color: red;",
                    tags$b("Warning:"), " Your dataset contains missing values.
    For better results, you should impute them first"
                )
            }
        })










        output$Agregation_nbPeptides_ui <- renderUI({
            req(isTRUE(rv.widgets$Agregation_filterProtAfterAgregation))

            widget <- numericInput(ns("Agregation_nbPeptides"),
                "Nb of peptides defining a protein",
                value = 0, min = 0, step = 1,
                width = "250px"
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Agregation"])
        })



        output$displayNbPeptides <- renderUI({
            req(rv$widgets$aggregation$filterProtAfterAgregation)

            if (rv$widgets$aggregation$filterProtAfterAgregation) {
                numericInput("nbPeptides", "Nb of peptides defining a protein",
                    value = 0, min = 0, step = 1,
                    width = "250px"
                )
            }
        })

        output$Agregation_btn_validate_ui <- renderUI({
            widget <- actionButton(ns("Agregation_btn_validate"),
                "Perform",
                class = btn_success_color
            )
            MagellanNTK::toggleWidget(widget, rv$steps.enabled["Agregation"])
        })


        observeEvent(input$Agregation_btn_validate, {
            rv.custom$temp.agregate <- aggregateFeatures4Prostar(
                object = rv$dataIn,
                i = length(rv$dataIn),
                name = "aggregated",
                fcol = "adjacencyMatrix",
                fun = rv.widgets$Agregation_method
            )


            if (is.null(rv.custom$temp.agregate$issues)) {
                rv$dataIn <- rv.custom$temp.agregate
                dataOut$trigger <- MagellanNTK::Timestamp()
                dataOut$value <- rv$dataIn
                rv$steps.status["Agregation"] <- global$VALIDATED
            } else {
                shinyjs::toggle("downloadAggregationIssues",
                    condition = length(rv.custom$temp.agregate$issues) > 0
                )
            }
        })



        output$downloadAggregationIssues <- downloadHandler(
            filename = "aggregation_issues.txt",
            content = function(file) {
                tmp.peptides <- lapply(rv.custom$temp.agregate$issues, 
                    function(x) paste0(x, collapse = ","))
                df <- data.frame(
                    Proteins = names(rv.custom$temp.agregate$issues),
                    Peptides = as.data.frame(do.call(rbind, tmp.peptides))
                )
                colnames(df) <- c("Proteins", "Peptides")
                write.table(df, 
                    file = file, 
                    row.names = FALSE, 
                    quote = FALSE, 
                    sep = "\t")
            }
        )





        ## --------------------------------------------------------
        ##
        ##
        ## ------------------------------------------------------
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
        })


        #----------------------------------------------------
        # Insert necessary code which is hosted by MagellanNTK
        # DO NOT MODIFY THIS LINE
        eval(parse(text = MagellanNTK::Module_Return_Func()))
    })
}

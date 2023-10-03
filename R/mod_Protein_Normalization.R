

#' @export
#'
#' @rdname pipeline-protein
#'
mod_Protein_Normalization_ui <- function(id) {
    ns <- NS(id)
}


#' @title xxx
#'
#' @description
#' This module contains the configuration informations for the corresponding 
#' pipeline. It is called by the nav_pipeline module of the package MagellanNTK
#'
#' @param id xxx
#'
#' @param dataIn The dataset
#'
#' @param steps.enabled A vector of boolean which has the same length of the 
#' steps of the pipeline. This information is used to enable/disable the 
#' widgets. It is not a communication variable between the caller and this 
#' module, thus there is no corresponding output variable
#'
#' @param remoteReset It is a remote command to reset the module. A boolean that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#'
#' @author Samuel Wieczorek
#'
#' @importFrom shinyjs disabled
#' @importFrom MagellanNTK Timestamp GlobalSettings
#'
#' @export
#'
#' @rdname pipeline-protein

mod_Protein_Normalization_server <- function(id,
    dataIn = reactive({NULL}),
    steps.enabled = reactive({NULL}),
    remoteReset = reactive({FALSE})) {

    #' @field config xxxx
    config <- list(
        name = "Normalization",
        parent = "Protein",
        steps = c("Description", "Normalization", "Save"),
        mandatory = c(TRUE, FALSE, TRUE)
    )

    # Define default selected values for widgets
    widgets.default.values <- list(
        norm_method = "None",
        norm_type = "overall",
        norm_spanLOESS = 0.7,
        norm_quantile = 0.15,
        norm_varReduction = FALSE,
        norm_sync = FALSE
    )


    ### -------------------------------------------------------------###
    ###                                                             ###
    ### ------------------- MODULE SERVER --------------------------###
    ###                                                             ###
    ### -------------------------------------------------------------###
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        rv.widgets <- reactiveValues(
            norm_method = widgets.default.values$norm_method,
            norm_type = widgets.default.values$norm_type,
            norm_spanLOESS = widgets.default.values$norm_spanLOESS,
            norm_quantile = widgets.default.values$norm_quantile,
            norm_varReduction = widgets.default.values$norm_varReduction,
            norm_sync = widgets.default.values$norm_sync
        )


        # Reactive values during the run of the process
        rv <- reactiveValues(
            dataIn = NULL,
            steps.status = NULL,
            reset = NULL,
            steps.enabled = NULL
        )

        # Returned value of the process
        dataOut <- reactiveValues(
            trigger = NULL,
            value = NULL
        )

        rv.norm <- reactiveValues(
            trackFromBoxplot = NULL,
            selectProt = NULL,
            resetTracking = FALSE
        )


        # Initialization of the module
        observeEvent(steps.enabled(), ignoreNULL = TRUE, {
            if (is.null(steps.enabled())) {
                rv$steps.enabled <- stats::setNames(rep(FALSE, length(config@steps)),
                    nm = config@steps
                )
            } else {
                rv$steps.enabled <- steps.enabled()
            }
        })

        # Set all the widgets to their default value after the remote Reset()
        observeEvent(remoteReset(), {
            lapply(names(rv.widgets), function(x) {
                rv.widgets[[x]] <- widgets.default.values[[x]]
            })

            rv.norm$resetTracking <- TRUE
        })



        # Call the modules server functions

        mod_plots_density_server(
            id = "densityPlot_Norm",
            obj = reactive({obj}),
            conds = reactive({conds}),
            legend = reactive({legend}),
            base_palette = reactive({Example_Palette()})
        )

        mod_popover_for_help_server("modulePopover_normQuanti",
            data = list(
                title = "Normalization quantile",
                content = "lower limit/noise (quantile = 0.15), 
                median (quantile = 0.5). Min value=0, max value=1"
            )
        )

        rv.norm$selectProt <- mod_plots_tracking_server(
            "master_tracking",
            obj = reactive({rv$current.obj}),
            params = reactive({NULL}),
            keyId = reactive({rv$dataIn@experimentData@other$proteinId}),
            reset = reactive({rv.norm$resetTracking}),
            slave = reactive({FALSE})
        )

        params <- list(
            typeSelect = "ProteinList",
            randSelect = "",
            colSelect = NULL,
            rand.indices = "",
            col.indices = NULL,
            listSelect = NULL,
            list.indices = ""
        )

        rv.norm$selectProt <- mod_plots_tracking_server("master_tracking",
            obj = reactive({rv$dataIn[[length(rv$dataIn)]]}),
            params = reactive({params}),
            keyId = reactive({keyId}),
            reset = reactive({FALSE}),
            slave = reactive({FALSE})
            )


        rv.norm$trackFromBoxplot <- callModule(mod_plots_intensity_server,
            "boxPlot_Norm",
            dataIn = reactive({rv$current.obj}),
            meta = reactive({fData(rv$current.obj)}),
            keyId = reactive({rv$current.obj@experimentData@other$proteinId}),
            conds = reactive({pData(rv$current.obj)$Condition}),
            pal = reactive({rv$PlotParams$paletteForConditions}),
            params = reactive({
                if (rv.widgets$norm_sync) {rv.norm$selectProt()} 
                else {NULL}
            }),
            reset = reactive({rv.norm$resetTracking}),
            slave = reactive({rv.widgets$norm_sync})
        )


        callModule(module_Not_a_numeric, "test_normQuant", reactive({
            rv$widgets$normalization$quantile
        }))


        ###### --------- Code for Description (step 0) -------#####
        output$Description <- renderUI({
            tagList(
                includeMarkdown(
                    paste0("md/", 
                        paste0(config@parent, "_", config@name, ".md"))),
                uiOutput(ns("datasetDescription")),
                uiOutput(ns("validationBtn_ui"))
            )
        })


        observeEvent(input$btn_validate_Description, 
            ignoreInit = TRUE, ignoreNULL = TRUE, {
            rv$dataIn <- dataIn()
            # rv$steps.status['Description'] <- global$VALIDATED
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
        })


        ###### ------------------- Code for step 1 ---------------#####


        # ObserveEvent of the widgets
        observeEvent(input$norm_method, {
            rv.widgets$norm_method <- input$norm_method
        })
        observeEvent(input$norm_type, {
            rv.widgets$norm_type <- input$norm_type
        })
        observeEvent(input$norm_spanLOESS, {
            rv.widgets$norm_spanLOESS <- input$norm_spanLOESS
        })
        observeEvent(input$norm_quantile, {
            rv.widgets$norm_quantile <- input$norm_quantile
        })
        observeEvent(input$norm_varReduction, {
            rv.widgets$norm_varReduction <- input$norm_varReduction
        })
        observeEvent(input$norm_sync, {
            rv.widgets$norm_sync <- input$norm_sync
        })




        output$norm_method_ui <- renderUI({
            # rv$steps.enabled
            rv.widgets$norm_method
            if (rv$steps.enabled["Normalization"]) {
                selectInput("norm_method", "Normalization method",
                    choices = normMethods,
                    selected = rv.widgets$norm_method,
                    width = "200px"
                )
            } else {
                shinyjs::disabled(
                    selectInput("norm_method", "Normalization method",
                        choices = normMethods,
                        selected = rv.widgets$norm_method,
                        width = "200px"
                    )
                )
            }
        })



        output$norm_type_ui <- renderUI({
            rv.widgets$norm_type
            if (rv$steps.enabled["Normalization"]) {
                selectInput("norm_type", "Normalization type",
                    choices = stats::setNames(
                        nm = c("overall", "within conditions")),
                    selected = rv.widgets$norm_type,
                    width = "150px"
                )
            } else {
                shinyjs::disabled(
                    selectInput("norm_type", "Normalization type",
                        choices = stats::setNames(nm = c("overall", 
                            "within conditions")),
                        selected = rv.widgets$norm_type,
                        width = "150px"
                    )
                )
            }
        })

        output$norm_spanLOESS_ui <- renderUI({
            rv.widgets$norm_type
            if (rv$steps.enabled["Normalization"]) {
                textInput("norm_spanLOESS", "Span",
                    value = rv.widgets$norm_spanLOESS,
                    width = "100px"
                )
            } else {
                shinyjs::disabled(
                    textInput("norm_spanLOESS", "Span",
                        value = rv.widgets$norm_spanLOESS,
                        width = "100px"
                    )
                )
            }
        })


        output$norm_quantile_ui <- renderUI({
            req(rv.widgets$norm_method == "QuantileCentering")
            if (rv$steps.enabled["Normalization"]) {
                tagList(
                    mod_popover_for_help_ui("modulePopover_normQuanti"),
                    textInput(ns("normalization.quantile"), NULL,
                        value = rv.widgets$norm_quantile,
                        width = "150px"
                    ),
                    module_Not_a_numericUI(ns("test_normQuant"))
                )
            } else {
                shinyjs::disabled(
                    tagList(
                        mod_popover_for_help_ui(ns("modulePopover_normQuanti")),
                        textInput(ns("normalization.quantile"), NULL,
                            value = rv.widgets$norm_quantile,
                            width = "150px"
                        ),
                        module_Not_a_numericUI(ns("test_normQuant"))
                    )
                )
            }
        })



        output$norm_scaling_ui <- renderUI({
            req(rv.widgets$norm_method)

            if (rv.widgets$norm_method == "MeanCentering") {
                # check if the normalisation has already been performed
                if (rv$steps.enabled["Normalization"]) {
                    checkboxInput(ns("norm_varReduction"), 
                        "Include variance reduction",
                        value = rv.widgets$norm_varReduction
                    )
                } else {
                    shinyjs::disabled(
                        checkboxInput(ns("norm_varReduction"), 
                            "Include variance reduction",
                            value = rv.widgets$norm_varReduction
                        )
                    )
                }
            }
        })




        outp$norm_masterTracking_ui <- renderUI({
            if (rv$steps.enabled["Normalization"]) {
                tagList(
                    mod_plots_tracking_ui(ns("mod_master_tracking")),
                    checkboxInput(ns("norm_sync"),
                        "Synchronise with selection above",
                        value = rv.widgets$norm_sync
                    )
                )
            } else {
                tagList(
                    mod_plots_tracking_ui(ns("mod_master_tracking")),
                    checkboxInput(ns("norm_sync"),
                        "Synchronise with selection above",
                        value = rv.widgets$norm_sync
                    )
                )
            }
        })

        # Buttons must be explicitly enabled/disabled with a full code
        # Otherwise, they do not disable
        output$norm_perform_btn_ui <- renderUI({
            if (rv$steps.enabled["Normalization"]) {
                actionButton(ns("perform.normalization"),
                    "Perform normalization",
                    class = MagellanNTK::GlobalSettings$actionBtnClass,
                    width = "170px"
                )
            } else {
                shinyjs::disabled(
                    actionButton(ns("perform.normalization"),
                        "Perform normalization",
                        class = GlobalSettings$actionBtnClass,
                        width = "170px"
                    )
                )
            }
        })



        output$helpForNormalizationMethods <- renderUI({
            req(rv.widgets$norm_method != "None")


            switch(rv.widgets$norm_method,
                GlobalQuantileAlignment = txt <- "This method proposes a 
                normalization of important magnitude that should be cautiously 
                used. It proposes to align the quantiles of all the replicates 
                as described in [Other ref. 1]; practically it amounts to 
                replace abundances by order statistics.",
                QuantileCentering = txt <- "These methods propose to shift the 
                sample distributions (either all of them at once, or within 
                each condition at a time) to align a specific quantile: the 
                median (under the assumption that up-regulations and 
                down-regulations are equally frequent), the 15% quantile 
                (under the assumption that the signal/noise ratio is roughly 
                the same in all the samples), or any other user's choice.",
                MeanCentering = txt <- "These methods propose to shift the 
                sample distributions (either all of them at once, or within 
                each condition at a time) to align their means. It is also 
                possible to force unit variance (or not).",
                SumByColumns = txt <- "These methods propose normalizations of 
                important magnitude that should be cautiously used. It operates
                on the original scale (not the log2 one) and propose to 
                normalize each abundance by the total abundance of the sample 
                (so as to focus on the analyte proportions among each sample).",
                LOESS = txt <- "This method proposes to apply a cyclic LOESS 
                [Other ref. 4, 5] normalization to the data (either all of them 
                at once, or on each condition independently). It relates to a
                combination of multiple regression models. The user can tune 
                the regression span (an higher span smooths the fit more, while 
                a lower span captures more trends).",
                vsn = txt <- "This method proposes to apply the Variance 
                Stabilization Normalization method [Other ref. 6] to the data 
                (either all of them at once, or on each condition 
                independently). No specific parameters required."
                )

            tags$p(txt)
        })


        observeEvent(rv.widgets$norm_method, {

            # shinyjs::toggle("perform.normalization", 
            # condition = rv.widgets$norm_method != "None")
            # shinyjs::toggle("spanLOESS", 
            # condition = rv.widgets$norm_method == "LOESS")
            #
            # shinyjs::toggle("normalization.type",
            # condition=( rv.widgets$norm_method %in% c("QuantileCentering", 
            # "MeanCentering", "SumByColumns", "LOESS", "vsn")))
            #
            # cond <- rv$current.obj@experimentData@other$typeOfData == 'protein'
            # trackAvailable <- rv.widgets$norm_method %in% normalizeMethods.dapar(withTracking=TRUE)
            # shinyjs::toggle('DivMasterProtSelection', condition= cond && trackAvailable)
        })



        output$viewComparisonNorm_hc <- renderHighchart({
            req(rv$dataIn)

            # ind <- grep("Normalized.", names(rv$dataset))
            # if (length(ind) > 0){
            #   obj1 <- rv$dataset[[ind - 1]]
            #   obj2 <- rv$dataset[[ind]]
            # }
            # else {
            #   obj1 <- rv$dataset[[input$datasets]]
            #   obj2 <- rv$current.obj
            # }
            #
            # if (is.null(obj1) || is.null(obj2))
            #   return(NULL)
            #
            # compareNormalizationD_HC(qDataBefore = Biobase::exprs(obj1),
            # qDataAfter = Biobase::exprs(obj2),
            # keyId = fData(rv$current.obj)[,rv$current.obj@experimentData@other$proteinId],
            # conds = Biobase::pData(obj1)$Condition,
            # pal = rv$PlotParams$paletteForConditions,
            # subset.view =   if(rv.norm$sync)
            # GetIndicesOfSelectedProteins_ForNorm()
            # else
            # GetIndicesOfSelectedProteins()
            # )
        })



        # GetIndicesOfSelectedProteins_ForNorm <- reactive({
        #   req(rv.norm$selectProt())
        #
        #   ind <- NULL
        #   ll <- fData(rv$current.obj)[ ,rv$current.obj@experimentData@other$proteinId]
        #   tt <- rv.norm$selectProt()$type
        #   switch(tt,
        #          ProteinList = ind <- rv.norm$selectProt()$list.indices,
        #          Random = ind <- rv.norm$selectProt()$rand.indices,
        #          Column = ind <- rv.norm$selectProt()$col.indices
        #   )
        #   if (length(ind)==0)
        #     ind <- NULL
        #   ind
        # })

        # GetIndicesOfSelectedProteins <- reactive({
        #   req(rv.norm$trackFromBoxplot())
        #
        #
        #   #print('in GetIndicesOfSelectedProteins')
        #   #print(rv.norm$trackFromBoxplot())
        #   ind <- NULL
        #   ll <- fData(rv$current.obj)[,rv$current.obj@experimentData@other$proteinId]
        #   tt <- rv.norm$trackFromBoxplot()$type
        #   switch(tt,
        #          ProteinList = ind <- rv.norm$trackFromBoxplot()$list.indices,
        #          Random = ind <- rv.norm$trackFromBoxplot()$rand.indices,
        #          Column = ind <- rv.norm$trackFromBoxplot()$col.indices
        #   )
        #   if (length(ind)==0)
        #     ind <- NULL
        #
        #   #print('ind = ')
        #   #print(ind)
        #   ind
        # })




        # ------------------------ STEP 1 : UI ------------------
        output$Normalization <- renderUI({
            name <- "Normalization"
            wellPanel(
                id = ns("toto"),
                tagList(
                    div(
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 20px;",
                            uiOutput(ns("norm_method_ui"))
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 20px;",
                            hidden(uiOutput(ns("norm_type_ui")))
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 20px;",
                            hidden(uiOutput(ns("norm_spanLOESS_ui"))),
                            module_Not_a_numericUI(ns("test_spanLOESS")),
                            uiOutput(ns("norm_quantile_ui")),
                            uiOutput(ns("norm_scaling_ui"))
                        ),
                        hidden(
                            div(
                                id = "DivMasterProtSelection",
                                style = "display:inline-block; 
                                vertical-align: middle; padding-right: 20px;",
                                uiOutput(ns("norm_masterTracking_ui"))
                            )
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 20px;",
                            hidden(uiOutput(ns("norm_perform_btn_ui")))
                        )
                    ),
                    uiOutput(ns("helpForNormalizationMethods")),
                    tags$hr(),
                    fluidRow(
                        column(width = 4, 
                            moduleDensityplotUI(ns("densityPlot_Norm"))),
                        column(
                            width = 4,
                            withProgress(message = "Building plot", 
                                detail = "", value = 0, {
                                mod_plots_intensity_ui(ns("boxPlot_Norm"))
                            })
                        ),
                        column(
                            width = 4,
                            withProgress(message = "Building plot", 
                                detail = "", value = 0, {
                                highchartOutput(ns("viewComparisonNorm_hc"))
                            })
                        )
                    )
                )
            )
        })






        observeEvent(input$btn_validate_Step1, ignoreInit = TRUE, {
            # Add your stuff code here
            # rv$widgets$normalization$method
            # rv$dataset[[input$datasets]]
            # # isolate({
            #
            # switch(rv$widgets$normalization$method,
            #        G_noneStr = rv$current.obj <- rv$dataset[[input$datasets]],
            #        GlobalQuantileAlignment = {
            #          rv$current.obj <- wrapper.normalizeD(rv$dataset[[input$datasets]], rv$widgets$normalization$method)
            #        },
            #        QuantileCentering = {
            #          quant <-NA
            #          if (!is.null(rv$widgets$normalization$quantile))
            #          {quant <- as.numeric(rv$widgets$normalization$quantile)}
            #
            #          rv$current.obj <- wrapper.normalizeD(obj = rv$dataset[[input$datasets]],
            #                                               method = rv$widgets$normalization$method,
            #                                               type = rv$widgets$normalization$type,
            #                                               cond = pData(rv$dataset[[input$datasets]])$Condition,
            #                                               quantile = quant,
            #                                               subset.norm = GetIndicesOfSelectedProteins_ForNorm())
            #
            #        } ,
            #        MeanCentering = {
            #          rv$current.obj <- wrapper.normalizeD(obj = rv$dataset[[input$datasets]],
            #                                               method = rv$widgets$normalization$method,
            #                                               conds = pData(rv$dataset[[input$datasets]])$Condition,
            #                                               type = rv$widgets$normalization$type,
            #                                               scaling=  rv$widgets$normalization$varReduction,
            #                                               subset.norm = GetIndicesOfSelectedProteins_ForNorm())
            #        },
            #        SumByColumns = {
            #          rv$current.obj <- wrapper.normalizeD(obj = rv$dataset[[input$datasets]],
            #                                               method = rv$widgets$normalization$method,
            #                                               conds = pData(rv$dataset[[input$datasets]])$Condition,
            #                                               type = rv$widgets$normalization$type,
            #                                               subset.norm = GetIndicesOfSelectedProteins_ForNorm())
            #
            #        },
            #        LOESS = { rv$current.obj <- wrapper.normalizeD(obj = rv$dataset[[input$datasets]],
            #                                                       method = rv$widgets$normalization$method,
            #                                                       conds = pData(rv$dataset[[input$datasets]])$Condition,
            #                                                       type = rv$widgets$normalization$type,
            #                                                       span = as.numeric(rv$widgets$normalization$spanLOESS))
            #        },
            #        vsn = {
            #          rv$current.obj <- wrapper.normalizeD(obj = rv$dataset[[input$datasets]],
            #                                               method = rv$widgets$normalization$method,
            #                                               conds = pData(rv$dataset[[input$datasets]])$Condition,
            #                                               type = rv$widgets$normalization$type)
            #        }
            # )
            # # })
            # rvModProcess$moduleNormalizationDone[1] <- TRUE
            # shinyjs::toggle("valid.normalization", condition=input$perform.normalization >= 1)

            # DO NOT MODIFY THE 2 FOLLOWING LINES
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
        })

        #-------------------------- Code for step 2 ------------------




        output$Save <- renderUI({
            rv$steps.enabled
            name <- "Save"
            wellPanel(
                if (rv$steps.enabled["Save"]) {
                    actionButton(ns(paste0("valid.normalization", name)),
                        "Save normalization",
                        class = MagellanNTK::GlobalSettings$btn_success_color,
                        width = "170px"
                    )
                } else {
                    shinyjs::disabled(
                        actionButton(ns(paste0("valid.normalization", name)),
                            "Save normalization",
                            class = MagellanNTK::GlobalSettings$btn_success_color,
                            width = "170px"
                        )
                    )
                }
            )
        })

        observeEvent(input$valid.normalization, ignoreInit = TRUE, {
            # Add your stuff code here
            # req(input$perform.normalization)
            # req(rv$current.obj)
            #
            # isolate({
            #   if (rv$widgets$normalization$method != G_noneStr) {
            #     rv$typeOfDataset <- rv$current.obj@experimentData@other$typeOfData
            #     name <- paste0("Normalized", ".", rv$typeOfDataset)
            #     rv$current.obj <- saveParameters(rv$current.obj,
            #                                      name,
            #                                      "Normalization",
            #                                      build_ParamsList_Normalization()
            #     )
            #
            #     rvModProcess$moduleNormalizationDone[2] <- TRUE
            #     UpdateDatasetWidget(rv$current.obj, name)
            #
            #   }
            #
            # } )
            #


            # DO NOT MODIFY THE FOLLOWING LINES
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn
        })


        ## -------------- END FOR UI DEFINITIONS  -------------------



        # Return value of module
        # DO NOT MODIFY THIS PART
        list(
            config = reactive({
                config@ll.UI <- stats::setNames(
                    lapply(
                        config@steps,
                        function(x) {
                            do.call("uiOutput", list(ns(x)))
                        }
                    ),
                    paste0("screen_", config@steps)
                )
                config
            }),
            dataOut = reactive({
                dataOut
            })
            # steps.status = reactive({rv$steps.status})
        )
    })
}

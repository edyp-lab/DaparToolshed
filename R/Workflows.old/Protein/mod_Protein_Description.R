
btn_style <- "display:inline-block; vertical-align: middle; padding: 7px"

mod_Protein_Description_ui <- function(id) {
    ns <- NS(id)
}



#' @title xxx
#'
#' @description
#' This module contains the configuration informations for the corresponding pipeline.
#' It is called by the nav_pipeline module of the package MagellanNTK
#'
#' @param id xxx
#'
#' @param dataIn The dataset
#'
#' @param steps.enabled A vector of boolen which has the same length of the steps
#' of the pipeline. This information is used to enable/disable the widgets. It is not
#' a communication variable between the caller and this module, thus there is no
#' corresponding output variable
#'
#' @param remoteReset It is a remote command to reset the module. A boolen that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#'
#' @author Samuel Wieczorek
#'
#' @export
#' @importFrom shinyjs disabled
#'
mod_Protein_Description_server <- function(id,
                                           dataIn = NULL,
                                           steps.enabled = reactive({
                                               NULL
                                           }),
                                           remoteReset = reactive({
                                               FALSE
                                           })) {

    #' @field global xxxx
    global <- list(
        VALIDATED = 1,
        UNDONE = 0,
        SKIPPED = -1
    )

    config <- list(
        name = "Protein_Description",
        steps = c("Description"),
        mandatory = c(TRUE)
    )

    # Define default selected values for widgets
    widgets.default.values <- list()

    ### -------------------------------------------------------------###
    ###                                                             ###
    ### ------------------- MODULE SERVER --------------------------###
    ###                                                             ###
    ### -------------------------------------------------------------###
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        rv.widgets <- reactiveValues()

        rv <- reactiveValues(
            dataIn = NULL,
            dataOut = NULL,
            status = NULL,
            remoteReset = NULL,
            steps.enabled = NULL
        )

        dataOut <- reactiveValues(
            trigger = NULL,
            value = NULL
        )





        # Initialization of the module
        observeEvent(steps.enabled(), ignoreNULL = TRUE, {
            if (is.null(steps.enabled())) {
                rv$steps.enabled <- stats::setNames(
                    rep(FALSE, length(config@steps)), 
                    config@steps)
            } else {
                rv$steps.enabled <- steps.enabled()
            }
        })


        observeEvent(remoteReset(), {
            lapply(names(rv.widgets), function(x) {
                rv.widgets[[x]] <- widgets.default.values[[x]]
            })
        })

        ### ----------------------------------------------------------------------------------------------------

        ###### ------------------- Code for Description (step 0) -------------------------    #####
        output$Description <- renderUI({
            rv$steps.enabled
            tagList(
                includeMarkdown(system.file("app/md", paste0(config@name, ".md"), package = "MagellanNTK")),
                uiOutput(ns("datasetDescription")),
                if (isTRUE(rv$steps.enabled["Description"])) {
                    actionButton(ns("btn_validate_Description"),
                        paste0("Start ", config@name),
                        class = "btn-success"
                    )
                } else {
                    shinyjs::disabled(
                        actionButton(ns("btn_validate_Description"),
                            paste0("Start ", config@name),
                            class = "btn-success"
                        )
                    )
                }
            )
        })


        Timestamp <- function() {
            if (verbose) cat(paste0("::Timestamp() from - ", id, "\n\n"))
            as.numeric(Sys.time())
        }

        observeEvent(input$btn_validate_Description, ignoreInit = TRUE, ignoreNULL = TRUE, {
            rv$dataIn <- dataIn()
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn

            # rv$status['Description'] <- stepStatus$VALIDATED
        })

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
            # status = reactive({rv$status})
        )
    })
}

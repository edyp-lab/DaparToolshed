btn_style <- "display:inline-block; vertical-align: middle; padding: 7px"

# source(file.path('.', 'mod_Protein_Description.R'), local=TRUE)$value
# source(file.path('.', 'mod_Protein_Normalization.R'), local=TRUE)$value


#' @title Module of Protein pipeline
#'
#' @description
#' This module contains the configuration information for the corresponding 
#' pipeline. It is called by the nav_pipeline module of the package MagellanNTK
#'
#' @name pipeline-protein
#' @author Samuel Wieczorek
#' 
#' @param id xxx
#'
#' @param dataIn The dataset
#'
#' @param steps.enabled A vector of boolean which has the same length of the 
#' steps of the pipeline. xxx
#'
#' @param remoteReset It is a remote command to reset the module. A boolean that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#'
#' @param status xxx
#' 
#' 
#' @examples
#' if(interactive()){
#' 
#' 
#' }
NULL


#' @export
#' @return NA
#' @rdname pipeline-protein
#'
mod_Protein_ui <- function(id) {
    ns <- NS(id)
}




#'
#' @return NA
#'
#'
#' @rdname pipeline-protein
#'
#' @export
#'
#' @import shiny
#'
mod_Protein_server <- function(id,
    dataIn = reactive({NULL}),
    steps.enabled = reactive({NULL}),
    remoteReset = reactive({FALSE}),
    status = reactive({NULL})) {

    #' @field global xxxx
    global <- list(
        VALIDATED = 1,
        UNDONE = 0,
        SKIPPED = -1
    )


    config <- reactiveValues(
        name = "Protein",
        steps = c("Description", "Normalization", "Filtering"),
        mandatory = c(TRUE, FALSE, FALSE)
    )

    ### -------------------------------------------------------------###
    ###                                                             ###
    ### ------------------- MODULE SERVER --------------------------###
    ###                                                             ###
    ### -------------------------------------------------------------###
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        rv <- reactiveValues(
            dataIn = NULL,
            dataOut = NULL,
            status = NULL,
            reset = NULL,
            steps.enabled = NULL
        )

        observeEvent(status(), {
            rv$status <- status()
        })

        # Initialization of the module
        observeEvent(steps.enabled(), ignoreNULL = TRUE, {
            if (is.null(steps.enabled())) {
                rv$steps.enabled <- stats::setNames(
                    rep(FALSE, length(config@steps)),
                    config@steps
                )
            } else {
                rv$steps.enabled <- steps.enabled()
            }
        })


        # Return value of module
        # DO NOT MODIFY THIS PART
        list(
            config = reactive({
                config
            }),
            dataOut = reactive({
                rv$dataOut
            }),
            status = reactive({
                rv$status
            })
        )
    })
}

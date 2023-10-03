#' @title  Tracking of entities within plots
#'
#' @description
#'
#' This shiny module offers a UI to select a subset of a dataset and 
#' superimpose quantitative values of this selection on the complete plot
#' Three modes of selection are implemented:
#'
#' - 'Protein list': xxx,
#' - 'Random': xxx,
#' - 'Specific column': xxx
#'
#' @name tracking
#'
#' @examples
#' if (interactive()) {
#'     ui <- tagList(
#'         mod_tracker_ui("track"),
#'         uiOutput("show")
#'     )
#'
#'     server <- function(input, output, session) {
#'         rv <- reactiveValues(
#'             tmp = NULL
#'         )
#'
#'         rv$tmp <- mod_tracker_server(
#'             id = "track",
#'             object = reactive({
#'                 ft[[1]]
#'             })
#'         )
#'
#'         output$show <- renderUI({
#'             p(paste0(rv$tmp(), collapse = " "))
#'         })
#'     }
#'     shinyApp(ui = ui, server = server)
#' }
NULL


#' @param id shiny id
#' @export
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs useShinyjs hidden
#'
#' @rdname tracking
#'
#' @return NA
#'
mod_tracker_ui <- function(id) {
    ns <- NS(id)

    tagList(
        useShinyjs(),
        actionButton(ns("rst_btn"), "Reset"),
        uiOutput(ns("typeSelect_ui")),
        uiOutput(ns("listSelect_ui")),
        uiOutput(ns("randSelect_ui")),
        uiOutput(ns("colSelect_ui"))
    )
}

#' @param id xxx
#' @param object A instance of the class `SummarizedExperiment`
#'
#' @rdname tracking
#'
#' @export
#'
#' @importFrom shinyjs toggle hidden show hide
#'
#' @return A `list()` of integers
#'
mod_tracker_server <- function(id,
    object) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        rv.track <- reactiveValues(
            typeSelect = "None",
            listSelect = NULL,
            randSelect = "",
            colSelect = NULL,
            indices = NULL
        )


        output$typeSelect_ui <- renderUI({
            tmp <- c("None", "ProteinList", "Random", "Column")
            nm <- c("None", "Protein list", "Random", "Specific Column")

            selectInput(ns("typeSelect"),
                "Type of selection",
                choices = stats::setNames(tmp, nm),
                selected = rv.track$typeSelect,
                width = "130px"
            )
        })


        output$listSelect_ui <- renderUI({
            widget <- selectInput(ns("listSelect"),
                "Select protein",
                choices = c("None", rowData(object())[, idcol(object())]),
                multiple = TRUE,
                selected = rv.track$listSelect,
                width = "200px",
                # size = 10,
                selectize = TRUE
            )

            if (rv.track$typeSelect == "ProteinList") {
                widget
            } else {
                hidden(widget)
            }
        })

        output$colSelect_ui <- renderUI({
            widget <- selectInput(ns("colSelect"),
                "Column of rowData",
                choices = c("", colnames(rowData(object()))),
                selected = rv.track$colSelect
            )
            if (rv.track$typeSelect == "Column") {
                widget
            } else {
                hidden(widget)
            }
        })

        output$randSelect_ui <- renderUI({
            widget <- textInput(ns("randSelect"),
                "Random",
                value = rv.track$randSelect,
                width = ("120px")
            )
            if (rv.track$typeSelect == "Random") {
                widget
            } else {
                hidden(widget)
            }
        })

        observeEvent(req(input$typeSelect), {
            rv.track$typeSelect <- input$typeSelect
        })


        observeEvent(input$rst_btn, {
            rv.track$typeSelect <- "None"
            rv.track$listSelect <- NULL
            rv.track$randSelect <- ""
            rv.track$colSelect <- NULL
            rv.track$indices <- NULL
        })

        # Catch event on the list selection
        observeEvent(input$listSelect, {
            rv.track$listSelect <- input$listSelect
            rv.track$randSelect <- ""
            rv.track$colSelect <- ""

            if (is.null(rv.track$listSelect)) {
                rv.track$indices <- NULL
            } else {
                rv.track$indices <- match(rv.track$listSelect, 
                    rowData(object())[[idcol(object())]])
            }
        })





        observeEvent(input$randSelect, {
            rv.track$randSelect <- input$randSelect
            rv.track$listSelect <- NULL
            rv.track$colSelect <- NULL
            cond <- is.null(rv.track$randSelect)
            cond <- cond || rv.track$randSelect == ""
            cond <- cond || (as.numeric(rv.track$randSelect) < 0)
            cond <- cond || (as.numeric(rv.track$randSelect) > nrow(object()))
            if (!cond) {
                rv.track$indices <- sample(seq_len(nrow(object())),
                    as.numeric(rv.track$randSelect),
                    replace = FALSE
                )
            }
        })

        observeEvent(input$colSelect, {
            rv.track$colSelect <- input$colSelect
            rv.track$listSelect <- NULL
            rv.track$randSelect <- ""

            if (rv.track$colSelect != "") {
                .op1 <- rowData(object())[, rv.track$colSelect]
                rv.track$indices <- which( .op1 == 1)
            }
        })

        return(reactive({
            rv.track$indices
        }))
    })
}

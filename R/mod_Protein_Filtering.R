btn_style <- "display:inline-block; vertical-align: middle; padding: 7px"



#' @title Module of filtering protein process
#'
#' @export
#'
#' @return NA
#'
#' @rdname pipeline-protein
#'
mod_Protein_Filtering_ui <- function(id) {
    ns <- NS(id)
}


#' @title Module of filtering protein process
#'
#' @param id xxx
#' @param dataIn xxx
#' @param steps.enabled xxx
#' @param remoteReset xxx
#'
#'
#' @importFrom shinyjs disabled
#' @export
#'
#' @return NA
#'
#' @rdname pipeline-protein
#'
mod_Protein_Filtering_server <- function(id,
    dataIn = NULL,
    steps.enabled = reactive({NULL}),
    remoteReset = reactive({FALSE})) {

    #' @field global xxxx
    global <- list(
        VALIDATED = 1,
        UNDONE = 0,
        SKIPPED = -1
    )



    # Define default selected values for widgets
    widgets.default.values <- list(
        select1 = 1,
        select2 = NULL,
        select3 = 1,
        select2_1 = 1,
        select2_2 = 1
    )


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
            reset = NULL,
            steps.enabled = NULL
        )


        dataOut <- reactiveValues(
            trigger = NULL,
            value = NULL
        )

        #' @field config xxxx
        config <- reactiveValues(
            name = "Protein_Filtering",
            steps = c("Description", "Step1", "Step2", "Step3"),
            mandatory = c(TRUE, FALSE, TRUE, TRUE)
        )

        rv.widgets <- reactiveValues(
            select1 = widgets.default.values$select1,
            select2 = widgets.default.values$select2,
            select3 = widgets.default.values$select3,
            select2_1 = widgets.default.values$select2_1,
            select2_2 = widgets.default.values$select2_2
        )



        # Initialization of the module
        observeEvent(steps.enabled(), ignoreNULL = TRUE, {
            if (is.null(steps.enabled())) {
                rv$steps.enabled <- stats::setNames(
                    rep(FALSE, length(config@steps)), config@steps)
            } else {
                rv$steps.enabled <- steps.enabled()
            }
        })


        observeEvent(remoteReset(), {
            lapply(names(rv.widgets), function(x) {
                rv.widgets[[x]] <- widgets.default.values[[x]]
            })
        })



        ###### ----------- Code for Description (step 0) ---------------#####
        output$Description <- renderUI({
            rv$steps.enabled
            # browser()
            wellPanel(
                tagList(
                    includeMarkdown(
                        system.file("app/md", paste0(config@name, ".md"), 
                            package = "MagellanNTK")),
                    uiOutput(ns("datasetDescription")),
                    if (isTRUE(rv$steps.enabled["Description"])) {
                        actionButton(ns("btn_validate_Description"),
                            paste0("Start ", config@name),
                            class = MagellanNTK::GlobalSettings$btn_success_color
                        )
                    } else {
                        shinyjs::disabled(
                            actionButton(ns("btn_validate_Description"),
                                paste0("Start ", config@name),
                                class = MagellanNTK::GlobalSettings$btn_success_color
                            )
                        )
                    }
                )
            )
        })



        observeEvent(input$btn_validate_Description, 
            ignoreInit = TRUE, ignoreNULL = TRUE, {
            rv$dataIn <- dataIn()
            dataOut$trigger <- Timestamp()
            dataOut$value <- rv$dataIn
            # rv$status['Description'] <- global$VALIDATED
        })


        ###### ------------------- Code for step 1 ----------------#####




        observeEvent(input$select1, {
            rv.widgets$select1 <- input$select1
        })
        observeEvent(input$select2, {
            rv.widgets$select2 <- input$select2
        })
        observeEvent(input$select3, {
            rv.widgets$select3 <- input$select3
        })
        observeEvent(input$select2_1, {
            rv.widgets$select2_1 <- input$select2_1
        })
        observeEvent(input$select2_2, {
            rv.widgets$select2_2 <- input$select2_2
        })


        output$test1 <- renderUI({
            # rv$steps.enabled
            rv.widgets$select1
            if (rv$steps.enabled["Step1"]) {
                selectInput(ns("select1"), "Select 1 in renderUI",
                    choices = seq_len(4),
                    selected = rv.widgets$select1,
                    width = "150px"
                )
            } else {
                shinyjs::disabled(
                    selectInput(ns("select1"), "Select 1 in renderUI",
                        choices = seq_len(4),
                        selected = rv.widgets$select1,
                        width = "150px"
                    )
                )
            }
        })



        output$test2 <- renderUI({
            rv$steps.enabled
            if (rv$steps.enabled["Step1"]) {
                selectInput(ns("select2"), "Select 2 in renderUI",
                    choices = seq_len(3),
                    selected = rv.widgets$select2,
                    width = "150px"
                )
            } else {
                shinyjs::disabled(
                    selectInput(ns("select2"), "Select 2 in renderUI",
                        choices = seq_len(4),
                        selected = rv.widgets$select2,
                        width = "150px"
                    )
                )
            }
        })




        # ------------------------ STEP 1 : UI ---------------------
        output$Step1 <- renderUI({
            # rv$steps.enabled
            name <- "Step1"
            wellPanel(
                id = ns("toto"),
                actionButton(ns("btn1"), "Btn 1"),
                tagList(
                    div(
                        id = ns("Step1a"),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle;padding-right: 20px;",
                            uiOutput(ns("test1"))
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle;padding-right: 20px;",
                            uiOutput(ns("test2"))
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 40px;",
                            if (rv$steps.enabled["Step1"]) {
                                selectInput(ns("select3"), "Select step 3",
                                    choices = seq_len(3),
                                    selected = rv.widgets$select3,
                                    width = "150px"
                                )
                            } else {
                                shinyjs::disabled(
                                    selectInput(ns("select3"), "Select step 3",
                                        choices = seq_len(5),
                                        selected = rv.widgets$select3,
                                        width = "150px"
                                    )
                                )
                            }
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle;padding-right: 20px;",
                            if (rv$steps.enabled["Step1"]) {
                                actionButton(ns(paste0("btn_validate_", name)),
                                    "Perform",
                                    class = MagellanNTK::GlobalSettings$btn_success_color
                                )
                            } else {
                                shinyjs::disabled(
                                    actionButton(
                                        ns(paste0("btn_validate_", name)),
                                        "Perform",
                                        class = MagellanNTK::GlobalSettings$btn_success_color
                                    )
                                )
                            }
                        )
                    )
                )
            )
        })


        observeEvent(input$btn_validate_Step1, ignoreInit = TRUE, {
            # Add your stuff code here
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn

            # rv$status['Step1'] <- global$VALIDATED
        })

        #-------------------------- Code for step 2 -------



        output$select2_1_UI <- renderUI({
            rv$steps.enabled
            if (rv$steps.enabled["Step2"]) {
                selectInput(ns("select2_1"), "Select 2_1 in renderUI",
                    choices = seq_len(3),
                    selected = rv.widgets$select2_1,
                    width = "150px"
                )
            } else {
                shinyjs::disabled(
                    selectInput(ns("select2_1"), "Select 2_1 in renderUI",
                        choices = seq_len(3),
                        selected = rv.widgets$select2_1,
                        width = "150px"
                    )
                )
            }
        })

        output$Step2 <- renderUI({
            rv$steps.enabled
            name <- "Step2"
            wellPanel(
                tagList(
                    div(
                        id = ns(name),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle;padding-right: 20px;",
                            uiOutput(ns("select2_1_UI"))
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle; padding-right: 40px;",
                            if (rv$steps.enabled["Step2"]) {
                                selectInput(ns("select2_2"), "Select 2_2",
                                    choices = seq_len(5),
                                    selected = rv.widgets$select2_2,
                                    width = "150px"
                                )
                            } else {
                                shinyjs::disabled(
                                    selectInput(ns("select2_2"),
                                        "Select 2_2",
                                        choices = seq_len(5),
                                        selected = rv.widgets$select2_2,
                                        width = "150px"
                                    )
                                )
                            }
                        ),
                        div(
                            style = "display:inline-block; 
                            vertical-align: middle;padding-right: 20px;",
                            if (rv$steps.enabled["Step2"]) {
                                actionButton(
                                    ns(paste0("btn_validate_", name)),
                                    "Perform",
                                    class = MagellanNTK::GlobalSettings$btn_success_color
                                )
                            } else {
                                shinyjs::disabled(
                                    actionButton(
                                        ns(paste0("btn_validate_", name)),
                                        "Perform",
                                        class = MagellanNTK::GlobalSettings$btn_success_color
                                    )
                                )
                            }
                        )
                    )
                )
            )
        })

        observeEvent(input$btn_validate_Step2, ignoreInit = TRUE, {
            # Add your stuff code here
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn

            # rv$status['Step2'] <- global$VALIDATED
        })


        #------------- Code for step 3 ---------------

        output$Step3 <- renderUI({
            rv$steps.enabled
            tagList(
                h3("Step 3"),
                if (rv$steps.enabled["Step3"]) {
                    actionButton(ns("btn_validate_Step3"),
                        "Perform",
                        class = MagellanNTK::GlobalSettings$btn_success_color
                    )
                } else {
                    shinyjs::disabled(
                        actionButton(ns("btn_validate_Step3"),
                            "Perform",
                            class = MagellanNTK::GlobalSettings$btn_success_color
                        )
                    )
                }
            )
        })


        observeEvent(input$btn_validate_Step3, ignoreInit = TRUE, {
            # Add your stuff code here
            rv$dataIn <- AddItemToDataset(rv$dataIn, config@name)
            dataOut$trigger <- MagellanNTK::Timestamp()
            dataOut$value <- rv$dataIn

            # rv$status['Step3'] <- global$VALIDATED
        })



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
            # status = reactive({rv$status})
        )
    })
}

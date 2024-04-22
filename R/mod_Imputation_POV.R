#' @title Filtering Shiny module
#'
#' @description
#' This function is a shiny module to xxx
#' This function is written with specifications of the package `MagellanNTK` so
#' as to be easily integrated into workflfow compliant with `MagellanNTK`.
#'
#' @name mod_Imputation_POV
#'
#' @return As for all modules used with `MagellanNTK`, the return value is a
#' `list()` of two items:
#' - trigger : xxx
#' - value: In this case, it contains a list() of three slots:
#'   - ll.var: a list() of instances of the class `Filtering`,
#'   - ll.query: a list of `character()` which describe the queries in natural
#'   language,
#'   - ll.widgets.value: a list of the values of widgets.
#'
#' @examplesIf interactive()
#' data(ft_na)
#' shiny::runApp(mod_Imputation_POV(ft_na[[1]]))
#' 
NULL





#' @export
#'
#' @rdname mod_Imputation_POV
#'
mod_Imputation_POV_ui <- function(id) {
  ns <- NS(id)
  .localStyle <- "display:inline-block; vertical-align: top;
                  padding-right: 20px;"
  wellPanel(
    # uiOutput for all widgets in this UI
    # This part is mandatory
    # The renderUI() function of each widget is managed by MagellanNTK
    # The dev only have to define a reactive() function for each
    # widget he want to insert
    # Be aware of the naming convention for ids in uiOutput()
    # For more details, please refer to the dev document.
    
    tags$div(
      tags$div(style = .localStyle, uiOutput(ns("POV_algorithm_ui"))),
      tags$div(style = .localStyle, uiOutput(ns("POV_Params"))),
      tags$div(style = .localStyle, uiOutput(ns("POV_showDetQuantValues")))
    ),
    # Insert validation button
    uiOutput(ns("mod_Imputation_POV_btn_validate_ui")),
    htmlOutput("helpForImputation"),
    tags$hr(),
    moduleMVPlotsUI("mvImputationPlots_MV")
  )
}




#' @param id xxx
#' @param obj An instance of the class `QFeatures`
#' @param reset A `integer(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
#' @rdname mod_Imputation_POV
#'
#' @export
#'
mod_Imputation_POV_server <- function(id,
  obj = reactive({NULL}),
  i = reactive({1}),
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    POV_algorithm = NULL
  )
  
  rv.custom.default.values <- list(
    
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
    
    
    observeEvent(obj(), ignoreNULL = TRUE,{
      req(obj())
      #browser()
      stopifnot(inherits(obj(), 'QFeatures'))
      rv$dataIn <- obj()
      
      
    }, priority = 1000)
    
    
    output$POV_algorithm_ui <- renderUI({
      widget <- selectInput(ns("POV_algorithm"), "Algorithm for POV",
        choices = imputationAlgorithmsProteins_POV,
        selected = rv.widgets$POV_algorithm,
        width = "150px"
      )
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$POV_showDetQuantValues <- renderUI({
      req(rv.widgets$POV_algorithm == "detQuantile")
       tagList(
          h5("The POV will be imputed by the following values :"),
          mod_DetQuantImpValues_ui("POV_DetQuantValues_DT")
        )
    })
    
    
    
    output$POV_Params_ui <- renderUI({
      req(rv.widgets$POV_algorithm)
      
      isolate({
        switch(rv$widgets$POV_algorithm,
          detQuantile = {
            tagList(
              tags$div(style = .localStyle,
                numericInput(ns("POV_detQuant_quantile"), "Quantile",
                  value = rv.widgets$POV_detQuant_quantile,
                  step = 0.5, min = 0, max = 100, width = "100px"
                )
              ),
              tags$div(style = .localStyle,
                numericInput(ns("POV_detQuant_factor"), "Factor",
                  value = rv.widgets$POV_detQuant_factor,
                  step = 0.1, min = 0, max = 10, width = "100px"
                )
              )
            )
          },
          KNN = {
            numericInput(ns("KNN_nbNeighbors"), "Neighbors",
              value = rv.widgets$POV_KNN_n, step = 1, min = 0,
              max = max(nrow(rv$current.obj), rv.widgets$POV_KNN_n),
              width = "100px"
            )
          }
        )
      })
    })
    
    
    
    output$mod_Imputation_POV_btn_validate_ui <- renderUI({
      #browser()
      req(xxx)
      
      widget <- actionButton(ns("mod_Imputation_POV_btn_validate"),
        "Perform POV imputation", class = "btn-success")
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$mod_Imputation_POV_btn_validate, {
      
      rv$current.obj
      m <- match.metacell(DAPAR::GetMetacell(rv$current.obj),
        pattern = "Missing POV",
        level = DAPAR::GetTypeofData(rv$current.obj)
      )
      nbPOVBefore <- length(which(m))
      
      withProgress(message = "", detail = "", value = 0, {
        incProgress(0.25, detail = "Find MEC blocks")
        
        .tmp <- NULL
        .tmp <- try({
          switch(rv$widgets$proteinImput$POV_algorithm,
            None = .tmp <- rv$dataset[[input$datasets]],
            slsa = {
              incProgress(0.5, detail = "slsa Imputation")
              .tmp <- wrapper.impute.slsa(
                rv$dataset[[input$datasets]])
            },
            detQuantile = {
              incProgress(0.5, detail = "det quantile Imputation")
              .tmp <- wrapper.impute.detQuant(
                obj = rv$dataset[[input$datasets]],
                qval = rv$widgets$proteinImput$POV_detQuant_quantile / 100,
                factor = rv$widgets$proteinImput$POV_detQuant_factor,
                na.type = 'Missing POV')
            },
            KNN = {
              incProgress(0.5, detail = "KNN Imputation")
              .tmp <- wrapper.impute.KNN(obj = rv$dataset[[input$datasets]],
                K = rv$widgets$proteinImput$POV_KNN_n)
            }
          )
        })
        
        if(inherits(.tmp, "try-error")) {
          mod_SweetAlert_server(id = 'sweetalert_perform_POVimputation_button',
            text = .tmp[[1]],
            type = 'error' )
        } else {
          # sendSweetAlert(
          #   session = session,
          #   title = "Success",
          #   type = "success"
          # )
          rv$current.obj <- .tmp
          # incProgress(0.75, detail = 'Reintroduce MEC blocks')
          incProgress(1, detail = "Finalize POV imputation")
          m <- match.metacell(DAPAR::GetMetacell(rv$current.obj),
            pattern = "Missing POV",
            level = DAPAR::GetTypeofData(rv$current.obj)
          )
          nbPOVAfter <- length(which(m))
          rv$nbPOVimputed <- nbPOVBefore - nbPOVAfter
          
          rv$impute_Step <- 1
          rv$imputePlotsSteps[["step1"]] <- rv$current.obj
          rvModProcess$moduleProtImputationDone[1] <- TRUE
          shinyjs::enable("perform.imputationMEC.button")
          shinyjs::enable("ValidImputation")
        }
      par <- rv.custom$funFilter()$value$ll.widgets.value
      params(tmp[[i()]], length(tmp[[i()]])) <- par
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- tmp
    })
      
    })
    
    return(reactive({dataOut}))
  })
  }




#' @export
#' @rdname mod_Imputation_POV
#' 
mod_Imputation_POV <- function(obj, i){
  ui <- mod_Imputation_POV_ui('pov')
  
  server <- function(input, output, session){
    
    res <- mod_Imputation_POV_server('pov',
      obj = reactive({obj}),
      i = reactive({i}))
    
    observeEvent(res()$trigger, {
      print(res()$value)
    })
  }
  
  app <- shiny::shinyApp(ui, server)
  
}
#' @title Filtering Shiny module
#'
#' @description
#' This function is a shiny module to xxx
#' This function is written with specifications of the package `MagellanNTK` so
#' as to be easily integrated into workflfow compliant with `MagellanNTK`.
#'
#' @name mod_Prot_Imputation_MEC
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
#' data(Exp1_R25_prot, package="DaparToolshedData")
#' obj <- Exp1_R25_prot[seq_len(100)]
#' shiny::runApp(mod_Prot_Imputation_MEC(obj, length(obj)))
#' 
NULL





#' @export
#'
#' @rdname mod_Prot_Imputation_MEC
#'
mod_Prot_Imputation_MEC_ui <- function(id) {
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
 
    #----------------------------------------------------
    
    
      uiOutput(ns('mv_impute_MEC'))
  )
}




#' @param id xxx
#' @param obj An instance of the class `QFeatures`
#' @param reset A `integer(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
#' @rdname mod_Prot_Imputation_MEC
#'
#' @export
#'
mod_Prot_Imputation_MEC_server <- function(id,
  obj = reactive({NULL}),
  i = reactive({1}),
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    MEC_algorithm = NULL,
    MEC_KNN_n = 10,
    MEC_detQuant_quantile = 2.5,
    MEC_detQuant_factor = 1,
    MEC_fixedValue = 0
  )
  
  rv.custom.default.values <- list(
    mv.present = FALSE
  )
  
  imputationAlgorithmsProteins_MEC <- list(
    "None" = "None",
    "Det quantile" = "detQuantile",
    "Fixed value" = "fixedValue"
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
      stopifnot(inherits(obj(), 'QFeatures'))
      rv$dataIn <- obj()
      
      qdata <- SummarizedExperiment::assay(rv$dataIn[[i()]])
      rv.custom$mv.present <- sum(is.na(qdata)) > 0
      
      # If there is no MV in the dataset, return it
      if (!rv.custom$mv.present){
        dataOut$trigger <- MagellanNTK::Timestamp()
        dataOut$value <- rv$dataIn
      }
      
    }, priority = 1000)
    
    
    
    
    output$mv_impute_MEC <- renderUI({
      req(rv$dataIn)
      widget <- NULL
      
      if (!rv.custom$mv.present) {
        widget <- tags$p("Your dataset does not contain any missing values.")
      } else {
        
      widget <- tagList(
        uiOutput(ns("warningMECImputation")),
        uiOutput(ns("MEC_chooseImputationMethod_ui")),
        uiOutput(ns("MEC_Params_ui")),
        uiOutput(ns("MEC_showDetQuantValues_ui")),
        uiOutput(ns("mod_Prot_Imputation_MEC_btn_validate_ui")),
        tags$hr(),
        withProgress(message = "", detail = "", value = 0, {
          incProgress(0.5, detail = "Building plots...")
          uiOutput(ns('mvplots_ui'))
        })
      )
      }
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$mvplots_ui <- renderUI({
      widget <- mod_mv_plots_ui(ns("mvplots"))
      toggleWidget(widget, is.enabled())
    })
    
    
    observe({
      req(rv$dataIn)
      mod_mv_plots_server("mvplots",
        data = reactive({rv$dataIn[[i()]]}),
        grp = reactive({get_group(rv$dataIn)}),
        mytitle = "MEC imputation",
        pal = reactive({NULL}),
        pattern = c("Missing", "Missing POV", "Missing MEC")
      )
    })
    
    
    
    output$MEC_chooseImputationMethod_ui <- renderUI({
      widget <- selectInput(ns("MEC_algorithm"), "Algorithm for MEC",
        choices = imputationAlgorithmsProteins_MEC,
        selected = rv.widgets$MEC_algorithm, 
        width = "150px"
      )
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$MEC_showDetQuantValues_ui <- renderUI({
      req(rv.widgets$MEC_algorithm == "detQuantile")
      tagList(
        h5("The MEC will be imputed by the following values :"),
        mod_DetQuantImpValues_ui(ns("MEC_DetQuantValues_DT"))
      )
    })
    
    
    
    output$MEC_Params_ui <- renderUI({
      req(rv.widgets$MEC_algorithm)
      
      widget <- switch(rv.widgets$MEC_algorithm,
          detQuantile = {
            tagList(
              numericInput(ns("MEC_detQuant_quantile"), "Quantile",
                  value = rv.widgets$MEC_detQuant_quantile,
                  step = 0.5,
                  min = 0,
                  max = 100,
                  width = "100px"
                ),
              numericInput(ns("MEC_detQuant_factor"), "Factor",
                  value = rv.widgets$MEC_detQuant_factor,
                  step = 0.1, min = 0, max = 10,
                  width = "100px"
                )
              )
          },
          fixedValue = {
            numericInput(ns("MEC_fixedValue"), "Fixed value",
              value = rv.widgets$MEC_fixedValue,
              step = 0.1, min = 0, max = 100,
              width = "100px"
            )
          }
        )
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    
    output$mod_Prot_Imputation_MEC_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("mod_Prot_Imputation_MEC_btn_validate"),
        "Perform MEC imputation", class = "btn-success")
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
   
    
    observeEvent(input$mod_Prot_Imputation_MEC_btn_validate, {
      
      req(rv$dataIn)
      withProgress(message = "", detail = "", value = 0, {
        incProgress(0.25, detail = "Reintroduce MEC")
        
      m <- match.metacell(
        get_metacell(rv$dataIn[[i()]]),
        pattern = "Missing MEC",
        level = get_type(rv$dataIn[[i()]])
      )
      nbMECBefore <- length(which(m))
      incProgress(0.75, detail = "MEC Imputation")
      
      #browser()
      withProgress(message = "", detail = "", value = 0, {
        incProgress(0.25, detail = "Find MEC blocks")
        
        .tmp <- NULL
        .param <- list()
        try({
          switch(rv.widgets$MEC_algorithm,
            None = .tmp <- rv$dataIn[[i()]],
            
            detQuantile = {
              incProgress(0.5, detail = "det quantile Imputation")
              .tmp <- wrapper.impute.detQuant(
                obj = rv$dataIn[[i()]],
                qval = rv.widgets$MEC_detQuant_quantile / 100,
                factor = rv.widgets$MEC_detQuant_factor,
                na.type = 'Missing MEC')
              .param <- list(
                MEC_algorithm = rv.widgets$MEC_algorithm,
                qval = rv.widgets$MEC_detQuant_quantile / 100,
                factor = rv.widgets$MEC_detQuant_factor,
                na.type = 'Missing MEC')
            },
            fixedValue = {
              .tmp <- wrapper.impute.fixedValue(
                obj = rv$dataIn[[i()]],
                fixVal = rv.widgets$MEC_fixedValue,
                na.type = "Missing MEC"
              )
              .param <- list(
                MEC_algorithm = rv.widgets$MEC_algorithm,
                fixVal = rv.widgets$MEC_fixedValue,
                na.type = 'Missing MEC')
            }
          )
        })
        
        if(inherits(.tmp, "try-error")) {
          mod_SweetAlert_server(id = 'sweetalert_perform_MECimputation_button',
            text = .tmp,
            type = 'error' )
        } else {
          # sendSweetAlert(
          #   session = session,
          #   title = "Success",
          #   type = "success"
          # )
          #rv$dataIn[[i()]] <- .tmp
          # incProgress(0.75, detail = 'Reintroduce MEC blocks')
          incProgress(1, detail = "Finalize MEC imputation")
          
          
          m <- match.metacell(omXplore::get_metacell(.tmp),
            pattern = "Missing MEC",
            level = omXplore::get_type(.tmp)
          )
          nbMECAfter <- length(which(m))
          rv$nbMECimputed <- nbMECBefore - nbMECAfter
          
        }
        
        params(.tmp) <- .param
        dataOut$trigger <- MagellanNTK::Timestamp()
        dataOut$value <- .tmp
      })
      })
    })
    
    return(reactive({dataOut}))
  })
}




#' @export
#' @rdname mod_Prot_Imputation_MEC
#' 
mod_Prot_Imputation_MEC <- function(obj, i){
  ui <- mod_Prot_Imputation_MEC_ui('mec')
  
  server <- function(input, output, session){
    
    res <- mod_Prot_Imputation_MEC_server('mec',
      obj = reactive({obj}),
      i = reactive({i}))
    
    observeEvent(res()$trigger, {
      print(res()$value)
    })
  }
  
  app <- shiny::shinyApp(ui, server)
  
}
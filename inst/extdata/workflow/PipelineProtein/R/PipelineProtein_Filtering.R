#' @title Shiny example process module.
#'
#' @description
#' This module contains the configuration informations for the corresponding pipeline.
#' It is called by the nav_pipeline module of the package MagellanNTK
#' 
#' The name of the server and ui functions are formatted with keywords separated by '_', as follows:
#' * first string `mod`: indicates that it is a Shiny module
#' * `pipeline name` is the name of the pipeline to which the process belongs
#' * `process name` is the name of the process itself
#' 
#' This convention is important because MagellanNTK call the different
#' server and ui functions by building dynamically their name.
#' 
#' In this example, `PipelineProtein_Filtering_ui()` and `PipelineProtein_Filtering_server()` define
#' the code for the process `PipelineProtein` which is part of the pipeline called `PipelineProtein`.
#'
#' @name example_module_process1
#' 
#' @param id xxx
#' @param dataIn The dataset
#' @param steps.enabled A vector of boolean which has the same length of the steps
#' of the pipeline. This information is used to enable/disable the widgets. It is not
#' a communication variable between the caller and this module, thus there is no
#' corresponding output variable
#' @param remoteReset It is a remote command to reset the module. A boolean that
#' indicates is the pipeline has been reseted by a program of higher level
#' Basically, it is the program which has called this module
#' @param steps.status xxx
#' @param current.pos xxx
#' @param path xxx
#' 
#' 
#' 
#' 
#' 
#' @author Samuel Wieczorek
#' 
#' 
NULL

#' @rdname example_module_process1
#' @export
#' 
PipelineProtein_Filtering_conf <- function(){
  Config(
    fullname = 'PipelineProtein_Filtering',
    mode = 'process',
    steps = c("Quanti metadata filtering", "Variable filtering"),
    mandatory = c(FALSE, FALSE)
  )
}


#' @rdname example_module_process1
#' 
#' @export
#'
PipelineProtein_Filtering_ui <- function(id){
  ns <- NS(id)
}



#' @rdname example_module_process1
#' 
#' @importFrom stats setNames rnorm
#' 
#' @export
#' 
PipelineProtein_Filtering_server <- function(id,
  dataIn = reactive({NULL}),
  steps.enabled = reactive({NULL}),
  remoteReset = reactive({FALSE}),
  steps.status = reactive({NULL}),
  current.pos = reactive({1}),
  path = NULL
  ){
  
  
  
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
    deleted.stringBased = NULL,
    deleted.metacell = NULL,
    deleted.numeric = NULL,
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
  
  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    core.code <- Get_Workflow_Core_Code(
      mode = 'process',
      name = id,
      w.names = names(widgets.default.values),
      rv.custom.names = names(rv.custom.default.values)
    )
    
    eval(str2expression(core.code))
    
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      file <- normalizePath(file.path(session$userData$workflow.path, 
        'md', paste0(id, '.md')))
      tagList(
        ### In this example, the md file is found in the extdata/module_examples 
        ### directory but with a real app, it should be provided by the package 
        ### which contains the UI for the different steps of the process module.
        ### system.file(xxx)
        
        if (file.exists(file))
          includeMarkdown(file)
        else
          p('No Description available'),
        
        
        # Used to show some information about the dataset which is loaded
        # This function must be provided by the package of the process module
        uiOutput(ns('datasetDescription_ui')),
        
        # Insert validation button
        uiOutput(ns('Description_btn_validate_ui'))
      )
    })
    
    output$datasetDescription_ui <- renderUI({
      # Insert your own code to visualize some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"),
        "Start",
        class = btn_success_color)
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- stepStatus$VALIDATED
    })
    
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$Quantimetadatafiltering <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        DT::dataTableOutput(ns("qMetacell_Filter_DT")),
        uiOutput(ns("Quantimetadatafiltering_buildQuery_ui")),
        uiOutput(ns("example_ui")),
        mod_ds_qMetacell_ui(ns("plots")),
        # Insert validation button
        uiOutput(ns("Quantimetadatafiltering_btn_validate_ui"))
      )
    })
    
    
    # >>> START: Definition of the widgets
    
    
    
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
      se = reactive({ mainAssay(rv$dataIn)}),
      init.pattern = "missing",
      conds = design.qf(rv$dataIn)$Condition
    )
    
    output$qMetacell_Filter_DT <- DT::renderDataTable(
      server = TRUE,{
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
      obj = reactive({mainAssay(rv$dataIn)}),
      conds = reactive({design.qf(rv$dataIn)$Condition}),
      list_tags = reactive({
        req(rv$dataIn)
        c(
          "None" = "None",
          qMetacell.def(typeDataset(mainAssay(rv$dataIn)))$node
        )
      }),
      keep_vs_remove = reactive({stats::setNames(nm = c("delete", "keep"))}),
      val_vs_percent = reactive({stats::setNames(nm = c("Count", "Percentage"))}),
      operator = reactive({stats::setNames(nm = SymFilteringOperators())})
    )
    
    
    output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
        "Perform qMetacell filtering",
        class = btn_success_color
      )
      
      cond <- length(rv.custom$funFilter$ll.fun()) > 0
      cond <- cond && rv$steps.enabled["Quantimetadatafiltering"]
      MagellanNTK::toggleWidget(widget, cond)
    })
    # >>> END: Definition of the widgets
    
    
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
      rv$steps.status["Quantimetadatafiltering"] <- stepStatus$VALIDATED
    })
    
    
    output$showPlot <- renderPlot({
      req(rv$dataIn)
      plot(as.matrix(rv$dataIn[[1]][,1]))
    })
    # <<< END ------------- Code for step 1 UI---------------
    
    
    # >>> START ------------- Code for step 2 UI---------------
    
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
      rv$steps.status["Variablefiltering"] <- stepStatus$VALIDATED
    })
    
    # <<< END ------------- Code for step 2 UI---------------
    
    
    # >>> START ------------- Code for step 'Save' UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Save_btn_validate_ui')),
        uiOutput(ns('dl_ui'))
      )
    })
    
    output$dl_ui <- renderUI({
      req(config@mode == 'process')
      req(rv$steps.status['Save'] == stepStatus$VALIDATED)
      dl_ui(ns('createQuickLink'))
    })
    
    output$Save_btn_validate_ui <- renderUI({
      toggleWidget(
        actionButton(ns("Save_btn_validate"), "Save",
                     class = btn_success_color),
        rv$steps.enabled['Save']
        )
    })
    observeEvent(input$Save_btn_validate, {
      # Do some stuff
      rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
                                          dataset = new.dataset,
                                          name = id)
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- stepStatus$VALIDATED
      dl_server('createQuickLink', dataIn = reactive({rv$dataIn}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

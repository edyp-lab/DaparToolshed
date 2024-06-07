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
#' @name PipelineProtein
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
#' @examplesIf interactive()
#' library(MagellanNTK)
#' data(ft_na)
#' path <- system.file('workflow/PipelineProtein', package = 'DaparToolshed')
#' shiny::runApp(workflowApp("PipelineProtein_Filtering", path, dataIn = ft_na))
#' 
#' 
#' 
#' @author Samuel Wieczorek
#' 
#' 
NULL

#' @rdname PipelineProtein
#' @export
#' 
PipelineProtein_Filtering_conf <- function(){
  Config(
    fullname = 'PipelineProtein_Filtering',
    mode = 'process',
    steps = c("Cell metadata filtering", "Variable filtering"),
    mandatory = c(FALSE, FALSE)
  )
}


#' @rdname PipelineProtein
#' 
#' @export
#'
PipelineProtein_Filtering_ui <- function(id){
  ns <- NS(id)
}



#' @rdname PipelineProtein
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
  widgets.default.values <- list()
  
  
  rv.custom.default.values <- list(
    deleted.stringBased = NULL,
    deleted.metacell = NULL,
    deleted.numeric = NULL,
    tmp1 = reactive({NULL}),
    tmp2 = reactive({NULL}),
    tmp = reactive({NULL})
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
        class = "btn-success")
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      req(dataIn())
      rv$dataIn <- rv.custom$tmp <- dataIn()
      
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- stepStatus$VALIDATED
    })
    
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$Cellmetadatafiltering <- renderUI({
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        uiOutput(ns("mod_metacell_filtering_ui")),
        # Insert validation button
        uiOutput(ns("Cellmetadatafiltering_btn_validate_ui"))
      )
    })
    
    
    observe({
     # req(rv$dataIn)
    rv.custom$tmp1 <- DaparToolshed::mod_Metacell_Filtering_server(
      id = "metaFiltering",
      obj = reactive({rv$dataIn}),
      i = reactive({length(rv$dataIn)}),
      is.enabled = reactive({rv$steps.enabled["Cellmetadatafiltering"]})
    )
    })
    
    # >>> START: Definition of the widgets
    output$mod_metacell_filtering_ui <- renderUI({

      widget <- DaparToolshed::mod_Metacell_Filtering_ui(ns("metaFiltering"))
      MagellanNTK::toggleWidget(widget, 
        rv$steps.enabled["Cellmetadatafiltering"])
    })
    

    
    output$Cellmetadatafiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Cellmetadatafiltering_btn_validate"),
        "Validate step",
        class = "btn-success"
      )
       MagellanNTK::toggleWidget(widget, 
        rv$steps.enabled["Cellmetadatafiltering"])
    })
    # >>> END: Definition of the widgets
    
    
    
    observeEvent(input$Cellmetadatafiltering_btn_validate, {
      rv.custom$tmp <- rv.custom$tmp1()$value
      
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- NULL
      rv$steps.status["Cellmetadatafiltering"] <- stepStatus$VALIDATED
    })
    

    # <<< END ------------- Code for step 1 UI---------------
    
    
    # >>> START ------------- Code for step 2 UI---------------
    
    output$Variablefiltering <- renderUI({
      wellPanel(
        uiOutput(ns("mod_variable_filtering_ui")),
        # Process the queries
        uiOutput(ns("Variablefiltering_btn_validate_ui"))
      )
    })
    

    observe({
      req(rv$dataIn)
      dataIn <- rv$dataIn
      # # If the previous step has been run and validated,
      # # Update dataIn to its result
      #if (rv$steps.status["Cellmetadatafiltering"] == stepStatus$VALIDATED)
      if (!is.null(rv.custom$tmp1()$value))
        dataIn <- rv.custom$tmp1()$value
      
      rv.custom$tmp2 <- DaparToolshed::mod_Variable_Filtering_server(
        id = "varFiltering",
        obj = reactive({dataIn}),
        i = reactive({length(dataIn)}),
        is.enabled = reactive({rv$steps.enabled["Variablefiltering"]})
      )
      
    })
    
    
    output$mod_variable_filtering_ui <- renderUI({
    widget <- DaparToolshed::mod_Variable_Filtering_ui(ns("varFiltering"))
    MagellanNTK::toggleWidget(widget, 
      rv$steps.enabled["Variablefiltering"])
  })
  

    output$Variablefiltering_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Variablefiltering_btn_validate"),
        "Validate step",
        class = "btn-success"
      )
      MagellanNTK::toggleWidget( widget, rv$steps.enabled["Variablefiltering"])
    })

    
    
    observeEvent(input$Variablefiltering_btn_validate, {
      
      rv.custom$tmp <- rv.custom$tmp2()$value
      
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- NULL
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
      req(rv$steps.status['Save'] == stepStatus$VALIDATED)
      req(config@mode == 'process')
      
      MagellanNTK::download_dataset_ui(ns('createQuickLink'))
    })
    
    output$Save_btn_validate_ui <- renderUI({
      toggleWidget(
        actionButton(ns("Save_btn_validate"), "Save",
                     class = "btn-success"),
        rv$steps.enabled['Save']
        )
    })
    
    observeEvent(input$Save_btn_validate, {
      # Do some stuff
      # Clean the result
      nTotal <- length(rv.custom$tmp)
      nOriginal <- length(rv$dataIn)
      if (nTotal- nOriginal > 1)
        rv.custom$tmp <- QFeatures::removeAssay(
          rv.custom$tmp, 
          (nOriginal+1):(nTotal-1))
      
      # Rename the new dataset with the name of the process
      names(rv.custom$tmp)[length(rv.custom$tmp)] <- 'Filtering'
      
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv.custom$tmp
      rv$steps.status['Save'] <- stepStatus$VALIDATED
      
      
      MagellanNTK::download_dataset_server('createQuickLink', 
        dataIn = reactive({rv.custom$tmp}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

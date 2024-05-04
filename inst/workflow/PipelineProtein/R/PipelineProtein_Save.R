

#' @export
#' 
PipelineProtein_Save_conf <- function(){
  Config(
    fullname = 'PipelineProtein_Save',
    mode = 'process'
  )
}



#' @export
PipelineProtein_Save_ui <- function(id){
  ns <- NS(id)
}


#' @export
PipelineProtein_Save_server <- function(id,
  dataIn = reactive({NULL}),
  steps.enabled = reactive({NULL}),
  remoteReset = reactive({FALSE}),
  steps.status = reactive({NULL}),
  current.pos = reactive({1})
){
  
  
  
  # Define default selected values for widgets
  # By default, this list is empty for the Save module
  # but it can be customized
  widgets.default.values <- NULL
  rv.custom.default.values <- NULL
  
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
    
    ###### ------------------- Code for Save (step 0) -------------------------    #####
    output$Save <- renderUI({
      file <- normalizePath(file.path(session$userData$workflow.path, 
        'md', paste0(id, '.md')))
      tagList(
        if (file.exists(file))
          includeMarkdown(file)
        else
          p('No Save available'),
        
        uiOutput(ns('datasetSave_ui')),
        
        # Insert validation button
        uiOutput(ns('Save_btn_validate_ui'))
      )
    })
    
    
    
    output$datasetSave_ui <- renderUI({
      # Insert your own code to visualize some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Save_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Save_btn_validate"),
        "Start",
        class = btn_success_color)
      toggleWidget(widget, rv$steps.enabled['Save'])
    })
    
    
    observeEvent(input$Save_btn_validate, {
      # In this process, there is no dataset resend to the server
      # This is why the dataOut$value is set to NULL. This triggers the 
      # validation of the step but without rebuilds the vector of datasets 
      # to send
      rv$dataIn <- dataIn()
      dataOut$trigger <- Timestamp()
      dataOut$value <- NULL
      rv$steps.status['Save'] <- stepStatus$VALIDATED
    })
    
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
    
  }
  )
}

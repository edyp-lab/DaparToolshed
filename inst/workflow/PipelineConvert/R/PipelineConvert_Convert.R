

#' @title xxx
#' @description xxx
#' @name mod_convert
#' @author Samuel Wieczorek 
#' @examplesIf interactive()
#' source("~/GitHub/DaparToolshed/inst/extdata/workflow/WorkflowConvert/R/Workflowconvert_Convert.R")
#' path <- system.file('extdata/workflow/WorkflowConvert', package = 'DaparToolshed')
#' shiny::runApp(MagellanNTK::workflowApp("Convert")
#' 
#' 
NULL

redBtnClass <- "btn-danger"
PrevNextBtnClass <- "btn-info"
btn_success_color <- 'info'

optionsBtnClass <- "info"
options(shiny.fullstacktrace = TRUE)

#' @rdname mod_convert
#' @export
#' 
PipelineConvert_Convert_conf <- function(){
  # This list contains the basic configuration of the process
  Config(
    fullname = 'PipelineConvert_Convert',
    # Define the type of module
    mode = 'process',
    # List of all steps of the process
    steps = c('Select File', 'Data Id', 'Exp and Feat Data', 'Design'),
    # A vector of boolean indicating if the steps are mandatory or not.
    mandatory = c(TRUE, TRUE, TRUE, TRUE)
    
  )
}



#' @title   mod_choose_pipeline_ui and mod_choose_pipeline_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#'
#' @rdname mod_convert
#'
#' @keywords internal
#' @export
#'
#' @importFrom shiny NS tagList
#' @import sos
#'
#' @return NA
#'
PipelineConvert_Convert_ui <- function(id) {
  ns <- NS(id)
}




#' Convert Server Function
#'
#' @param id xxx
#' @param dataIn xxx
#' @param steps.enabled xxx
#' @param remoteReset xxx
#'
#' @import QFeatures
#' @importFrom shinyalert shinyalert
#' @importFrom shinyjs disabled
#'
#' @export
#'
#' @rdname mod_convert
#'
#' @return NA
#'
PipelineConvert_Convert_server <- function(id,
                           dataIn = reactive({NULL}),
                           steps.enabled = reactive({NULL}),
                           remoteReset = reactive({FALSE}),
                           steps.status = reactive({NULL}),
                           current.pos = reactive({1})
) {
  
  pkgs.require("openxlsx")
  
  
  
  widgets.default.values <- list(
    SelectFile_software = '',
    SelectFile_file = NULL,
    SelectFile_typeOfData = NULL,
    SelectFile_checkDataLogged = "no",
    SelectFile_replaceAllZeros = TRUE,
    SelectFile_XLSsheets = '',
    
    DataId_datasetId = NULL,
    DataId_parentProteinID = NULL,
    
    ExpandFeatData_idMethod = NULL,
    ExpandFeatData_quantCols = NULL,
    ExpandFeatData_inputGroup = NULL
  )
  
  rv.custom.default.values <- list(
    tab = NULL,
    design = NULL
  )
  
  ### -------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ### -------------------------------------------------------------###
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
    
    # #Define local reactive variables
    rv.convert <- reactiveValues(
      tab = NULL,
      design = NULL,
      inputGroup = NULL
    )
    
    # >>>
    # >>> START ------------- Code for Description UI---------------
    # >>> 
    
    
    output$Description <- renderUI({
      file <- normalizePath(file.path(session$userData$workflow.path, 
        'md', paste0(id, '.md')))
      
      
      tagList(
        # In this example, the md file is found in the extdata/module_examples directory
        # but with a real app, it should be provided by the package which
        # contains the UI for the different steps of the process module.
        uiOutput(ns('Description_btn_validate_ui')),
        
        # Used to show some information about the dataset which is loaded
        # This function must be provided by the package of the process module
        uiOutput(ns('datasetDescription_ui')),
        
        if (file.exists(file))
          includeMarkdown(file)
        else
          p('No Description available')
      )
    })
    
    output$datasetDescription_ui <- renderUI({
      # Insert your own code to visualize some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"), "Start",
                             class = 'info')
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- stepStatus$VALIDATED
    })
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$SelectFile <- renderUI({
      fluidPage(
        wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        uiOutput(ns('SelectFile_software_ui')),
        uiOutput(ns('SelectFile_file_ui')),
        uiOutput(ns('SelectFile_ManageXlsFiles_ui')),
        uiOutput(ns('SelectFile_typeOfData_ui')),
        uiOutput(ns('SelectFile_checkDataLogged_ui')),
        uiOutput(ns('SelectFile_replaceAllZeros_ui')),
        
        # Insert validation button
        uiOutput(ns('SelectFile_btn_validate_ui')),
        )
      )
    })
    
    
    # >>> START: Definition of the widgets
    
    
    
    output$SelectFile_software_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_software"), 
                             "Software to import from",
                             choices = setNames(nm = GetSoftAvailables()),
                             selected = rv.widgets$SelectFile_software)
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    # This part must be customized by the developer of a new module
    output$SelectFile_file_ui <- renderUI({
      req(rv.widgets$SelectFile_software)
      fluidRow(
        column(width = 2,
          MagellanNTK::mod_popover_for_help_server(
                 "help_chooseFile",
                 title = "Data file",
                 content = "Select one (.txt, .csv, .tsv, .xls, .xlsx) file."
               ),
          MagellanNTK::mod_popover_for_help_ui(ns("help_chooseFile"))
        ),
        column(width = 10,
               widget <- fileInput(
                 ns("SelectFile_file"), "",
                 multiple = FALSE,
                 accept = c(".txt", ".tsv", ".csv", ".xls", ".xlsx")
               )
        )
      )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    fileExt.ok <- reactive({
      req(rv.widgets$SelectFile_file$name)
      authorizedExts <- c("txt", "csv", "tsv", "xls", "xlsx")
      ext <- GetExtension(rv.widgets$SelectFile_file$name)
      !is.na(match(ext, authorizedExts))
    })
    
    ############ Read text file to be imported ######################
    observeEvent(req(rv.widgets$SelectFile_file), {
      rv.widgets$SelectFile_XLSsheets
      
      ext <- GetExtension(rv.widgets$SelectFile_file$name)
      
      if (((ext %in% c("xls", "xlsx"))) && 
          is.null(rv.widgets$SelectFile_XLSsheets))
        return(NULL)
      
      
      authorizedExts <- c("txt", "csv", "tsv", "xls", "xlsx")
      
      if (!fileExt.ok()) {
        shinyjs::info("Warning : this file is not a text nor an Excel file !
     Please choose another one.")
      } else {
        tryCatch({
          
          shinyjs::disable("SelectFile_file")
          f.path <- rv.widgets$SelectFile_file$datapath
          rv.convert$tab <- switch(ext,
                                   txt = read.csv(f.path, header = TRUE, sep = "\t", as.is = T),
                                   csv = read.csv(f.path, header = TRUE, sep = ";", as.is = T),
                                   tsv = read.csv(f.path, header = TRUE, sep = "\t", as.is = T),
                                   xls = readExcel(f.path, ext, sheet = rv.widgets$SelectFile_XLSsheets),
                                   xlsx = readExcel(f.path, ext, sheet = rv.widgets$SelectFile_XLSsheets)
          )
          
          colnames(rv.convert$tab) <- gsub(".", "_", colnames(rv.convert$tab), fixed = TRUE)
          colnames(rv.convert$tab) <- gsub(" ", "_", colnames(rv.convert$tab), fixed = TRUE)
        },
        warning = function(w) {
          shinyjs::info(conditionMessage(w))
          return(NULL)
        },
        error = function(e) {
          shinyjs::info(conditionMessage(e))
          return(NULL)
        },
        finally = {
          # cleanup-code
        })
      }
    })
    
    output$SelectFile_ManageXlsFiles_ui <- renderUI({
      req(rv.widgets$SelectFile_software)
      req(rv.widgets$SelectFile_file)
      
      req(GetExtension(rv.widgets$SelectFile_file$name) %in% c("xls", "xlsx"))
      
      tryCatch({   
        sheets <- listSheets(rv.widgets$SelectFile_file$datapath)
        widget <- selectInput(ns("SelectFile_XLSsheets"), 
                              "sheets", 
                              choices = as.list(sheets), 
                              width = "200px",
                              selected = rv.widgets$SelectFile_XLSsheets)
      },
      warning = function(w) {
        shinyjs::info(conditionMessage(w))
        return(NULL)
      },
      error = function(e) {
        shinyjs::info(conditionMessage(e))
        return(NULL)
      },
      finally = {
        # cleanup-code
      }
      )
      toggleWidget(widget, rv$steps.enabled['SelectFile'])
    })
    
    
    
    output$SelectFile_typeOfData_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_typeOfData"), 
                             "Is it a peptide or protein dataset ?",
                             choices = c("peptide dataset" = "peptide",
                                         "protein dataset" = "protein"
                             ),
                             selected = rv.widgets$SelectFile_typeOfData)
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    output$SelectFile_checkDataLogged_ui <- renderUI({
      widget <- radioButtons(ns("SelectFile_checkDataLogged"), 
                             "Are your data already log-transformed ?",
                             choices = c("yes (they stay unchanged)" = "yes",
                                         "no (they wil be automatically transformed)" = "no"),
                             selected = rv.widgets$SelectFile_checkDataLogged
      )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    output$SelectFile_replaceAllZeros_ui <- renderUI({
      widget <- checkboxInput(ns("SelectFile_replaceAllZeros"), 
                              "Replace all 0 and NaN by NA",
                              value = rv.widgets$SelectFile_replaceAllZeros
      )
      
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
    })
    
    
    
    output$SelectFile_btn_validate_ui <- renderUI({
      widget <-  actionButton(ns("SelectFile_btn_validate"), "Perform",
                              class = 'info')
      toggleWidget(widget, rv$steps.enabled['SelectFile'] )
      
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$SelectFile_btn_validate, {
      # Do some stuff
      # rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
      #                                     dataset = rnorm(1:5),
      #                                     name = paste0('temp_',id))
      # 
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['SelectFile'] <- stepStatus$VALIDATED
      
    })
    
    
    # <<< END ------------- Code for step 1 UI---------------
    
    
    
    
    # >>> START ------------- Code for step 2 UI---------------
    
    output$DataId <- renderUI({
      wellPanel(
        br(), br(),
        tags$div(
          tags$div(
            # Two examples of widgets in a renderUI() function
            uiOutput(ns('DataId_datasetId_ui')),
            uiOutput(ns('DataId_warningNonUniqueID_ui'))
          ),
          tags$div(
            style = "display:inline-block; vertical-align: top;",
            uiOutput("DataId_parentProteinID_ui"),
            uiOutput("DataId_previewProteinID_ui")
          ),
          
          # Insert validation button
          # This line is necessary. DO NOT MODIFY
          uiOutput(ns('DataId_btn_validate_ui'))
        )
      )
    })
    
    
    output$DataId_datasetId_ui <- renderUI({
      req(rv.convert$tab)
      
      #.choices <- setNames(nm = c("AutoID", colnames(rv.convert$tab)))
      #names(.choices) <- c("Auto ID", colnames(rv.convert$tab))
      
      MagellanNTK::mod_popover_for_help_server("help_convertIdType",
                             title = "ID definition",
                             content = "If you choose the automatic ID, 
                            Prostar will build an index.")
      
      tagList(
        MagellanNTK::mod_popover_for_help_ui(ns("help_convertIdType")),
        widget <- selectInput(ns("DataId_datasetId"), 
                              label = "", 
                              choices = setNames(nm = c("AutoID", colnames(rv.convert$tab))),
                              selected = rv.widgets$DataId_datasetId
        )
      )
      
      toggleWidget(widget, rv$steps.enabled['DataId'] )
    })
    
    
    datasetID_Ok <- reactive({
      req(rv.widgets$DataId_datasetId)
      req(rv.convert$tab)
      if (rv.widgets$DataId_datasetId == "AutoID") {
        t <- TRUE
      } else {
        t <- (length(as.data.frame(rv.convert$tab)[, rv.widgets$DataId_datasetId])
              == length(unique(as.data.frame(rv.convert$tab)[, rv.widgets$DataId_datasetId])))
      }
      t
    })
    
    
    output$DataId_warningNonUniqueID_ui <- renderUI({
      # req(rv.widgets$DataId_datasetId != "AutoID")
      # req(rv.convert$tab)
      # 
      # df <- as.data.frame(rv.convert$tab)
      # t <- (length(df[, rv.widgets$DataId_datasetId]) == length(unique(df[, rv.widgets$DataId_datasetId])))
      # 
      if (!datasetID_Ok()) {
        text <- "<img src=\"images/Problem.png\" height=\"24\"></img><font color=\"red\">
        Warning ! Your ID contains duplicate data. Please choose another one."
      } else {
        text <- "<img src=\"images/Ok.png\" height=\"24\"></img>"
      }
      HTML(text)
    })
    
    
    output$DataId_parentProteinId_ui <- renderUI({
      req(rv.convert$tab)
      req(rv.widgets$SelectFile_typeOfData != "protein")
      
      MagellanNTK::mod_popover_for_help_server("help_ProteinId",
                             title = "Select protein IDs",
                             content = "Select the column containing the parent protein IDs."
      )
      
      tagList(
        MagellanNTK::mod_popover_for_help_ui(ns("help_ProteinId")),
        widget <- selectInput(ns("DataId_parentProteinId"),
                              "",
                              choices = setNames(nm =c("", colnames(rv.convert$tab))),
                              selected = rv.widgets$DataId_parentProteinId)
      )
      toggleWidget(widget, rv$steps.enabled['DataId'] )
    })
    
    
    output$helpTextDataID <- renderUI({
      req(rv.widgets$SelectFile_typeOfData)
      
      t <- switch(rv.widgets$SelectFile_typeOfData,
                  protein = "proteins",
                  peptide = "peptides",
                  default = "")
      
      txt <- paste("Please select among the columns of your data the one that
                corresponds to a unique ID of the ", t, ".", sep = " ")
      helpText(txt)
    })
    
    
    
    
    output$previewProteinID_UI <- renderUI({
      req(rv.widgets$DataId_parentProteinId)
      
      tagList(
        p(style = "color: black;", "Preview"),
        tableOutput("previewProtID")
      )
    })
    
    output$previewProtID <- renderTable(
      head(rv.convert$tab[, rv.widgets$DataId_parentProteinId]), colnames = FALSE
    )
    
    
    
    output$DataId_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("DataId_btn_validate"), "Perform",
                             class = 'info')
      toggleWidget(widget, rv$steps.enabled['DataId'] )
    })
    
    observeEvent(input$DataId_btn_validate, {
      # Do some stuff
      # new.dataset <- 10*rv$dataIn[[length(rv$dataIn)]]
      # rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
      #                                     dataset = new.dataset,
      #                                     name = paste0('temp_',id))
      # 
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['DataId'] <- stepStatus$VALIDATED
    })
    
    # <<< END ------------- Code for step 2 UI---------------
    
    
    # >>> START ------------- Code for step 3 UI---------------
    
    output$ExpandFeatData <- renderUI({
      wellPanel(
        shinyjs::useShinyjs(),
        fluidRow(
          column(width = 4, uiOutput(ns('ExpandFeatData_idMethod_ui'))),
          #column(width = 4, uiOutput(ns("ExpandFeatData_checkIdTab_ui")),
          column(width = 4, shinyjs::hidden(
            div(id = "warning_neg_values",
                p("Warning : Your original dataset may contain negative values",
                  "so that they cannot be logged. Please check back the dataset or",
                  "the log option in the first tab.")
            )
          ))
        ),
        fluidRow(
          column(width = 4, uiOutput(ns("ExpandFeatData_quantCols_ui"), width = "400px")),
          column(width = 8, shinyjs::hidden(
            uiOutput(ns("ExpandFeatData_inputGroup_ui")
                     #, width = "600px"
            )
          ))
        ),
        
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('ExpandFeatData_btn_validate_ui'))
      )
    })
    
    
    output$ExpandFeatData_idMethod_ui <- renderUI({
      widget <- radioButtons(ns("ExpandFeatData_idMethod"), 
                             "Provide identification method",
                             choices = list(
                               "No (default values will be computed)" = FALSE,
                               "Yes" = TRUE),
                             selected = rv.widgets$ExpandFeatData_idMethod)
      toggleWidget(widget, rv$steps.enabled['ExpandFeatData'] )
    })
    
    
    output$ExpandFeatData_inputGroup_ui <- renderUI({
      req(as.logical(rv.widgets$ExpandFeatData_idMethod))
      rv.widgets$ExpandFeatData_quantCols
      
      rv.convert$inputGroup <- DaparToolshed::mod_inputGroup_server('inputGroup',
                            df = rv.convert$tab,
                            quantCols = rv.widgets$ExpandFeatData_quantCols)
      mod_inputGroup_ui(ns('inputGroup'))
    })
    
    output$ExpandFeatData_quantCols_ui <- renderUI({
      req(rv.convert$tab)
      
      MagellanNTK::mod_popover_for_help_server("help_ExpandFeatData_quantCols",
                             title = "Quantitative data",
                             content = "Select the columns that are quantitation values
            by clicking in the field below.")
      
      tagList(
        MagellanNTK::mod_popover_for_help_ui(ns("help_ExpandFeatData_quantCols")),
        widget <- selectInput(ns("ExpandFeatData_quantCols"),
                              label = "",
                              choices = setNames(nm=colnames(rv.convert$tab)),
                              multiple = TRUE, selectize = FALSE ,
                              width = "200px", size = 20,
                              selected = rv.widgets$ExpandFeatData_quantCols
        )
      )
      
      toggleWidget(widget, rv$steps.enabled['ExpandFeatData'] )
    })
    
    
    
    
    observe({
      shinyjs::toggle("warning_neg_values",
                      condition = !is.null(rv.widgets$ExpandFeatData_quantCols) &&
                        length(which(rv.convert$tab[, rv.widgets$ExpandFeatData_quantCols] < 0)) > 0
      )
      shinyjs::toggle("ExpandFeatData_idMethod", condition = !is.null(rv.widgets$ExpandFeatData_quantCols))
      shinyjs::toggle("ExpandFeatData_inputGroup_ui", condition = as.logical(rv.widgets$idMethod) == TRUE)
    })
    
    
    
    output$ExpandFeatData_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("ExpandFeatData_btn_validate"),
                             "Perform",
                             class = 'info')
      toggleWidget(widget, rv$steps.enabled['ExpandFeatData'] )
    })
    
    observeEvent(input$ExpandFeatData_btn_validate, {
      # Do some stuff
      # new.dataset <- 10*rv$dataIn[[length(rv$dataIn)]]
      # rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
      #                                     dataset = new.dataset,
      #                                     name = paste0('temp_',id))
      # 
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['ExpandFeatData'] <- stepStatus$VALIDATED
    })
    
    # <<< END ------------- Code for Design UI---------------
    
    
    ############# STEP 4 ######################
    output$Design <- renderUI({
      
      rv.convert$design <- DaparToolshed::mod_buildDesign_server("designEx", 
        rv.widgets$ExpandFeatData_quantCols)
      
      wellPanel(
        mod_buildDesign_ui(ns("designEx")),
        uiOutput(ns('Design_btn_validate_ui'))
        
      )
    })
    
    
    output$Design_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("Design_btn_validate"), "Validate design",
                             class = 'info')
      toggleWidget(widget, rv$steps.enabled['Design'] )
    })
    
    observeEvent(input$Design_btn_validate, {
      
      req(rv.convert$design)
      # Do some stuff
      # new.dataset <- 10*rv$dataIn[[length(rv$dataIn)]]
      # rv$dataIn <- Add_Datasets_to_Object(object = rv$dataIn,
      #                                     dataset = new.dataset,
      #                                     name = paste0('temp_',id))
      # 
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Design'] <- stepStatus$VALIDATED
    })
    
    # >>> START ------------- Code for Save UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Save_btn_validate_ui')),
        uiOutput(ns('mod_dl_ui'))
      )
    })
    
    output$mod_dl_ui <- renderUI({
      req(config@mode == 'process')
      req(rv$steps.status['Save'] == stepStatus$VALIDATED)
      MagellanNTK::mod_download_dataset_ui(ns('createQuickLink'))
    })
    
    output$Save_btn_validate_ui <- renderUI({
      toggleWidget(actionButton(ns("Save_btn_validate"), 
        "Create QFeatures dataset", 
        class = 'info'),
                   rv$steps.enabled['Save']
      )
    })
    
    
    observeEvent(input$Save_btn_validate, {
      # Create QFeatures dataset file
      rv$dataIn <- createQFeatures(data = rv.convert$tab, 
                                   sample = as.data.frame(rv.convert$design()$design),
                                   indQData = rv.widgets$ExpandFeatData_quantCols,
                                   keyId = rv.widgets$DataId_datasetId,
                                   analysis = "analysis",
                                   indexForMetacell = rv.widgets$ExpandFeatData_inputGroup,
                                   typeDataset = rv.widgets$SelectFile_typeOfData,
                                   parentProtId = rv.widgets$DataId_parentProteinID,
                                   force.na = rv.widgets$SelectFile_replaceAllZeros,
                                   software = rv.widgets$SelectFile_software)
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- stepStatus$VALIDATED
      MagellanNTK::mod_download_dataset_server('createQuickLink', 
                dataIn = reactive({rv$dataIn}),
                extension = c('csv', 'xlsx', 'RData'))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  })
}


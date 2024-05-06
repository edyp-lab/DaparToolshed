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
#' In this example, `PipelineProtein_DA_UI()` and `PipelineProtein_DA_server()` define
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
#' data(Exp1_R25_prot, package = "DaparToolshedData")
#' obj <- Exp1_R25_prot
#' # Simulate imputation of missing values
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' qData <- as.matrix(assay(obj[[2]]))
#' sTab <- MultiAssayExperiment::colData(obj)
#' limma <- limmaCompleteTest(qData, sTab)
#' 
#' new.dataset <- obj[[length(obj)]]
#' HypothesisTest(new.dataset) <- limma
#' obj <- DaparToolshed::addDatasets(obj, new.dataset, 'HypothesisTest')
#' 
#' 
#' path <- system.file('workflow/PipelineProtein', package = 'DaparToolshed')
#' shiny::runApp(workflowApp("PipelineProtein_Filtering", path, dataIn = obj))
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
PipelineProtein_DA_conf <- function(){
  Config(
    fullname = 'PipelineProtein_DA',
    mode = 'process',
    steps = c("Step1", "Step2", "Step3"),
    mandatory = c(FALSE, FALSE, FALSE)
  )
}


#' @rdname PipelineProtein
#' 
#' @export
#'
PipelineProtein_DA_UI <- function(id){
  ns <- NS(id)
}



#' @rdname PipelineProtein
#' 
#' @importFrom stats setNames rnorm
#' 
#' @export
#' 
PipelineProtein_DA_server <- function(id,
  dataIn = reactive({NULL}),
  thlogfc = reactive({0}),
  steps.enabled = reactive({NULL}),
  remoteReset = reactive({FALSE}),
  steps.status = reactive({NULL}),
  current.pos = reactive({1})
  ){
  
  
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    Step1_Comparison = "None",
    DA_Condition1 = "",
    DA_Condition2 = "",
    DA_val_vs_percent = "Value",
    DA_ChooseFilters = "None",
    DA_seuilNA_percent = 0,
    DA_seuilNA = 0,
    DA_filter_th_NA = 0,
    Step2_numericValCalibration = "None",
    Step2_calibrationMethod = NULL,
    DA_numValCalibMethod = 0,
    DA_th_pval = 0,
    DA_type_pval = '-log10()',
    DA_FDR = 0,
    DA_NbSelected = 0,
    DA_nBinsHistpval = 80,
    Step3_viewAdjPval = FALSE
  )
  
  
  rv.custom.default.values <- list(
    comps = NULL,
    nbTotalAnaDiff = NULL,
    nbSelectedAnaDiff = NULL,
    nbSelectedTotal_Step3 = NULL,
    nbSelected_Step3 = NULL,
    conditions = list(cond1 = NULL, cond2 = NULL),
    calibrationRes = NULL,
    errMsgcalibrationPlot = NULL,
    errMsgcalibrationPlotALL = NULL,
    pi0 = NULL,
    filename = NULL
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
        uiOutput(ns('datasetDescription_UI')),
        
        # Insert validation button
        uiOutput(ns('Description_btn_validate_UI'))
      )
    })
    
    output$datasetDescription_UI <- renderUI({
      # Insert your own code to visualize some information
      # about your dataset. It will appear once the 'Start' button
      # has been clicked
      
    })
    
    output$Description_btn_validate_UI <- renderUI({
      widget <- actionButton(ns("Description_btn_validate"),
        "Start",
        class = "btn-success")
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      
      # Find the assay containing the hypothesis tests comparisons
      
      rv.custom$res_AllPairwiseComparisons <- unlist(lapply(seq(length(dataIn())), function(x){
        if(!is.null(HypothesisTest(dataIn()[[x]])))
          HypothesisTest(dataIn()[[x]])
      }))
      
      
      
      rv$dataIn <- dataIn()[[length(dataIn())]]
      
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- stepStatus$VALIDATED
    })
    
    
    
    # >>>
    # >>> START ------------- Code for step 1 UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$Step1 <- renderUI({
      .style <- "display:inline-block; vertical-align: top; padding-right: 60px"
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        tagList(
            tags$div(
              tags$div(style = .style, uiOutput(ns('Step1_Comparison_UI')))
            ),
            #uiOutput(ns("pushpval_UI")),
            tags$hr(),
            tags$div(
              tags$div(style = .style, uiOutput(ns("Step1_volcano_UI"))),
              tags$div(style = .style, uiOutput(ns("Step1_tooltipInfo_UI")))
              )
          ),
        # Insert validation button
        uiOutput(ns("Step1_btn_validate_UI"))
      )
    })
    
    output$Step1_Comparison_UI <- renderUI({
      req(rv.custom$res_AllPairwiseComparisons$logFC)
      ll <- unlist(strsplit(
        colnames(rv.custom$res_AllPairwiseComparisons$logFC), "_logFC"))

      
    widget <- selectInput(ns("Step1_Comparison"), "Select a comparison",
      choices = ll,
      selected = rv.widgets$Step1_Comparison,
      width = "300px"
    )
    MagellanNTK::toggleWidget(widget, 
      rv$steps.status['Step1'] == stepStatus$VALIDATED)
    })
    
    
    
    
    mod_volcanoplot_server(
      id = ns("Step1_volcano"),
      dataIn = reactive({rv$dataIn}),
      data = reactive({xxx}),
      comp = reactive({omXplore::get_group(rv^dataIn)}),
      thlogfc = reactive({0}),
      thpval = reactive({0}),
      tooltip = reactive({NULL}),
      reset = reactive({NULL}),
      is.enabled = reactive({rv$steps.enabled["Step1"]})
    )
    
    output$Step1_volcano_UI <- renderUI({
      widget <- mod_volcanoplot_ui(ns("Step1_volcano"))
      MagellanNTK::toggleWidget(widget, rv$steps.enabled["Step1"])
    })
    
    
    output$Step1_tooltipInfo_UI <- renderUI({
      # req(rv.widgets$Step1_Comparison != "None")
      
      widget <- tagList(
        MagellanNTK::mod_popover_for_help_ui(ns("modulePopover_volcanoTooltip")),
        selectInput(ns("Step1_tooltipInfo"),
          label = NULL,
          choices = colnames(SummarizedExperiment::rowData(rv$dataIn)),
          selected = rv.widgets$Step1_tooltipInfo,
          multiple = TRUE,
          selectize = FALSE,
          width = "300px", size = 5
        ),
        actionButton(ns("Step1_validTooltipInfo"),  "Validate tooltip choice", 
          class = actionBtnClass)
      )
      
      MagellanNTK::toggleWidget(widget, 
        rv$steps.status['Step1'] == stepStatus$VALIDATED)
    })
    
    
    
    MagellanNTK::mod_popover_for_help_server("modulePopover_volcanoTooltip",
      title = "Tooltip",
      content = "Infos to be displayed in the tooltip of volcanoplot"
    )
    
    
    # Fill the variable 'rv.custom$resAnaDiff' with informations relatives to
    # the comparison choosen by the user.
    # Concertely, it extracts data from the variable rv.custom$res_AllPairwiseComparisons
    # which contains all info (logFC and pValue) for all comparisons.
    UpdateCompList <- reactive({
      req(rv.widgets$Step1_Comparison)
      
          index <- which(
            paste(
              as.character(rv.widgets$Step1_Comparison),
              "_logFC",
              sep = ""
            ) == colnames(rv.custom$res_AllPairwiseComparisons$logFC)
          )
          
          # Update of the list rv.custom$resAnaDiff
          rv.custom$Condition1 <- strsplit(as.character(rv.widgets$Step1_Comparison), "_vs_")[[1]][1]
          rv.custom$Condition2 <- strsplit(as.character(rv.widgets$Step1_Comparison), "_vs_")[[1]][2]
          
          rv.custom$resAnaDiff <- list(
            logFC = (rv.custom$res_AllPairwiseComparisons$logFC)[, index],
            P_Value = (rv.custom$res_AllPairwiseComparisons$P_Value)[, index],
            condition1 = rv.custom$Condition1,
            condition2 = rv.custom$Condition2
          )
        }
      )
    
    
    # By default, the tooltip for volcanoplot is set to the proteinId
    observe({
      req(rv$dataIn)
      if (is.null(rv.widgets$Step1_tooltipInfo)) {
        rv.widgets$Step1_tooltipInfo <- omXplore::get_proteinID(rv$dataIn)
      }
    })
    
    
    
    
    
    
    
    
    
    output$pushpval_UI <- renderUI({
      req(rv.widgets$Step1_Comparison != "None")
      
      MagellanNTK::mod_popover_for_help_server("modulePopover_pushPVal",
        title = h3("Push p-value"),
        content = "This functionality is useful in case of multiple pairwise comparisons 
              (more than 2 conditions): At the filtering step, a given analyte X
              (either peptide or protein) may have been kept because it contains
              very few missing values in a given condition (say Cond. A), even
              though it contains (too) many of them in all other conditions
              (say Cond B and C only contains 'MEC' type missing values).
              Thanks to the imputation step, these missing values are no
              longer an issue for the differential analysis, at least from
              the computational viewpoint. However, statistically speaking,
              when performing B vs C, the test will rely on too many imputed
              missing values to derive a meaningful p-value: It may be wiser
              to consider analyte X as non-differentially abundant, regardless
              the test result (and thus, to push its p-value to 1). This is just
              the role of the P-value push parameter. It makes it possible to
              introduce a new filtering step that only applies to each pairwise
              comparison, and which assigns a p-value of 1 to analytes that, for
              the considered comparison are assumed meaningless due to too many
              missing values (before imputation)."
      )
      
      wellPanel(
        MagellanNTK::mod_popover_for_help_ui(ns("modulePopover_pushPVal")),
        mod_query_metacell_UI("AnaDiff_query"))
    })
    
    #---------------------------
    
    # Extract conditions of the current dataset to represent the
    # selected comparison. It returns a subset of the current dataset that will 
    # be used to filter the data within the 'Push p-value' feature
    Get_Dataset_to_Analyze <- reactive({
      
      datasetToAnalyze <- NULL
      if (rv.widgets$Step1_Comparison == "None" || is.null(rv$dataIn))
        return(NULL)
      
      if (length(grep("all-", rv.widgets$Step1_Comparison)) == 1) {
        .conds <- omXplore::get_group(rv$dataIn)$Condition
        condition1 <- strsplit(as.character(rv.widgets$Step1_Comparison), "_vs_")[[1]][1]
        ind_virtual_cond2 <- which(.conds != condition1)
        datasetToAnalyze <- rv$dataIn[[length(rv^dataIn)]]
        Biobase::pData(datasetToAnalyze)$Condition[ind_virtual_cond2] <- "virtual_cond_2"
      } else {
        condition1 <- strsplit(as.character(rv.widgets$Step1_Comparison), "_vs_")[[1]][1]
        condition2 <- strsplit(as.character(rv.widgets$Step1_Comparison), "_vs_")[[1]][2]
        
        if (substr(condition1, 1, 1) == "(" &&
            substr(condition1, nchar(condition1), nchar(condition1)) == ")") {
          condition1 <- sub("^.(.*).$", "\\1", condition1)
        }
        
        if (substr(condition2, 1, 1) == "(" &&
            substr(condition2, nchar(condition2), nchar(condition2)) == ")") {
          condition2 <- sub("^.(.*).$", "\\1", condition2)
        }
        
        
        
        ind <- c(
          which(omXplore::get_group(rv$dataIn) == condition1),
          which(omXplore::get_group(rv$dataIn) == condition2)
        )
        
        datasetToAnalyze <- rv$dataIn[[length(rv$dataIn)]][, ind]
        datasetToAnalyze@experimentData@other$names_metacell <-
          rv$dataIn@experimentData@other$names_metacell[ind]
      }
      
      datasetToAnalyze
    }) %>% bindCache(rv$dataIn, rv.widgets$Step1_Comparison)
    
    
    
    GetPairwiseCompChoice <- reactive({
      req(rv.custom$res_AllPairwiseComparisons$logFC)
      ll <- unlist(strsplit(colnames(rv.custom$res_AllPairwiseComparisons$logFC), "_logFC"))
      ll
    })
    
    GetFiltersScope <- function()
      c("Whole Line" = "WholeLine",
        "Whole matrix" = "WholeMatrix",
        "For every condition" = "AllCond",
        "At least one condition" = "AtLeastOneCond"
      )
    
    
    
    
    observe({
      rv.custom$AnaDiff_indices <- mod_query_metacell_server(id = "AnaDiff_query",
        obj = reactive({req(Get_Dataset_to_Analyze())}),
        reset = reactive({rv_anaDiff$local.reset}),
        op_names = reactive({c('Push p-value', 'Keep original p-value')})
      )
      
    })
    
    
    observeEvent(req(rv.custom$AnaDiff_indices()$indices),{
      # shinyjs::toggleState("AnaDiff_performFilteringMV",
      #                 condition = length(AnaDiff_indices()$indices > 0))
      #  
      UpdateCompList()
      .ind <- rv.custom$AnaDiff_indices()
      params <- rv.custom$AnaDiff_indices()$params
      .protId <- omXplore::get_proteinID(rv$dataIn)

      rv.widgets$Step1_tooltipInfo <- .protId
      
      #--------------------------------
      
      if (!is.null(rv.custom$AnaDiff_indices()$indices) && 
          length(.ind$indices) < nrow(Get_Dataset_to_Analyze())) {
        
        if (rv$widgets$anaDiff$KeepRemove == 'delete')
          indices_to_push <- .ind$indices
        else if (rv$widgets$anaDiff$KeepRemove == 'keep')
          indices_to_push <- seq_len(nrow(Get_Dataset_to_Analyze()))[-(.ind$indices)]
        
        rv.custom$resAnaDiff$P_Value[indices_to_push] <- 1
        n <- length(rv.custom$resAnaDiff$P_Value)
        rv.custom$resAnaDiff$pushed <- seq(n)[indices_to_push]
        
        #shinyjs::toggleState("div_AnaDiff_query", condition = FALSE)
      }
      #shinyjs::toggleState("div_AnaDiff_query", condition = TRUE)
    })
    
    
    
    observeEvent(rv.widgets$Step1_Comparison, ignoreInit = TRUE, {
      rv.widgets$Step1_tooltipInfo <- omXplore::get_proteinID(rv$dataIn[[length(rv$dataIn)]])
      UpdateCompList()
      
      req(rv.widgets$Step1_Comparison)
      cond1 <- rv.custom$Condition1
      cond2 <- rv.custom$Condition2
      
      rv.custom$filename <- paste0("anaDiff_", cond1, "_vs_", cond2, ".xlsx")
    })
    
    
    
    
    
    ## --------------------------------------------------------
    ## ---------------------------------------------------------
    
    
    
    
    
    MagellanNTK::mod_popover_for_help_server("modulePopover_keepLines",
      title = "n values",
      content = "Keep the lines which have at least n intensity values."
    )
    
    
    
    GetBackToCurrentResAnaDiff <- reactive({
      req(rv.custom$res_AllPairwiseComparisons)
      
      index <- which(
        paste(
          as.character(rv.widgets$Step1_Comparison), "_logFC",
          sep = ""
        ) ==
          colnames(rv.custom$res_AllPairwiseComparisons$logFC)
      )
      rv.custom$resAnaDiff <- list(
        logFC = (rv.custom$res_AllPairwiseComparisons$logFC)[, index],
        P_Value = (rv.custom$res_AllPairwiseComparisons$P_Value)[, index],
        condition1 = strsplit(
          as.character(rv.widgets$Step1_Comparison), "_vs_"
        )[[1]][1],
        condition2 = strsplit(
          as.character(rv.widgets$Step1_Comparison), "_vs_"
        )[[1]][2]
      )
      rv.custom$resAnaDiff
    })
    
    
    
    
    
    not_a_numeric <- function(input) {
      if (is.na(as.numeric(input))) {
        "Please input a number"
      } else {
        NULL
      }
    }
    
    
    
    
    
    
    
    
    
    output$Step1_btn_validate_UI <- renderUI({
      widget <- actionButton(ns("Step1_btn_validate"),
        "Validate step",
        class = "btn-success"
      )
       MagellanNTK::toggleWidget(widget,  rv$steps.enabled["Step1"])
    })
    # >>> END: Definition of the widgets
    
    
    
    observeEvent(input$Step1_btn_validate, {
      
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- NULL
      rv$steps.status["Step1"] <- stepStatus$VALIDATED
    })
    

    # <<< END ------------- Code for step 1 UI---------------
    
    
    # >>> START ------------- Code for step 2 UI---------------
    output$Step2 <- renderUI({
      
      .style <- "display:inline-block; 
        vertical-align: middle; 
        padding-right: 40px;"
      
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        
        tagList(
          tags$div(
            tags$div(style = .style,
              uiOutput(ns('Step2_calibrationMethod_UI'))
            ),
            tags$div(style = .style,
              uiOutput(ns("Step2_numericValCalibration_UI"))
            ),
            tags$div(style = .style,
              uiOutput(ns("nBins_UI"))
            ),
            tags$div(style = .style,
              p(tags$strong(
                paste0("value of pi0: ", round(as.numeric(rv.custom$pi0), digits = 2))
                ))
            )
          ),
          tags$hr(),
          fluidRow(
            column(width = 6, fluidRow(style = "height:800px;",
              imageOutput("calibrationPlotAll", height = "800px")
            )),
            column(width = 6, fluidRow(style = "height:400px;",
                imageOutput("calibrationPlot", height = "400px")
              ),
              fluidRow(style = "height:400px;", 
                highchartOutput("histPValue"))
            )
          )
        ),

        # Insert validation button
        uiOutput(ns("Step2_btn_validate_UI"))
      )
    })
    
    
    
    output$Step2_calibrationMethod_UI <- renderUI({
      widget <- selectInput("Step2_calibrationMethod", "Calibration method",
        choices = c("None" = "None", calibMethod_Choices),
        selected = rv.widgets$Step2_calibrationMethod,
        width = "200px"
      )
    MagellanNTK::toggleWidget(widget,  rv$steps.enabled["Step2"])
    })
    
    
    
    output$Step2_numericValCalibration_UI <- renderUI({
      req(rv.widgets$Step2_calibrationMethod == "numeric value")
      widget <- numericInput(ns("numericValCalibration"),
        "Proportion of TRUE null hypohtesis",
        value = rv.widgets$Step2_numericValCalibration,
        step = 0.05
      )
      MagellanNTK::toggleWidget(widget,  rv$steps.enabled["Step2"])
  })
    

    
    observeEvent(rv.widgets$Step2_calibrationMethod, {
      shinyjs::toggle("numericValCalibration",
        condition = rv.widgets$Step2_calibrationMethod == "numeric value"
      )
    })
    

    
    output$Step2_nBins_UI <- renderUI({
      req(rv.custom$resAnaDiff)
      req(rv.custom$pi0)
      req(rv.widgets$Step2_nBinsHistpval)
      
      widget <- selectInput(
        ns("Step2_nBinsHistpval"), 
        "n bins of p-value histogram",
        choices = c(1, seq(from = 0, to = 100, by = 10)[-1]),
        selected = rv.widgets$Step2_nBinsHistpval, 
        width = "80px")
      MagellanNTK::toggleWidget(widget, rv$steps.enabled["Step2"])
    })
    
    
    histPValue <- reactive({
      req(rv.custom$resAnaDiff)
      req(rv.cutoms$pi0)
      req(rv.widgets$Step2_nBinsHistpval)
      req(data()$logFC)
      req(!is.(data()$logFC))
      req(length(data()$logFC) > 0)
      
     
      m <- match.metacell(DAPAR::GetMetacell(rv$dataIn),
        pattern = c("Missing", "Missing POV", "Missing MEC"),
        level = "peptide"
      )
      if (length(which(m)) > 0) {
        return()
      }
      
      # isolate({
      t <- NULL
      method <- NULL
      t <- data()$P_Value
      t <- t[which(abs(data()$logFC) >= thlogfc())]
      toDelete <- which(t == 1)
      if (length(toDelete) > 0) {
        t <- t[-toDelete]
      }
      histPValue_HC(t,
        bins = as.numeric(rv.widgets$Step2_nBinsHistpval),
        pi0 = rv.custom$pi0
      )
      # })
    })
    
    output$histPValue <- renderHighchart({
      histPValue()
    })
    
    
    
    
    
    
    output$calibrationResults <- renderUI({
      req(rv.custom$calibrationRes)
      rv$dataIn
      
      
      txt <- paste("Non-DA protein proportion = ",
        round(100 * rv.cutom$calibrationRes$pi0, digits = 2), "%<br>",
        "DA protein concentration = ",
        round(100 * rv.cutom$calibrationRes$h1.concentration, digits = 2),
        "%<br>",
        "Uniformity underestimation = ",
        rv.cutom$calibrationRes$unif.under, "<br><br><hr>",
        sep = ""
      )
      HTML(txt)
    })
    
    
    
    
    calibrationPlot <- reactive({
      req(rv.widgets$Step2_calibrationMethod != "None")
      rv.custom$resAnaDiff
      req(rv$dataIn)
      
      if (length(rv.custom$resAnaDiff$logFC) == 0) {
        return()
      }
      
      m <- match.metacell(omXplore::get_metacell(rv$dataIn),
        pattern = c("Missing", "Missing POV", "Missing MEC"),
        level = "peptide")
      if (length(which(m)) > 0) {
        return()
      }
      
      #cond <- c(rv.custom$resAnaDiff$condition1, rv.custom$resAnaDiff$condition2)
      
      
      t <- NULL
      method <- NULL
      t <- rv.custom$resAnaDiff$P_Value
      t <- t[which(abs(rv.custom$resAnaDiff$logFC) >= logfc())]
      toDelete <- which(t == 1)
      if (length(toDelete) > 0) {
        t <- t[-toDelete]
      }
      
      l <- NULL
      ll <- NULL
      result <- tryCatch(
        {
          if ((rv.widgets$Step2_calibrationMethod == "numeric value") &&
              !is.null(rv.widgets$Step2_numValCalibMethod)) {
            ll <- catchToList(
              wrapperCalibrationPlot(
                t,
                rv.widgets$Step2_numValCalibMethod
              )
            )
            .warns <- ll$warnings[grep("Warning:", ll$warnings)]
            rv.custom$errMsgCalibrationPlot <- .warns
          } else if (rv.widgets$Step2_calibrationMethod == "Benjamini-Hochberg") {
            ll <- catchToList(wrapperCalibrationPlot(t, 1))
            .warns <- ll$warnings[grep("Warning:", ll$warnings)]
            rv.custom$errMsgCalibrationPlot <- .warns
          } else {
            ll <- catchToList(
              wrapperCalibrationPlot(t, rv.widgets$Step2_calibrationMethod)
            )
            .warns <- ll$warnings[grep("Warning:", ll$warnings)]
            rv.custom$errMsgCalibrationPlot <- .warns
          }
          rv.custom$pi0 <- ll$value$pi0
          
          rvModProcess$moduleAnaDiffDone[2] <- !is.null(rv.custom$pi0)
        },
        warning = function(w) {
          shinyjs::info(paste("Calibration plot", ":",
            conditionMessage(w),
            sep = " "
          ))
        },
        error = function(e) {
          shinyjs::info(paste("Calibration plot", ":",
            conditionMessage(e),
            sep = " "
          ))
        },
        finally = {
          # cleanup-code
        }
      )
    })
    
    output$calibrationPlot <- renderImage(
      {
        outfile <- tempfile(fileext = ".png")
        
        # Generate a png
        png(outfile, width = 600, height = 500)
        calibrationPlot()
        dev.off()
        
        # Return a list
        list(
          src = outfile,
          alt = "This is alternate text"
        )
      },
      deleteFile = TRUE
    )
    
    
    
    
    output$errMsgCalibrationPlot <- renderUI({
      req(rv.custom$errMsgCalibrationPlot)
      req(rv$dataIn)
      
      txt <- NULL
      
      for (i in 1:length(rv.custom$errMsgCalibrationPlot)) {
        txt <- paste(txt, "errMsgCalibrationPlot: ",
          rv.custom$errMsgCalibrationPlot[i], "<br>",
          sep = ""
        )
      }
      
      div(HTML(txt), style = "color:red")
    })
    
    
    output$errMsgCalibrationPlotAll <- renderUI({
      rv.custom$errMsgCalibrationPlotAll
      req(rv$dataIn)
      if (is.null(rv.custom$errMsgCalibrationPlotAll)) {
        return()
      }
      
      txt <- NULL
      for (i in 1:length(rv.custom$errMsgCalibrationPlotAll)) {
        txt <- paste(txt, "errMsgCalibrationPlotAll:",
          rv.custom$errMsgCalibrationPlotAll[i], "<br>",
          sep = ""
        )
      }
      
      div(HTML(txt), style = "color:red")
    })
    
    
    
    calibrationPlotAll <- reactive({
      rv.custom$resAnaDiff
      req(rv$dataIn)
      req(is.na(logfc()))
      req((length(rv.custom$resAnaDiff$logFC) == 0)) 
      
      m <- match.metacell(omXplore::get_metacell(rv$dataIn),
        pattern = c("Missing", "Missing POV", "Missing MEC"),
        level = "peptide")
      req(length(which(m)) == 0)
      
      cond <- c(rv.custom$resAnaDiff$condition1, rv.custom$resAnaDiff$condition2)
      
      t <- NULL
      method <- NULL
      t <- rv.custom$resAnaDiff$P_Value
      t <- t[which(abs(rv.custom$resAnaDiff$logFC) >= thlogfc())]
      toDelete <- which(t == 1)
      if (length(toDelete) > 0) {
        t <- t[-toDelete]
      }
      
      l <- NULL
      result <- tryCatch(
        {
          l <- catchToList(wrapperCalibrationPlot(t, "ALL"))
          .warns <- l$warnings[grep("Warning:", l$warnings)]
          rv.custom$errMsgCalibrationPlotAll <- .warns
          rvModProcess$moduleAnaDiffDone[2] <- !is.null(rv.custom$pi0)
        },
        warning = function(w) {
          shinyjs::info(paste("Calibration Plot All methods", ":",
            conditionMessage(w),
            sep = " "
          ))
        },
        error = function(e) {
          shinyjs::info(paste("Calibration Plot All methods", ":",
            conditionMessage(e),
            sep = " "
          ))
        },
        finally = {
          # cleanup-code
        }
      )
    })
    
    
    
    #--------------------------------------------------
    output$calibrationPlotAll <- renderImage(
      {
        outfile <- tempfile(fileext = ".png")
        
        # Generate a png
        png(outfile, width = 600, height = 500)
        calibrationPlotAll()
        dev.off()
        
        # Return a list
        list(
          src = outfile,
          alt = "This is alternate text"
        )
      },
      deleteFile = TRUE
    )
    
    
    
    
    
    
    
    output$Step2_btn_validate_UI <- renderUI({
      widget <- actionButton(ns("Step2_btn_validate"),
        "Validate step",
        class = "btn-success"
      )
      MagellanNTK::toggleWidget(widget, rv$steps.enabled["Step2"])
    })
    # >>> END: Definition of the widgets
    
    
    
    observeEvent(input$Step2_btn_validate, {
      
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- NULL
      rv$steps.status["Step2"] <- stepStatus$VALIDATED
    })
    
    # <<< END ------------- Code for step 2 UI---------------
    
    
    # >>> START ------------- Code for step 2 UI---------------
    output$Step3 <- renderUI({
      widget <- wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUI() function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        tagList(
          fluidRow(
            column(width = 5,
              mod_set_pval_threshold_UI("Title"),
              uiOutput(ns("nbSelectedItems")),
              actionButton(ns('validate_pval'), "Validate threshold", class = actionBtnClass)
            ),
            column(width = 7,
              withProgress(message = "", detail = "", value = 1, {
                uiOutput(ns('Step3_volcanoplot_UI'))
              })
            )
          ),
          tags$hr(),
          fluidRow(
            column(width = 4,
              checkboxInput(ns('Step3_viewAdjPval'), 
                'View adjusted p-value', 
                value = rv.widgets$Step3_viewAdjPval)
            ),
            column(width = 4,
              downloadButton(ns("Step3_download_SelectedItems_UI"), 
                "Download (Excel file)", class = actionBtnClass)
            )
          ),
          DT::DTOutput(ns("Step3_anaDiff_selectedItems"))
        ),
        # Insert validation button
        uiOutput(ns("Step3_btn_validate_UI"))
      )
      
      MagellanNTK::toggleWidget(widget, 
        rv$steps.enabled["Step3"] && 
          as.character(rv.widgets$Step1_Comparison) != "None")
    })
    
    
    
    #-------------------------------------------------------------------
    #
    mod_volcanoplot_server(
      id = ns("Step3_volcano"),
      dataIn = reactive({rv$dataIn}),
      data = reactive({xxx}),
      comp = reactive({omXplore::get_group(rv^dataIn)}),
      thlogfc = reactive({0}),
      thpval = reactive({0}),
      tooltip = reactive({NULL}),
      reset = reactive({NULL}),
      is.enabled = reactive({rv$steps.enabled["Step3"]})
    )
    
    output$Step3_volcanoplot_UI <- renderUI({
      widget <- mod_volcanoplot_ui(ns("Step3_volcano"))
      MagellanNTK::toggleWidget(widget, rv$steps.enabled["Step3"])
    })
    
    
    output$nbSelectedItems <- renderUI({
      rv.widgets$Step1_thpval
      rv$dataIn
      req(Build_pval_table())
      
      
      m <- match.metacell(omXplore::get_metacell(rv$dataIn),
        pattern = c("Missing", "Missing POV", "Missing MEC"),
        level = "peptide"
      )
      #req(length(which(m)) > 0)
      
      p <- Build_pval_table()
      upItemsPVal <- NULL
      upItemsLogFC <- NULL
      
      
      upItemsLogFC <- which(abs(p$logFC) >= as.numeric(logfc()))
      upItemsPVal <- which(-log10(p$P_Value) >= as.numeric(
        rv.widgets$Step1_thpval
      ))
      
      rv.custom$nbTotalAnaDiff <- nrow(SummarizedExperiment::assay(rv$dataIn))
      rv.custom$nbSelectedAnaDiff <- NULL
      t <- NULL
      
      if (!is.null(rv.widgets$Step1_thpval) && !is.null(logfc())) {
        t <- intersect(upItemsPVal, upItemsLogFC)
      } else if (!is.null(rv.widgets$Step1_thpval) && is.null(logfc())) {
        t <- upItemsPVal
      } else if (is.null(rv.widgets$Step1_thpval) && !is.null(logfc())) {
        t <- upItemsLogFC
      }
      rv.custom$nbSelectedAnaDiff <- length(t)
      
      
      ##
      ## Condition: A = C + D
      ##
      A <- rv.custom$nbTotalAnaDiff
      B <- A - length(rv.custom$resAnaDiff$pushed)
      C <- rv.custom$nbSelectedAnaDiff
      D <- ( A - C)
      # 
      # txt <- paste("Total number of ", rv$typeOfDataset, "(s) = ", A , "<br>",
      #   "\t <em>Total remaining after push p-values = ", B , "</em><br>",
      #     paste("Number of selected ", rv$typeOfDataset, "(s) = ", C, sep=''),
      #     paste("Number of non selected ", rv$typeOfDataset, "(s) = ", D, sep = ''),
      #     sep = ""
      # )
      
      div(id="bloc_page",
        style = "background-color: lightgrey; width: 300px",
        p(paste("Total number of ", omXplore::get_type(rv.dataIn), "(s) = ", A, sep = '' )),
        tags$em(p(style = "padding:0 0 0 20px;", paste("Total remaining after push p-values = ", B, sep=''))),
        p(paste("Number of selected ", omXplore::get_type(rv.dataIn), "(s) = ", C, sep = '')),
        p(paste("Number of non selected ", omXplore::get_type(rv.dataIn), "(s) = ", D, sep = ''))
      )
      #HTML(txt)
    })
    
    
    
    #################################################################
    ###### Set code for widgets managment
    ################################################################
    
    logpval <- mod_set_pval_threshold_server(id = "Title",
      pval_init = reactive({10^(-rv.widgets$Step1_thpval)}),
      fdr = reactive({Get_FDR()}))
    
    
    observeEvent(logpval(), {
      req(logpval())
      tmp <- gsub(",", ".", logpval(), fixed = TRUE)
      
      rv.widgets$Step1_thpval <- as.numeric(tmp)
      
    })
    
    
    observeEvent(input$Step1_validTooltipInfo, {
      
      .tmp <- c(omXplore__get_proteinID(proteinId), 
        rv.widgets$Step1_tooltipInfo)
      rv.widgets$Step1_tooltipInfo <- unique(.tmp)
      
    })
    
    
    output$Step3_selectedItems_UI <- DT::renderDT({
      df <- Build_pval_table()
      
      if (rv.widgets$Step3_viewAdjPval){
        df <- df[order(df$Adjusted_PValue, decreasing=FALSE), ]
        .coldefs <- list(list(width = "200px", targets = "_all"))
      } else {
        name <- paste0(c('Log_PValue (', 'Adjusted_PValue ('), 
          as.character(rv.widgets$Step1_Comparison), ")")
        .coldefs <- list(
          list(width = "200px", targets = "_all"),
          list(targets = (match(name, colnames(df)) - 1), visible = FALSE))
      }
      
      
      DT::datatable(df,
        escape = FALSE,
        rownames = FALSE,
        selection = 'none',
        options = list(initComplete = initComplete(),
          dom = "frtip",
          pageLength = 100,
          scrollY = 500,
          scroller = TRUE,
          server = FALSE,
          columnDefs = .coldefs,
          ordering = !rv.widgets$Step3_viewAdjPval
        )
      ) %>%
        DT::formatStyle(
          paste0("isDifferential (",
            as.character(rv.widgets$Step1_Comparison), ")"),
          target = "row",
          backgroundColor = DT::styleEqual(c(0, 1), c("white", orangeProstar))
        )
      
    })
    
    
    output$Step3_download_SelectedItems_UI <- downloadHandler(
      
      
      filename = function() {rv.custom$filename},
      content = function(fname) {
        DA_Style <- openxlsx::createStyle(fgFill = orangeProstar)
        hs1 <- openxlsx::createStyle(fgFill = "#DCE6F1",
          halign = "CENTER",
          textDecoration = "italic",
          border = "Bottom")
        
        wb <- openxlsx::createWorkbook() # Create wb in R
        openxlsx::addWorksheet(wb, sheetName = "DA result") # create sheet
        openxlsx::writeData(wb,
          sheet = 1,
          as.character(rv.widgets$Step1_Comparison),
          colNames = TRUE,
          headerStyle = hs1
        )
        openxlsx::writeData(wb,
          sheet = 1,
          startRow = 3,
          Build_pval_table(),
        )
        
        .txt <- paste0("isDifferential (",
          as.character(rv.widgets$Step1_Comparison),
          ")")
        
        ll.DA.row <- which(Build_pval_table()[, .txt] == 1)
        ll.DA.col <- rep(which(colnames(Build_pval_table()) == .txt),
          length(ll.DA.row) )
        
        openxlsx::addStyle(wb,
          sheet = 1, 
          cols = ll.DA.col,
          rows = 3 + ll.DA.row, 
          style = DA_Style
        )
        
        openxlsx::saveWorkbook(wb, file = fname, overwrite = TRUE)
      }
    )

    
    
    Get_FDR <- reactive({
      req(rv.widgets$Step1_thpval)
      req(Build_pval_table())
      
      adj.pval <- Build_pval_table()$Adjusted_PValue
      logpval <- Build_pval_table()$Log_PValue
      upitems_logpval <- which(logpval >= rv.widgets$Step1_thpval)
      
      fdr <- max(adj.pval[upitems_logpval], na.rm = TRUE)
      rv.custom$FDR <- as.numeric(fdr)
      as.numeric(fdr)
    })
    
    observeEvent(input$validate_pval,{
      rvModProcess$moduleAnaDiffDone[3] <- TRUE
    })
    
    Get_Nb_Significant <- reactive({
      nb <- length(
        which(
          Build_pval_table()[paste0(
            "isDifferential (",
            as.character(rv.widgets$Step1_Comparison), ")"
          )] == 1
        )
      )
      rv$widgets$anaDiff$NbSelected <- nb
      nb
    })
    
    observe({
      req(Get_FDR())
      req(Get_Nb_Significant())
      
      th <- Get_FDR() * Get_Nb_Significant()
      
      if (th < 1) {
        warntxt <- paste0("With such a dataset size (",
          Get_Nb_Significant(), " selected discoveries), an FDR of ",
          round(100 * Get_FDR(), digits = 2),
          "% should be cautiously interpreted as strictly less than one
        discovery (", round(th, digits = 2), ") is expected to be false"
        )
        mod_errorModal_server('warn_FDR',
          title = 'Warning',
          text = warntxt)
      }
      
    })
    
    
    
    
    GetCalibrationMethod <- reactive({
      req(rv.widgets$Step2_numValCalibMethod)
      req(rv.widgets$Step2_calibrationMethod != 'None')
      .calibMethod <- NULL
      if (rv.widgets$Step2_calibrationMethod == "Benjamini-Hochberg") {
        .calibMethod <- 1
      } else if (rv.widgets$Step2_calibrationMethod == "numeric value") {
        .calibMethod <- as.numeric(rv.widgets$Step2_numValCalibMethod)
      } else {
        .calibMethod <- rv.widgets$Step2_calibrationMethod
      }
      .calibMethod
      
    })
    
    
    Build_pval_table <- reactive({
      req(rv.custom$resAnaDiff$logFC)
      req(rv.custom$resAnaDiff$P_Value)
      req(rv$dataIn)
      req(GetCalibrationMethod())
      rv.widgets$Step1_thpval
      
      .digits <- 3
      
      pval_table <- data.frame(
        id = rownames(SummarizedExperiment::assay(rv$dataIn)),
        logFC = round(rv.custom$resAnaDiff$logFC, digits = .digits),
        P_Value = rv.custom$resAnaDiff$P_Value,
        Log_PValue = -log10(rv.custom$resAnaDiff$P_Value),
        Adjusted_PValue = rep(NA, length(rv.custom$resAnaDiff$logFC)),
        isDifferential = rep(0, length(rv.custom$resAnaDiff$logFC))
      )
      
      thpval <- rv.widgets$Step1_thpval
      
      #
      # Determine significant proteins
      signifItems <- intersect(which(pval_table$Log_PValue >= thpval),
        which(abs(pval_table$logFC) >= thlogfc())
      )
      pval_table[signifItems,'isDifferential'] <- 1
      
      upItems_pval <- which(-log10(rv.custom$resAnaDiff$P_Value) >= thpval)
      upItems_logFC <- which(abs(rv.custom$resAnaDiff$logFC) >= thlogfc())
      rv.custom$adjusted_pvalues <- diffAnaComputeAdjustedPValues(
        rv.custom$resAnaDiff$P_Value[upItems_logFC],
        GetCalibrationMethod())
      pval_table[upItems_logFC, 'Adjusted_PValue'] <- rv.custom$adjusted_pvalues
      
      
      
      # Set only significant values
      pval_table$logFC <- signif(pval_table$logFC, digits = 4)
      pval_table$P_Value <- signif(pval_table$P_Value, digits = 4)
      pval_table$Adjusted_PValue <- signif(pval_table$Adjusted_PValue, digits = 4)
      pval_table$Log_PValue <- signif(pval_table$Log_PValue, digits = 4)
      
      
      
      tmp <- as.data.frame(
        SummarizedExperiment::rowData(rv$dataIn)[, rv.widgets$Step1_tooltipInfo]
      )
      names(tmp) <- rv.widgets$Step1_tooltipInfo
      pval_table <- cbind(pval_table, tmp)
      
      colnames(pval_table)[2:6] <- paste0(colnames(pval_table)[2:6], " (", as.character(rv.widgets$Step1_Comparison), ")")
      
      pval_table
    })
    
    
    
    isContainedIn <- function(strA, strB) {
      return(all(strA %in% strB))
    }
    
    #-------------------------------------------------------------------
    
    output$Step3_btn_validate_UI <- renderUI({
      widget <- actionButton(ns("Step3_btn_validate"),
        "Validate step",
        class = "btn-success"
      )
      MagellanNTK::toggleWidget(widget, 
        rv$steps.enabled["Step3"])
    })
    # >>> END: Definition of the widgets
    
    
    
    observeEvent(input$Step3_btn_validate, {
      
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- NULL
      rv$steps.status["Step3"] <- stepStatus$VALIDATED
    })
    
    # <<< END ------------- Code for step 2 UI---------------
    
    
    
    # >>> START ------------- Code for step 'Save' UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Save_btn_validate_UI')),
        uiOutput(ns('dl_UI'))
      )
    })
    
    output$dl_UI <- renderUI({
      req(rv$steps.status['Save'] == stepStatus$VALIDATED)
      req(config@mode == 'process')
      
      MagellanNTK::mod_download_dataset_UI(ns('createQuickLink'))
    })
    
    output$Save_btn_validate_UI <- renderUI({
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
      
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv.custom$tmp
      rv$steps.status['Save'] <- stepStatus$VALIDATED
      
      
      MagellanNTK::mod_download_dataset_server('createQuickLink', 
        dataIn = reactive({rv.custom$tmp}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

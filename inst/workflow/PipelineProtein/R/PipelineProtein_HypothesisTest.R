#' @title Shiny example process module.
#'
#' @description
#' This module contains the configuration information for the corresponding pipeline.
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
#' In this example, `PipelineProtein_HypothesisTest_ulength(rv$dataIn)` and `PipelineProtein_HypothesisTest_server()` define
#' the code for the process `ProcessProtein` which is part of the pipeline called `PipelineProtein`.


#' @rdname PipelineProtein
#' @export
#' 
PipelineProtein_HypothesisTest_conf <- function(){
  Config(
    fullname = 'PipelineProtein_HypothesisTest',
    mode = 'process',
    steps = c('HypothesisTest'),
    mandatory = c(FALSE)
  )
}


#' @param id xxx
#' 
#' @rdname PipelineProtein
#' 
#' @author Samuel Wieczorek
#' 
#' @export
#'
PipelineProtein_HypothesisTest_ui <- function(id){
  ns <- NS(id)
}


#' @param id xxx
#'
#' @param dataIn An instance of the class 
#'
#' @param steps.enabled A vector of boolean which has the same length of the steps
#' of the pipeline. This information is used to enable/disable the widgets. It is not
#' a communication variable between the caller and this module, thus there is no
#' corresponding output variable
#'
#' @param remoteReset It is a remote command to reset the module. A boolean that
#' indicates if the pipeline has been reset by a program of higher level
#' Basically, it is the program which has called this module
#' 
#' @param steps.status xxx
#' 
#' @param current.pos xxx
#'
#' @rdname PipelineProtein
#' 
#' @importFrom stats setNames rnorm
#' 
#' @export
#' 
PipelineProtein_HypothesisTest_server <- function(id,
  dataIn = reactive({NULL}),
  steps.enabled = reactive({NULL}),
  remoteReset = reactive({FALSE}),
  steps.status = reactive({NULL}),
  current.pos = reactive({1})
){
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    HypothesisTest_design = "None",
    HypothesisTest_method = "None",
    HypothesisTest_ttestOptions = 'Student',
    HypothesisTest_thlogFC = 0
  )
  
  rv.custom.default.values <- list(
    listNomsComparaison = NULL,
    n = NULL,
    swap.history = NULL,
    AllPairwiseComp = NULL
  )
  
  ###-------------------------------------------------------------###
  ###                                                             ###
  ### ------------------- MODULE SERVER --------------------------###
  ###                                                             ###
  ###-------------------------------------------------------------###
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    requireNamespace('DaparToolshed')
    
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
        # In this example, the md file is found in the extdata/module_examples directory
        # but with a real app, it should be provided by the package which
        # contains the UI for the different steps of the process module.
        # Insert validation button
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
      widget <- actionButton(ns("Description_btn_validate"),
        "Start",
        class = "btn-success")
      toggleWidget(widget, rv$steps.enabled['Description'])
    })
    
    
    observeEvent(input$Description_btn_validate, {
      rv$dataIn <- dataIn()
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Description'] <- stepStatus$VALIDATED
    })
    
    
    
    # >>>
    # >>> START ------------- Code for HypothesisTest UI---------------
    # >>> 
    
    # >>>> -------------------- STEP 1 : Global UI ------------------------------------
    output$HypothesisTest <- renderUI({
      shinyjs::useShinyjs()

      .style <- "display:inline-block; vertical-align: middle; 
      padding-right: 20px;"
      m <- DaparToolshed::match.metacell(
        omXplore::get_metacell(rv$dataIn[[length(rv$dataIn)]]),
        pattern = c("Missing", "Missing POV", "Missing MEC"),
        level = omXplore::get_type(rv$dataIn[[length(rv$dataIn)]])
      )
      NA.count <- length(which(m))
      
      wellPanel(
        # uiOutput for all widgets in this UI
        # This part is mandatory
        # The renderUlength(rv$dataIn) function of each widget is managed by MagellanNTK
        # The dev only have to define a reactive() function for each
        # widget he want to insert
        # Be aware of the naming convention for ids in uiOutput()
        # For more details, please refer to the dev document.
        
        if (NA.count > 0) {
          tags$p("Your dataset contains missing values. Before using the
      Hypothesis test, you must filter/impute them.")
        } else {
        tagList(
          div(
            div(style = .style,
              uiOutput(ns('HypothesisTest_warning_conditions_ui'))
            ),
            div(style = .style,
              uiOutput(ns('HypothesisTest_design_ui'))
            ),
            
            div(style = .style,
              uiOutput(ns('HypothesisTest_method_ui'))
            ),
            
            div(style = .style,
              uiOutput(ns('HypothesisTest_ttestOptions_ui'))
            ),
            
            div(style = .style,
              uiOutput(ns('HypothesisTest_thlogFC_ui'))
            ),
            
            tags$div(style = .style,
              uiOutput(ns("HypothesisTest_correspondingRatio_ui"))
            ),
            
            div(style = .style,
              uiOutput(ns('HypothesisTest_info_Limma_disabled_ui')),
            ),
            tags$hr(),
            div(style = .style,
              uiOutput(ns("HypothesisTest_swapConds_ui"))
              )
            ),
          div(style = .style,
            uiOutput(ns("HypothesisTest_perform_btn_ui"))
          ),
          
          
          
            div(style = .style,
              uiOutput(ns("HypothesisTest_btn_validate_ui"))
            ),
          tags$hr()
          ,div(style = .style,
            highchartOutput(ns("FoldChangePlot"), height = "100%")
          )

        )
        }
      )
    })
    
    
    
    output$HypothesisTest_warning_conditions_ui <- renderUI({
      req(length(unique(omXplore::get_group(rv$dataIn))) > 26)
      req(getDesignLevel(MultiAssayExperiment::colData(rv$dataIn)) > 1)
      h3('Limma with this version of Prostar does not handle datasets with 
      more than 26 conditions. Such, the Limma option is desactivated for the 
        current dataset')
    })
    # >>> START: Definition of the widgets
    
    
    output$HypothesisTest_design_ui <- renderUI({
      widget <- selectInput(ns("HypothesisTest_design"), "Contrast",
        choices = c("None" = "None", 
          "One vs One" = "OnevsOne",
          "One vs All" = "OnevsAll"),
        selected = rv.widgets$HypothesisTest_design,
        width = "150px"
      )
      MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'] )
    })
    
    output$HypothesisTest_method_ui <- renderUI({
      .methods <- c("None" = "None", "t-tests" = "ttests")
      if(enable_Limma())
        .methods <- c(.methods, "Limma" = "Limma")
      
      widget <- selectInput(ns("HypothesisTest_method"), "Statistical test",
          choices = .methods,
          selected = rv.widgets$HypothesisTest_method,
          width = "150px"
        )
      MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'])
    })
    
    
    
    output$HypothesisTest_ttestOptions_ui <- renderUI({
      req(rv.widgets$HypothesisTest_method == "ttests")
      widget <- radioButtons(ns("HypothesisTest_ttestOptions"), 
        "t-tests options",
      choices = c("Student", "Welch"),
      selected = rv.widgets$HypothesisTest_ttestOptions,
      width = "150px"
    )
      MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'])
    })
    
    
    
    output$HypothesisTest_thlogFC_ui <- renderUI({
      widget <- textInput(ns("HypothesisTest_thlogFC"),
        "log(FC) threshold",
        value = rv.widgets$HypothesisTest_thlogFC,
        width = "150px"
      )
      MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'])
    })
    
    
    
    output$HypothesisTest_correspondingRatio_ui <- renderUI({
      ratio <- as.numeric(rv.widgets$HypothesisTest_thlogFC)
      p("(FC = ", 2^(ratio), ")")
    })
    
 
    
    output$HypothesisTest_perform_btn_ui <- renderUI({
        widget <- actionButton(ns("HypothesisTest_perform_btn"),
            "Compute log FC plot",
            class = actionBtnClass
          )
    MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'])
    })
    
    
    output$FoldChangePlot <- renderHighchart({
      req(rv.custom$AllPairwiseComp)
      withProgress(message = "Computing plot...", detail = "", value = 0.5, {
          tmp.df <- as.data.frame(rv.custom$AllPairwiseComp$logFC)
          th <- as.numeric(rv.widgets$HypothesisTest_thlogFC)
          DaparToolshed::hc_logFC_DensityPlot(tmp.df, th)
       })
    })
    

    
    output$showConds <- renderUI({
      req(rv.custom$AllPairwiseComp)
      .style <- "align: center; display:inline-block; vertical-align: middle;
      padding-right: 50px; padding-bottom: 50px;"
      
      
      widget <- lapply(seq_len(rv.custom$n), function(i) {
        ll.conds <- unlist(
          strsplit(
            colnames(rv.custom$AllPairwiseComp$logFC)[i],
            split = "_"
          )
        )
        
        div(
          div(style = .style, p(gsub("[()]", "", ll.conds[1]))),
          div(style = .style, p(gsub("[()]", "", ll.conds[3]))),
          div(style = .style,
            actionButton(ns(paste0("compswap", i)), "",
              icon("sync", lib = "font-awesome"),
              style = "border-width: 0px; padding: 0px",
              width = "30px", height = "30px",
              class = actionBtnClass
            )
          )
        )
      })
      
      do.call(tagList, widget)
      MagellanNTK::toggleWidget(widget, rv$steps.enabled['HypothesisTest'])
      
    })
    
    
    observeEvent(req(sum(GetSwapShinyValue()) > 0), {
      req(rv.custom$AllPairwiseComp)
      swap <- GetSwapShinyValue()
      
      isolate({
        ind.swap <- which(swap != rv.custom$swap.history)
        rv.custom$swap.history <- swap
        if (length(ind.swap) > 0) {
          for (i in ind.swap) {
            current.comp <- colnames(rv.custom$AllPairwiseComp$logFC)[i]
            
            # Swap comparisons names
            ll <- unlist(strsplit(current.comp, split = "_"))
            tmp.cond1 <- gsub("[( )]", "", ll[1])
            tmp.cond2 <- gsub("[( )]", "", ll[3])
            tmp.logFC <- paste0("(", tmp.cond2, ")_vs_(", tmp.cond1, ")_logFC" )
            tmp.pval <- paste0( "(",  tmp.cond2, ")_vs_(", tmp.cond1, ")_pval" )
            colnames(rv.custom$AllPairwiseComp$logFC)[i] <- tmp.logFC
            colnames(rv.custom$AllPairwiseComp$P_Value)[i] <- tmp.pval
            
            # Swap logFC values
            .logFC <- rv.custom$AllPairwiseComp$logFC
            rv.custom$AllPairwiseComp$logFC[, i] <- -.logFC[, i]
          }
        }
      })
    })
    
    GetSwapShinyValue <- reactive({
      req(rv.custom$n)

      unlist(lapply(seq(rv.custom$n), function(x) 
        input[[paste0("compswap", x)]]
      ))
    })
    
    
    ### Computation of comparisons selected in the variable
    # 'rv$widgets$hypothesisTest$design'
    ComputeComparisons <- reactive({
      req(rv.widgets$HypothesisTest_method != "None")
      req(rv.widgets$HypothesisTest_design != "None")
      
      #browser()
      m <- DaparToolshed::match.metacell(omXplore::get_metacell(rv$dataIn[[length(rv$dataIn)]]),
        pattern = "Missing",
        level = omXplore::get_type(rv$dataIn[[length(rv$dataIn)]]))
      
      req(length(which(m)) == 0)
      
      rv.custom$AllPairwiseComp <- NULL
      rv.custom$AllPairwiseComp <- switch(rv.widgets$HypothesisTest_method,
        Limma = {
          DaparToolshed::limmaCompleteTest(
            qData = SummarizedExperiment::assay(rv$dataIn, length(rv$dataIn)),
            sTab = MultiAssayExperiment::colData(rv$dataIn),
            comp.type = rv.widgets$HypothesisTest_design
          )
        },
        ttests = {
          DaparToolshed::compute_t_tests(
            obj = rv$dataIn,
            i = length(rv$dataIn),
            contrast = rv.widgets$HypothesisTest_design,
            type = rv.widgets$HypothesisTest_ttestOptions
          )
        }
      )
      
      rv.custom$listNomsComparaison <- colnames(rv.custom$AllPairwiseComp$logFC)

      rv.custom$AllPairwiseComp
    }) %>% bindCache(
      rv$dataIn,
      rv.widgets$HypothesisTest_method,
      rv.widgets$HypothesisTest_design,
      rv.widgets$HypothesisTest_ttestOptions
    )
    
    
    
    
    
    
    
    output$HypothesisTest_btn_validate_ui <- renderUI({
      widget <- actionButton(ns("HypothesisTest_btn_validate"),
        "Run HypothesisTest",
        class = "btn-success")
      toggleWidget(widget, rv$steps.enabled['HypothesisTest'] )
      
    })
    
    observeEvent(input$HypothesisTest_btn_validate, {
      # Do some stuff
      
      # DO NOT MODIFY THE THREE FOLLOWINF LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- NULL
      rv$steps.status['HypothesisTest'] <- stepStatus$VALIDATED
    })
    
    enable_Limma <- reactive({
      req(rv$dataIn)
      
      enable <- TRUE
      
      nConds <-length(unique(omXplore::get_group(rv$dataIn)))
      design <- MultiAssayExperiment::colData(rv$dataIn)
      nLevel <- DaparToolshed::getDesignLevel(design)   
      enable <- (nConds <= 26 && nLevel == 1) ||
        (nConds < 10 && (nLevel%in% c(2,3)))
      
      enable
    })
    
    output$HypothesisTest_info_Limma_disabled_ui <- renderUI({
      req(!enable_Limma())
      tagList(
        tags$p('Info: Limma has been disabled because the design of your dataset:'),
        tags$ul(
          tags$li(p('is of level 1 and contains more than 26 conditions,')),
          tags$li('is of level 2 or 3 and contains more than 9 conditions.')
        ),
        tags$p('Prostar does not currently handle these cases.')
      )
    })
    
    
    
    output$HypothesisTest_swapConds_ui <- renderUI({
      req(rv.custom$AllPairwiseComp)
      tagList(
        h3("Swap conditions"),
          uiOutput(ns("showConds")),
      )
    })
    
    # >>> END: Definition of the widgets
    
    observeEvent(input$HypothesisTest_perform_btn, {
      # Do some stuff 
      
      req(rv.widgets$HypothesisTest_method)
      req(rv$dataIn)

      rv.widgets$HypothesisTest_thlogFC <- as.numeric(
        rv.widgets$HypothesisTest_thlogFC)
      #browser()
      rv.custom$AllPairwiseComp <- ComputeComparisons()
      
      if(is.null(rv.custom$AllPairwiseComp)){} 
      else if(inherits(rv.custom$AllPairwiseComp, "try-error")) {
        
        MagellanNTK::mod_SweetAlert_server(id = 'sweetalert_PerformLogFCPlot',
          text = rv.custom$AllPairwiseComp[[1]],
          type = 'error' )
      } else {
        # sendSweetAlert(
        #   session = session,
        #   title = "Success",
        #   type = "success"
        # )
        
        rv.custom$n <- ncol(rv.custom$AllPairwiseComp$logFC)
        rv.custom$swap.history <- rep(0, rv.custom$n)
      }
      
  })
    
    # <<< END ------------- Code for step 1 UI---------------
    
    
    # >>> START ------------- Code for step 3 UI---------------
    output$Save <- renderUI({
      tagList(
        # Insert validation button
        # This line is necessary. DO NOT MODIFY
        uiOutput(ns('Save_btn_validate_ui'))
      )
    })
    
    output$Save_btn_validate_ui <- renderUI({
      tagList(
        toggleWidget( 
          actionButton(ns("Save_btn_validate"), "Validate step",
            class = "btn-success"),
          rv$steps.enabled['Save']
        ),
        if (config@mode == 'process' && 
            rv$steps.status['Save'] == stepStatus$VALIDATED) {
          download_dataset_ui(ns('createQuickLink'))
        }
      )
      
    })
    observeEvent(input$Save_btn_validate, {
      # Do some stuff

      new.dataset <- rv$dataIn[[length(rv$dataIn)]]
      df <- cbind(rv.custom$AllPairwiseComp$logFC, 
        rv.custom$AllPairwiseComp$P_Value)
      DaparToolshed::HypothesisTest(new.dataset) <- DataFrame(df)
      rv$dataIn <- QFeatures::addAssay(rv$dataIn, new.dataset, 'HypothesisTest')
      
      # DO NOT MODIFY THE THREE FOLLOWING LINES
      dataOut$trigger <- Timestamp()
      dataOut$value <- rv$dataIn
      rv$steps.status['Save'] <- stepStatus$VALIDATED
      download_dataset_server('createQuickLink', data = reactive({rv$dataIn}))
      
    })
    # <<< END ------------- Code for step 3 UI---------------
    
    
    
    # Insert necessary code which is hosted by MagellanNTK
    # DO NOT MODIFY THIS LINE
    eval(parse(text = Module_Return_Func()))
  }
  )
}

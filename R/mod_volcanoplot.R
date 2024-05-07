

#' @title   popover_for_help_ui and popover_for_help_server
#' @description  A shiny Module.
#'
#' @param id xxx
#' @param dataIn An instance of the class `SummarizedExperiment` containing
#' results of hypothesis testing stored as a data.frame in its metadata
#' @param thlogfc xxxx
#' @param thpval xxx
#' @param tooltip xxx
#' @param reset xxx
#' @param is.enabled xxx
#' 
#' 
#' @name volcanoplot
#' 
#' @examplesIf interactive()
#' library(highcharter)
#' data(Exp1_R25_prot, package="DaparToolshedData")
#' obj <- Exp1_R25_prot
#' # Simulate imputation of missing values
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' qData <- as.matrix(assay(obj[[2]]))
#' sTab <- MultiAssayExperiment::colData(obj)
#' limma <- limmaCompleteTest(qData, sTab)
#' 
#' df <- data.frame(
#'     x = limma$logFC[['25fmol_vs_10fmol_logFC']], 
#'     y = -log10(limma$P_Value[['25fmol_vs_10fmol_pval']]),
#'     index = as.character(rownames(obj[[2]])))
#' colnames(df) <- c("x", "y", "index")
#' tooltipSlot <- c("Fasta_headers", "Sequence_length")
#' df <- cbind(df, Biobase::fData(obj)[, tooltipSlot])
#' colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
#' if (ncol(df) > 3) {
#'     colnames(df)[seq.int(from = 4, to = ncol(df))] <-
#'         paste("tooltip_", colnames(df)[seq.int(from = 4, to = ncol(df))],
#'          sep = "")
#' }
#' hc_clickFunction <- JS("function(event) {
#' Shiny.onInputChange('eventPointClicked',
#' [this.index]+'_'+ [this.series.name]);}")
#' cond <- c("25fmol", "10fmol")
#' diffAnaVolcanoplot_rCharts(
#' df, 
#' th_pval = 2.5, 
#' th_logfc = 1, 
#' conditions = cond, 
#' clickFunction = hc_clickFunction)
#' shiny::runApp(volcanoplot(xxxx))
#'
NULL




#' @importFrom shiny NS tagList
#' @importFrom shinyjs inlineCSS useShinyjs
#' @export
#' @rdname volcanoplot
#'
mod_volcanoplot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    highcharter::highchartOutput(ns("volcanoPlot"), 
      width = "600px", height = "600px"),
    uiOutput(ns("quantiDT"))
  )
}



#' @export
#' @rdname volcanoplot
#'
mod_volcanoplot_server <- function(
    id,
    dataIn = reactive({NULL}),
  comparison = reactive({NULL}),
    group = reactive({NULL}),
  thlogfc = reactive({0}),
    thpval = reactive({0}),
    tooltip = reactive({NULL}),
    reset = reactive({NULL}),
    is.enabled = reactive({TRUE})
    ) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    DifferentialAnalysis_thpval = NULL,
    widget2 = NULL,
    widget3 = NULL
  )
  
  
  rv.custom.default.values <- list(
    colorsVolcanoplot = NULL,
    data = NULL,
    conditions = NULL
  )
  
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # DO NOT MODIFY THIS FUNCTION CALL
    eval(
      str2expression(
        Get_AdditionalModule_Core_Code(
          w.names = names(widgets.default.values),
          rv.custom.names = names(rv.custom.default.values)
        )
      )
    )
    
    observeEvent(dataIn(), {
      rv$dataIn <- dataIn()
       rv.custom$data <- HypothesisTest(dataIn())
       print('toto')
    })
    
    
    # output$quantiDT <- renderUI({
    #   req(input$eventPointClicked)
    #   
    #   if (omXplore::get_type(rv$dataIn) == "protein") {
    #     if (is.null(omXplore::get_adjacencyMatrix(rv$dataIn))) {
    #       shinyBS::bsCollapse(id = ns("collapseVolcanoInfos"),
    #         open = "Protein",
    #         multiple = TRUE,
    #         shinyBS::bsCollapsePanel("Protein",
    #           tagList(
    #             uiOutput(ns("Warning_Infos")),
    #             DT::dataTableOutput(ns("Infos"))
    #           ), style = "info")
    #       )
    #     } else {
    #       shinyBS::bsCollapse(id = ns("collapseVolcanoInfos"),
    #         open = c("Protein", "Specific peptides", "Shared peptides"),
    #         multiple = TRUE,
    #         shinyBS::bsCollapsePanel("Protein",
    #           tagList(
    #             uiOutput(ns("Warning_Infos")),
    #             DT::dataTableOutput(ns("Infos"))
    #           ), style = "info"),
    #         shinyBS::bsCollapsePanel("Specific peptides",
    #           tagList(
    #             uiOutput(ns("Warning_specificPeptidesInfos")),
    #             DT::dataTableOutput(ns("specificPeptidesInfos"))
    #           ), style = "primary"),
    #         shinyBS::bsCollapsePanel("Shared peptides",
    #           tagList(
    #             uiOutput(ns("Warning_sharedPeptidesInfos")),
    #             DT::dataTableOutput(ns("sharedPeptidesInfos"))
    #           ), style = "primary")
    #       )
    #     }
    #   } else if (omXplore::get_type(rv$dataIn) == "peptide") {
    #     shinyBS::bsCollapse(id = ns("collapseVolcanoInfos"),
    #       open = "Peptide",
    #       multiple = TRUE,
    #       shinyBS::bsCollapsePanel("Peptide",
    #         tagList(
    #           uiOutput(ns("Warning_Infos")),
    #           DT::dataTableOutput(ns("Infos"))
    #         ),
    #         style = "info"
    #       )
    #     )
    #   }
    # })
    
    
    
    GetSortingIndices <- reactive({
      req(comparison())
      
      condition1 <- comparison()[1]
      condition2 <- comparison()[2]
      .ind <- which(group() == comparison()[1])
      
      if (length(grep("all", comparison()[2])) == 0) {
        ind <- c(
          which(group() == comparison()[1]),
          which(group() == comparison()[2])
        )
      } else {
        ind <- c(which(group() == comparison()[1]),
          c(seq(length(group())))[-.ind])
      }
      ind
    })
    
    GetBorderIndices <- reactive({
      conds <- group()[GetSortingIndices()]
      ## build index for border-formatting
      borders_index <- unlist(lapply(
        unique(conds),
        function(x) {
          which.max(x == conds)
        }
      ))
      borders_index
    })
    
    
    
    output$Warning_sharedPeptidesInfos <- renderUI({
      GetDataFor_sharedPeptidesInfos()
      if (nrow(GetDataFor_sharedPeptidesInfos()) > 153) {
        p(MSG_WARNING_SIZE_DT)
      }
    })
    
    
    
    
    GetDataFor_sharedPeptidesInfos <- reactive({
      req(comparison())
      
      ind <- GetSortingIndices()
      borders_index <- GetBorderIndices()
      
      #.ind <- last(grep(pattern = "peptide", names(rv$dataIn)))
      #prev.dataset <- rv$dataIn[[names(rv$dataIn)[.ind]]]
      
      prot <- GetExprsClickedProtein()
      prot.indice <- rownames(prot)
      
      data <- getDataForExprs(rv$dataIn, 3)
      data <- data[, c(ind, (ind + ncol(data) / 2))]
      
      Xspec <- ExtractUniquePeptides(omXplore::get_adjacencyMatrix(rv$dataIn))
      Xshared <- omXplore::get_adjacencyMatrix(rv$dataIn)
      
      i <- which(colnames(Xspec) == prot.indice)
      specificPeptidesIndices <- which(Xspec[, i] == 1)
      allPeptidesIndices <- which(Xshared[, i] == 1)
      peptidesIndices <- setdiff(allPeptidesIndices, specificPeptidesIndices)
      data <- data[peptidesIndices, ]
      data
    })
    
    output$sharedPeptidesInfos <- DT::renderDataTable(server = TRUE, {
      data <- GetDataFor_sharedPeptidesInfos()
      c.tags <- BuildColorStyles(rv$dataIn)$tags
      c.colors <- BuildColorStyles(rv$dataIn)$colors
      
      dt <- DT::datatable(
        data,
        extensions = c("Scroller"),
        options = list(
          initComplete = initComplete(),
          dom = "frtip",
          blengthChange = FALSE,
          displayLength = 20,
          ordering = FALSE,
          server = FALSE,
          columnDefs = list(
            list(
              targets = c(((ncol(data) / 2) + 1):(ncol(data))),
              visible = FALSE
            )
          )
        )
      ) %>%
        DT::formatStyle(
          colnames(data)[1:(ncol(data) / 2)],
          colnames(data)[((ncol(data) / 2) + 1):(ncol(data))],
          backgroundColor = DT::styleEqual(c.tags, c.colors)
        ) %>%
        DT::formatStyle(GetBorderIndices(), borderLeft = "3px solid #000000")
      
      dt
    })
    
    output$Warning_specificPeptidesInfos <- renderUI({
      GetDataFor_specificPeptidesInfos()
      if (nrow(GetDataFor_specificPeptidesInfos()) > 153) {
        p(MSG_WARNING_SIZE_DT)
      }
    })
    
    GetDataFor_specificPeptidesInfos <- reactive({
      
      ind <- GetSortingIndices()
      borders_index <- GetBorderIndices()
      #.ind <- last(grep(pattern = "peptide", names(rv$dataIn)))
      #prev.dataset <- rv$dataIn[[names(rv$dataIn)[.ind]]]
      
      prot <- GetExprsClickedProtein()
      prot.indice <- rownames(prot)
      
      data <- getDataForExprs(rv$dataIn, 3)
      data <- data[, c(ind, (ind + ncol(data) / 2))]
      
      
      Xspec <- ExtractUniquePeptides(omXplore::get_adjacencyMatrix(rv$dataIn))
      
      i <- which(colnames(Xspec) == prot.indice)
      peptidesIndices <- which(Xspec[, i] == 1)
      data <- data[peptidesIndices, ]
      data
    })
    
    
    output$specificPeptidesInfos <- DT::renderDataTable(server = TRUE, {
      data <- GetDataFor_specificPeptidesInfos()
      c.tags <- BuildColorStyles(rv$dataIn)$tags
      c.colors <- BuildColorStyles(rv$dataIn)$colors
      
      dt <- DT::datatable(
        data,
        extensions = c("Scroller"),
        options = list(
          initComplete = initComplete(),
          dom = "frtip",
          blengthChange = FALSE,
          displayLength = 20,
          ordering = FALSE,
          columnDefs = list(
            list(
              targets = c(
                ((ncol(data) / 2) + 1):(ncol(data))
              ),
              visible = FALSE
            )
          )
        )
      ) %>%
        DT::formatStyle(
          colnames(data)[1:(ncol(data) / 2)],
          colnames(data)[((ncol(data) / 2) + 1):(ncol(data))],
          backgroundColor = DT::styleEqual(c.tags, c.colors)
        ) %>%
        DT::formatStyle(GetBorderIndices(), borderLeft = "3px solid #000000")
      
      dt
    })
    
    
    ## ---------------------------------------------------------
    GetExprsClickedProtein <- reactive({
      req(rv$dataIn)
      req(input$eventPointClicked)
      
      
      ind <- GetSortingIndices()
      # browser()
      this.index <- as.integer(strsplit(input$eventPointClicked, "_")[[1]][1])
      this.series.name <- strsplit(input$eventPointClicked, "_")[[1]][2]
      
      .digits <- 3
      data <- getDataForExprs(rv$dataIn, .digits)
      data <- data[, c(ind, (ind + ncol(data) / 2))]
      
      index.g1 <- which((-log10(rv.custom$data$P_Value) >= thpval()
      ) & (abs(rv.custom$data$logFC) >= as.numeric(thlogfc())))
      
      data.g1 <- data[index.g1, ]
      data.g2 <- data[-index.g1, ]
      
      switch(this.series.name,
        g1 = data <- data.g1[this.index + 1, ],
        g2 = data <- data.g2[this.index + 1, ]
      )
      data
    })
    
    
    
    output$Warning_Infos <- renderUI({
      GetDataFor_Infos()
      if (nrow(GetDataFor_Infos()) > 153) {
        p(MSG_WARNING_SIZE_DT)
      }
    })
    
    
    
    
    GetDataFor_Infos <- reactive({
      data <- GetExprsClickedProtein()
      data
    })
    
    ## -------------------------------------------------------------
    output$Infos <- DT::renderDataTable(server = TRUE, {
      
      borders_index <- GetBorderIndices()
      data <- GetExprsClickedProtein()
      c.tags <- BuildColorStyles(rv$dataIn)$tags
      c.colors <- BuildColorStyles(rv$dataIn)$colors
      
      dt <- DT::datatable(
        data,
        extensions = c("Scroller"),
        options = list(
          initComplete = initComplete(),
          dom = "frtip",
          blengthChange = FALSE,
          displayLength = 20,
          ordering = FALSE,
          header = FALSE,
          columnDefs = list(
            list(
              targets = c(((ncol(data) / 2) + 1):(ncol(data))),
              visible = FALSE
            )
          )
        )
      ) %>%
        DT::formatStyle(
          colnames(data)[1:(ncol(data) / 2)],
          colnames(data)[((ncol(data) / 2) + 1):(ncol(data))],
          backgroundColor = DT::styleEqual(c.tags, c.colors)
        ) %>%
        DT::formatStyle(borders_index, borderLeft = "3px solid #000000")
      dt
    })
    
    ## ---------------------------------------------------------------------
    output$volcanoPlot <- highcharter::renderHighchart2({
      req(rv$dataIn)
      comparison()
      thlogfc()
      thpval()
      
      print('IN output$volcanoPlot <- renderHighchart')
        isolate({
          withProgress(message = "Building plot...", detail = "", value = 0, {
          m <- match.metacell(omXplore::get_metacell(rv$dataIn),
            pattern = c("Missing", "Missing POV", "Missing MEC"),
            level = omXplore::get_type(rv$dataIn)
          )
          if (length(which(m)) > 0) {
            return(NULL)
          }
          

          ht <- HypothesisTest(rv$dataIn)
          prefix <- paste0(comparison()[1], '_vs_', comparison()[2])
          df <- data.frame(
            x = ht[, paste0(prefix, '_logFC')],
            y = -log10(ht[, paste0(prefix, '_pval')]),
            index = seq(nrow(ht))
          )
          
          if (length(tooltip()) > 0 && sum(is.na(tooltip())) == 0) {
            df <- cbind(df, 
              SummarizedExperiment::rowData(rv$dataIn)[tooltip()])
          }
          
          colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
          if (ncol(df) > 3) {
            .range <- seq.int(from = 4, to = ncol(df), by = 1)
            colnames(df)[.range] <-
              paste("tooltip_", colnames(df)[.range], sep = "")
          }
          
          clickFun <-
            JS(paste0(
              "function(event) {Shiny.onInputChange('",
              ns("eventPointClicked"),
              "', [this.index]+'_'+ [this.series.name]);}"
            ))
          
          widget <- diffAnaVolcanoplot_rCharts(
            df,
            th_logfc = as.numeric(thlogfc()),
            th_pval = as.numeric(thpval()),
            conditions = group(),
            clickFunction = clickFun,
            pal = rv.custom$colorsVolcanoplot
          )
})
        })
        widget
        #MagellanNTK::toggleWidget(widget, is.enabled())
      })
  })
}



#' @export
#' @rdname volcanoplot
#' 
volcanoplot <- function(
    dataIn, 
    comparison,
  group, 
    thlogfc, 
    thpval,
  tooltip){
  
  ui <- fluidPage(
  mod_volcanoplot_ui("volcano"))
  
  server <- function(input, output) {
    logpval <- mod_volcanoplot_server(
      id = "volcano",
      dataIn = reactive({dataIn}),
      comparison = reactive({comparison}),
      group = reactive({group}),
      thlogfc = reactive({thlogfc}),
      thpval = reactive({thpval}),
      tooltip = reactive({tooltip})
     # fdr = reactive({3.8})
      )
  
    observeEvent(logpval(), {
      print(logpval())
    })
  }

  app <- shiny::shinyApp(ui, server)
}

#' @title Filtering Shiny module
#'
#' @description
#' This function is a shiny module to create a list of queries (instances of 
#' the class `Filtering` to filter the quantitative metadata of an instance
#' of the class `SummarizedExperiment`).
#' This function is written with specifications of the package `MagellanNTK` so
#' as to be easily integrated into workflfow compliant with `MagellanNTK`.
#'
#' @name mod_Metacell_Filtering
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
#' shiny::runApp(mod_Metacell_Filtering(ft_na[[1]]))
#' 
NULL





#' @export
#'
#' @rdname mod_Metacell_Filtering
#'
mod_Metacell_Filtering_ui <- function(id) {
  ns <- NS(id)
  wellPanel(
    # uiOutput for all widgets in this UI
    # This part is mandatory
    # The renderUI() function of each widget is managed by MagellanNTK
    # The dev only have to define a reactive() function for each
    # widget he want to insert
    # Be aware of the naming convention for ids in uiOutput()
    # For more details, please refer to the dev document.
    
    uiOutput(ns("Quantimetadatafiltering_buildQuery_ui")),

    DT::dataTableOutput(ns("qMetacell_Filter_DT"))
    ,uiOutput(ns('plots_ui')),
    # Insert validation button
    uiOutput(ns("Quantimetadatafiltering_btn_validate_ui"))
  )
}




#' @param id xxx
#' @param obj An instance of the class `QFeatures`
#' @param keep_vs_remove A character(1) indicating whether to keep or delete 
#' items. Default value is "delete"
#' @param operator xxx
#' @param reset A `ìnteger(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
#' @rdname mod_Metacell_Filtering
#'
#' @export
#'
mod_Metacell_Filtering_server <- function(id,
  obj = reactive({NULL}),
  i = reactive({1}),
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list()
  
  rv.custom.default.values <- list(
    indices = NULL,
    Filtering = NULL,
    query = list(),
    fun.list = list(),
    widgets.value = list(),
    funFilter = reactive({NULL}),
    qMetacell_Filter_SummaryDT = NULL
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
    
    
    observeEvent(req(obj()), ignoreNULL = TRUE,{

      stopifnot(inherits(obj(), 'QFeatures'))
      rv$dataIn <- obj()
      rv.custom$qMetacell_Filter_SummaryDT <- data.frame(
        query = "-",
        nbDeleted = "-",
        TotalMainAssay = "-",
        stringsAsFactors = FALSE
      )
      
      
      mod_ds_metacell_Histos_server(
        id = "plots",
        obj = reactive({rv$dataIn[[i()]]}),
        pattern = reactive({"Missing"}),
        group = reactive({omXplore::get_group(rv$dataIn)})
      )
      
      
      
    }, priority = 1000)

    
    
    output$plots_ui <- renderUI({

      widget <- mod_ds_metacell_Histos_ui(ns("plots"))
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    
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
    
    output$qMetacell_Filter_DT <- DT::renderDataTable(server = TRUE,{
        req(rv.custom$qMetacell_Filter_SummaryDT)
        req(rv.custom$funFilter()$value$ll.query)
        df <- rv.custom$qMetacell_Filter_SummaryDT
        query <- rv.custom$funFilter()$value$ll.query
        df[, "query"] <- ConvertListToHtml(query)
        showDT(df)
      })
    
    
    observe({
      req(is.enabled())
      req(rv$dataIn)
      rv.custom$funFilter <- mod_qMetacell_FunctionFilter_Generator_server(
        id = "query",
        obj = reactive({rv$dataIn[[i()]]}),
        conds = reactive({omXplore::get_group(rv$dataIn)}),
        keep_vs_remove = reactive({stats::setNames(c('Push p-value', 'Keep original p-value'), nm = c("delete", "keep"))}),
        val_vs_percent = reactive({stats::setNames(nm = c("Count", "Percentage"))}),
        operator = reactive({stats::setNames(nm = SymFilteringOperators())})
      )
    })
    
    output$Quantimetadatafiltering_buildQuery_ui <- renderUI({
      widget <- mod_qMetacell_FunctionFilter_Generator_ui(ns("query"))
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    

    output$Quantimetadatafiltering_btn_validate_ui <- renderUI({
      #browser()
      req(length(rv.custom$funFilter()$value$ll.fun) > 0)
      
      widget <- actionButton(ns("Quantimetadatafiltering_btn_validate"),
        "Perform qMetacell filtering",
        class = "btn-success"
      )
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$Quantimetadatafiltering_btn_validate, {

      tmp <- filterFeaturesOneSE(
        object = rv$dataIn,
        i = i(),
        name = "qMetacellFiltered",
        filters = rv.custom$funFilter()$value$ll.fun
      )
      
      indices <- GetIndices(rv.custom$funFilter()$value$ll.fun)
      
      # Add infos
      nBefore <- nrow(tmp[[length(tmp) - 1]])
      nAfter <- nrow(tmp[[length(tmp)]])
      
      rv.custom$qMetacell_Filter_SummaryDT[, "nbDeleted"] <- nBefore - nAfter
      rv.custom$qMetacell_Filter_SummaryDT[, "TotalMainAssay"] <- nrow(assay(tmp[[length(tmp)]]))
      
      par <- rv.custom$funFilter()$value$ll.widgets.value
      params(tmp[[i()]], length(tmp[[i()]])) <- par
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- tmp
    })
    
    return(reactive({dataOut}))
  })
}



#' @export
#' @rdname mod_Metacell_Filtering
#' 
mod_Metacell_Filtering <- function(obj = NULL, i = 1){
  ui <- mod_Metacell_Filtering_ui('query')
  
  server <- function(input, output, session){
    
    res <- mod_Metacell_Filtering_server('query',
      obj = reactive({obj}),
      i = reactive({i}))
    
    observeEvent(res()$trigger, {
      print(res()$value)
    })
  }
  
  app <- shiny::shinyApp(ui, server)
  
}
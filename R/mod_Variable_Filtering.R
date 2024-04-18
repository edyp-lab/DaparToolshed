#' @title Filtering Shiny module
#'
#' @description
#' This function is a shiny module to create a list of queries (instances of 
#' the class `Filtering` to filter the quantitative metadata of an instance
#' of the class `SummarizedExperiment`).
#' This function is written with specifications of the package `MagellanNTK` so
#' as to be easily integrated into workflfow compliant with `MagellanNTK`.
#'
#' @name mod_Variable_Filtering
#' 
#' @param id xxx
#' @param obj An instance of the class `QFeatures`
#' @param keep_vs_remove A character(1) indicating whether to keep or delete 
#' items. Default value is "delete"
#' @param operator xxx
#' @param reset A `ìnteger(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
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
#' shiny::runApp(mod_Variable_Filtering(ft_na[[1]]))
#' 
NULL





#' @export
#'
#' @rdname mod_Variable_Filtering
#'
mod_Variable_Filtering_ui <- function(id) {
  ns <- NS(id)
  wellPanel(
    # uiOutput for all widgets in this UI
    # This part is mandatory
    # The renderUI() function of each widget is managed by MagellanNTK
    # The dev only have to define a reactive() function for each
    # widget he want to insert
    # Be aware of the naming convention for ids in uiOutput()
    # For more details, please refer to the dev document.
    
    uiOutput(ns("variable_buildQuery_ui")),
    shinydashboardPlus::box(
      DT::dataTableOutput(ns("variable_Filter_DT"))
      ,uiOutput(ns('plots_ui'))
    ),
    # Insert validation button
    uiOutput(ns("variable_btn_validate_ui"))
  )
}





#' @rdname mod_Variable_Filtering
#'
#' @export
#'
mod_Variable_Filtering_server <- function(id,
  obj,
  i,
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
    variable_Filter_SummaryDT = data.frame(
      query = "-",
      nbDeleted = "-",
      TotalMainAssay = "-",
      stringsAsFactors = FALSE
    )
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
    
    
    observe({
      stopifnot(inherits(obj(), 'QFeatures'))
    })
    
    
    output$plots_ui <- renderUI({
      
      mod_ds_metacell_Histos_server(
        id = "plots",
        obj = reactive({obj()[[i()]]}),
        pattern = reactive({"Missing"}),
        group = reactive({omXplore::get_group(obj())})
      )
      mod_ds_metacell_Histos_ui(ns("plots"))
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
    
    output$variable_Filter_DT <- DT::renderDataTable(server = TRUE,{
      req(rv.custom$funFilter()$value$ll.query)

      df <- rv.custom$variable_Filter_SummaryDT
      query <- rv.custom$funFilter()$value$ll.query
      df[, "query"] <- ConvertListToHtml(query)
      showDT(df)
    })
    
    output$variable_buildQuery_ui <- renderUI({
      widget <- mod_VariableFilter_Generator_ui(ns("query"))
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    rv.custom$funFilter <- mod_VariableFilter_Generator_server(
      id = "query",
      obj = reactive({obj()[[i()]]})
    )
    
    
    
    observeEvent(rv.custom$funFilter()$trigger, {
    })
    
    
    output$variable_btn_validate_ui <- renderUI({
      #browser()
      req(length(rv.custom$funFilter()$value$ll.var) > 0)
      
      widget <- actionButton(ns("variable_btn_validate"),
        "Perform filtering", class = "btn-success")
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    # >>> END: Definition of the widgets
    
    
    observeEvent(input$variable_btn_validate, {
      browser()
      tmp <- filterFeaturesOneSE(
        object = obj(),
        i = i(),
        name = "variableFiltered",
        filters = rv.custom$funFilter()$value$ll.var
      )
      
      # Add infos
      #browser()
      nBefore <- nrow(tmp[[length(tmp) - 1]])
      nAfter <- nrow(tmp[[length(tmp)]])
      
      rv.custom$variable_Filter_SummaryDT[, "nbDeleted"] <- nBefore - nAfter
      rv.custom$variable_Filter_SummaryDT[, "TotalMainAssay"] <- nrow(assay(tmp[[length(tmp)]]))
      
      par <- rv.custom$funFilter()$value$ll.widgets.value
      params(tmp[[i()]], length(tmp[[i()]])) <- par
      dataOut$trigger <- MagellanNTK::Timestamp()
      dataOut$value <- tmp
    })
    
    return(reactive({dataOut}))
  })
}



#' @export
#' @rdname mod_Variable_Filtering
#' 
mod_Variable_Filtering <- function(obj, i){
  ui <- mod_Variable_Filtering_ui('query')
  
  server <- function(input, output, session){
    
    res <- mod_Variable_Filtering_server('query',
      obj = reactive({obj}),
      i = reactive({i}))
    
    observeEvent(res()$trigger, {
      print(res()$value)
    })
  }
  
  app <- shiny::shinyApp(ui, server)
  
}
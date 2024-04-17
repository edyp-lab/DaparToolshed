#' @title Build Variable Filter function
#'
#' @description
#' This function is a shiny module to create a list of queries (instances of 
#' the class `VariableFilter` to filter the quantitative metadata of an instance
#' of the class `SummarizedExperiment`).
#' This function is written with specifications of the package `MagellanNTK` so
#' as to be easily integrated into workflfow compliant with `MagellanNTK`.
#'
#' @name mod_VariableFilter_Generator
#'
#' @return As for all modules used with `MagellanNTK`, the return value is a
#' `list()` of two items:
#' - trigger : xxx
#' - value: In this case, it contains a list() of three slots:
#'   - ll.var: a list() of instances of the class `VariableFilter`,
#'   - ll.query: a list of `character()` which describe the queries in natural
#'   language,
#'   - ll.widgets.value: a list of the values of widgets.
#'
#' @examplesIf interactive()
#' data(ft_na)
#' obj <- ft_na[[1]]
#' operator = setNames(nm = SymFilteringOperators())
#' keep_vs_remove <- setNames(nm = c("delete", "keep"))
#' value = 3
#' shiny::runApp(
#' mod_VariableFilter_Generator(obj, keep_vs_remove, value, operator))
#' 
NULL





#' @export
#'
#' @rdname mod_VariableFilter_Generator
#'
mod_VariableFilter_Generator_ui <- function(id) {
  ns <- NS(id)
  wellPanel(
    DT::dataTableOutput(ns("VarFilter_DT")),
    # Build queries
    uiOutput(ns("chooseKeepRemove_ui")),
    uiOutput(ns("cname_ui")),
    uiOutput(ns("value_ui")),
    uiOutput(ns("operator_ui")),
    uiOutput(ns("addFilter_btn_ui")),
    # Show example
    uiOutput(ns("example_ui")),
    # Process the queries
    uiOutput(ns("btn_validate_ui"))
  )
}




#' @param id xxx
#' @param obj An instance of the class `SummarizedExperiment`
#' @param keep_vs_remove xxx
#' @param operator xxx
#' @param reset A `ìnteger(1)` xxxx
#' @param is.enabled A `logical(1)` that indicates whether the module is
#' enabled or disabled. This is a remote command.
#'
#' @rdname mod_VariableFilter_Generator
#'
#' @export
#'
mod_VariableFilter_Generator_server <- function(id,
  obj,
  cname = reactive({NULL}),
  value = reactive({NULL}),
  operator = reactive({NULL}),
  keep_vs_remove = reactive({'delete'}),
  reset = reactive({NULL}),
  is.enabled = reactive({TRUE})) {
  
  # Define default selected values for widgets
  # This is only for simple workflows
  widgets.default.values <- list(
    cname = "None",
    value = "None",
    keep_vs_remove = "delete",
    operator = "<="
  )
  
  rv.custom.default.values <- list(
    indices = NULL,
    variableFilter = NULL,
    query = list(),
    fun.list = list(),
    widgets.value = list()
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
    
    
    # output$example_ui <- renderUI({
    #   req(length(rv.custom$varFilters) > 0)
    #   req(rv$steps.status["Variablefiltering"] == 0)
    #   
    #   temp <- filterFeaturesOneSE(
    #     object = mainAssay(rv$dataIn),
    #     filters = rv.custom$varFilters
    #   )
    #   
    #   mod_filterExample_server(
    #     id = "varFilterExample",
    #     objBefore = reactive({mainAssay(rv$dataIn)}),
    #     objAfter = reactive({temp}),
    #     query = reactive({rv.custom$varQueries})
    #   )
    #   
    #   widget <- mod_filterExample_ui(ns("varFilterExample"))
    #   MagellanNTK::toggleWidget(widget, 
    #     rv$steps.enabled["Variablefiltering"])
    # })
    # 
    
    
    
    output$chooseKeepRemove_ui <- renderUI({
 
      widget <- radioButtons(ns("keepRemove"),
        "Type of filter operation",
        choices = setNames(nm = c('delete', 'keep')),
        selected = rv.widgets$keep_vs_remove
      )
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$VarFilter_DT <- DT::renderDataTable(server = TRUE,{
      req(rv.custom$varQueries)
        rv.custom$varFilter_DT[, "query"] <- ConvertListToHtml(rv.custom$varQueries)
        showDT(rv.custom$varFilter_DT)
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
    
    
    output$cname_ui <- renderUI({
      .choices <- c("None", colnames(rowData(obj())))
      
      widget <- selectInput(ns("cname"),
        "Column name",
        choices = stats::setNames(.choices, nm = .choices),
        selected = rv.widgets$cname,
        width = "300px"
      )
      
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    output$operator_ui <- renderUI({
      req(rv.widgets$value)
      if (is.na(as.numeric(rv.widgets$value))) {
        .operator <- c("==", "!=", "startsWith", "endsWith", "contains")
      } else {
        .operator <- DaparToolshed::SymFilteringOperators()
      }
      
      
      widget <- selectInput(ns("operator"),
        "operator",
        choices = stats::setNames(nm = .operator),
        selected = rv.widgets$operator,
        width = "100px"
      )
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    output$value_ui <- renderUI({
      widget <- textInput(ns("value"),
        "value",
        value = rv.widgets$value,
        width = "100px"
      )
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    

    output$addFilter_btn_ui <- renderUI({
      widget <- actionButton(ns("addFilter_btn"), "Add filter")
      MagellanNTK::toggleWidget(widget, is.enabled())
    })
    
    
    Extract_Value <- function(value){
      val <- NULL
      val <- tryCatch({
        as.numeric(value)
        },
        warning = function(w){value},
        error = function(e) NA
)
      return(val)
    }
    
    BuildVariableFilter <- reactive({
      req(rv.widgets$value)
      
      

        VariableFilter(
          field = rv.widgets$cname,
          value = Extract_Value(rv.widgets$value),
          condition = rv.widgets$operator
        )
    })
    
    
    WriteQuery <- reactive({
      #req()
      
      value <- Extract_Value(rv.widgets$value)
      query <- paste0(rv.widgets$cname, " ", rv.widgets$operator, " ", value)
      query
      
    })
    
    observeEvent(input$addFilter_btn, {
      rv.custom$ll.var <- append(rv.custom$ll.var, BuildVariableFilter())
      rv.custom$ll.query <- append(rv.custom$ll.query, WriteQuery())
      rv.custom$ll.widgets.value <- append(rv.custom$ll.widgets.value,
        list(reactiveValuesToList(rv.widgets)))
      
      
      # Append a new FunctionFilter to the list
      dataOut$trigger <- as.numeric(Sys.time())
      dataOut$value <- list(
        ll.var = rv.custom$ll.var,
        ll.query = rv.custom$ll.query,
        ll.widgets.value = rv.custom$ll.widgets.value
      )
    })
    
    return(reactive({dataOut}))
  })
}



#' @export
#' @rdname mod_VariableFilter_Generator
#' 
mod_VariableFilter_Generator <- function(
    obj,
  cname = NULL,
  value = NULL,
  operator = NULL,
  keep_vs_remove = 'delete'){
  
  
  ui <- mod_VariableFilter_Generator_ui('query')
  
  server <- function(input, output, session){
    
    res <- mod_VariableFilter_Generator_server('query',
      obj = reactive({obj}),
      cname = reactive({cname}),
      value = reactive({value}),
      operator = reactive({operator}),
      keep_vs_remove = reactive({keep_vs_remove})
      )
    
    observeEvent(res()$trigger, {
      print(res()$value)
    })
  }
  
  app <- shiny::shinyApp(ui, server)
  
}
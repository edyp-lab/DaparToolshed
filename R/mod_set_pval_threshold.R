

#' @title Module set pval threshold
#' @description  A shiny Module.
#' 
#' @param id xxx
#' @param pval_init xxx
#' @param fdr xxx
#' @param options xxx
#' 
#' 
#' 
#' @name mod_set_pval_threshold
#' 
#' @examplesIf interactive()
#' shiny::runApp(mod_set_pval_threshold())
NULL




#'
#' @export
#' @importFrom shiny NS tagList
#' @importFrom shinyjs inlineCSS useShinyjs
#' @rdname mod_set_pval_threshold
#'
mod_set_pval_threshold_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style("#pvalPanel {width: 100%;}"),
    wellPanel(id = 'pvalPanel',
      useShinyjs(),
      uiOutput(ns('warn_NULL_fdr_UI')),
      MagellanNTK::mod_popover_for_help_ui(ns("modulePopover_pValThreshold")),
      br(),
      tags$div(style = "align: center;display:inline-block; vertical-align: top;",
        tags$style(HTML(".form-control {height: 20px; font-size: 13px;}")),
        tags$style(HTML(".radio {padding-right: 0px; padding-bottom: 10px;}")),
        uiOutput(ns('thresholdType_UI'))),
      tags$div(style = "align: center;display:inline-block; vertical-align: center;",
        uiOutput(ns('text1_UI')),
        disabled(uiOutput(ns('text2_UI')))),
      tags$div(style = "align: center;display:inline-block; vertical-align: center;",
        uiOutput(ns('warn_text1_UI')),
        uiOutput(ns('warn_text2_UI'))),
      br(),
      tags$div(style = "align: center;display:inline-block; vertical-align: center; ",
        actionButton(ns('ApplyThreshold'), 'Apply threshold', class = actionBtnClass)
      ),
      tags$div(style = "align: center;display:inline-block; vertical-align: center; padding-left: 20px;",
        uiOutput(ns('showFDR_UI'))
      )
      
    )
  )
}

#' @rdname mod_set_pval_threshold
#'
#' @export
#'
mod_set_pval_threshold_server <- function(id,
  pval_init = reactive({1}),
  fdr = reactive({0}),
  threshold.type = reactive({'logpval'}),
  reset = reactive({FALSE}),
  is.enabled = reactive({TRUE})) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    dataOut <- reactiveVal()
    threshold_type <- reactiveVal('pval')
    
    .head <- "To perform the selection using a FDR threshold of x% : "
    .pt1 <- "Display in the table below the adjusted p-values. The proteins are then automatically sorted by increasing adjusted p-values"
    .pt2 <- "Spot the protein P which has the largest adjusted p-value below x%"
    .pt3 <- "Tune the p-value (or log p-value) threshold using a value between the p-value (or log p-value) of P and of the next protein below in the list."
    MagellanNTK::mod_popover_for_help_server("modulePopover_pValThreshold",
      title = h4("Significant threshold"),
      content = HTML(paste0(.head, "<br>", 
        "<ul>", 
        "<li>", .pt1, "</li>", 
        "<li>", .pt2, "</li>",
        "<li>", .pt3, "</li>",
        "</ul>"))
    )
    
    
    
    output$warn_NULL_fdr_UI <- renderUI({
      req(is.null(fdr()))
      
      mod_SweetAlert_server('sweetAlert_fdr',
        text = "FDR cannot be computed which this threshold.",
        showClipBtn = TRUE,
        type = 'warning')
    })
    
    output$thresholdType_UI <- renderUI({

      radioButtons(ns('thresholdType'), NULL, 
        choices = c('-log10(p-value)' = 'logpval', 'p-value' = 'pval'),
        selected = threshold.type())
    })
    
    
    observeEvent(input$thresholdType, {threshold_type(input$thresholdType)})
    
    output$showFDR_UI <- renderUI({
      req(fdr())
      txt <- "FDR = NA"
      if (!is.infinite(fdr())) {
        txt <- paste0("FDR = ", round(100 * fdr(), digits = 2), ' %' )
      }
      p(txt)
    })
    
    output$warn_text1_UI <- renderUI({
      req(0 > as.numeric(input$text2))
      p(style='color: red;', 'Must be greater than 0.')
    })
    
    output$warn_text2_UI <- renderUI({
      req(0 > as.numeric(input$text1) || as.numeric(input$text1) > 1)
      p(style='color: red;', 'Must be between 0 and 1.')
    })
    
    output$text1_UI <- renderUI({
      pval_init()
      textInput(ns('text2'), NULL, 
        value = -log10(pval_init()), 
        width = '100px')
    })
    
    output$text2_UI <- renderUI({
      pval_init()
      textInput(ns('text1'), NULL, 
        value = pval_init(), 
        width = '100px')
    })
    
    observe({
      shinyjs::toggleState('text2', condition = input$thresholdType == 'logpval')
      shinyjs::toggleState('text1', condition = input$thresholdType == 'pval')
    })
    
    observeEvent(input$text1, ignoreInit = TRUE, {
      req(input$thresholdType == 'pval')
      updateTextInput(session, 'text2', value = -log10(as.numeric(input$text1)))
      
    })
    
    observeEvent(input$text2, ignoreInit = TRUE, {
      req(input$thresholdType == 'logpval')
      updateTextInput(session, 'text1', value = 10^(-as.numeric((input$text2))))
    })
    
    
    observeEvent(input$ApplyThreshold, {
      dataOut(as.numeric(input$text2))
    })
    
    return(reactive({dataOut()}))
  }
  )
}



#' @rdname mod_set_pval_threshold
#' @export
#' 
mod_set_pval_threshold <- function(
    pval_init = 1,
  fdr = 0,
  threshold.type = 'logpval'){
  
  ui <- fluidPage(
    mod_set_pval_threshold_ui("Title")
    )
  
server <- function(input, output) {
  
  rv <- reactiveValues(
    logpval = reactive({NULL}))
  
  #output$test <- renderUI({
    rv$logpval <- mod_set_pval_threshold_server(id = "Title",
      pval_init = reactive({pval_init}),
      fdr = reactive({fdr}),
      threshold.type = reactive({threshold.type})
    )
    
  #})
  
  observeEvent(req(rv$logpval()), {
    print(rv$logpval())
  })
}

app <- shiny::shinyApp(ui, server)
}

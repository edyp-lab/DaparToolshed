#' @title xxx
#' @description xxx
#' 
#' 
#' @param data xxx
#' @param title xxx
#' @param pal xxx
#' @param pattern xxx
#' @param ... xxx
#' 
#' @name mod_mv_plots
#' 
#' 
#' @examplesIf interactive()
#' data(Exp1_R25_prot, package = 'DaparToolshedData')
#' shiny::runApp(mod_mv_plots(Exp1_R25_prot[[1]]))
#' 
NULL



#' @export
#' @rdname mod_mv_plots
#' 
mod_mv_plots_ui <- function(id) {
  ns <- NS(id)
  .style <- "display:inline-block; vertical-align: top; padding-right: 20px;"
  tagList(
    tags$div(
      tags$div( style = .style,
        highchartOutput(ns("plot_viewNAbyMean"), width = "600px")
      ),
      tags$div(style = .style,
        tagList(
          uiOutput(ns("WarnForImageNA")),
          imageOutput(ns("plot_showImageNA"), width = "600px")
        )
      )
    )
  )
}





#' @export
#' @rdname mod_mv_plots
#' 
mod_mv_plots_server <- function(id, 
  data = reactive({NULL}), 
  grp = reactive({NULL}),
  mytitle = NULL, 
  pal = reactive({NULL}), 
  pattern,
  is.enabled = reactive({TRUE})) {
  
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      grp = grp(),
      title = mytitle,
      palette = pal()
    )
    
    observeEvent(req(data()), {
      stopifnot(inherits(data(), 'SummarizedExperiment'))
      
      if(is.null(rv$grp))
        rv$grp <- paste0('grp_', seq(ncol(data())))
      
      if(is.null(rv$palette))
        rv$title <- 'my title'
      
      if(is.null(rv$palette))
        rv$palette <- GetColorsForConditions(rv$grp)
    })
    
  output$plot_viewNAbyMean <- renderHighchart({
    req(data())
    
    hc_mvTypePlot2( obj = data(),
      group = rv$grp,
      pattern = pattern,
      title = rv$title,
      pal = rv$palette
    )
  })
  
  
  
  
  output$WarnForImageNA <- renderUI({
    req(!is.null(wrapper.mvImage(data(), rv$grp)))
    wellPanel(
      p(style = 'color: red;',
        "The 'MEC plot' cannot be showed as the dataset contains empty lines.")
    )
    # tryCatch(
    #     {
    #         wrapper.mvImage(data())
    #     },
    #     warning = function(w) {
    #         #p(conditionMessage(w))
    #         p('toto')
    #     },
    #     error = function(e) {
    #         #p(conditionMessage(e))
    #         p('toto')
    #     },
    #     finally = {
    #         # cleanup-code
    #     }
    # )
  })
  
  output$plot_showImageNA <- renderImage({
      # req(wrapper.mvImage(data()))
      
      # A temp file to save the output. It will be deleted after
      # renderImage
      # sends it, because deleteFile=TRUE.
      outfile <- tempfile(fileext = ".png")
      #browser()
      png(outfile)
      wrapper.mvImage(obj = data(), group = rv$grp)
      dev.off()
      
      # Return a list
      list(
        src = outfile,
        alt = "This is alternate text"
      )
    }, deleteFile = TRUE
  )
})

 } 



#' @export
#' @rdname mod_mv_plots
#' 
mod_mv_plots <- function(
    data, 
  title = NULL,
  pal = NULL,
  pattern = NULL){
  
  ui <- fluidPage(
    mod_mv_plots_ui("mvImputationPlots_MV")
  )
    
    
  server <- function(input, output, session){
    
    #observe({
      res.imp <- mod_mv_plots_server("mvImputationPlots_MV",
      data = reactive({data}),
      title = reactive(title),
      pal = reactive(pal),
      pattern = c("Missing", "Missing POV", "Missing MEC")
      )
    #})
  }
  
  app <- shiny::shinyApp(ui, server)
}


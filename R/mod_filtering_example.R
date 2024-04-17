#' @title Filtering example Shiny module
#'
#' @description
#' xxxx
#'
#' @name filtering-example
#' 
#' @param id xxx
#' @param obj An instance of the class `SummarizedExperiment`
#' @param indices xxx
#' @param operation A character(1) that indicates whether to keep or remove
#' lines identified by indices. Available values are 'keep' (default) 
#' or 'delete'
#' @param title xxx
#' 
#' @return NA
#'
#' @examplesIf interactive()
#' data(Exp1_R25_prot, package = 'DaparToolshedData')
#' obj <- Exp1_R25_prot[[1]]
#' indices <- 1:5
#' operation <- "delete"
#' shiny::runApp(mod_filtering_example(obj, indices, operation))
#' 
NULL


#' @rdname filtering-example
#' @importFrom shinyBS bsModal
#' @importFrom shiny NS actionLink tagList radioButtons uiOutput 
#' @importFrom DT dataTableOutput
#' @import shiny
#' @export
#' 
mod_filtering_example_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    actionLink(ns("show_filtering_example"), "Preview filtering"),
    shinyBS::bsModal(ns("example_modal"),
      title = "Example preview of the filtering result.",
      size = "large",
      trigger = ns("show_filtering_example"),
      tagList(
        uiOutput(ns("show_title")),
        radioButtons(ns("run_btn"), "Example dataset",
          choices = setNames(
            nm = c("original dataset", "simulate filtered dataset")
          )
        ),
        DT::dataTableOutput(ns("example_tab_filtered"))
      ),
      tags$head(tags$style(paste0("#", ns("example_modal"), " .modal-footer{ display:none}"))),
      tags$head(tags$style(paste0("#", ns("example_modal"), " .modal-dialog{ width:1000px}"))),
      tags$head(tags$style(paste0("#", ns("example_modal"), " .modal-body{ min-height:700px}"))
      )
    )
  )
}




#' @rdname filtering-example
#' @importFrom shinyBS bsModal
#' @importFrom shiny NS renderUI moduleServer 
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
#' @import shiny
#' @export
#' 
mod_filtering_example_server <- function(id, 
  obj, 
  indices = NULL, 
  operation = 'keep', 
  title = 'myTitle') {
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    output$show_title <- renderUI({
      h3(title())
    })
    
    
    # ###############
    # # options modal
    # jqui_draggable(paste0("#","example_modal"," .modal-content"),
    #                options = list(revert=FALSE)
    # )
    # ###############
    
    # colorsTypeMV = list(MEC = 'orange',
    #                     POV = 'lightblue',
    #                     identified = 'white',
    #                     recovered = 'lightgrey',
    #                     combined = 'red')
    
    legendTypeMV <- list(
      MEC = "Missing in Entire Condition (MEC)",
      POV = "Partially Observed Value (POV)",
      identified = "Quant. by direct id",
      recovered = "Quant. by recovery",
      combined = "Combined tags"
    )
    
    
    rgb2col <- function(rgbmat) {
      ProcessColumn <- function(col) {
        rgb(rgbmat[1, col],
          rgbmat[2, col],
          rgbmat[3, col],
          maxColorValue = 255
        )
      }
      sapply(1:ncol(rgbmat), ProcessColumn)
    }
    
    
    
    DarkenColors <- function(ColorsHex) {
      # Convert to rgb
      # This is the step where we get the matrix
      ColorsRGB <- col2rgb(ColorsHex)
      
      # Darken colors by lowering values of RGB
      ColorsRGBDark <- round(ColorsRGB * 0.5)
      
      # Convert back to hex
      ColorsHexDark <- rgb2col(ColorsRGBDark)
      
      return(ColorsHexDark)
    }
    
    output$example_tab_filtered <- DT::renderDataTable({
      df <- Build_enriched_qdata(obj())
      .colDef <- range.invisible <- NULL
      is.enriched <- !isTRUE(all.equal(df, obj()))
      index2darken <- NULL
      # Darken lines that will be filtered
      if (!is.null(indices()) &&
          input$run_btn == "simulate filtered dataset") {
        if (operation() == "keep") {
          index2darken <- (1:nrow(obj()))[-indices()]
        } else if (operation() == "delete") {
          index2darken <- indices()
        }
      }   
      
      
      if (is.enriched){
        .style <- BuildColorStyles(omXplore::get_type(obj()))
        c.tags <- names(.style)
        c.colors <- unlist(.style, use.names = FALSE)
        
        range.invisible <- (((ncol(df)-1) / 2) + 2):ncol(df)
        
        for (i in index2darken) {
            df[i, range.invisible] <- paste0("darken_", df[i, range.invisible])
          }
        c.tags <- c(c.tags, paste0("darken_", c.tags))
        c.colors <- c(c.colors, DarkenColors(c.colors))
          
        .colDef <- list(
            list(
              targets = range.invisible,
              visible = FALSE
            )
          )
        }

      dt <- DT::datatable(df,
        extensions = c("Scroller"),
        options = list(
          dom = "Brtip",
          pageLength = 15,
          orderClasses = TRUE,
          autoWidth = TRUE,
          deferRender = TRUE,
          bLengthChange = FALSE,
          scrollX = 200,
          scrollY = 500,
          scroller = TRUE,
          server = FALSE
          ,columnDefs = .colDef
        )
      )


      if (is.enriched){
        dt <- dt %>%
        DT::formatStyle(
          colnames(df)[2:(((ncol(df)-1) / 2) + 1)],
          colnames(df)[range.invisible],
          backgroundColor = DT::styleEqual(c.tags, c.colors)
        )
      }

      dt
    })
  })
}





#' @rdname filtering-example
#' @importFrom shinyBS bsModal
#' @importFrom shiny fluidPage 
#' @importFrom DT renderDataTable datatable formatStyle styleEqual
#' @import shiny
#' @export
#' 
mod_filtering_example <- function(
    obj,
    indices = NULL,
    operation = 'keep',
    title = 'myTitle'){
  
ui <- fluidPage(
  mod_filtering_example_ui('tree')
)

server <- function(input, output) {
  #utils::data('Exp1_R25_prot', package='DAPARdata')
  #obj <- Exp1_R25_prot[1:20]
  # filtering.query <- list(
  #   MetacellTag = c('Missing POV', 'Missing MEC'),
  #   MetacellFilters = "WholeMatrix",
  #   KeepRemove = "delete",
  #   metacell_value_th = 1,
  #   metacell_percent_th = 0,
  #   val_vs_percent = "Count",
  #   metacellFilter_operator = ">="
  # )
  
  
  # indices <- GetIndices_MetacellFiltering(
  #   obj = obj,
  #   level = omXplore::get_type(obj),
  #   pattern = filtering.query$MetacellTag,
  #   type = filtering.query$MetacellFilters,
  #   percent = filtering.query$val_vs_percent == "Percentage",
  #   op = filtering.query$metacellFilter_operator,
  #   th = filtering.query$metacell_value_th)
  
  mod_filtering_example_server('tree',
    obj = reactive({obj}),
    indices = reactive({indices}),
    operation = reactive({operation}),
    title = reactive({title}))
}

app <- shiny::shinyApp(ui = ui, server = server)

}

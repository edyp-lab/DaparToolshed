# Module UI

#' @title   open_dataset_ui and open_dataset_server
#' 
#' @description  A shiny Module.
#' 
#' @param id xxx
#' 
#' @name open_dataset
#'
#' @keywords internal
#' 
#' @examples 
#' if (interactive()){
#' ui <- fluidPage(
#' tagList(
#'   open_dataset_ui("qf_file"),
#'   textOutput('res')
#' )
#' )
#' 
#' server <- function(input, output, session) {
#'   rv <- reactiveValues(
#'     obj = NULL
#'   )
#'   rv$obj <- open_dataset_server("qf_file")
#'   
#'   output$res <- renderText({
#'     rv$obj()
#'     paste0('Names of the datasets: ', names(rv$obj()))
#'   })
#' }
#' 
#' shinyApp(ui, server)
#' }
#' 
NULL




#' @export 
#' @rdname open_dataset
#' @importFrom shiny NS tagList 
#' @import shinyjs
#' 
open_dataset_ui <- function(id){
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "Open file", multiple = FALSE),
    actionButton(ns('load_qf_btn'), 'Load file'),
    mod_infos_dataset_ui(ns("overview"))
    #uiOutput(ns("infoAboutAggregationTool"))
  )
}


#' @rdname open_dataset
#' @export
#' @importFrom shinyjs info
#' @import QFeatures
#' 
open_dataset_server <- function(id){
  
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    
    rv.openqf <- reactiveValues(
      dataRead = NULL,
      dataOut = NULL
    )
    
    # output$infoAboutAggregationTool <- renderUI({
    #   req(rv.openqf$dataRead)
    #   if (NeedsUpdate()) {
    #     tags$div(
    #       tags$div(
    #         style = "display:inline-block; vertical-align: top;",
    #         tags$img(src = "images/Problem.png", height = 25)
    #       ),
    #       tags$div(
    #         style = "display:inline-block; vertical-align: top;",
    #         HTML("The dataset was created with a former version of ProStaR,
    #            which experimental design is not compliant with the current
    #              software functionalities. Please update the design below")
    #       )
    #     )
    #   } else {
    #     m <- match.metacell(GetMetacell(rv$current.obj),
    #       pattern = "Missing",
    #       level = GetTypeofData(rv$current.obj)
    #     )
    #     NA.count <- length(which(m))
    #     nb.empty.lines <- sum(apply(m, 1, all))
    #     
    #     tagList(
    #       tags$h3("Info"),
    #       if (rv$typeOfDataset == "protein") {
    #         tags$p("Note: the aggregation tool
    #            has been disabled because the dataset contains
    #            protein quantitative data.")
    #       },
    #       if (NA.count > 0) {
    #         tags$p("As your dataset contains missing values, you should
    #            impute them prior to proceed", br(), "
    #            to the differential analysis.")
    #       },
    #       if (nb.empty.lines > 0) {
    #         tags$p("As your dataset contains lines with no values, you
    #            should remove them with the filter", br(), " tool
    #            prior to proceed to the analysis of the data.")
    #       }
    #     )
    #   }
    # })
    
    
    
    
    ## -- Open a MSnset File --------------------------------------------
    observeEvent(input$load_qf_btn, ignoreInit = TRUE, {
      input$file
      #ClearMemory()
      #ClearUI()
      
      authorizedExtension <- "qf"
      
      warn.wrong.file <- "Warning : this file is not a QFeatures file !
       Please choose another one."
      tryCatch(
        {
          if (length(grep(GetExtension(input$file$datapath),
                          authorizedExtension,ignore.case = TRUE)) == 0) {
            warning(warn.wrong.file)
          } else {
            rv.openqf$dataRead <- readRDS(input$file$datapath)
          }
          
          if (!inherits(rv.openqf$dataRead,'QFeatures')) {
            warning(warn.wrong.file)
          }
          
          
          # rv$current.obj.name <- DeleteFileExtension(input$file$name)
          # rv$typeOfDataset <- rv$current.obj@experimentData@other$typeOfData
          # rv$indexNA <- which(is.na(Biobase::exprs(rv$current.obj)))
          # rv$updateDesign_designChecked <- check.design(Biobase::pData(rv$current.obj))
          # colnames(Biobase::fData(rv$current.obj)) <- gsub(".", "_",
          #   colnames(Biobase::fData(rv$current.obj)),
          #   fixed = TRUE
          # )
          # names(rv$current.obj@experimentData@other) <- gsub(".", "_",
          #   names(rv$current.obj@experimentData@other),
          #   fixed = TRUE
          # )
          # 
          # .protId <- rv$current.obj@experimentData@other$proteinId
          # rv$widgets$aggregation$proteinId <- .protId
          # rv$proteinId <- rv$current.obj@experimentData@other$proteinId
          # if (is.null(rv$current.obj@experimentData@other$RawPValues)) {
          #   rv$current.obj@experimentData@other$RawPValues <- FALSE
          # } else if (isTRUE(rv$current.obj@experimentData@other$RawPValues)) {
          #   .params <- rv$current.obj@experimentData@other$Params
          #   rv$method <- .params[["HypothesisTest"]]$method
          # }
          # nn <- names(rv$current.obj@experimentData@other$Params)
          # ind <- grep("HypothesisTest", nn)
          # if (length(ind) > 0) {
          #   rv$res_AllPairwiseComparisons <- Get_AllComparisons(rv$current.obj)
          #   rv$listNomsComparaison <- colnames(rv$res_AllPairwiseComparisons$logFC)
          # }
          # 
          # l.params <- list(filename = rv$current.obj.name)
          # retroCompatibility()
          # loadObjectInMemoryFromConverter()
          # 
          rv.openqf$dataOut <- rv.openqf$dataRead
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
      
    })
    
    
    
    #observeEvent(rv.openqf$dataRead, {
    mod_infos_dataset_server("overview", reactive({rv.openqf$dataRead}))
    #})
    
    reactive({rv.openqf$dataOut})
  })
  
}






#----------------------------------------------------

ui <- fluidPage(
  tagList(
    open_dataset_ui("qf_file"),
    textOutput('res')
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    obj = NULL,
    result = NULL
  )
  
 
    rv$result <- open_dataset_server("qf_file")
  
  observeEvent(req(rv$result()), {
    rv$obj <- rv$result()
  })
  
  
  output$res <- renderText({
    rv$obj
    print(rv$obj)
    paste0('Names of the datasets: ', names(rv$obj))
  })
}

shinyApp(ui, server)

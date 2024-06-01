#' @title  Help popover windows
#'
#' @description  A shiny Module.
#' 
#' @param id A `character(1)` xxx
#' @param df xxx
#' @param quanCols A vector of  
#'
#' @name mod_inputGroup
#' 
#' @return A shiny app
#'
#' @examplesIf interactive()
#' shiny::runApp(mod_inputGroup())
#' 
NULL



#' @rdname mod_inputGroup
#' @export
#' @importFrom shiny NS tagList
#' @importFrom shinyjs inlineCSS useShinyjs
#'
mod_inputGroup_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("inputGroup")),
    uiOutput(ns("checkIdentificationTab"))
  )
}


#'
#' @rdname mod_inputGroup
#'
#' @export
#'
mod_inputGroup_server <- function(id, df, quantCols) {
  pkgs.require("shinyBS")

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      dataOut = NULL
    )
    # Extra functions
    #
    #
    shinyOutput <- function(FUN, id, num, ...) {
      inputs <- character(num)
      for (i in seq_len(num)) 
        inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
      inputs
    }
    
    
    # function for dynamic inputs in DT
    shinyInput <- function(FUN, id, num, ...) {
      inputs <- character(num)
      for (i in seq_len(num))
        inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
      
      inputs
    }
    
    
    # function to read DT inputs
    shinyValue <- function(id, num) {
      unlist(lapply(seq_len(num), function(i) {
        value <- input[[paste0(id, i)]]
        if (is.null(value)) NA else value
      }))
    }
    
    #
    #
    # End of extra functions
    
    output$inputGroup <- renderUI({
      # if (is.null(input$choose_quantitative_columns) || is.null(df))
      #  return(NULL)
      
      n <- length(quantCols)
      
      input_list <- lapply(seq_len(n), function(i) {
        inputName <- ns(paste("input_", i, sep = ""))
        div(
          div(
            style = "align: center;display:inline-block; vertical-align:
          middle;padding-right: 10px;",
            p(tags$strong(paste0("Identification col. for ", quantCols[i])))
          ),
          div(
            style = "align: center;display:inline-block; vertical-align: middle;",
            selectInput(inputName, "", choices = c("None", colnames(df))
            )
          )
        )
      })
      do.call(tagList, input_list)
    })
    
    
    observeEvent(input[["input_1"]], ignoreInit = T, ignoreNULL = F, {
      n <- length(quantCols)
      lapply(seq(2, n), function(i) {
        inputName <- paste("input_", i, sep = "")
        start <- which(colnames(df) == input[["input_1"]])
        
        if (input[["input_1"]] == "None") {
          .select <- "None"
        } else {
          .select <- colnames(df)[(i - 1) + start]
        }
        updateSelectInput(session, inputName, selected = .select)
      })
    })
    
    
    
    isOk <- reactive({
      quantCols
      shinyValue("input_", length(quantCols))
      temp <- shinyValue( "input_", length(quantCols))

      res <- NULL
      if (length(which(temp == "None")) > 0) {
        txt <- "The identification method is not appropriately defined for
      each sample."
        res <- list(trigger = Timestamp(), ok = FALSE, temp = temp, txt = txt)
        rv$dataOut <- NULL
      } else {
        if (length(temp) != length(unique(temp))) {
          txt <- "There are duplicates in identification columns."
          res <- list(trigger = Timestamp(), ok = FALSE, temp = temp, txt = txt)
          rv$dataOut <- NULL
        } else {
          txt <- "Correct"
          res <- list(trigger = Timestamp(), ok = TRUE, temp = temp, txt = txt)
          rv$dataOut <- temp
        }
      }

      res
    })
    
   #   observeEvent(isOk()$trigger, ignoreInit=TRUE, {
   #     print(isOk()$temp)
   #     if (isOk()$ok)
   #       rv$dataOut <- isOk()$temp
   # })
    
    output$checkIdentificationTab <- renderUI({
      
      if (isOk()$ok)
            img <- "images/Ok.png"
          else
            img <- "images/Problem.png"

      tags$div(
        tags$div(
          tags$div(style = "display:inline-block;",
            tags$img( src = img,height = 25)
          ),
          tags$div(style = "display:inline-block;", tags$p(isOk()$txt))
        )
      )
    })
    
   reactive({rv$dataOut})
    
  })
}



#' @export
#' @rdname mod_inputGroup
#' 
mod_inputGroup <- function(){

ui <- mod_inputGroup_ui("inputGroup")

server <- function(input, output, session) {
  
  file <- system.file('extdata/Exp1_R25_prot.txt', package='DaparToolshedData')
  df <- read.csv(file, header = TRUE, sep = "\t", as.is = T)
  
  toto <- DaparToolshed::mod_inputGroup_server("inputGroup", df, colnames(df)[49:54])
    
  observeEvent(toto(), ignoreNULL=FALSE,{
    print(toto())
  })
  
}

app <- shinyApp(ui = ui, server = server)
}

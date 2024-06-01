#' @title xxxx
#'
#' @description
#' xxxx
#'
#' @name build-design
#' 
#' @param id xxx
#' @param quantCols xxx
#' 
#' @return NA
#'
#' @exampleIf interactive()
#' shiny::runApp(mod_buildDesign(letters[seq(6)]))
#'
NULL

options(shiny.reactlog=TRUE) 

#' @rdname build-design
#' @import shiny
#' @export
#' 
mod_buildDesign_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    tags$p("If you do not know how to fill the experimental design, you can
            click on the '?' next to each design in the list that appear
            once the conditions are checked or got to the ",
      actionLink(ns("linkToFaq1"), "FAQ", style = "background-color: white"),
      " page."),
    fluidRow(
      column(width = 6,
        tags$b("1 - Fill the \"Condition\" column to identify
                the conditions to compare.")
      ),
      column(width = 6, uiOutput(ns("UI_checkConditions")))
    ),
    fluidRow(
      column(width = 6, uiOutput(ns("UI_hierarchicalExp"))),
      column(width = 6, uiOutput(ns("checkDesign")))
    ),
    hr(),
    selectInput(ns("convert_reorder"), "Order by conditions ?",
                choices = setNames(nm = c("No", "Yes")),
                width = "100px"
    ),
    tags$div(
      tags$div(style = "display:inline-block; vertical-align: top;",
        uiOutput(ns("viewDesign"), width = "100%")
      ),
      tags$div(style = "display:inline-block; vertical-align: top;",
        shinyjs::hidden(div(id = "showExamples", uiOutput(ns("designExamples"))))
      )
    )
    #shinyjs::disabled(actionButton(ns('validateDesign'), 'Validate design'))
  )
}


#' @rdname build-design
#' @importFrom magrittr "%>%"
#' @export
mod_buildDesign_server <- function(id,
                                   quantCols) {
  pkgs.require("shinyBS")
  
  requireNamespace('magrittr')
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    
    rv <- reactiveValues(
      hot = data.frame(
        Sample.name = as.character(quantCols),
        Condition = rep("", length(quantCols)),
        stringsAsFactors = FALSE
      ),
      conditionsChecked = NULL,
      newOrder = NULL
    )
    
    
    dataOut <- reactiveValues(
      design = NULL
    )
    
    color_renderer <- reactive({
      rv$hot$Condition
      conds <- rv$hot$Condition
      if (length(which(conds == "")) == 0) {
        uniqueConds <- unique(conds)
      } else {
        uniqueConds <- unique(conds[-which(conds == "")])
      }
      
      nUniqueConds <- length(uniqueConds)
      pal <- ExtendPalette(nUniqueConds)
      
      txt <- "function (instance, td, row, col, prop, value, cellProperties) {
  Handsontable.renderers.TextRenderer.apply(this, arguments);"
      c <- 1
      for (i in 1:length(conds)) {
        if (conds[i] != "") {
          txt <- paste0(txt, "if(row==", (i - 1), " && col==",
                        c, ") {td.style.background = '",
                        pal[which(conds[i] == uniqueConds)], "';}"
          )
        }
      }
      txt <- paste0(txt, "}")
      
      return(txt)
    })
    
    
    
    #----------------------------------------------------------
    observeEvent(input$btn_checkConds, {
      input$convert_reorder
      
      # if (length(grep("Bio.Rep", colnames(rv$hot))) > 0) {
      #   return(NULL)
      # }
      req(!("Bio.Rep" %in% colnames(rv$hot)))
      
      if (input$convert_reorder == "Yes") {
        rv$newOrder <- order(rv$hot["Condition"])
        rv$hot <- rv$hot[rv$newOrder, ]
      }
      
      rv$conditionsChecked <- check.conditions(rv$hot$Condition)
    })
    
    
    
    
    
    observeEvent(req(input$linkToFaq1), {
      updateTabsetPanel(session, "navPage", "faqTab")
    })
    
    
    
    #-------------------------------------------------------------
    output$hot <- rhandsontable::renderRHandsontable({
      rv$hot
      input$chooseExpDesign
      
      if (is.null(rv$hot)) {
        rv$hot <- data.frame(
          Sample.name = as.character(input$choose_quantitative_columns),
          Condition = rep("", length(input$choose_quantitative_columns)),
          stringsAsFactors = FALSE
        )
      }
      
      hot <- rhandsontable::rhandsontable(
        rv$hot,
        rowHeaders = NULL,
        fillHandle = list(
          direction = "vertical",
          autoInsertRow = FALSE,
          maxRows = nrow(rv$hot)
        )
      ) %>%
        rhandsontable::hot_rows(rowHeights = 30) %>%
        rhandsontable::hot_context_menu(
          allowRowEdit = TRUE,
          allowColEdit = FALSE,
          allowInsertRow = FALSE,
          allowInsertColumn = FALSE,
          allowRemoveRow = TRUE,
          allowRemoveColumn = FALSE,
          autoInsertRow = FALSE
        ) %>%
        rhandsontable::hot_cols(renderer = color_renderer()) %>%
        rhandsontable::hot_col(col = "Sample.name", readOnly = TRUE)
      
      if (!is.null(input$chooseExpDesign)) {
        switch(input$chooseExpDesign,
               FlatDesign = {
                 if ("Bio.Rep" %in% colnames(rv$hot)) {
                   hot <- hot %>%
                     rhandsontable::hot_col(
                       col = "Bio.Rep",
                       readOnly = TRUE
                     )
                 }
               },
               twoLevelsDesign = {
                 if ("Tech.Rep" %in% colnames(rv$hot)) {
                   hot <- hot %>%
                     rhandsontable::hot_col(
                       col = "Tech.Rep",
                       readOnly = TRUE
                     )
                 }
               },
               threeLevelsDesign = {
                 if ("Analyt.Rep" %in% colnames(rv$hot)) {
                   hot <- hot %>%
                     rhandsontable::hot_col(col = "Analyt.Rep", readOnly = TRUE)
                 }
               }
        )
      }
      hot
    })
    
    
    
    #--------------------------------------------------------------------------
    observeEvent(input$hot, {rv$hot <- rhandsontable::hot_to_r(input$hot)})
    
    #----------------------------------------------------------
    output$UI_checkConditions <- renderUI({
      req(rv$hot)
      rv$conditionsChecked
      input$convert_reorder
      
      if ((sum(rv$hot$Condition == "") == 0) && (input$convert_reorder != "None")) {
        tags$div(
          tags$div(style = "display:inline-block;",
            actionButton(ns("btn_checkConds"), "Check conditions", 
              class = "btn btn-primary")
          ),
          tags$div(style = "display:inline-block;",
            if (!is.null(rv$conditionsChecked)) {
              if (isTRUE(rv$conditionsChecked$valid)) {
                txt <- "<img src=\"images/Ok.png\" height=\"24\"></img>Correct conditions."
              } else {
                txt <- "<img src=\"images/Problem.png\" height=\"24\"></img><font color=\"red\">Invalid conditions."
              }
              tagList(
                tags$div(style = "display:inline-block;", HTML(txt)),
                if (!isTRUE(rv$conditionsChecked$valid))
                  tags$p(rv$conditionsChecked$warn)

              )
            }
          )
        )
      } else {
        tagList(
          br(), br(), br(), br()
        )
      }
    })
    
    
    
    #--------------------------------------------------------------------------
    output$UI_hierarchicalExp <- renderUI({
      req(rv$conditionsChecked)
      req(rv$conditionsChecked$valid)
      
      tagList(
          div(
            div(style = "display:inline-block; vertical-align: middle;",
              tags$b("2 - Choose the type of experimental design and complete it accordingly")
            ),
            div(style = "display:inline-block; vertical-align: middle;",
              tags$button(id = "btn_helpDesign", tags$sup("[?]"),
                class = "Prostar_tooltip"
              )
            )
          ),
          radioButtons(ns("chooseExpDesign"), "",
                       choices = c(
                         "Flat design (automatic)" = "FlatDesign",
                         "2 levels design (complete Bio.Rep column)" = "twoLevelsDesign",
                         "3 levels design (complete Bio.Rep and Tech.Rep columns)" = "threeLevelsDesign"
                       ),
                       selected = character(0)
          )
        )

    })
    
    
    
    
    
    
    #------------------------------------------------------------------------------
    output$viewDesign <- renderUI({
      #req(!(rv$designSaved))
      
      tagList(
        h4("Design"),
        rhandsontable::rHandsontableOutput(ns("hot"))
      )
    })
    
    #------------------------------------------------------------------------------
    output$designExamples <- renderUI({
      req(input$chooseExpDesign)
      
      switch(input$chooseExpDesign,
             FlatDesign = {
               tags$p("There is nothing to do for the flat design: the 'Bio.Rep'
           column is already filled.")
             },
             twoLevelsDesign = {
               tagList(
                 h4("Example for a 2-levels design"),
                 DaparToolshed::mod_designExample_server("buildDesignExampleTwo", 2),
                 DaparToolshed::mod_designExample_ui(ns("buildDesignExampleTwo"))
               )
             },
             threeLevelsDesign = {
               tagList(
                 h4("Example for a 3-levels design"),
                 DaparToolshed::mod_designExample_server("buildDesignExampleThree", 3),
                 DaparToolshed::mod_designExample_ui(ns("buildDesignExampleThree"))
               )
             }
      )
    })
    
    
    #------------------------------------------------------------------------------
    observe({
      shinyjs::onclick("btn_helpDesign", {
        shinyjs::toggle(id = "showExamples", anim = TRUE)
      })
    })
    
    #------------------------------------------------------------------------------
    observeEvent(input$chooseExpDesign, {
      rv$hot
      rv$designChecked <- NULL
      switch(input$chooseExpDesign,
             FlatDesign = {
               rv$hot <- data.frame(rv$hot[, 1:2],
                                    Bio.Rep = seq(1:nrow(rv$hot)),
                                    stringsAsFactors = FALSE
               )
             },
             twoLevelsDesign = {
               rv$hot <- data.frame(rv$hot[, 1:2],
                                    Bio.Rep = rep("", nrow(rv$hot)),
                                    Tech.Rep = seq(1:nrow(rv$hot)),
                                    stringsAsFactors = FALSE
               )
             },
             threeLevelsDesign = {
               rv$hot <- data.frame(rv$hot[, 1:2],
                                    Bio.Rep = rep("", nrow(rv$hot)),
                                    Tech.Rep = rep("", nrow(rv$hot)),
                                    Analyt.Rep = seq(1:nrow(rv$hot)),
                                    stringsAsFactors = FALSE
               )
             }
      )
    })
    
    
    
    
    
    
    
    #--------------------------------------------------------------------------
    observeEvent(input$btn_checkDesign, {
      rv$designChecked <- check.design(rv$hot)
    })
    
    #--------------------------------------------------------------------------
    output$checkDesign <- renderUI({
      req(input$chooseExpDesign)
      rv$designChecked
      req(rv$conditionsChecked)
      
      req(rv$conditionsChecked$valid)
      
      switch(isolate({
        input$chooseExpDesign
      }),
      FlatDesign = {},
      twoLevelsDesign = {
        if (sum(rv$hot$Bio.Rep == "") > 0) {
          return(NULL)
        }
      },
      threeLevelsDesign = {
        if ((sum(rv$hot$Bio.Rep == "") + sum(rv$hot$Tech.Rep == "")) > 0) {
          return(NULL)
        }
      }
      )
      
      
      tags$div(
        tags$div(style = "display:inline-block;",
          actionButton(ns("btn_checkDesign"), "Check design", 
            class = "btn btn-primary")
        ),
        tags$div(
          style = "display:inline-block;",
          if (!is.null(rv$designChecked)) {
            if (isTRUE(rv$designChecked$valid)) {
              shinyjs::enable("validateDesign")
              img <- "images/Ok.png"
              txt <- "Correct design"
            } else {
              img <- "images/Problem.png"
              txt <- "Invalid design"
            }
            
            
            tagList(
              tags$div(
                tags$div(style = "display:inline-block;",
                  tags$img(src = img, height = 25)
                ),
                tags$div(style = "display:inline-block;",
                  tags$p(txt)
                )
              ),
              if (!isTRUE(rv$designChecked$valid)) {
                shinyjs::disable("validateDesign")
                warn.txt <- unique(rv$designChecked$warn)
                tags$ul(
                  lapply(
                    warn.txt,
                    function(x) {
                      tags$li(x)
                    }
                  )
                )
              } else {
                shinyjs::enable("validateDesign")
              }
            )
          } else {
            shinyjs::disable("validateDesign")
          }
        )
      )
    })
    
    
    
    
    observeEvent(rv$designChecked$valid, {
      dataOut$trigger <- MagellanNTK::Timestamp()
      if (isTRUE(rv$designChecked$valid))
        dataOut$design <- rv$hot
      else
        dataOut$design <- NULL
    })
    
    
    reactive({dataOut})
  })
} # end of 











#' @rdname build-design
#' @export
mod_buildDesign <- function(quantCols) {
  
  ui <- mod_buildDesign_ui('buildDesign')
  server <- function(input, output, session) {
    mod_buildDesign_server('buildDesign', quantCols)
  }
  
  app <- shinyApp(ui, server)
}












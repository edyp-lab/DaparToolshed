# Module UI

#' @title   mod_settings_ui and mod_settings_server
#'
#' @description  A shiny Module.
#'
#' @param id shiny id
#'
#' @param input internal
#'
#' @param output internal
#'
#' @param session internal
#'
#' @return NA
#'
#' @rdname mod_settings
#'
#' @keywords internal
#'
#' @export
#'
#' @importFrom shiny NS tagList
#'
mod_settings_ui <- function(id) {
    ns <- NS(id)
    tagList(
        tabsetPanel(
            tabPanel(
                "Miscallenous",
                div(
                    div(
                        style = "display:inline-block; 
                        vertical-align: middle; padding-right: 20px;",
                        mod_popover_for_help_ui(ns("modulePopover_numPrecision"))
                    ),
                    div(
                        style = "display:inline-block; 
                        vertical-align: middle;",
                        uiOutput(ns("settings_nDigits_UI"))
                    )
                ),
                br(), hr(),
                tags$p(style = "font-size: 18px;", 
                    tags$b("Figure export options")),
                tagList(
                    tags$div(
                        style = "display:inline-block; 
                        vertical-align: middle; padding-right: 40px;",
                        selectInput(ns("sizePNGplots"), 
                            "Size of images (PNG)", 
                            choices = c("1200 * 800"), 
                            width = "150px")
                    ),
                    tags$div(
                        style = "display:inline-block; 
                        vertical-align: middle; padding-right: 40px;",
                        selectInput(ns("resoPNGplots"), 
                            "Resolution", 
                            choices = c(150), 
                            width = "100px")
                    )
                )
            ),
            tabPanel(
                "Colors",
                div(id = "showInfoColorOptions", 
                    tags$p("Color customization is available after data 
                        loading only.")),
                uiOutput(ns("defineColorsUI"))
            ),
            tabPanel("Plots", uiOutput(ns("gradient_ui")))
        )
    )
    # )
}



# Module Server

#' @rdname mod_settings
#'
#' @param id xxx
#' @param obj xxx
#'
#' @return NA
#'
#' @export
#'
#' @keywords internal
#'
#' @import shiny
#'
#' @import highcharter
#' @importFrom shinyjs onclick toggle
#'
#'
mod_settings_server <- function(id, obj) {
    
  
  pkgs.require(c("colourpicker", "shinyBS"))    
  
    moduleServer(id, function(input, output, session) {
        ns <- session$ns


        grey <- "#FFFFFF"
        orangeProstar <- "#E97D5E"

        observe({
            req(obj())

            if (!is(obj(), "QFeatures")) {
                warning("mod_settings: 'obj()' is not of class 'QFeatures'.")
                return(NULL)
            }
            rv.settings$conditions <- colData(obj())[["Condition"]]
        })


        rv.settings <- reactiveValues(
            nDigits = 10,
            conditions = c("A", "B", "A", "B", "B", "B"),
            colorsVolcanoplot = reactiveValues(
                In = orangeProstar,
                Out = "lightgrey"
            ),
            colorsTypeMV = reactiveValues(
                MEC = orangeProstar,
                POV = "lightblue"
            ),
            choosePalette = "Dark2",
            typeOfPalette = "predefined",
            examplePalette = NULL,
            basePalette = NULL,
            defaultGradientRate = 0.9,
            legDS = NULL,
            corrMatrixGradient = 0.9,
            legDS_Violinplot = NULL,
            legendForSamples = NULL
        )


        listBrewerPalettes <- c(
            "Dark2 (qualit.)" = "Dark2",
            "Accent (qualit.)" = "Accent",
            "Paired (qualit.)" = "Paired",
            "Pastel1 (qualit.)" = "Pastel1",
            "Pastel2 (qualit.)" = "Pastel2",
            "Set1 (qualit.)" = "Set1",
            "Set2 (qualit.)" = "Set2",
            "Set3 (qualit.)" = "Set3",
            "BrBG (diverging)" = "BrBG",
            "PiYG (diverging)" = "PiYG",
            "PRGn (diverging)" = "PRGn",
            "PuOr (diverging)" = "PuOr",
            "RdBu (diverging)" = "RdBu",
            "RdGy (diverging)" = "RdGy",
            "RdYlBu (diverging)" = "RdYlBu",
            "RdYlGn (diverging)" = "RdYlGn",
            "Spectral (diverging)" = "Spectral"
        )


        output$gradient_ui <- renderUI({
            sliderInput(ns("defaultGradientRate"),
                "Default color gradient for correlation matrix",
                min = 0, 
                max = 1, 
                value = rv.settings$defaultGradientRate, 
                step = 0.01
            )
        })



        mod_popover_for_help_server("modulePopover_numPrecision",
            data = list(
                title = HTML(
                    paste0("<strong><font size=\"4\">
                        Numerical precisions</font></strong>")),
                content = "Set the number of decimals to display for 
                numerical values."
            )
        )

        output$settings_nDigits_UI <- renderUI({
            numericInput(ns("settings_nDigits"), "", 
                value = rv.settings$nDigits, 
                min = 0, 
                width = "100px")
        })

        observeEvent(input$settings_nDigits, {
            rv.settings$nDigits <- input$settings_nDigits
        })




        ########################
        observeEvent(c(rv.settings$typeOfPalette, rv.settings$choosePalette), {
            rv.settings$choosePalette
            nbConds <- length(unique(rv.settings$conditions))
            nbColors <- max(3, nbConds)
            rv.settings$basePalette <- NULL

            switch(rv.settings$typeOfPalette,
                predefined = {
                    rv.settings$basePalette <- RColorBrewer::brewer.pal(nbColors, 
                        rv.settings$choosePalette)[seq_len(nbConds)]
                },
                custom = {
                    for (i in seq_len(nbConds)) {
                        rv.settings$basePalette <- c(
                            rv.settings$basePalette,
                            input[[ns(paste0("customColorCondition_", i))]]
                        )
                        if (is.null(rv.settings$basePalette)) {
                            return(NULL)
                        }
                    }
                }
            )

            SetExamplePalette()
        })



        # observeEvent(input$shinythemeSelector,{
        #   tags$script("$('#shinythemeSelector')\n  .on('change', 
        #   function(el) {\n      curThemePath = 'shinythemes/css/' + 
        #   curTheme + '.min.css';\n    }\n\n    // Find the <link> 
        #   element with that has the bootstrap.css\n 
        #   var $link = $('link').filter(function() {\n
        #   var theme = $(this).attr('href');\n 
        #   theme = theme.replace(/^.*\\//, '').replace(/(\\.min)?\\.css$/, '');\n
        #   return $.inArray(theme, allThemes) !== -1;\n    });\n\n
        #   // Set it to the correct path\n
        #   $link.attr('href', curThemePath);\n  });")
        #   #theme = shinytheme(input$shinythemeSelector)
        # })

        observe({
            shinyjs::onclick("btn_configConditionsColors", {
                shinyjs::toggle(id = "defineColorsForConditionsUI", anim = TRUE)
            })

            shinyjs::onclick("btn_configMVColors", {
                shinyjs::toggle(id = "defineColorsForMVUI", anim = TRUE)
            })


            shinyjs::onclick("btn_configVolcanoColors", {
                shinyjs::toggle(id = "defineColorsForVolcanoUI", anim = TRUE)
            })


            shinyjs::onclick("btn_configFCColors", {
                shinyjs::toggle(id = "defineColorsForFCUI", anim = TRUE)
            })
        })



        ##########
        output$defineColorsUI <- renderUI({
            shinyBS::bsCollapse(
                id = "collapseExample", open = "",
                shinyBS::bsCollapsePanel("Colors for conditions",
                    uiOutput(ns("defineColorsForConditionsUI")),
                    style = "primary"
                ),
                shinyBS::bsCollapsePanel("Colors for missing values", tagList(
                    colourpicker::colourInput(ns("colMEC"), 
                        "Select colour for MEC", 
                        orangeProstar, 
                        showColour = "background"),
                    colourpicker::colourInput(ns("colPOV"), 
                        "Select colour for POV", "lightblue", 
                        showColour = "background")
                ), style = "primary"),
                shinyBS::bsCollapsePanel("Colors for volcanoplots",
                    colourpicker::colourInput(ns("colVolcanoIn"), 
                        "Select colour for selected entities",
                        rv.settings$colorsVolcanoplot$In,
                        showColour = "background"
                    ),
                    colourpicker::colourInput(ns("colVolcanoOut"), 
                        "Select colour for filtered out entities", 
                        rv.settings$colorsVolcanoplot$Out, 
                        showColour = "background"),
                    style = "primary"
                ),
                shinyBS::bsCollapsePanel("logFC distribution", 
                    "todo", 
                    style = "primary")
            )
        })




        output$defineColorsForConditionsUI <- renderUI({
            tagList(
                fluidRow(
                    column(width = 3, radioButtons(ns("typeOfPalette"),
                        "Type of palette for conditions",
                        choices = c(
                            "predefined" = "predefined",
                            "custom" = "custom"
                        ),
                        selected = rv.settings$typeOfPalette
                    )),
                    column(
                        width = 6, uiOutput(ns("predefinedPaletteUI")),
                        uiOutput(ns("customPaletteUI"), width = "200px")
                    )
                ),
                highcharter::highchartOutput(ns("displayPalette")),
                hr()
            )
        })

        observeEvent(input$choosePalette, {
            rv.settings$choosePalette <- input$choosePalette
        })


        observeEvent(input$typeOfPalette, {
            rv.settings$typeOfPalette <- input$typeOfPalette
        })


        output$predefinedPaletteUI <- renderUI({
            rv.settings$typeOfPalette

            if (rv.settings$typeOfPalette == "custom") {
                return(NULL)
            }

            selectInput(ns("choosePalette"), "Predefined palettes",
                choices = listBrewerPalettes,
                selected = rv.settings$choosePalette, width = "200px"
            )
        })

        output$customPaletteUI <- renderUI({
            rv.settings$typeOfPalette
            if (rv.settings$typeOfPalette == "predefined") {
                return(NULL)
            }

            mod_observe_dynamic_colourPicker_input_ui(ns("colourPickerInputs"))
        })


        rv.settings$dynColors <- mod_observe_dynamic_colourPicker_input_server(
            "colourPickerInputs",
            n = reactive({length(unique(rv.settings$conditions))}),
            label = reactive({unique(rv.settings$conditions)})
            )


        SetExamplePalette <- reactive({
            .conds <- rv.settings$conditions
            nbConds <- length(unique(.conds))
            for (i in seq_len(nbConds)) {
                .ind <- which(.conds == unique(.conds)[i])
                rv.settings$examplePalette[.ind] <- rv.settings$basePalette[i]
            }
        })

        observe({
            req(rv.settings$conditions)
            req(rv.settings$dynColors())
            rv.settings$typeOfPalette
            if (rv.settings$typeOfPalette == "predefined") {
                return(NULL)
            }

            rv.settings$basePalette <- rv.settings$dynColors()
            SetExamplePalette()
        })

        observeEvent(input$colMEC, {
            rv.settings$colorsTypeMV$MEC <- input$colMEC
        })

        observeEvent(input$colPOV, {
            rv.settings$colorsTypeMV$POV <- input$colPOV
        })

        observeEvent(input$colVolcanoIn, {
            rv.settings$colorsVolcanoplot$In <- input$colVolcanoIn
        })

        observeEvent(input$colVolcanoOut, {
            rv.settings$colorsVolcanoplot$Out <- input$colVolcanoOut
        })



        output$displayPalette <- highcharter::renderHighchart({
            rv.settings$examplePalette

            .conds <- rv.settings$conditions
            nbConds <- length(unique(.conds))
            df <- data.frame(y = abs(10 + stats::rnorm(length(.conds))))

            highcharter::highchart() %>%
                dapar_hc_chart(chartType = "column") %>%
                highcharter::hc_add_series(data = df, 
                    type = "column", 
                    colorByPoint = TRUE) %>%
                highcharter::hc_colors(rv.settings$examplePalette) %>%
                highcharter::hc_plotOptions(
                    column = list(
                        stacking = "normal"), 
                    animation = list(duration = 1)
                    ) %>%
                highcharter::hc_legend(enabled = FALSE) %>%
                highcharter::hc_yAxis(labels = FALSE, 
                    title = list(text = "")) %>%
                highcharter::hc_xAxis(categories = seq_len(nbConds), 
                    title = list(text = ""))
        })

        return(reactive({
            rv.settings
        }))
    })
}

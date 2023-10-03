#' @title  Help popover windows
#'
#' @description  A shiny Module.
#' 
#' @param id A `character(1)` xxx
#' @param title A `character()` xxx
#' @param content A `character()` xxx 
#'
#' @name mod_helpPopover
#' 
#' @return NA
#'
#' @example inst/extdata/examples/ex_mod_helpPopover.R
NULL



#' @rdname mod_helpPopover
#' @export
#' @importFrom shiny NS tagList
#' @importFrom shinyjs inlineCSS useShinyjs
#'
mod_helpPopover_ui <- function(id) {
    ns <- NS(id)
    fluidPage(
        shinyjs::useShinyjs(),
        shinyjs::inlineCSS(pop_css),
        div(
            div(
                # edit1
                style = "display:inline-block; vertical-align: middle; 
                padding-bottom: 5px;",
                uiOutput(ns("write_title_ui"))
            ),
            div(
                style = "display:inline-block; vertical-align: middle;
                padding-bottom: 5px;",
                uiOutput(ns("dot")),
                uiOutput(ns("show_Pop"))
            )
        )
    )
}


#'
#' @rdname mod_helpPopover
#'
#' @export
#'
mod_helpPopover_server <- function(id, title, content) {
  pkgs.require("shinyBS")
  
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        output$write_title_ui <- renderUI({
            HTML(paste0("<strong><font size=\"4\">", title, "</font></strong>"))
        })

        output$dot <- renderUI({
            tags$button(tags$sup("?"), class = "custom_tooltip")
        })

        output$show_Pop <- renderUI({
            shinyBS::bsTooltip(ns("dot"),
                content,
                trigger = "hover"
            )
        })
    })
}


pop_css <- "button.custom_tooltip {
    background:none;
    color: #2EA8B1;
    border:none;
    padding-left:1ch;
    font: inherit;
    /*border is optional*/
        /*border-bottom:1px solid #444;*/
    cursor: pointer;
    font-weight: bold;
    display: inline-block;
    padding:0;
}

button.custom_tooltip_white {
    background:none;
    color: white;
    border:none;
    padding-left:1ch;
    font: inherit;
    /*border is optional*/
        /*border-bottom:1px solid #444;*/
    cursor: pointer;
    font-weight: bold;
    display: inline-block;
    padding:0;
}

.input-color {
    position: relative;
}
.input-color input {
    padding-left: 15px;
    border: 0px;
    background: transparent;
}
.input-color .color-box {
    width: 15px;
    height: 15px;
    display: inline-block;
    background-color: #ccc;
    position: absolute;
    left: 5px;
    top: 5px;

}"

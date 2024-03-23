#' @title xxx
#'
#' @description
#' xxxxx
#'
#' @name mod_filterExample
#' @return NA
#'
#' @examples
#' if (interactive()) {
#'     data(ft_na, package='DaparToolshed')
#'     ui <- mod_filterExample_ui("example")
#'
#'     server <- function(input, output, session) {
#'         mod_filterExample_server(
#'             id = "example",
#'             objBefore = reactive({ft_na[[1]]}),
#'             objAfter = reactive({ft_na[[1]][-c(2, 6)]}),
#'             query = reactive({"query"})
#'         )
#'     }
#'
#'     shinyApp(ui = ui, server = server)
#' }
NULL


#' @param id
#'
#' @rdname mod_filterExample
#'
#' @export
#'
mod_filterExample_ui <- function(id) {
  pkgs.require(c("shinyBS", "DT"))
  
    ns <- NS(id)

    fluidPage(
        actionButton(ns("show_filtering_example"),
            "Preview filtering",
            class = actionBtnClass
        ),
        shinyBS::bsModal(ns("example_modal"),
            title = "Example preview of the filtering result.",
            size = "large",
            trigger = ns("show_filtering_example"),
            tagList(
                uiOutput(ns("show_text")),
                radioButtons(ns("run_btn"), "Example dataset",
                    choices = stats::setNames(nm = c("original dataset", 
                        "simulate filtered dataset"))
                ),
                dataTableOutput(ns("example_tab_filtered"))
            ),
            tags$head(
                tags$style(
                    paste0("#", ns("example_modal"), 
                        " .modal-footer{ display:none}"))),
            tags$head(
                tags$style(
                    paste0("#", ns("example_modal"), 
                        " .modal-dialog{ width:1000px}"))),
            tags$head(
                tags$style(
                    paste0("#", ns("example_modal"), 
                        " .modal-body{ min-height:700px}")))
        )
    )
}


#' @param id A `character(1)`
#' @param objBefore An instance of the class `SummarizedExperiment`
#' @param objAfter xxx
#' @param query A `character()`
#'
#' @rdname mod_filterExample
#'
#' @export
#'
mod_filterExample_server <- function(id,
    objBefore,
    objAfter,
    query) {
    
    
  pkgs.require(c("grDevices", "DT"))

    moduleServer(
        id,
        function(input, output, session) {
            indices <- reactiveVal()

            observe({
                req(objBefore(), objAfter())
                id.Before <- rowData(objBefore())[, idcol(objBefore())]
                id.After <- rowData(objAfter())[, idcol(objAfter())]
                indices(which(is.na(match(id.Before, id.After))))
            })

            output$show_text <- renderUI({
                h3(query())
            })


            # ###############
            # # options modal
            # jqui_draggable(paste0("#","example_modal"," .modal-content"),
            #                options = list(revert=FALSE)
            # )

            legendTypeMV <- list(
                MEC = "Missing in Entire Condition (MEC)",
                POV = "Partially Observed Value (POV)",
                identified = "Identified",
                recovered = "Recovered",
                combined = "Combined"
            )


            rgb2col <- function(rgbmat) {
                ProcessColumn <- function(col) {
                    grDevices::rgb(rgbmat[1, col],
                        rgbmat[2, col],
                        rgbmat[3, col],
                        maxColorValue = 255
                    )
                }
                sapply(seq(ncol(rgbmat)), ProcessColumn)
            }



            DarkenColors <- function(ColorsHex) {
                # Convert to rgb
                # This is the step where we get the matrix
                ColorsRGB <- grDevices::col2rgb(ColorsHex)

                # Darken colors by lowering values of RGB
                ColorsRGBDark <- round(ColorsRGB * 0.5)

                # Convert back to hex
                ColorsHexDark <- rgb2col(ColorsRGBDark)

                return(ColorsHexDark)
            }

            output$example_tab_filtered <- DT::renderDataTable({

                # Build df to integrate qMetacell values
                .id <- idcol(objBefore())
                df <- cbind(
                    keyId = rowData(objBefore())[, .id],
                    round(assay(objBefore()), digits = 3),
                    qMetacell(objBefore())
                )

                # browser()
                colors <- custom_qMetacell_colors()
                c.tags <- names(colors)
                c.colors <- unname(unlist(colors))
                range.invisible <- c(((2 + (ncol(df) - 1) / 2)):ncol(df))


                .cond1 <- input$run_btn == "simulate filtered dataset"
                if (!is.null(indices()) && .cond1) {
                    index2darken <- indices()

                    for (i in index2darken) {
                        df[i, range.invisible] <- paste0("darken_", 
                            df[i, range.invisible])
                    }

                    c.tags <- c(c.tags, paste0("darken_", c.tags))
                    c.colors <- c(c.colors, DarkenColors(c.colors))
                }

                DT::datatable(df,
                    # rownames = FALSE,
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
                        server = FALSE,
                        columnDefs = list(
                            list(
                                targets = range.invisible,
                                visible = FALSE
                            )
                        )
                    )
                ) %>%
                    DT::formatStyle(
                        colnames(df)[2:(1 + (ncol(df) - 1) / 2)],
                        colnames(df)[range.invisible],
                        backgroundColor = DT::styleEqual(c.tags, c.colors)
                    )
            })
        }
    )
}

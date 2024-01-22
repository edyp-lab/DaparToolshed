#' @title Distribution of Observed values with respect to intensity values
#'
#' @description
#'
#' This method shows density plots which represents the repartition of
#' Partial Observed Values for each replicate in the dataset.
#' The colors correspond to the different conditions (Condition in colData of
#' the object of class `QFeatures`).
#' The x-axis represent the mean of intensity for one condition and one
#' entity in the dataset (i. e. a protein)
#' whereas the y-axis count the number of observed values for this entity
#' and the considered condition.
#'
#' @section xxx:
#'
#' xxxx
#'
#' @section missing value image:
#' Plots a heatmap of the quantitative data. Each column represent one of
#' the conditions in the object of class \code{MSnSet} and
#' the color is proportional to the mean of intensity for each line of
#' the dataset.
#' The lines have been sorted in order to visualize easily the different
#' number of missing values. A white square is plotted for missing values.
#'
# #' @seealso The package `DaparToolshed` for the 'qMetadata' metadata.
#'
#' @name plot-mv
#' @return A plot
#'
#' @param id A `character(1)` which is the 'id' of the shiny module.
#' @param vData An instance of a class `DaparViz`.
#' @param pal.name  A `character(1)` which is the name of the palette (from
#' the package [RColorBrewer] to use.
#' @author Samuel Wieczorek, Enora Fremy
#'
#' @example inst/extdata/examples/ex_mv_imputation.R
#' 
NULL


#' @rdname plot-mv
#' @export
#' @importFrom shiny NS tagList plotOutput
#' @importFrom highcharter highchartOutput
#'
mv_imputation_ui <- function(id) {
    ns <- NS(id)
    tagList(
        tags$div(
            tags$div(
                style = "display:inline-block; vertical-align: top; padding-right: 20px;",
                highcharter::highchartOutput(ns("nabyMean_ui"), width = "600px")
            ),
            tags$div(
                style = "display:inline-block; vertical-align: top; padding-right: 20px;",
                plotOutput(ns("imageNA_ui"), width = "600px")
            )
        )
    )
}

#'
#' @rdname plot-mv
#' @export
#' @importFrom highcharter renderHighchart
#'
mv_imputation_server <- function(id,
                                 vData = reactive({NULL}),
                                 pal.name = reactive({NULL})) {
    
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Initialisation of the module
        observe({
            stopifnot(inherits(vData(), "DaparViz"))
        })


        output$nabyMean_ui <- renderHighchart({
            req(vData())

            withProgress(message = "Making MV Intensity plot", value = 100, {
              hc_mvTypePlot2(DaparViz = vData(), pattern = "Missing POV")
            })
        })


        output$imageNA_ui <- renderPlot({
            req(vData())
            withProgress(message = "Making MV Heatmap plot", value = 100, {
                mv.image(DaparViz = vData())
            })
        })
    })
}

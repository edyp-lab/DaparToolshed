#' library(QFeatures)
#' data(subR25prot, package='omXplore')
#' data <- assay(subR25prot, 1)
#' conds <- colData(subR25prot)$Condition
#'
#' #-----------------------------
#' # xxx
#' #-----------------------------
#'
#' mv.density(subR25prot[[1]], conds, pattern = "Missing MEC")
#'
#' mv.mec.heatmap(subR25prot[[1]], conds)
#'
#' #-----------------------------
#' # xxx
#' #-----------------------------
#'
#' if (interactive()) {
#'     data(subR25prot, package='omXplore')
#'
#'     ui <- mod_mv_imputation_ui("plot")
#'
#'     server <- function(input, output, session) {
#'         mod_mv_imputation_server("plot",
#'             se = reactive({
#'                 subR25prot[[1]]
#'             }),
#'             conds = colData(subR25prot)$Condition
#'         )
#'     }
#'
#'     shinyApp(ui = ui, server = server)
#' }
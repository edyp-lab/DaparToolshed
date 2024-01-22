#' data(ft_na, package='DaparViz')
#' data <- assay(ft_na, 1)
#' conds <- colData(ft_na)$Condition
#'
#' #-----------------------------
#' # xxx
#' #-----------------------------
#'
#' mv.density(ft_na[[1]], conds, pattern = "missing MEC")
#'
#' mv.mec.heatmap(ft_na[[1]], conds)
#'
#' #-----------------------------
#' # xxx
#' #-----------------------------
#'
#' if (interactive()) {
#'     data(ft_na, package='DaparViz')
#'
#'     ui <- mod_mv_imputation_ui("plot")
#'
#'     server <- function(input, output, session) {
#'         mod_mv_imputation_server("plot",
#'             se = reactive({
#'                 ft_na[[1]]
#'             }),
#'             conds = colData(ft_na)$Condition
#'         )
#'     }
#'
#'     shinyApp(ui = ui, server = server)
#' }
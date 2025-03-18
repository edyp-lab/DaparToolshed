#' @title Filter features of one SE based on their rowData
#'
#' @description
#'
#' The `filterFeaturesOneSE` methods enables users to filter features
#' based on a variable in their `rowData`. It is directly inspired of the
#' function `filterFeature` of the package `QFeatures`.
#' The first difference is that the filter only applies to one 
#' `SummarizedExperiment` contained in the object rather than applying on 
#' all the SE.
#' This method generates a new `SummarizedExperiment` object which is added 
#' to the `QFeatures` object. If the SE on which the filter applies is the 
#' last one of the object, then a new xxxx. If it is not the last one, the 
#' new SE is added and all the further SE are deleted. The features matching 
#' the. The filters can be provided as instances of class `AnnotationFilter` 
#' (see the package `QFeatures`) or of class `FunctionFilter` (see below).
#'
#' @section Function filters:
#'
#' The function filters are filters as defined in the
#' [DaparToolshed] package. Each filter is defined by a name (which is the 
#' name of a function) and a list which contains the parameters passed to the 
#' function. Those filters can be created with the `FunctionFilter` 
#' constructor.
#'
#' Those functions are divided into two main categories:
#'  - the one that filter on one rowData feature,
#'  - the one based on a two-dimensional information such as the adjacency 
#'  matrix
#'
#'  for the first category, all filters of class [AnnotationFilter] can be 
#'  used as they are used in `QFeatures`
#'
#'  For the second category, the package `DaparToolshed` provides filter 
#'  functions based either on the adjacency matrix:
#'  - [DaparToolshed::topnPeptides()]: xxxx
#'  - [DaparToolshed::sharedPeptides()]: xxx
#'  - [DaparToolshed::specPeptides()]: xxx
#'
#' Or based on the quantitative metadata (identification):
#'  - [DaparToolshed::qMetacellWholeMatrix()]: xxx
#'  - [DaparToolshed::qMetacellWholeLine()]: xxx
#'  - [DaparToolshed::qMetacellOnConditions()]: xxx
#'
#' @return A filtered `QFeature` object
#'
#' @author Samuel Wieczorek
#'
#' @name QFeatures-filtering-oneSE
#'
#' @rdname QFeatures-filtering-oneSE
#'
#' @aliases filterFeaturesOneSE filterFeaturesOneSE, DaparToolsehd, 
#' FunctionFilter, VariableFilter
#'
#' @examples
#' data(feat1, package = 'QFeatures')
#' ## ----------------------------------------
#' ## Creating function filters
#' ## ----------------------------------------
#'
#' #FunctionFilter('FUN',
#' #               param1 = 'value_of_param1',
#' #               param2 = 'value_of_param2')
#'
#' FunctionFilter('qMetacellWholeLine',
#'                cmd = 'delete',
#'                pattern = 'imputed POV')
#'
#' ## ----------------------------------------------------------------
#' ## Filter the last assay to keep only specific peptides. This filter
#' ## only applies on peptide dataset.
#' ## ----------------------------------------------------------------
#'
#' spec.filter <- FunctionFilter('specPeptides', list())
#' ## using a user-defined character filter
#' filterFeaturesOneSE(feat1, list(FunctionFilter('specPeptides', list())))
#'
#'
#' ## ----------------------------------------------------------------
#' ## Filter the last assay to keep only specific peptides and topn 
#' ## peptides. The two filters are run sequentially.
#' ## ----------------------------------------------------------------
#'
#' lst.filters <- list(FunctionFilter('specPeptides', list()))
#' lst.filters <- append(lst.filters,
#' FunctionFilter('topnPeptides',
#' fun = 'rowSums',
#' top = 2))
#' filterFeaturesOneSE(feat1, lst.filters)
#'
#' ## ----------------------------------------------------------------
#' ## Filter the last assay to delete peptides where, in at least one 
#' ## condition, there is less than 80% of samples marked as 'imputed POV'
#' ## ----------------------------------------------------------------
#'
#' filter <- FunctionFilter('qMetacellOnConditions',
#' cmd = 'delete',
#' mode = 'AtLeastOneCond',
#' pattern = 'imputed POV',
#' conds = colData(feat1)$Condition,
#' percent = TRUE,
#' th = 0.8,
#' operator = '<')
#'
#'  filterFeaturesOneSE(feat1, filter)
#'
#'
NULL


#' @exportClass FunctionFilter
#' @rdname QFeatures-filtering-oneSE
setClass("FunctionFilter",
    slots = c(
        name = "character",
        params = "list"
    ),
    prototype = list(
        name = character(1),
        params = list()
    )
)


#' @param name `character(1)` refering to the name of the function 
#' to apply the filter on.
#' @param ... Additional arguments
#'
#' @export FunctionFilter
#' @rdname QFeatures-filtering-oneSE
#'
#' @importFrom methods new
#' 
FunctionFilter <- function(name, ...) {
    new("FunctionFilter",
        name = name,
        params = list(...)
    )
}


#' @param object An instance of class `QFeatures` or `SummarizedExperiment`.
#'
#' @param i The index or name of the assay which features will be
#'     filtered the create the new assay.
#'
#' @param name A `character(1)` naming the new assay. Default is
#'     `newAssay`. Note that the function will fail if there's
#'     already an assay with `name`.
#'
#' @param filters A `list()` containing instances of class `AnnotationFilter` 
#' or `FunctionFilter`
#'
#' @exportMethod filterFeaturesOneSE
#' @importFrom S4Vectors metadata
#' @import AnnotationFilter
#' @rdname QFeatures-filtering-oneSE
setMethod(
    "filterFeaturesOneSE", "QFeatures",
    function(object, i, name = "newAssay", filters) {
        if (length(object) == 0) {
            return(object)
        }
        if (name %in% names(object)) {
            stop("There's already an assay named '", name, "'.")
        }
        if (missing(i)) {
            i <- length(object)
        }

        if (missing(filters)) {
            return(object)
        }
        ## Create the filtered assay
        
        new.se <- object[[i]]
   
        for (x in filters) {
          if (inherits(x, "AnnotationFilter")) {
            new.se <- filterFeatures(object, x, i)[[i]]
          } else if (inherits(x, "FunctionFilter")) {
            new.se <- do.call(x@name, append(list(object = object[[i]]), x@params))
          }
        }

        ## Add the assay to the QFeatures object
        object <- addAssay(object,
            new.se,
            name = name
        )

        if (nrow(new.se) > 0) {
            idcol <- S4Vectors::metadata(object[[i]])$idcol
            if (is.null(idcol)) {
                warning("xxx")
                S4Vectors::metadata(object[[i]])$idcol <- "_temp_ID_"
                idcol <- "_temp_ID_"
            }


            ## Link the input assay to the aggregated assay
            object <- addAssayLink(object,
                from = names(object)[i],
                to = name,
                varFrom = idcol,
                varTo = idcol
            )
        }

        return(object)
    }
)



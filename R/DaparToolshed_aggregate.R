
##' @title Aggregate an assay's quantitative features
##'
##' @description
##'
##' This function aggregates the quantitative features of an assay,
##' applying a summarization function (`fun`) to sets of features.
##' The `fcol` variable name points to a rowData column that defines
##' how to group the features during aggregate. This variable can
##' either be a vector (we then refer to an *aggregation by vector*)
##' or an adjacency matrix (*aggregation by matrix*).
##'
##' The quantitative metadata are aggregated with a function (`fun.qmeta`).
##'
##' The list of agregation methods can be obtained with the function
##' aggregateMethods()]. This function compiles both methods from the 
##' packages `DaparToolshed` and `QFeatures`.
##'
##' @param object An instance of class `QFeatures` or `SummarizedExperiment`
##'
##' @param i The index or name of the assay which features will be
##'     aggregated the create the new assay.
##'
##' @param fcol A `character(1)` naming a rowdata variable (of assay
##'     `i` in case of a `QFeatures`) defining how to aggregate the
##'     features of the assay. This variable is either a `character`
##'     or a (possibly sparse) matrix. See below for details.
##'
##' @param name A `character(1)` naming the new assay. Default is
##'     `newAssay`. Note that the function will fail if there's
##'     already an assay with `name`.
##'
##' @param fun A function used for quantitative feature
##'     aggregation. See Details for examples.
##'
##' @param ... Additional parameters passed the `fun` and `fun.qmeta`.
##'
##' @return A `QFeatures` object with an additional assay or a
##'  `SummarizedExperiment` object (or subclass thereof).
##'
##' @details xxxxxxx
##'
##' @section Iterative aggregation function:
##'
##' xxxxxx
##' xxxxx
##'
##' @section Quantitative metadata aggregation:
##'
##' xxxxxx
##' xxxx
##'
##' The function to aggregate the quantitative metadata is
##'
##' - `aggQmetadat()` xxxxx
##'
##' @seealso The *QFeatures* vignette provides an extended example and
##'     the *Aggregation* vignette, for a complete quantitative
##'     proteomics data processing pipeline.
##'
##' @name DaparToolshed-aggregate
##' 
##' @return NA
##'
##' @examples
##'
##' ## ---------------------------------------
##' ## An example QFeatures with PSM-level data
##' ## ---------------------------------------
##' data(ft, package='DaparToolshed')
##' ft
##'
##' ## Aggregate peptides into proteins
##' ## using the adjacency matrix
##' feat1 <- aggregateFeatures4Prostar(object = ft,
##' i = 1,
##' name = 'aggregated',
##' fcol = 'adjacencyMatrix',
##' fun = colSumsMat)
##' feat1
##'
##' assay(feat1[[1]])
##' assay(feat1[[2]])
##' aggcounts(feat1[[2]])
##' assay(feat1[[3]])
##' aggcounts(feat1[[3]])
##' rowData(ft[[2]])
NULL

##' @exportMethod aggregateFeatures4Prostar
##' @rdname DaparToolshed-aggregate
##' @importFrom MsCoreUtils robustSummary
setMethod(
    "aggregateFeatures4Prostar", "QFeatures",
    function(object, i, fcol, name = "newAssay",
             fun = MsCoreUtils::robustSummary, ...) {
        if (length(object) == 0) {
            return(object)
        }
        if (name %in% names(object)) {
            stop("There's already an assay named '", name, "'.")
        }
        if (missing(i)) {
            i <- length(object)
        }


        # Add stats on agregation
        aggAssay <- aggregateFeatures4Prostar(
            object = object[[i]],
            fcol = fcol,
            fun = fun,
            conds = design.qf(object)$Condition,
            ...
        )
        # browser()
        ## Add the assay to the QFeatures object
        object <- addAssay(object,
            aggAssay,
            name = name
        )
        ## Link the input assay to the aggregated assay
        addAssayLink(object,
            from = i,
            to = name,
            varFrom = fcol,
            varTo = fcol
        )
    }
)


##' @exportMethod aggregateFeatures4Prostar
##' @rdname DaparToolshed-aggregate
##' @importFrom MsCoreUtils robustSummary
setMethod(
    "aggregateFeatures4Prostar", "SummarizedExperiment",
    function(object, fcol, fun = MsCoreUtils::robustSummary, conds, ...) {
        .aggregateFeatures4Prostar(object, fcol, fun, conds, ...)
    }
)


.aggregateFeatures4Prostar <- function(object, fcol, fun, conds, ...) {
    X <- adjacencyMatrix(object)

    # add agregation of qMetacell
    # Aggregate the quantitative metdata
    #if(is.null(fun.qmeta) || is.missing(fun.qmeta))
    aggQ <- aggQmetacell(
        qMeta = qMetacell(object),
        X = adjacencyMatrix(object),
        level = typeDataset(object),
        conds = conds
    )

    ## Remove the qMetacell that should be dropped anyway and not be aggregated
    ## within QFeatures::aggregateFeatures
    rowData(object)[["qMetacell"]] <- NULL

    ## Create the aggregated assay
    aggAssay <- aggregateFeatures(object, fcol, fun, ...)

    ## Add the qMetacell to the new assay
    qMetacell(aggAssay) <- aggQ
    X.spec <- X.shared <- X

    X.spec[which(rowSums(as.matrix(X.spec)) > 1), ] <- 0
    X.shared[which(rowSums(as.matrix(X.shared)) == 1), ] <- 0
    
    .allPep <- t(as.matrix(X)) %*% !is.na(assay(object))
    .specPep <- t(as.matrix(X.spec)) %*% !is.na(assay(object))
    .sharedPep <- t(as.matrix(X.shared)) %*% !is.na(assay(object))
    rowData(aggAssay)[["allPeptidesUsed"]] <-.allPep
    rowData(aggAssay)[["specPeptidesUsed"]] <- .specPep
    rowData(aggAssay)[["sharedPeptidesUsed"]] <- .sharedPep


    ## Enrich the new assay
    typeDataset(aggAssay) <- "proteins"
    idcol(aggAssay) <- NULL


    return(aggAssay)
}




#' @title Get the type of dataset
#' @description xxx
#'
#' @param qMeta  An object of class 'SummarizedExperiment'
#' @param X xxxx
#' @param level A `character(1)` which is the type of dataset
#' @param conds A `character()` vector which is the names of conditions
#' for each sample in the dataset.
#'
#' @return xxxxx
#'
#' @examples
#' data(ft, package='DaparToolshed')
#' qMeta <- qMetacell(ft, 1)
#' X <- adjacencyMatrix(ft, 1)
#' level <- typeDataset(ft, 1)
#' conds <- colData(ft)$Condition
#' aggQmeta <- aggQmetacell(qMeta, X, level, conds)
#'
#' @rdname DaparToolshed-aggregate
#'
aggQmetacell <- function(qMeta, X, level, conds) {
    # stopifnot(inherits(object, "SummarizedExperiment"))

    rowcol <- function(meta.col, X.col) {
        meta.col[X.col > 0]
    }

    df <- data.frame(stringsAsFactors = TRUE)
    for (j in seq(ncol(qMeta))) {
        for (i in seq(ncol(X))) {
            df[i, j] <- metacombine(
                rowcol(qMeta[, j], X[, i]),
                level
            )
        }
    }

    df[df == "NA"] <- NA
    dimnames(df) <- list(colnames(X), colnames(qMeta))
    # Delete protein with only NA

    # Post processing of metacell to discover 'imputed POV', 'imputed MEC'
    df <- Set_POV_MEC_tags(conds, df, level)

    df
}


#' @export
#' @rdname DaparToolshed-aggregate
aggregateMethods <- function() {
    stats::setNames(
        c("medianPolish",
            "robustSummary",
            "colMeansMat",
            "colSumsMat"
            ),
        nm = c(
            "median Polish",
            "robust Summary",
            "col Means Mat",
            "col Sums Mat"
            )
        )
}

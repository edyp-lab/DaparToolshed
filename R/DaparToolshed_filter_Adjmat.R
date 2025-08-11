#' @title Filter a peptide assay on the basis of its adjacency matrix.
#'
#' @description
#'
#' These functions filters (delete) peptides of an assay, applying a function
#' on peptides and proteins. They can be used alone but the usual usage is to
#' create an instance of a class [FunctionFilter] and to pass it to the function
#' [filterFeaturesOneSE] in order to create a new assay, embedded into the 
#' [QFeatures] object.
#'
#' @param object An object of class `SummarizedExperiment`
#'
#' @param fun A `list()` of additional parameters
#'
#' @param top A `integer(1)` which is the number of xxx
#' 
#' @param ... Additional arguments
#' @param X xxx
#' @param qData xxx
#'
#' @details
#'
#' This function builds an intermediate matrix with scores for each peptide
#' based on 'fun' parameter. Once this matrix is built, one select the 
#' 'n' peptides which have the higher score
#'
#' The list of filter functions is given by `adjMatFilters()`:
#'
#' - `specPeptides()`: returns a new assay of class `SummazizedExperiment` 
#' with only specific peptides;
#'
#' - `sharedpeptides()`: returns a new assay of class `SummazizedExperiment` 
#' with only shared peptides;
#'
#' - `opnPeptides()`: returns a new assay of class `SummazizedExperiment` with 
#' only the 'n' peptides which best satisfies the condition. The condition is 
#' represented by functions which calculates a score for each peptide among 
#' all samples. The list of these functions is given by `topnFunctions()`:
#'
#' - `rowMedians()`: xxx;
#'
#' - `rowMeans()`: xxx;
#'
#' - `rowSums()`: xxx;
#'
#'
#' @seealso The [QFeatures-filtering-oneSE] man page for the 
#' class `FunctionFilter`.
#'
#' @name adjacency-matrix-filter
#'
#' @author Samuel Wieczorek
#' 
#' @return NA
#'
#' @examples
#' library(Matrix)
#' library(QFeatures)
#' #------------------------------------------------
#' # This function will keep only specific peptides
#' #------------------------------------------------
#'
#' f1 <- FunctionFilter("specPeptides", list())
#'
#' #------------------------------------------------
#' # This function will keep only shared peptides
#' #------------------------------------------------
#'
#' f2 <- FunctionFilter("sharedPeptides", list())
#'
#' #------------------------------------------------
#' # This function will keep only the 'n' best peptides
#' # w.r.t the quantitative sum of each peptides among
#' # all samples
#' #------------------------------------------------
#'
#' f3 <- FunctionFilter("topnPeptides", fun = "rowSums", top = 2)
#'
#' #------------------------------------------------------
#' # To run the filter(s) on the dataset, use [xxx()]
#' # IF several filters must be used, store them in a list
#' #------------------------------------------------------
#'
#' data(ft, package='DaparToolshed')
#' lst.filters <- list()
#' lst.filters <- append(lst.filters, f1)
#' lst.filters <- append(lst.filters, f3)
#'
#' ft <- filterFeaturesOneSE(
#'     object = ft,
#'     i = 1,
#'     name = "filtered",
#'     filters = lst.filters
#' )
#'
NULL


#' @export
#' @rdname adjacency-matrix-filter
#' @importFrom stats setNames
#' 
AdjMatFilters <- function() {
    stats::setNames(c(
        "allPeptides",
        "specPeptides",
        "sharedPeptides",
        "topnPeptides"
    ),
    nm = c(
        "all peptides",
        "Only specific peptides",
        "Only shared peptides",
        "N most abundant"
    )
    )
}


#' @export
#' @rdname adjacency-matrix-filter
#' @import SummarizedExperiment
#' 
allPeptides <- function(object, ...) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    stopifnot("adjacencyMatrix" %in% names(SummarizedExperiment::rowData(object)))

    return(object)
}


#' @export
#' @rdname adjacency-matrix-filter
#' @importFrom QFeatures adjacencyMatrix
specPeptides <- function(object, ...) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    stopifnot("adjacencyMatrix" %in% names(SummarizedExperiment::rowData(object)))

    X <- QFeatures::adjacencyMatrix(object)
    X.specific <- subAdjMat_specificPeptides(X)
    object <- .UpdateSEBasedOnAdjmat(object, X.specific)

    return(object)
}



#' @export
#' @rdname adjacency-matrix-filter
subAdjMat_specificPeptides <- function(X) {
    # Mask for only shared peptides
    x.spec <- as.matrix(X) != 0
    ind <- which(rowSums(x.spec) > 1)
    if (length(ind) > 0) {
        X[ind, ] <- 0
    }
    X
}


#' @export
#' @rdname adjacency-matrix-filter
sharedPeptides <- function(object, ...) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    stopifnot("adjacencyMatrix" %in% names(SummarizedExperiment::rowData(object)))

    X <- QFeatures::adjacencyMatrix(object)
    X.shared <- subAdjMat_sharedPeptides(X)
    object <- .UpdateSEBasedOnAdjmat(object, X.shared)
    return(object)
}

#' @export
#' @rdname adjacency-matrix-filter
subAdjMat_sharedPeptides <- function(X) {
    # Mask for only shared peptides
    x.spec <- as.matrix(X) != 0
    ind <- which(rowSums(x.spec) == 1)
    if (length(ind) > 0) {
        X[ind, ] <- 0
    }
    X
}

#'
#' @noRd
#' @importFrom MagellanNTK pkgs.require
.UpdateSEBasedOnAdjmat <- function(object, X) {
  
  MagellanNTK::pkgs.require('PSMatch')
  
  SummarizedExperiment::rowData(object)$adjacencyMatrix <- X
    # Identify and delete the empty lines in the dataset
    emptyLines <- which(rowSums(as.matrix(X)) == 0)
    if (length(emptyLines) > 0) {
        object <- object[-emptyLines]
    }

    # Reload the adjacency matrix after lines deletion
    # Identify empty columns in the adjacency matrix
    X <- QFeatures::adjacencyMatrix(object)
    emptyCols <- which(colSums(as.matrix(X)) == 0)

    if (length(emptyCols) > 0) {
        X <- X[, -emptyCols]
        SummarizedExperiment::rowData(object)$adjacencyMatrix <- X
        idcol <- S4Vectors::metadata(object)$parentProtId
        .val <- PSMatch::makePeptideProteinVector(X)
        SummarizedExperiment::rowData(object)[, idcol] <- .val
    }

    return(object)
}


#' @export
#' @rdname adjacency-matrix-filter
topnFunctions <- function() {
    c("rowMedians", "rowMeans", "rowSums")
}


#' @export
#' @rdname adjacency-matrix-filter
topnPeptides <- function(object, fun, top) {
    stopifnot(inherits(object, "SummarizedExperiment"))
    .names <- names(SummarizedExperiment::rowData(object))
    stopifnot("adjacencyMatrix" %in% .names)

    X <- QFeatures::adjacencyMatrix(object)
    qData <- SummarizedExperiment::assay(object)
    X.topn <- subAdjMat_topnPeptides(X, qData, fun, top)
    object <- .UpdateSEBasedOnAdjmat(object, X.topn)
    return(object)
}



#' @export
#' @rdname adjacency-matrix-filter
#' @importFrom methods as
subAdjMat_topnPeptides <- function(X, qData, fun, top) {
    if (!(fun %in% topnFunctions())) {
        warning("'fun' must be one of the following: 'rowMedians', 
            'rowMeans', 'rowSums'")
        return(NULL)
    }

    temp.X <- as(X * do.call(fun, list(qData)), "dgCMatrix")

    # Get the 'n' entities with the best score for each column
    for (c in seq_len(ncol(temp.X))) {
        v <- order(temp.X[, c], decreasing = TRUE)[seq_len(top)]
        l <- v[which((temp.X[, c])[v] != 0)]
        if (length(l) > 0) {
            X[-l, c] <- 0
        }
    }
    X
}

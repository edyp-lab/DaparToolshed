#' @title Utility funcitons to dela with QFeatures objects.
#'
#' @description
#'
#' xxxxx
#'
#' @name QFeatures-utils
#' 
#' @return NA
#'
#' @examples
#'
#' NULL
#'
NULL


#' @param object An instance of the class `QFeatures`
#' @rdname QFeatures-utils
#' @export
last_assay <- function(object) {
    object[[length(object)]]
}

#' @param object An instance of the class `QFeatures`
#' @rdname QFeatures-utils
#' @export
n_assays_in_qf <- function(object) {
    length(object)
}


#' @param obj.se An instance of the class `QFeatures`
#' @param colData A data.frame() which contains data for the future object,
#' @param metadata.qf A data.frame() which contains the metadata for the future object,
#' @param name The name of the assay in the QFeatures object
#' @rdname QFeatures-utils
#' @return An instance of QFeatures class
#' 
#' @examples
#' # example code
#' 
#' @export
#'
QFeaturesFromSE <- function(obj.se, 
  colData = data.frame(), 
  metadata.qf = data.frame(),
  name = 'myname'){
  
  stopifnot(inherits(obj.se, "SummarizedExperiment"))
  new.obj <- QFeatures::readQFeatures(
    assayData = assay(obj.se),
    colData = colData,
    metadata = metadata.qf
  )
  
  rowData(new.obj[[1]]) <- rowData(obj.se)
  metadata(new.obj[[1]]) <- metadata(obj.se)
  
  names(new.obj)[1] <- name
  return(new.obj)
}
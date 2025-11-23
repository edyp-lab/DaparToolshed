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


#' @param obj An instance of the class `QFeatures`
#' @param range A vector of integers
#' @rdname QFeatures-utils
#' @return description An instance of QFeatures class
#' @export
SubQFeatures <- function(obj, range){
  
  stopifnot(inherist(obj, "QFeatures"))
  new.obj <- QFeatures::readQFeatures(
    assayData = assay(obj[[range]]),
    colData = colData(obj),
    metadata = metadata(obj)
  )
  
  for (i in range){
    rowData(new.obj[[i]]) <- rowData(obj[[i]])
    metadata(new.obj[[i]]) <- metadata(obj[[i]])
  }
  
  return(new.obj)
}
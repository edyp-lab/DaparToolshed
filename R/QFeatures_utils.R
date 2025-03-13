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

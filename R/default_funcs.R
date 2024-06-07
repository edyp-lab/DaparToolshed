
#' @title
#' Datasets processing
#'
#' @description
#' This manual page describes manipulation methods using [list] objects. In 
# 'the following functions, if `object` is of class `list`, and optional array
#' index or name `i` can be specified to define the array (by name of
#' index) on which to operate.
#'
#' The following functions are currently available:
#'
#' - `keepDatasets(object, range)` keep datasets in object which
#' are in range
#'
#' - `addDatasets(object, dataset, name)` add the 'dataset' to the 
#' object (of type list)
#'
#' - `Save(object, file)` stores the object to a .RData file
#'
#' @details
#' The object must be of type list. Thetwo functions are implemented here for 
# 'a simple list. For other dataset classes, their implementation must be part 
#' of the package which uses MagellanNTK
#'
#' @param object An object of class `list`.
#'
#' @param range A xxxx
#'
#' @param dataset `character(1)` providing the base with respect to which
#'     logarithms are computed. Default is log2.
#'
#' @param name A `character(1)` naming the new array name.
#'
#' @return An processed object of the same class as `object`.
#'
#' @aliases keepDatasets keepDatasets,list-method
#' @aliases addDatasets addDatasets,list-method
#'
#' @name dataset-processing
#'
#' @importFrom methods setMethod new
#' 
#'
NULL


#' @title Adds a dataset to the list
#' @description This function appends a dataset in the list with customization
#' if necessary
#' @param object An instance of the class `QFeatures` 
#' Must get TRUE to inherits(object, 'QFeatures')
#' @param dataset An instance of class `SummarizedExperiment`
#' @param name the name to associate to the dataset in the final object
#' 
#' @rdname addDatasets
#' 
#' @export
addDatasets <- function(object, dataset, name) {
    stopifnot(inherits(object, "QFeatures"))
    stopifnot(inherits(dataset, "SummarizedExperiment"))
  
    object <- QFeatures::addAssay(object, dataset, name)
    
    object
}




#' @title Get a subset of the object
#' @description This function deletes the items not included in the
#' range parameter
#' @param object An instance of type list. Must get TRUE to inherits(object, 'list')
#' @param range xxx
#' 
#' @rdname keepDatasets
#' 
#' @export
#'
keepDatasets <- function(object, range = seq(length(object))) {
  stopifnot(is.Magellan.compliant(object))
  if (!is.numeric(range))
    stop('Provide numeric range of array to be processed')
  
    if (min(range) < 1 || max(range) > length(object))
      stop('Provide numeric range with values consistent with the number of assays')
    
  if (is.null(object)) {
    return()
  }
  
  object[range]
}
#' @rdname QFeatures-accessors
setGeneric("paramshistory", 
  function(object, ...) standardGeneric("paramshistory"))
setGeneric("paramshistory<-", 
  function(object, ..., value) standardGeneric("paramshistory<-"))




#' @export
#' @exportMethod paramshistory
#' @rdname QFeatures-accessors
setMethod(
  "paramshistory", "QFeatures",
  function(object, i, slotName = 'paramshistory') {
    lapply(object[[i]], paramshistory, slotName )
  }
)


#' @export
#' @exportMethod paramshistory
#' @rdname QFeatures-accessors
#' @import SummarizedExperiment
#' 
setMethod(
  "paramshistory", "SummarizedExperiment",
  function(object, slotName = "paramshistory") {
    .GetMetadataSlot(object, slotName)
  }
)




#' @export
#' @rdname QFeatures-accessors
"paramshistory<-" <- function(object, i, slotName = "paramshistory", value) {
  if (inherits(object, "SummarizedExperiment")) {
    S4Vectors::metadata(object)[[slotName]] <- value
    return(object)
  } else {
    if(inherits(object, "QFeatures")){
      if (length(i) != 1) {
        stop("'i' must be of length one. Repeat the call to add a matrix to 
            multiple assays.")
      }
      if (is.numeric(i) && i > length(object)) {
        stop("Subscript is out of bounds.")
      }
      if (is.character(i) && !(i %in% names(object))) {
        stop("Assay '", i, "' not found.")
      }
      S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
  }
  return(object)
}


#' @exportMethod paramshistory
#' @rdname QFeatures-accessors
setMethod(
  "paramshistory", "QFeatures",
  function(object, i, slotName = "paramshistory") {
    lapply(object[[i]],
      paramshistory,
      slotName = slotName
    )
  }
)

#' @exportMethod paramshistory
#' @rdname QFeatures-accessors
setMethod(
  "paramshistory", "SummarizedExperiment",
  function(object, slotName = "paramshistory") {
    .GetMetadataSlot(object, slotName)
  }
)


#' @export
#' @rdname QFeatures-accessors
"paramshistory<-" <- function(object, i, value) {
  if (inherits(object, "SummarizedExperiment")) {
    S4Vectors::metadata(object)[["paramshistory"]] <- value
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
      se <- object[[i]]
      S4Vectors::metadata(object[[i]])[["paramshistory"]] <- value
    }
  }
  return(object)
}

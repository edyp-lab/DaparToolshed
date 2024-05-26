
#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
  "params", "QFeatures",
  function(object, i, slotName = "params") {
    lapply(object[[i]],
      params,
      slotName = slotName
    )
  }
)

#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
  "params", "SummarizedExperiment",
  function(object, slotName = "params") {
    .GetMetadataSlot(object, slotName)
  }
)


#' @export
#' @rdname QFeatures-accessors
"params<-" <- function(object, i, value) {
  if (inherits(object, "SummarizedExperiment")) {
    S4Vectors::metadata(object)[["params"]] <- value
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
      S4Vectors::metadata(object[[i]])[["params"]] <- value
    }
  }
  return(object)
}

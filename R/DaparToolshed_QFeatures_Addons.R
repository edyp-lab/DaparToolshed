#' @title Additional accessors for instances of class `Qfeatures`.
#'
#' @description
#'
#' These names are common to all assays contained in the object. This is why
#' they are stored in the global metadata. This function is used whenever it i
#' s necessary to (re)detect MEC and POV (new dataset or when post processing 
#' protein qMetacell after aggregation)
#' 
#' @name QFeatures-accessors
#'
#' @param object An instance of class `SummarizedExperiment` or `QFeatures`.
#'
#' @param i The index or name of the assays to extract the quantitative 
#' metadata from. All must have a rowdata variable named as `slotName`
#'
#' @param slotName xxx
#'
#' @param value xxx
#' 
#' @return NA
#'
#'
#' @details
#'
#' Additional slots for rowdata of a `SummarizedExperiment` object:
#'  - qMetacell: xxx
#'
#' Additional slots for Metadata for a `QFeatures` object:
#'  - xxx: xxxx
#'
#' Additional slots for Metadata for a `SummarizedExperiment` object:
#'  - qMetacell: xxxx
#'  - parentProtId: xxx
#'  - idcol: xxxx
#'  - typeDataset: xxx
#'
#'
#' @section Quantitative metadata:
#'
#' Default slotName is `"qMetacell"`.
#'  The value is an adjacency matrix with row and column names. The
#'     matrix will be coerced to compressed, column-oriented sparse
#'     matrix (class `dgCMatrix`) as defined in the `Matrix` package,
#'     as generaled by the [sparseMatrix()] constructor.
#'
#'
#' @examples
#' data(ft, package='DaparToolshed')
#' design.qf(ft)
NULL

#' @exportMethod qMetacell
#' @rdname QFeatures-accessors
#' @return NA
setMethod(
    "qMetacell", "QFeatures",
    function(object, i) {
        .GetRowdataSlot(object[[i]], "qMetacell")
    }
)




#' @export
#' @rdname QFeatures-accessors
#' @return NA
setMethod(
    "qMetacell", "SummarizedExperiment",
    function(object) {
        .GetRowdataSlot(object, "qMetacell")
    }
)



#' @export
#' @rdname QFeatures-accessors
#' @return NA
"qMetacell<-" <- function(object,
                          i,
                          slotName = "qMetacell",
                          value) {
    if (is.null(colnames(value)) | is.null(rownames(value))) {
        stop("The DataFrame must have row and column names.")
    }
    ## Coerse to a data.frame
    # value <- as(value, "data.frame")
    if (inherits(object, "SummarizedExperiment")) {
        if (!identical(rownames(value), rownames(object))) {
            stop("Row names of the SummarizedExperiment and the DataFrame must match.")
        }
        # if (slotName %in% colnames(rowData(object)))
        #  stop("Found an existing variable ", slotName, ".")
        rowData(object)[[slotName]] <- value
        return(object)
    } else {
    if (inherits(object, "QFeatures"))
      {
      if (length(i) != 1) {
        stop("'i' must be of length one. Repeat the call to add a matrix to multiple assays.")
    }
    if (is.numeric(i) && i > length(object)) {
        stop("Subscript is out of bounds.")
    }
    if (is.character(i) && !(i %in% names(object))) {
        stop("Assay '", i, "' not found.")
    }
    
    qMetacell(object[[i]], slotName) <- value
    }
}
    return(object)
}







#' @exportMethod GetUniqueTags
#' @rdname QFeatures-accessors
#' @return NA
setMethod(
  "GetUniqueTags", "QFeatures",
  function(object, i) {
    GetUniqueTags(object[[i]])
  }
)




#' @export
#' @rdname QFeatures-accessors
#' @return NA
setMethod(
  "GetUniqueTags", "SummarizedExperiment",
  function(object) {
    df <- qMetacell(object)
    tmp <- sapply(colnames(df), function(x) unique(df[,x]))
    ll <- unique(as.vector(tmp))
    return(ll)
  }
)







# 
# 
# #' @exportMethod adjacencyMatrix
# #' @rdname QFeatures-accessors
# #' @return NA
# setMethod(
#   "adjacencyMatrix", "QFeatures",
#   function(object, i, slotName = "adjacencyMatrix") {
#     lapply(
#       object[[i]],
#       function(x) {
#         .GetRowdataSlot(x, slotName = slotName)
#       }
#     )
#   }
# )
# 
# 
# 
# 
# #' @export
# #' @rdname QFeatures-accessors
# #' @return NA
# setMethod(
#   "adjacencyMatrix", "SummarizedExperiment",
#   function(object, slotName = "adjacencyMatrix") {
#     .GetRowdataSlot(object, slotName)
#   }
# )
# 
# 
# 
# #' @export
# #' @rdname QFeatures-accessors
# #' @return NA
# "adjacencyMatrix<-" <- function(object,
#                           i,
#                           slotName = "adjacencyMatrix",
#                           value) {
#   if (is.null(colnames(value)) | is.null(rownames(value))) {
#     stop("The DataFrame must have row and column names.")
#   }
#   ## Coerse to a data.frame
#   # value <- as(value, "data.frame")
#   if (inherits(object, "SummarizedExperiment")) {
#     if (!identical(rownames(value), rownames(object))) {
#       stop("Row names of the SummarizedExperiment and the DataFrame 
#                 must match.")
#     }
#     # if (slotName %in% colnames(rowData(object)))
#     #  stop("Found an existing variable ", slotName, ".")
#     rowData(object)[[slotName]] <- value
#     return(object)
#   }
#   stopifnot(inherits(object, "QFeatures"))
#   if (length(i) != 1) {
#    stop("'i' must be of length one. Repeat the call to add a matrix to 
#             multiple assays.")
#  }
#   if (is.numeric(i) && i > length(object)) {
#    stop("Subscript is out of bounds.")
#   }
#   if (is.character(i) && !(i %in% names(object))) {
#     stop("Assay '", i, "' not found.")
#   }
#   se <- object[[i]]
#   object[[i]] <- adjacencyMatrix(se, slotName) <- value
#   return(object)
# }






#' @importFrom S4Vectors metadata
#' @return NA
#' @rdname QFeatures-accessors
.GetMetadataSlot <- function(object, slotName = NULL) {
    S4Vectors::metadata(object)[[slotName]]
}


#' @return NA
#' @rdname QFeatures-accessors
.GetRowdataSlot <- function(object, slotName = NULL) {
    rowData(object)[[slotName]]
}


#---------------------------------------------------------------------------

#' @exportMethod typeDataset
#' @rdname QFeatures-accessors
setMethod("ConnectedComp", "QFeatures",
          function(object, i, slotName = "ConnectedComp") {
            lapply(object[[i]],
                   .GetMetadataSlot,
                   slotName = slotName
            )
          }
)


#' @export
#' @rdname QFeatures-accessors
setMethod("ConnectedComp", "SummarizedExperiment",
          function(object, slotName = "ConnectedComp") {
            .GetMetadataSlot(object, slotName)
          }
)


#' @export
#' @rdname QFeatures-accessors
"ConnectedComp<-" <- function(object, i, slotName = "ConnectedComp", value) {
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
  
  se <- object[[i]]
  S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
  }
  return(object)
}


#' @exportMethod typeDataset
#' @rdname QFeatures-accessors
setMethod("typeDataset", "QFeatures",
    function(object, i, slotName = "typeDataset") {
        lapply(object[[i]],
            .GetMetadataSlot,
            slotName = slotName
        )
    }
)
#' @export
#' @rdname QFeatures-accessors
setMethod("typeDataset", "SummarizedExperiment",
    function(object, slotName = "typeDataset") {
        .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @rdname QFeatures-accessors
"typeDataset<-" <- function(object, i, slotName = "typeDataset", value) {
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
          se <- object[[i]]
          S4Vectors::metadata(object[[i]])[[slotName]] <- value
      }
}
    return(object)
}



#' @exportMethod idcol
#' @rdname QFeatures-accessors
setMethod(
    "idcol", "QFeatures",
    function(object, i, slotName = "idcol") {
        stopifnot(!is.null(object))
        lapply(object[[i]],
            .GetMetadataSlot,
            slotName = slotName
        )
    }
)
#' @export
#' @rdname QFeatures-accessors
setMethod(
    "idcol", "SummarizedExperiment",
    function(object, slotName = "idcol") {
        .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @rdname QFeatures-accessors
"idcol<-" <- function(object, i, slotName = "idcol", value) {
    if (inherits(object, "SummarizedExperiment")) {
        S4Vectors::metadata(object)[[slotName]] <- value
        return(object)
    } else {
      if (inherits(object, "QFeatures")){
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
    S4Vectors::metadata(object[[i]])[[slotName]] <- value
      }
    }
    return(object)
}




#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
setMethod(
    "parentProtId", "QFeatures",
    function(object, i, slotName = "parentProtId") {
        lapply(object[[i]],
            parentProtId,
            slotName = slotName
        )
    }
)

#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
setMethod(
    "parentProtId", "SummarizedExperiment",
    function(object, slotName = "parentProtId") {
        if (typeDataset(object) == "peptide") {
            .GetMetadataSlot(object, slotName)
        }
    }
)


#' @export
#' @rdname QFeatures-accessors
"parentProtId<-" <- function(object, i, slotName = "parentProtId", value) {
    if (inherits(object, "SummarizedExperiment")) {
        if (typeDataset(object) != "peptide") {
            stop("The dataset must contain peptides.")
        }
        S4Vectors::metadata(object)[[slotName]] <- value
        return(object)
    } else {
    if(inherits(object, "QFeatures")){
    if (typeDataset(object[[i]]) != "peptide") {
        stop("The dataset must contain peptides.")
    }
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
    S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
    }
    return(object)
}





#' @exportMethod analysis
#' @rdname QFeatures-accessors
setMethod(
    "analysis", "QFeatures",
    function(object, i, slotName = "analysis") {
        analysis(object[[i]], slotName)
    }
)

#' @export
#' @rdname QFeatures-accessors
setMethod(
    "analysis", "SummarizedExperiment",
    function(object, slotName = "analysis") {
        .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @rdname QFeatures-accessors
"analysis<-" <- function(object, i, slotName = "analysis", value) {
    if (inherits(object, "SummarizedExperiment")) {
        S4Vectors::metadata(object)[[slotName]] <- value
        return(object)
    } else {
    if (inherits(object, "QFeatures")){
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
    S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
    }
    return(object)
}





#' @exportMethod version
#' @rdname QFeatures-accessors
setMethod(
    "version", "QFeatures",
    function(object, slotName = "version") {
        .GetMetadataSlot(object, slotName = slotName)
    }
)



#' @export
#' @rdname QFeatures-accessors
"version<-" <- function(object, slotName = "version", value) {
    stopifnot(inherits(object, "QFeatures"))
    S4Vectors::metadata(object)[[slotName]] <- value
    return(object)
}






#' @exportMethod design.qf
#' @rdname QFeatures-accessors
setMethod(
    "design.qf", "QFeatures",
    function(object, slotName = "design") {
        colData(object)
    }
)



#' @export
#' @rdname QFeatures-accessors
"design.qf<-" <- function(object, slotName = "design", value) {
    stopifnot(inherits(object, "QFeatures"))
    colData(object)@listData <- value
    return(object)
}

#' @export
#' @rdname QFeatures-accessors
mainAssay <- function(object) {
    object[[length(object)]]
}




#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
  "HypothesisTest", "QFeatures",
  function(object, i, slotName = "HypothesisTest") {
    SummarizedExperiment::rowData(object[[i]])[[slotName]]
  }
)

#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
  "HypothesisTest", "SummarizedExperiment",
  function(object, slotName = "HypothesisTest") {
    SummarizedExperiment::rowData(object)[[slotName]]
  }
)


#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
"HypothesisTest<-" <- function(object, i, slotName = "HypothesisTest", value) {
  if (inherits(object, "SummarizedExperiment")) {
    SummarizedExperiment::rowData(object)[[slotName]] <- value
    return(object)
  } else {
    if (inherits(object, "QFeatures")){
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
      # se <- object[[i]]
      SummarizedExperiment::rowData(object[[i]])[slotName] <- value
    }
  }
  return(object)
}




#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
    "params", "QFeatures",
    function(object, i, slotName = "params") {
        params(object[[i]], slotName)
    }
)

#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
setMethod(
    "params", "SummarizedExperiment",
    function(object, slotName = "params") {
        .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @exportMethod params
#' @rdname QFeatures-accessors
"params<-" <- function(object, i, slotName = "params", value) {
    if (inherits(object, "SummarizedExperiment")) {
        S4Vectors::metadata(object)[[slotName]] <- value
        return(object)
    } else {
    if (inherits(object, "QFeatures")){
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
    # se <- object[[i]]
    S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
    }
    return(object)
}









#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
setMethod(
  "names_metacell", "QFeatures",
  function(object, i, slotName = "names_metacell") {
    lapply(object[[i]], names_metacell, slotName = slotName)
  }
)

#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
setMethod(
  "names_metacell", "SummarizedExperiment",
  function(object, slotName = "names_metacell") {
    .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @rdname QFeatures-accessors
"names_metacell<-" <- function(object, i, slotName = "names_metacell", value) {
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
  se <- object[[i]]
  S4Vectors::metadata(object[[i]])[[slotName]] <- value
    }
  }
  return(object)
}



#' @title xxx
#' @description xxx
#' @param obj xxx
#' @param i xxx
#' @export
#' @return xxx
#' 
NAIsZero <- function(obj, i){
  
  assay(obj[[i]])[which(is.na(assay(obj[[i]])))] <- 0
  
  m <- match.metacell(
    omXplore::get_metacell(obj[[i]]),
    pattern = c("Missing", "Missing POV", "Missing MEC"),
    level = omXplore::get_type(obj[[i]])
  )
  
  qMetacell(obj[[i]])[m] <- "Quantified"
  
  return(obj)
  
}
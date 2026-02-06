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
#' @param i The index or name of the assays to extract the quantitative 
#' metadata from. All must have a rowdata variable named as `slotName`
#' @param slotName A character(0) which is the name of the slot in the metadata
#' @param value The content of the slot in the metadata
#' 
#' @return NA
#'
#'
#' @details
#'
#' Additional slots for Metadata for a `SummarizedExperiment` object:
#'  - qMetacell: A `data.frame()`
#'  - parentProtId: A `character()`
#'  - idcol: A `character()`
#'  - typeDataset: A `character()`
#'
#'
#' @section Quantitative metadata:
#'
#' Default slotName is "qMetacell".
#'  The value is an adjacency matrix with row and column names. The
#'     matrix will be coerced to compressed, column-oriented sparse
#'     matrix (class `dgCMatrix`) as defined in the `Matrix` package,
#'     as generaled by the [sparseMatrix()] constructor.
#'
#'
#' @examples
#' data(subR25pept)
#' design.qf(subR25pept)
NULL



#' @rdname QFeatures-accessors
setGeneric("qMetacell", 
  function(object, ...) standardGeneric("qMetacell"))
setGeneric("qMetacell<-", 
  function(object, ..., value) standardGeneric("qMetacell<-"))

#' @exportMethod qMetacell
#' @rdname QFeatures-accessors
setMethod(
    "qMetacell", "QFeatures",
    function(object, i) {
        .GetRowdataSlot(object[[i]], "qMetacell")
    }
)




#' @export
#' @rdname QFeatures-accessors
setMethod(
    "qMetacell", "SummarizedExperiment",
  #' @param object n instance of class `QFeatures` or `SummarizedExperiment`
    function(object) {
        .GetRowdataSlot(object, "qMetacell")
    }
)



#' @export
#' @rdname QFeatures-accessors
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
      SummarizedExperiment::rowData(object)[[slotName]] <- value
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




#' @rdname QFeatures-accessors
setGeneric("GetUniqueTags", 
  function(object, ...) standardGeneric("GetUniqueTags"))



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
#' @import SummarizedExperiment
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
#' @export
.GetMetadataSlot <- function(object, slotName = NULL) {
    S4Vectors::metadata(object)[[slotName]]
}


#' @return NA
#' @rdname QFeatures-accessors
.GetRowdataSlot <- function(object, slotName = NULL) {
  SummarizedExperiment::rowData(object)[[slotName]]
}


#---------------------------------------------------------------------------
#' @rdname QFeatures-accessors
setGeneric("ConnectedComp", 
  function(object, ...) standardGeneric("ConnectedComp"))
setGeneric("ConnectedComp<-", 
  function(object, ..., value) standardGeneric("ConnectedComp<-"))


#' @exportMethod ConnectedComp
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
#' @import SummarizedExperiment
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


#' @rdname QFeatures-accessors
setGeneric("typeDataset", 
  function(object, ...) standardGeneric("typeDataset"))
setGeneric("typeDataset<-", 
  function(object, ..., value) standardGeneric("typeDataset<-"))


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


#' @rdname QFeatures-accessors
setGeneric("idcol", 
  function(object, ...) standardGeneric("idcol"))
setGeneric("idcol<-", 
  function(object, ..., value) standardGeneric("idcol<-"))


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


#' @rdname QFeatures-accessors
setGeneric("parentProtId", 
  function(object, ...) standardGeneric("parentProtId"))
setGeneric("parentProtId<-", 
  function(object, ..., value) standardGeneric("parentProtId<-"))


#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
setMethod(
    "parentProtId", "QFeatures",
    function(object, i, slotName = "parentProtId") {
        lapply(object[[i]], parentProtId, slotName)
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



#' @rdname QFeatures-accessors
setGeneric("filename", 
  function(object, ...) standardGeneric("filename"))
setGeneric("filename<-", 
  function(object, ..., value) standardGeneric("filename<-"))



#' @exportMethod filename
#' @rdname QFeatures-accessors
setMethod(
  "filename", "QFeatures",
  function(object, slotName = "filename") {
    .GetMetadataSlot(object, slotName)
    }
)


#' @export
#' @rdname QFeatures-accessors
"filename<-" <- function(object, slotName = "filename", value) {
  if (inherits(object, "QFeatures")){
      S4Vectors::metadata(object)[[slotName]] <- value
    }
  return(object)
}






#' @rdname QFeatures-accessors
setGeneric("analysis", 
  function(object, ...) standardGeneric("analysis"))
setGeneric("analysis<-", 
  function(object, ..., value) standardGeneric("analysis<-"))


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




#' @rdname QFeatures-accessors
setGeneric("version", 
  function(object, ...) standardGeneric("version"))
setGeneric("version<-", 
  function(object, ..., value) standardGeneric("version<-"))


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



#' @rdname QFeatures-accessors
setGeneric("design.qf", 
  function(object, ...) standardGeneric("design.qf"))
setGeneric("design.qf<-", 
  function(object, ..., value) standardGeneric("design.qf<-"))

#' @exportMethod design.qf
#' @rdname QFeatures-accessors
setMethod(
    "design.qf", "QFeatures",
    function(object, slotName = "design") {
      SummarizedExperiment::colData(object)
    }
)



#' @export
#' @rdname QFeatures-accessors
"design.qf<-" <- function(object, slotName = "design", value) {
    stopifnot(inherits(object, "QFeatures"))
  SummarizedExperiment::colData(object)@listData <- value
    return(object)
}

#' @export
#' @rdname QFeatures-accessors
mainAssay <- function(object) {
    object[[length(object)]]
}



#' @rdname QFeatures-accessors
setGeneric("HypothesisTest", 
  function(object, ...) standardGeneric("HypothesisTest"))
setGeneric("HypothesisTest<-", 
  function(object, ..., value) standardGeneric("HypothesisTest<-"))

#' @export
#' @exportMethod HypothesisTest
#' @rdname QFeatures-accessors
setMethod(
  "HypothesisTest", "QFeatures",
  function(object, i, slotName = "HypothesisTest") {
    SummarizedExperiment::rowData(object[[i]])[[slotName]]
  }
)

#' @export
#' @exportMethod HypothesisTest
#' @rdname QFeatures-accessors
setMethod(
  "HypothesisTest", "SummarizedExperiment",
  function(object, slotName = "HypothesisTest") {
    SummarizedExperiment::rowData(object)[[slotName]]
  }
)


#' @export
#' @exportMethod HypothesisTest
#' @rdname QFeatures-accessors
#' @import SummarizedExperiment
#' 
"HypothesisTest<-" <- function(object, i, value) {
  if (inherits(object, "SummarizedExperiment")) {
    SummarizedExperiment::rowData(object)[["HypothesisTest"]] <- value
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
      SummarizedExperiment::rowData(object[[i]])["HypothesisTest"] <- value
    }
  }
  return(object)
}





#' @rdname QFeatures-accessors
setGeneric("DifferentialAnalysis", 
  function(object, ...) standardGeneric("DifferentialAnalysis"))
setGeneric("DifferentialAnalysis<-", 
  function(object, ..., value) standardGeneric("DifferentialAnalysis<-"))




#' @export
#' @exportMethod DifferentialAnalysis
#' @rdname QFeatures-accessors
setMethod(
  "DifferentialAnalysis", "QFeatures",
  function(object, i, slotName = "DifferentialAnalysis") {
    SummarizedExperiment::rowData(object[[i]])[[slotName]]
  }
)

#' @export
#' @exportMethod DifferentialAnalysis
#' @rdname QFeatures-accessors
setMethod(
  "DifferentialAnalysis", "SummarizedExperiment",
  function(object, slotName = "DifferentialAnalysis") {
    SummarizedExperiment::rowData(object)[[slotName]]
  }
)


#' @export
#' @exportMethod DifferentialAnalysis
#' @rdname QFeatures-accessors
#' @import SummarizedExperiment
#' 
"DifferentialAnalysis<-" <- function(object, i, value) {
  if (inherits(object, "SummarizedExperiment")) {
    SummarizedExperiment::rowData(object)[["DifferentialAnalysis"]] <- value
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
      SummarizedExperiment::rowData(object[[i]])["DifferentialAnalysis"] <- value
    }
  }
  return(object)
}




#' @rdname QFeatures-accessors
setGeneric("names_metacell", 
  function(object, ...) standardGeneric("names_metacell"))

#' @exportMethod names_metacell
#' @rdname QFeatures-accessors
setMethod(
  "names_metacell", "QFeatures",
  function(object, i, slotName = "names_metacell") {
    lapply(object[[i]], names_metacell, slotName = slotName)
  }
)

#' @exportMethod parentProtId
#' @rdname QFeatures-accessors
#' @import SummarizedExperiment
#' 
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



#' @title Set NA values to 0
#' @param obj An instance of QFeatures class
#' @param i An integer which is the index of the assay in the QFeatures object
#' @export
#' @return An instance of QFeatures class
#' 
#' @import SummarizedExperiment
#' 
NAIsZero <- function(obj, i){
  stopifnot(inherits(obj, 'QFeatures'))
  assay(obj[[i]])[which(is.na(assay(obj[[i]])))] <- 0
  
  m <- match.metacell(
    qMetacell(obj[[i]]),
    pattern = c("Missing", "Missing POV", "Missing MEC"),
    level = typeDataset(obj[[i]])
  )
  
  qMetacell(obj[[i]])[m] <- "Quantified"
  
  return(obj)
  
}
#' @title Creates an object of class `QFeatures` from text file.
#'
#' @description
#'
#' Creates an object of class `QFeatures` from a
#' single tabulated-like file for quantitative and meta-data and a dataframe
#' for the samples description.
#'
#' @param data The name of a tab-separated file that contains the data.
#'
#' @param file A `character(1)`. The name of a file xxx
#'
#' @param sample A dataframe describing the samples (in lines).
#'
#' @param indQData A vector of string where each element is the name
#' of a column in designTable that have to be integrated in
#' the `rowData()` table of the `QFeatures` object.
#'
#' @param keyId A `character(1)` or `numeric(1)` which is the indice of the 
#' column containing the ID of entities (peptides or proteins)
#'
#' @param indexForMetacell xxxxxxxxxxx
#'
#' @param force.na A `boolean` that indicates if the '0' and 'NaN' values of
#' quantitative values  must be replaced by 'NA' (Default is FALSE)
#'
#' @param typeDataset A string that indicates whether the dataset is about
#'
#' @param parentProtId A `character(1)` For peptide entities, a string which 
#' is the name of a column in rowData. It contains the id of parent proteins 
#' and is used to generate adjacency matrix and process to aggregation.
#'
#' @param analysis A `character(1)` which is the name of the MS study.
#'
#' @param processes A vector of A `character()` which contains the name of 
#' processes which has already been run on the data. Default is 'original'.
#'
#' @param typePipeline A `character(1)` The type of pipeline used with this 
#' dataset. The list of predefined pipelines in DaparToolshed can be obtained 
#' with the function `pipelines()`. Default value is NULL
#'
#' @param software A `character(1)`
#'
#' @param name A `character(1)` which is the name of the assay in the 
#' QFeatures object. Default is 'original'
#'
#' @return An instance of class `QFeatures`.
#'
#' @author Samuel Wieczorek
#'
#' @examples inst/extdata/examples/ex_createQFeatures.R
#'
#' @importFrom QFeatures readQFeatures
#' @importFrom utils installed.packages read.table
#'
#' @export
#'
#' @rdname import-export-QFeatures
#'
createQFeatures <- function(data = NULL,
  file = NULL,
  sample,
  indQData,
  keyId = "AutoID",
  indexForMetacell = NULL,
  logData = FALSE,
  force.na = TRUE,
  typeDataset,
  parentProtId = NULL,
  analysis = "foo",
  processes = NULL,
  typePipeline = NULL,
  software = NULL,
  name = "original") {


  pkgs.require('QFeatures')
  
    # Check parameters validity
    if (missing(data) && missing(file)) {
        stop("Either 'data' or 'file' is required")
    } else if (!missing(data) && !missing(file)) {
        stop("Only 'data' or 'file' is required at a time. Please choose 
            one of them.")
    }

    if (!missing(data) && !is(data, "data.frame")) {
        stop("'data' must be a data.frame")
    }

    if (!missing(file)) {
        if (!is(file, "xxx")) {
            stop("'file' must be a connection")
        } else {
            data <- read.table(file,
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE
            )
        }
    }



    if (missing(sample)) {
        stop("'sample' is required")
    } else if (!is(sample, "data.frame")) {
        stop("'sample' must be a data.frame")
    }


    if (missing(indQData)) {
        stop("'indQData' is required")
    }
  # else if (!is.numeric(indQData)) {
  #       stop("'indQData' must be a vector of integer")
  #   }

    if (missing(indexForMetacell)) {
        stop("'indexForMetacell' is required")
    }
  # else if (!is.numeric(indexForMetacell)) {
  #       stop("'indexForMetacell' must be a vector of integer")
  #   }

    if (!is.null(keyId) && !is.character(keyId)) {
        stop("'keyId' must be either NULL nor a string")
    }

    if (missing(typeDataset)) {
        stop("'typeDataset' is required")
    }

    # Standardize all colnames
    colnames(data) <- ReplaceSpecialChars(colnames(data))
    indexForMetacell  <- ReplaceSpecialChars(indexForMetacell)
    
    
    if (is.numeric(indQData))
        indQData <- colnames(data)[indQData]
    
    indQData <- ReplaceSpecialChars(indQData)
    
    if (is.numeric(indexForMetacell))
      indexForMetacell <- colnames(data)[indexForMetacell]


    # Standardizes names
    keyId <- ReplaceSpecialChars(keyId)
    typeDataset <- ReplaceSpecialChars(typeDataset)
    parentProtId <- ReplaceSpecialChars(parentProtId)
    analysis <- ReplaceSpecialChars(analysis)
    #processes <- ReplaceSpecialChars(processes)
    typePipeline <- ReplaceSpecialChars(typePipeline)
    software <- ReplaceSpecialChars(software)
    indexForMetacell  <- ReplaceSpecialChars(indexForMetacell)
    



    if (keyId == "AutoID") {
        auto <- rep(paste(typeDataset, "_", seq_len(nrow(data)), sep = ""))
        data <- cbind(data, AutoID = auto)
        rownames(data) <- auto
    } else {
        rownames(data) <- data[, keyId]
    }
    

    # Creates the QFeatures object
    obj <- QFeatures::readQFeatures(
      table = data,
      quantCols = indQData,
      ecol = indQData,
      name = "original",
      fnames = keyId
      )

    ## Encoding the sample data
    sample <- lapply(sample, function(x) {ReplaceSpecialChars(x)})
    design.qf(obj) <- sample

  
    # Get the metacell info
    tmp.qMetacell <- NULL
    if (!is.null(indexForMetacell)) {
      tmp.qMetacell <- data[, indexForMetacell]
      #tmp.qMetacell <- apply(tmp.qMetacell, 2, tolower)
      #tmp.qMetacell <- apply(tmp.qMetacell, 2, function(x) gsub("\\s", "", x))
      tmp.qMetacell <- as.data.frame(tmp.qMetacell, stringsAsFactors = FALSE)
      colnames(tmp.qMetacell) <- gsub(".", "_", colnames(tmp.qMetacell), fixed = TRUE)
      
      qMetacell <- BuildMetacell(
        from = software,
        level = typeDataset,
        qdata = assay(obj),
        conds = colData(obj)$Condition,
        df = tmp.qMetacell
        )
    
      colnames(qMetacell) <- gsub(".", "_", colnames(qMetacell), fixed = TRUE)
    
      # Add the quantitative cell metadata info
      qMetacell(obj[["original"]]) <- qMetacell
    
      # Remove the identification columns which became useless
      .ind <- -match(indexForMetacell, colnames(rowData(obj[[1]])))
      rowData(obj[[1]]) <- rowData(obj[[1]])[, .ind]
      }

    
    
    # Enrich the metadata for whole QFeatures object
    S4Vectors::metadata(obj)$versions <- ProstarVersions()
    S4Vectors::metadata(obj)$analysis <- list(
        analysis = analysis
        #typePipeline = typePipeline
        #processes = c("original", processes)
    )

    # Fill the metadata for the first assay
    typeDataset(obj[["original"]]) <- typeDataset
    idcol(obj[["original"]]) <- keyId

    if (tolower(typeDataset) == "peptide") {
      pkgs.require('PSMatch')
      parentProtId(obj[["original"]]) <- parentProtId
      
      # Create the adjacency matrix
      #X <- PSMatch::makeAdjacencyMatrix(rowData(obj[[1]])[, parentProtId])
      #rownames(X) <- rownames(rowData(obj[[1]]))
      #adjacencyMatrix(obj[[1]]) <- X
      
      # Create the connected components
      #ConnectedComp(obj[[1]]) <- PSMatch::ConnectedComponents(X)
    }

    
    if (force.na) {
      obj <- QFeatures::zeroIsNA(obj, seq_along(obj))
    }
    
    if (logData) {
      obj <- QFeatures::logTransform(obj, seq_along(obj))
    }
    
    return(obj)
}

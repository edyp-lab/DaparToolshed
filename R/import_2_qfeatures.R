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
#' @param indexForMetacell They must be in the same order as the samples in
#' the experimental design
#' @param logData xxx
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
#' @param description A text which describes the study.
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
#' @author Samuel Wieczorek, Manon Gaudin
#'
#' @examples inst/extdata/examples/ex_createQFeatures.R
#'
#' @import QFeatures
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
  description = NULL,
  processes = NULL,
  typePipeline = NULL,
  software = NULL,
  name = "original") {


    # Check parameters validity
    if (missing(data)) {
        stop("'data' is required")
    } else if (!is(data, "data.frame")) {
        stop("'data' must be a data.frame")
      }
    

  
  # Standardize all colnames
  colnames(data) <- ReplaceSpecialChars(colnames(data))
  
  
  
    if (missing(sample)) {
        stop("'sample' is required")
    } else if (!is(sample, "data.frame")) {
        stop("'sample' must be a data.frame")
    } else{
      sample <- cbind (sample, quantCols = names(data)[indQData])
      
    }


    if (missing(indQData)) {
        stop("'indQData' is required")
    }
  
  
  
  
  
  
  # else if (!is.numeric(indQData)) {
  #       stop("'indQData' must be a vector of integer")
  #   }

    # if (missing(indexForMetacell)) {
    #     stop("'indexForMetacell' is required")
    # }
  # else if (!is.numeric(indexForMetacell)) {
  #       stop("'indexForMetacell' must be a vector of integer")
  #   }

    if (!is.null(keyId) && !is.character(keyId)) {
        stop("'keyId' must be either NULL nor a string")
    }

    if (missing(typeDataset)) {
        stop("'typeDataset' is required")
    }

    
   
    
    #if (is.numeric(indQData))
     #   indQData <- colnames(data)[indQData]
    
    #indQData <- ReplaceSpecialChars(indQData)

 
    qdata <- data[,indQData]
    tmp.qMetacell <- NULL
    #ind2delete <- colnames(data)[indQData]
    if(!is.null(indexForMetacell)){
      if (is.numeric(indexForMetacell))
        indexForMetacell <- colnames(data)[indexForMetacell]
      else
        indexForMetacell  <- ReplaceSpecialChars(indexForMetacell)
      
      tmp.qMetacell <- data[, indexForMetacell]
      tmp.qMetacell <- as.data.frame(tmp.qMetacell, stringsAsFactors = FALSE)
      colnames(tmp.qMetacell) <- gsub(".", "_", colnames(tmp.qMetacell), fixed = TRUE)
      
      #ind2delete <- c(ind2delete, indexForMetacell)
      ind2delete <- indexForMetacell
    } 
    #browser()
    ind2delete <- match(ind2delete, colnames(data))
    data <- data[, -ind2delete]
    #data <- cbind(qdata, data)
    #indQData <- 1:length(indQData)
    #indexForMetacell <- length(indQData) + 1:length(indexForMetacell)
    # Applying new order for qdata columns
    

    # Standardizes names
    keyId <- ReplaceSpecialChars(keyId)
    typeDataset <- ReplaceSpecialChars(typeDataset)
    parentProtId <- ReplaceSpecialChars(parentProtId)
    analysis <- ReplaceSpecialChars(analysis)
    #processes <- ReplaceSpecialChars(processes)
    typePipeline <- ReplaceSpecialChars(typePipeline)
    software <- ReplaceSpecialChars(software)
    



    if (keyId == "AutoID") {
        auto <- rep(paste(typeDataset, "_", seq_len(nrow(data)), sep = ""))
        data <- cbind(data, AutoID = auto)
        rownames(data) <- auto
    } else {
        rownames(data) <- data[, keyId]
    }

    #
    # Creates the QFeatures object
    obj <- QFeatures::readQFeatures(
      assayData = data,
      #quantCols = indQData,
      colData = sample,
      name = "original",
      fnames = keyId
      )

    #browser()
    # The function readQFeatures changes non alphanumeric characters, like '+' converted to '.'
    # Force the use of original colnames
    #colnames(rowData(obj[[1]])) <- colnames(data)[-match(indQData, colnames(data))]

    ## Encoding the sample data
    
    design.qf(obj) <- lapply(sample, function(x) {ReplaceSpecialChars(x)})

 
    # Get the metacell info

      qMetacell <- BuildMetacell(
        from = software,
        level = typeDataset,
        qdata = assay(obj),
        conds = colData(obj)$Condition,
        df = tmp.qMetacell
        )
    
      colnames(qMetacell) <- gsub(".", "_", colnames(qMetacell), fixed = TRUE)
      rownames(qMetacell) <- rownames(data)
      

      # Add the quantitative cell metadata info
      qMetacell(obj[["original"]]) <- qMetacell
    
      # if (!is.null(indexForMetacell)) {
      # # Remove the identification columns which became useless
      # .ind <- -match(indexForMetacell, colnames(rowData(obj[[1]])))
      # rowData(obj[[1]]) <- rowData(obj[[1]])[, .ind]
      # }
      
    
    # Enrich the metadata for whole QFeatures object
    S4Vectors::metadata(obj)$versions <- ProstarVersions()
    S4Vectors::metadata(obj)$file <- file
    S4Vectors::metadata(obj)$analysis <- list(
        analysis = analysis,
        description = description
        #typePipeline = typePipeline
        #processes = c("original", processes)
    )

    # Fill the metadata for the first assay
    typeDataset(obj[["original"]]) <- typeDataset
    idcol(obj[["original"]]) <- keyId

    if (tolower(typeDataset) == "peptide" && !is.null(parentProtId)) {
      pkgs.require('PSMatch')
      parentProtId(obj[["original"]]) <- parentProtId
      
      # Create the adjacency matrix
      X <- PSMatch::makeAdjacencyMatrix(rowData(obj[[1]])[, parentProtId])
      rownames(X) <- rownames(rowData(obj[[1]]))
      adjacencyMatrix(obj[[1]]) <- X
      
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

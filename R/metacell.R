
#' @title Quantitative metadata vocabulary for entities
#' 
#' @description
#' This function gives the vocabulary used for the quantitative metadata of each entity in
#' each condition.
#' 
#' @section Glossary:
#' 
#' Peptide-level vocabulary
#'
#' |-- 'Any'
#' |    |
#' |    |-- 1.0 'Quantified'
#' |    |    |
#' |    |    |-- 1.1 "Quant. by direct id" (color 4, white)
#' |    |    |
#' |    |    |-- 1.2 "Quant. by recovery" (color 3, lightgrey)
#' |    |
#' |    |-- 2.0 "Missing" (no color)
#' |    |    |
#' |    |    |-- 2.1 "Missing POV" (color 1)
#' |    |    |
#' |    |    |-- 2.2 'Missing MEC' (color 2)
#' |    |
#' |    |-- 3.0 'Imputed'
#' |    |    |
#' |    |    |-- 3.1 'Imputed POV' (color 1)
#' |    |    |
#' |    |    |-- 3.2 'Imputed MEC' (color 2)
#'
#'
#'
#' Protein-level vocabulary:
#' |-- 'Any'
#' |    |
#' |    |-- 1.0 'Quantified'
#' |    |    |
#' |    |    |-- 1.1 "Quant. by direct id" (color 4, white)
#' |    |    |
#' |    |    |-- 1.2 "Quant. by recovery" (color 3, lightgrey)
#' |    |
#' |    |-- 2.0 "Missing"
#' |    |    |
#' |    |    |-- 2.1 "Missing POV" (color 1)
#' |    |    |
#' |    |    |-- 2.2 'Missing MEC' (color 2)
#' |    |
#' |    |-- 3.0 'Imputed'
#' |    |    |
#' |    |    |-- 3.1 'Imputed POV' (color 1)
#' |    |    |
#' |    |    |-- 3.2 'Imputed MEC' (color 2)
#' |    |
#' |    |-- 4.0 'Combined tags' (color 3bis, lightgrey)
#'
#' 
#' @section Conversion to the glossary:
#' 
#' A generic conversion
#' 
#' Conversion for Proline datasets
#' 
#' Conversion from Maxquant datasets
#' 
#' 
#' @name q_metacell
#' 
#' @examples 
#' 
#' metacell.def('protein')
#' metacell.def('peptide')
#' 
NULL



#' @param level A string designing the type of entity/pipeline. 
#' Available values are: `peptide`, `protein`
#' 
#' @return A data.frame containing the different tags and corresponding colors for the level
#' given in parameter
#' 
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @export
#'
#' @rdname q_metacell
#' 
metacell.def <- function(level){
  if(missing(level))
    stop("'level' is required.")
  
  def <- switch(level,
                peptide = {
                  node <- c(
                    "Any",
                    "Quantified",
                    "Quant. by direct id",
                    "Quant. by recovery",
                    "Missing",
                    "Missing POV",
                    "Missing MEC",
                    "Imputed",
                    "Imputed POV",
                    "Imputed MEC"
                  )
                  parent <- c(
                    "",
                    "Any",
                    "Quantified",
                    "Quantified",
                    "Any",
                    "Missing",
                    "Missing",
                    "Any",
                    "Imputed",
                    "Imputed"
                  )
                  data.frame(
                    node = node,
                    parent = parent
                  )
                },
                protein = {
                  node <- c(
                    "Any",
                    "Quantified",
                    "Quant. by direct id",
                    "Quant. by recovery",
                    "Missing",
                    "Missing POV",
                    "Missing MEC",
                    "Imputed",
                    "Imputed POV",
                    "Imputed MEC",
                    "Combined tags"
                  )
                  parent <- c(
                    "",
                    "Any",
                    "Quantified",
                    "Quantified",
                    "Any",
                    "Missing",
                    "Missing",
                    "Any",
                    "Imputed",
                    "Imputed",
                    "Any"
                  )
                  
                  data.frame(
                    node = node,
                    parent = parent
                  )
                }
  )
  
  
  colors <- custom_metacell_colors()
  
  def <- cbind(def, color = rep("white", nrow(def)))
  
  for (n in seq_len(nrow(def))) {
    def[n, "color"] <- colors[[def[n, "node"]]]
  }
  
  return(def)
}



#' @export
#' @rdname q_metacell
custom_metacell_colors <- function()
  list("Any" = "white",
       "Missing" = "#CF8205",
       "Missing POV" = "#E5A947",
       "Missing MEC" = "#F1CA8A",
       "Quantified" = "#0A31D0",
       "Quant. by recovery" = "#B9C4F2",
       "Quant. by direct id" = "#6178D9",
       "Combined tags" = "#1E8E05",
       "Imputed" = "#A40C0C",
       "Imputed POV" = "#E34343",
       "Imputed MEC" = "#F59898")



#' @title Number of each metacell tags
#' @param obj A instance of the class `SummarizedExperiment`
#' @examples
#' data(subR25prot)
#' GetNbTags(subR25prot[[1]])
#' 
#' @export
#' @import omXplore
#' 
GetNbTags <- function(obj){
  df <- omXplore::get_metacell(obj)
  level <- omXplore::get_type(obj)
  nodes <- metacell.def(level)$node
  
  nb <- sapply(nodes, function(x) length(which(df == x)))
  return(nb)
}



#' @title Parent name of a node
#' @description xxx
#' @param level xxx
#' @param node xxx
#' 
#' #' @examples 
#' Parent('protein', 'Missing')
#' Parent('protein', 'Missing POV')
#' Parent('protein', c('Missing POV', 'Missing MEC'))
#' Parent('protein', c('Missing', 'Missing POV', 'Missing MEC'))
#' 
#' 
#' @export
Parent <- function(level, node=NULL){
  parents <- NULL
  tags <- metacell.def(level)
  
  if (!is.null(node) && length(node) > 0){
    for (p in node){
      ind <- match(p, tags$node)
      if (length(ind) > 0)
        parents <- unique(c(parents, tags$parent[ind]))
    }
  }
  
  
  return(parents)
}

#' @title Names of all chidren of a node
#' @description xxx
#' @param level xxx
#' @param parent xxx
#' 
#' @examples 
#' Children('protein', 'Missing')
#' Children('protein', 'Missing POV')
#' Children('protein', c('Missing POV', 'Missing MEC'))
#' Children('protein', c('Missing', 'Missing POV', 'Missing MEC'))
#' @export
#' 
#' @return A vector
#' 
#' 
Children <- function(level, parent = NULL){
  childrens <- NULL
  tags <- metacell.def(level)
  if (!is.null(parent) && length(parent) > 0){
    for (p in parent){
      ind <- grepl(p, tags$parent)
      if (length(ind) > 0)
        childrens <- unique(c(childrens, tags$node[ind]))
    }
  }
  return(childrens)
}




#' @title List of metacell tags
#'
#' @description
#' This function gives the list of metacell tags available.
#' 
#' - onlyPresent: In this case, the function gives the tags found in a dataset.
#' In addition, and w.r.t to the hierarchy of tags, if all leaves of a node are
#' present, then the tag corresponding to this node is added.
#'
#' @param object An object of class `QFeatures`
#' @param i xxx
#' @param ... xxx
#'
#' @return A vector of tags..
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' GetMetacellTags(subR25pept, 1, level="peptide")
#' GetMetacellTags(subR25pept, 1, level="peptide", onlyPresent=TRUE)
#'
#' @export
#'
#'

#' @rdname QFeatures-accessors
setGeneric("GetMetacellTags", 
  function(object, ...) standardGeneric("GetMetacellTags"))



#' @exportMethod GetMetacellTags
#' @rdname QFeatures-accessors
#' @return NA
setMethod("GetMetacellTags", "QFeatures",
  function(object, i, ...) {
    GetMetacellTags(object[[i]], ...)
  }
)


#' @export
#' @rdname QFeatures-accessors
#' @return NA
setMethod("GetMetacellTags", "SummarizedExperiment",
  function(object, ...) {
    .GetMetacellTags(qMetacell(object), ...)
    }
)


#' @rdname QFeatures-accessors
setMethod("GetMetacellTags", "data.frame",
          function(object, ...) {
            .GetMetacellTags(object, ...)
          }
)


.GetMetacellTags <- function(object = NULL,
                            level = NULL,
                            onlyPresent = FALSE) {
  
  if (is.null(object))
    stop('object in NULL')
  
  if (is.null(level))
    stop('level in NULL')
  
  if (is.null(level) && all)
    stop("level must be defined is 'onlyPresent' equals to FALSE")

    
   ll <- NULL
  if(onlyPresent) {
    # Compute unique tags
    tmp <- sapply(colnames(object), function(x) unique(object[,x]))
    ll <- unique(as.vector(tmp))
    
    # Check if parent must be added
    test <- match (Children(level, 'Any'), ll)
    if (length(test) == length(Children(level, 'Any')) && !all(is.na(test)))
      ll <- c(ll, 'Any')
    
    test <- match (Children(level, 'Quantified'), ll)
    if (length(test) == length(Children(level, 'Quantified')) && !all(is.na(test)))
      ll <- c(ll, 'Quantified')
    
    test <- match (Children(level, 'Missing'), ll)
    if (length(test) == length(Children(level, 'Missing')) && !all(is.na(test)))
      ll <- c(ll, 'Missing')
    
    test <- match (Children(level, 'Imputed'), ll)
    if (length(test) == length(Children(level, 'Imputed')) && !all(is.na(test)))
      ll <- c(ll, 'Imputed')
    
    test <- match (Children(level, 'Combined tags'), ll)
    if (length(test) == length(Children(level, 'Combined tags')) && !all(is.na(test)))
      ll <- c(ll, 'Combined tags')

  } else {
    ll <- metacell.def(level)$node[-which(metacell.def(level)$node =='Any')]
  }
  
  return(ll)
  
}



#' @title Sets the MEC tag in the qMetacell
#' 
#' @description 
#' 
#' This function is based on the qMetacell dataframe to look for either missing
#' values (used to update an initial dataset) or imputed values (used when
#' post processing protein qMetacell after aggregation)
#' 
#' @param obj xxx
#' @param conds A 1-col dataframe with the condition associated to each sample.
#' @param df An object of class `SummarizedExperiment`
#' @param level Type of entity/pipeline
#' 
#' @return An instance of class `QFeatures`.
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' library(QFeatures)
#' data(subR25prot)
#' conds <- design.qf(subR25prot)$Condition
#' df <- Set_POV_MEC_tags(subR25prot[[2]], conds)
#'
#' @name q_metacell
#' 
NULL



#' @export
#' @rdname q_metacell
#' @import SummarizedExperiment
Set_POV_MEC_tags <- function(obj, conds){
  stopifnot(inherits(obj, "SummarizedExperiment"))
  u_conds <- unique(conds)
  
  df <- assay(obj)
  qMeta <- qMetacell(obj)
  level <- typeDataset(obj)
    
    
  for (i in seq_len(length(u_conds))) {
    ind.samples <- which(conds == u_conds[i])
    
    ind.imputed <- match.metacell(qMeta[, ind.samples], "Imputed", level)
    
    ind.missing <- match.metacell(qMeta[, ind.samples], "Missing", level)
    
    ind.missing.pov <- ind.missing & 
      rowSums(ind.missing) < length(ind.samples) & 
      rowSums(ind.missing) > 0
    
    ind.missing.mec <- ind.missing & 
      rowSums(ind.missing) == length(ind.samples)
    
    ind.imputed.pov <- ind.imputed & 
      rowSums(ind.imputed) < length(ind.samples) & 
      rowSums(ind.imputed) > 0
    
    ind.imputed.mec <- ind.imputed & 
      rowSums(ind.imputed) == length(ind.samples)
    
    qMeta[, ind.samples][ind.imputed.mec] <- "Imputed MEC"
    qMeta[, ind.samples][ind.missing.mec] <- "Missing MEC"
    qMeta[, ind.samples][ind.imputed.pov] <- "Imputed POV"
    qMeta[, ind.samples][ind.missing.pov] <- "Missing POV"
  }
  return(qMeta)
}


#' @export
#' @rdname q_metacell
Set_POV_MEC_tags2 <- function(conds, df, level){
  u_conds <- unique(conds)
  
  for (i in seq_len(length(u_conds))) {
    ind.samples <- which(conds == u_conds[i])
    
    ind.imputed <- match.metacell(df[, ind.samples], "Imputed", level)
    
    ind.missing <- match.metacell(df[, ind.samples], "Missing", level)
    
    ind.missing.pov <- ind.missing & 
      rowSums(ind.missing) < length(ind.samples) & 
      rowSums(ind.missing) > 0
    
    ind.missing.mec <- ind.missing & 
      rowSums(ind.missing) == length(ind.samples)
    
    ind.imputed.pov <- ind.imputed & 
      rowSums(ind.imputed) < length(ind.samples) & 
      rowSums(ind.imputed) > 0
    
    ind.imputed.mec <- ind.imputed & 
      rowSums(ind.imputed) == length(ind.samples)
    
    df[, ind.samples][ind.imputed.mec] <- "Imputed MEC"
    df[, ind.samples][ind.missing.mec] <- "Missing MEC"
    df[, ind.samples][ind.imputed.pov] <- "Imputed POV"
    df[, ind.samples][ind.missing.pov] <- "Missing POV"
  }
  return(df)
}

#' @title The set of available softwares to convert from
#' 
#' @examples 
#' GetSoftAvailables()
#' @export

GetSoftAvailables <- function(){
  
  
  # funcs <- ls('package:DaparToolshed')
  # funcs <- funcs[grep('Metacell_', funcs)]
  # funcs <- strsplit(funcs, 'Metacell_')
  # funcs <- unlist(lapply(funcs, function(x) x[[2]]))
  # funcs <- funcs[-which(funcs=='generic')]
  
  funcs <- c("DIA_NN", "maxquant", "proline")
  
  return(funcs)
}




#' @title metacell function which xxx
#' 
#' @description xxx
#' 
#' @param from xxx
#' @param level xxx
#' @param qdata A matrix of quantitative data
#' @param conds xxx
#' @param df A data.frame which contains the type of identification of the
#' entities. It must have the same dimensions as `qData`.
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
#' @example inst/extdata/examples/ex_BuildMetacell.R
#' 
#' @export
#' 
#' @rdname BuildMetacell 
#' 
#'  
BuildMetacell <- function(from = NULL, 
                          level,
                          qdata = NULL,
                          conds = NULL, 
                          df = NULL){
  if (missing(from)) {
    stop("'from' is required.")
  }
  if (!(from %in% GetSoftAvailables()))
    stop("'from' must be one of the following")
  if (missing(level)) {
    stop("'level' is required.")
  }
  if (is.null(qdata)) {
    stop("'qdata' is required.")
  }
  if (is.null(conds)) {
    stop("'conds' is required.")
  }
  
  if(!is.null(df))
    if (nrow(df) != nrow(qdata) || ncol(df) != ncol(qdata))
     stop("'df' and 'qdata' must be dataframes of the same dimensions.")
  
  if (is.null(df)) {
    df <- Metacell_generic(qdata, conds, level)
  } else {
    switch(from,
           maxquant = df <- Metacell_maxquant(qdata, conds, df, level),
           proline = df <- Metacell_proline(qdata, conds, df, level),
           'DIA-NN' = df <- Metacell_proline(qdata, conds, df, level)
    )
  }
  
  return(df)
}





#' @title Sets the qMetacell dataframe
#' 
#' @description
#' 
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0. 
#' Conversion rules
#' Quanti			Tag		
#' NA or 0		NA		
#'
#' 
#' @param qdata A matrix of quantitative data
#' @param conds xxx
#' @param level xxx
#' 
#' @return xxxxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples 
#' library(SummarizedExperiment)
#' data(subR25pept)
#' conds <- design.qf(subR25pept)$Condition
#' qdata <- assay(subR25pept[[2]])
#' df <- Metacell_generic(qdata, conds, 'peptide')
#' 
#' @export
#' 
#' @rdname q_metacell
#' 
#'  
Metacell_generic <- function(qdata, conds, level) {
  if (missing(qdata)) {
    stop("'qdata' is required")
  }
  if (missing(conds)) {
    stop("'conds' is required.")
  }
  if (missing(level)) {
    stop("'level' is required.")
  }
  
  df <- data.frame(
    matrix(rep("Quantified", nrow(qdata) * ncol(qdata)),
           nrow = nrow(qdata),
           ncol = ncol(qdata)
    ),
    stringsAsFactors = FALSE
  )
  
  # Rule 1
  qdata[qdata == 0] <- NA
  df[is.na(qdata)] <- "Missing"
  df <- Set_POV_MEC_tags2(conds, df, level)
  
  colnames(df) <- paste0("metacell_", colnames(qdata))
  colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
  
  return(df)
}



#' @title Sets the metacell dataframe for datasets which are from Dia-NN software
#'
#' @description
#' Actually, this function uses the generic function to generate metacell info
#' 
#' @param qdata An object of class `MsnSet`
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' conds <- design.qf(subR25pept)$Condition
#' qdata <- SummarizedExperiment::assay(subR25pept[[2]])
#' df <- Metacell_DIA_NN(qdata, conds, df, level = "peptide")
#' 
#'
#' @export
#'
#'
Metacell_DIA_NN <- function(qdata, conds, df, level = NULL) {
  if (missing(qdata)) {
    stop("'qdata' is required")
  }
  if (missing(conds)) {
    stop("'conds' is required.")
  }
  if (missing(level)) {
    stop("'level' is required.")
  }
  
  
  return(df)
}

#' @title Sets the metacell dataframe for datasets which are from Proline software
#'
#' @description
#' In the quantitative columns, a missing value is identified by no value rather
#' than a value equal to 0.
#' 
#' In these datasets, the metacell info is computed from the 'PSM count' columns.
#' 
#' Conversion rules
#' Initial conversion rules for proline
#' |--------------|-----------------|-----|
#' | Quanti       |    PSM count    | Tag |
#' |--------------|-----------------|-----|
#' |  == 0 | N.A. |   whatever      | 2.0 |
#' |  > 0         |    > 0          | 1.1 |
#' |  > 0         |    == 0         | 1.2 |
#' |  > 0         |   unknown col   | 1.0 |
#' |--------------|-----------------|-----|
#'
#' @param qdata An object of class `MsnSet`
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' conds <- design.qf(subR25pept)$Condition
#' qdata <- SummarizedExperiment::assay(subR25pept[[2]])
#' df <- Metacell_proline(qdata, conds, level = "peptide")
#' 
#'
#' @export
#'
#'
Metacell_proline <- function(qdata, conds, df = NULL, level = NULL) {
  if (missing(qdata)) {
    stop("'qdata' is required")
  }
  if (missing(conds)) {
    stop("'conds' is required.")
  }
  if (missing(level)) {
    stop("'level' is required.")
  }
  
  
  if (is.null(df)) {
    df <- data.frame(matrix(rep("Quantified", nrow(qdata) * ncol(qdata)),
                            nrow = nrow(qdata),
                            ncol = ncol(qdata)
    ),
    stringsAsFactors = FALSE
    )
  }
  
  # Rule 1
  df[is.na(qdata)] <- "Missing"
  df <- Set_POV_MEC_tags2(conds, df, level)
  
  # Rule 2
  df[df > 0 & qdata > 0] <- "Quant. by direct id"
  
  # Rule 3
  df[df == 0 & qdata > 0] <- "Quant. by recovery"
  
  colnames(df) <- paste0("metacell_", colnames(qdata))
  colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
  
  return(df)
}




#' @title Sets the metacell dataframe
#'
#' @description
#' Initial conversion rules for maxquant
#' |------------|-----------------------|--------|
#' | Quanti     |     Identification    |    Tag |
#' |------------|-----------------------|--------|
#' |  == 0      |       whatever        |    2.0 |
#' |  > 0       |       'By MS/MS'      |    1.1 |
#' |  > 0       |      'By matching'    |    1.2 |
#' |  > 0       |       unknown col     |    1.0 |
#' |------------|-----------------------|--------|
#'
#' @param qdata An object of class `MsnSet`
#'
#' @param conds xxx
#'
#' @param df A list of integer xxxxxxx
#'
#' @param level xxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' conds <- design.qf(subR25pept)$Condition
#' qdata <- SummarizedExperiment::assay(subR25pept[[2]])
#' df2 <- Metacell_maxquant(qdata, conds, level = "peptide")
#'
#' @export
#'
#'
Metacell_maxquant <- function(qdata, conds, df = NULL, level = NULL) {
  if (missing(qdata)) {
    stop("'qdata' is required")
  }
  if (missing(conds)) {
    stop("'conds' is required.")
  }
  if (missing(level)) {
    stop("'level' is required.")
  }
  
  
  if (is.null(df)) {
    df <- data.frame(matrix(rep(
      "Quantified",
      nrow(qdata) * ncol(qdata)
    ),
    nrow = nrow(qdata),
    ncol = ncol(qdata)
    ),
    stringsAsFactors = FALSE
    )
  }
  
  
  # Rule 1
  qdata[qdata == 0] <- NA
  
  # Rule 2
  df[df == "byms/ms"] <- "Quant. by direct id"
  
  # Rule 3
  df[df == "bymatching"] <- "Quant. by recovery"
  
  # Add details for NA values
  df[is.na(qdata)] <- "Missing"
  df <- Set_POV_MEC_tags2(conds, df, level)
  
  colnames(df) <- paste0("metacell_", colnames(qdata))
  colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
  
  return(df)
}




#' @title Similar to the function `is.na()` but focused on the equality
#' with the paramter 'type'.
#'
#' @param metadata A data.frame
#'
#' @param pattern The value to search in the dataframe
#'
#' @param level xxx
#'
#' @return A boolean dataframe
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' metadata <- qMetacell(subR25pept[[1]])
#' m <- match.metacell(metadata, pattern = "Missing", level = "peptide")
#' m <- match.metacell(metadata, pattern = 'Missing POV', level = "peptide")
#' m <- match.metacell(metadata, pattern = c('Missing', 'Missing POV'), level = "peptide")
#' @export
#'
match.metacell <- function(metadata, pattern = NULL, level) {
  if (missing(metadata))
    stop("'metadata' is required")
  
  if (missing(pattern) || is.null(pattern))
    stop("'pattern' is required.")

  if (missing(level))
    stop("'level' is required.")
  
  
  #is.subset <- pattern == intersect(pattern,  metacell.def(level)$node)
  if (sum(pattern == intersect(pattern,  metacell.def(level)$node)) !=  length(pattern)) {
    stop(paste0(
      "'pattern' is not correct. Available values are: ",
      paste0(metacell.def(level)$node, collapse = " ")
    ))
  }
  
  ll.res <- lapply(pattern, function(x) {metadata == x})
  
  res <- NULL
  for (i in seq_len(length(ll.res))) {
    if (i == 1) {
      res <- ll.res[[1]]
    } else {
      res <- res | ll.res[[i]]
    }
  }
  
  return(res)
}

#' @rdname q_metacell
setGeneric("UpdateMetacellAfterImputation", 
  function(object, ...) standardGeneric("UpdateMetacellAfterImputation"))




#' @title Update quantitative metadata after imputation
#' 
#' @description
#' Update the quantitative metadata information of missing values that were imputed
#' 
#' @param object xxx
#' @param from xxx
#' @param to xxx
#' @param ... xxx
#' 
#' @examples
#' data(subR25prot)
#' subR25prot[[2]] <- UpdateMetacellAfterImputation(subR25prot[[2]], 'Missing', 'Imputed')
#' 
#' @author Samuel Wieczorek
#'
#' 
#' @return NA
#' 
#' @rdname q_metacell
#' @exportMethod UpdateMetacellAfterImputation
#' @export
#' 
#' @import omXplore
#' @import SummarizedExperiment
#' 
setMethod("UpdateMetacellAfterImputation", "SummarizedExperiment",
          function(object,
                   from,
                   to,
                   ...) {
            
            if (missing(object))
              stop("'object' is required.")
            #level <- typeDataset(object)
            # if (missing(na.type)){
            #   values <- unname(search.qMetacell.tags('Missing', level))
            #   stop("'na.type' is required. Available values are: ", 
            #        paste0(values, collapse=' '))
            # }
            
            
            ind <- match.metacell(
              metadata = omXplore::get_metacell(object), 
              pattern = c('Missing', 'Missing POV', 'Missing MEC'), 
              level = omXplore::get_type(object)) & !is.na(assay(object))
            
            rowData(object)$qMetacell[ind] <- gsub(
              pattern = "Missing", 
              replacement = "Imputed",
              x = omXplore::get_metacell(object)[ind],
              fixed = TRUE)
            object
          })




#' @title
#' Search pattern in qMetacell vocabulary
#' 
#' @description
#' Gives all the tags of the metadata vocabulary containing the pattern 
#' (parent and all its children).
#' 
#' @param pattern The string to search.
#' 
#' @param level The available levels are : names()
#' 
#' @param depth xxx
#' 
#' @author Samuel Wieczorek
#' 
#' @examples
#' search.metacell.tags('Missing POV', 'peptide')
#' search.metacell.tags('Quantified', 'peptide')
#' 
#' @export
#' @return NA
#' 
#' @rdname q_metacell
#' 
#' 
search.metacell.tags <- function(pattern, level, depth = "1") {
  if (missing(pattern)) {
    stop("'pattern' is required.")
  } else if (!(pattern %in% metacell.def(level)$node)) {
    stop(paste0("'pattern' must be one of the following: ", 
                paste0(metacell.def(level)$node, collapse = " ")))
  }
  
  if (missing(level)) {
    stop("'level' is required.")
  }
  if (!(depth %in% c("0", "1", "*"))) {
    stop("'depth' must be one of the following: 0, 1 or *")
  }
  
  .ind <- which(metacell.def(level)$parent == pattern)
  tags <- NULL
  tags <- switch(depth,
                 "0" = pattern,
                 "1" = c(pattern, 
                         metacell.def(level)$node[.ind]),
                 "*" = {
                   if (length(metacell.def(level)$node[.ind]) == 0) {
                     search.metacell.tags(pattern, level, "0")
                   } else {
                     c(pattern, unlist(lapply(
                       metacell.def(level)$node[.ind],
                       function(x) {
                         search.metacell.tags(x, level, depth)
                       }
                     )))
                   }
                 }
  )
  
  return(tags)
}



#' @title Combine peptide metadata to build protein metadata
#' 
#' @description
#' 
#' Agregation rules for the cells quantitative metadata of peptides. 
#' Please refer to the qMetacell.def vocabulary in `qMetacell.def()`
#' 
#' # Basic agreagtion
#' Agregation of non imputed values (2.X) with quantitative values 
#' (1.0, 1.X, 3.0, 3.X)
#' |----------------------------
#' Not possible
#' |----------------------------
#' 
#' Agregation of different types of missing values (among 2.1, 2.2)
#' |----------------------------
#' * Agregation of 2.1 peptides between each other gives a missing value 
#'   non imputed (2.0)
#' * Agreagtion of 2.2 peptides between each other givesa missing value 
#'   non imputed (2.0)
#' * Agregation of a mix of 2.1 and 2.2 gives a missing value non imputed (2.0)
#' |----------------------------
#' 
#' 
#' Agregation of a mix of quantitative values (among 1.0, 1.1, 1.2, 3.0, 3.X)
#' |----------------------------
#' * if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, 
#'   then the final metadata is set the this tag
#' * if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, 
#'   then the final metadata is set to 1.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0, 
#'   then the final metadata is set to 3.0
#' * if the set of metacell to agregate is a mix of 3.X and 3.0 and other (X.X),
#'   then the final metadata is set to 4.0
#' |----------------------------
#' 
#' # Post processing
#' Update metacell with POV/MEC status for the categories 2.0 and 3.0
#' TODO
#' 
#' @param met xxx
#' 
#' @param level xxx
#' 
#' @examples
#' \donttest{
#' ll <- omXplore::metacell.def('peptide')$node
#' for (i in 1:length(ll))
#' test <- lapply(combn(ll, i, simplify = FALSE), 
#' function(x) tag <- metacombine(x, 'peptide'))
#' }
#' 
#' @export
#' @rdname q_metacell
#' @importFrom utils combn
#' 
#' 
metacombine <- function(met, level) {
  tag <- NULL
  if (length(met) == 0) {
    return("Missing")
  }
  
  u_met <- unique(met)
  
  # ComputeNbTags <- function(tag) {
  #     sum(
  #         unlist(
  #             lapply(search.metacell.tags(tag, level),
  #                 function(x) length(grep(x, u_met)))))
  # }
  # 
  # 
  # nb.tags <- lapply(metacell.def(level)$node, 
  #     function(x) as.numeric(x %in% u_met))
  # n.imputed <- ComputeNbTags("Imputed")
  # n.missing <- ComputeNbTags("Missing")
  # n.quanti <- ComputeNbTags("Quantified")
  
  .missing_exists <- sum(c('Missing', 'Missing POV', 'Missing MEC') %in% u_met) > 0
  .quanti_exists <- sum(c('Quantified', 'Quant. by direct id', 'Quant. by recovery') %in% u_met) > 0
  .imputed_exists <- sum(c('Imputed', 'Imputed POV', 'Imputed MEC') %in% u_met) > 0
  
  ###
  ### RULE 1
  ###
  #if (n.missing > 0 && (n.imputed > 0 || n.quanti > 0)) tag <- "STOP"
  if ( .missing_exists && (.quanti_exists ||.imputed_exists ))
    tag <- "STOP"
  
  
  ###
  ### RULE 2
  ###
  if (setequal(u_met,'Missing POV')) tag <- 'Missing'
  
  ###
  ### RULE 3
  ###
  if(setequal(u_met, 'Missing MEC')) tag <- 'Missing'
  
  ###
  ### RULE 4
  ###
  if(setequal(u_met, c('Missing POV', 'Missing MEC'))) tag <- 'Missing'
  
  
  # --- RULE 5 ---
  # if the type of all the peptides to agregate is 1.0, 1.1 or 1.2, then 
  # the final metadata is set the this tag
  if (length(u_met) == 1 && 
      (u_met %in%  c('Quantified', "Quant. by direct id", "Quant. by recovery")))
    tag <- u_met
  
  
  # --- RULE 5 bis---
  # if the type of all the peptides to agregate is 3.0, 3.1 or 3.2, then 
  # the final metadata is set the this tag
  if (length(u_met) == 1 &&
      (u_met %in%  c('Imputed', "Imputed POV", "Imputed MEC")))
    tag <- u_met
  
  
  # --- RULE 6 ---
  # if the set of metacell to agregate is a mix of 1.0, 1.1 or 1.2, then 
  # the final metadata is set to 1.0
  inter <- intersect(u_met, c('Quantified', "Quant. by direct id", "Quant. by recovery"))
  if ((length(inter) >=2 ) && length(u_met)==length(inter) && sum(u_met == inter)== length(u_met) )
    tag <- 'Quantified'
  
  
  # --- RULE 7 ---
  # if the set of metacell to agregate is a mix of 3.0, 3.1 or 3.2, then 
  # the final metadata is set to 3.0
  inter <- intersect(u_met, c('Imputed', "Imputed POV", "Imputed MEC"))
  if ((length(inter) >= 2) &&  length(u_met)==length(inter) && sum(u_met == inter)== length(u_met) )
    tag <- 'Imputed'
  
  
  
  # If the set of metacell to agregate is a mix of 3.X and 3.0 and quantitative values 
  # (1.X), then the final metadata is set to 4.0
  if (length(u_met) > 1 && !.missing_exists && .quanti_exists && .imputed_exists) 
    tag <- "Combined tags"
  
  
  return(tag)
}


#' 
#' #' @title
#' #' Symbolic product of matrices
#' #'
#' #' @description
#' #' Execute a product two matrices: the first is an adjacency one  while the
#' #' second if a simple dataframe
#' #'
#' #' @param X An adjacency matrix between peptides and proteins
#' #'
#' #' @param obj.pep A dataframe of the cell metadata for peptides
#' #'
#' #' @return xxxx
#' #'
#' #' @author Samuel Wieczorek
#' #'
#' #' @examples
#' #' \donttest{
#' #' data(subR25pept)
#' #' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' #' protID <- parentProtId(obj.pep[[length(obj.pep)]])
#' #' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]], protID)
#' #' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' #' agg.meta <- AggregateMetacell(obj.pep[[length(obj.pep)]])
#' #' }
#' #'
#' #' @export
#' #' 
#' #' @import omXplore
#' #' @import stats
#' #'
#' AggregateMetacell <- function(obj.pep) {
#'   
#'   issues <- NULL
#'   meta <- qMetacell(obj.pep)
#'   level <- TypeDataset(obj.pep)
#'   X <- adjacencyMatrix(obj.pep)
#'   rowcol <- function(meta.col, X.col) (meta.col)[X.col > 0]
#'   
#'   df <- data.frame(stringsAsFactors = TRUE)
#'   for (j in seq_len(ncol(meta))) {
#'     for (i in seq_len(ncol(X))) {
#'       df[i, j] <- metacombine(rowcol(meta[, j], X[, i]), level)
#'     }
#'   }
#'   
#'   df[df == "NA"] <- NA
#'   colnames(df) <- obj.pep@experimentData@other$names_metacell
#'   rownames(df) <- colnames(X)
#'   # Delete protein with only NA
#'   
#'   
#'   
#'   # Post processing of metacell to discover 'Imputed POV', 'Imputed MEC'
#'   conds <- design.qf(obj.pep)$Condition
#'   df <- Set_POV_MEC_tags2(conds, df, level)
#'   
#'   
#'   # Search for issues
#'   prot.ind <- unique(rownames(which(df == "STOP", arr.ind = TRUE)))
#'   if (!is.null(prot.ind)) {
#'     issues <- stats::setNames(
#'       lapply(
#'         prot.ind,
#'         function(x) {
#'           rownames(X)[which(X[, which(colnames(X) == x)] == 1)]
#'         }
#'       ),
#'       prot.ind
#'     )
#'   }
#'   
#'   list(
#'     metacell = df,
#'     issues = issues
#'   )
#' }
#' 


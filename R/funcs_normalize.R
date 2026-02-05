
#' @title Normalisation
#' 
#' @description 
#' Provides several methods to normalize quantitative data from
#' a \code{SummarizedExperiment} object.
#' They are organized in six main families : GlobalQuantileAlignement,
#' sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
#' For the first family, there is no type.
#' For the five other families, two type categories are available :
#' "Overall" which means that the value for each protein
#' (ie line in the expression data tab) is computed over all the samples ;
#' "within conditions" which means that the value for each protein
#' (ie line in the \code{SummarizedExperiment::assay()} data tab) is computed condition
#' by condition.
#'
#' @param target Category of normalization method to show. 
#' Either "all", "withTracking" or "withoutTracking".
#'
#' @param qData A data.frame with quantitative data to normalize. 
#'
#' @param conds A `character()` vector which is the names of conditions
#' for each sample in the dataset.
#'
#' @param type "overall" (shift all the sample distributions at once) or
#' "within conditions" (shift the sample distributions within each
#' condition at a time).
#'
#' @param subset.norm A vector of index indicating rows to be used for
#' normalization.
#'
#' @param quantile A float that corresponds to the quantile used to
#' align the data.
#' 
#' @param scaling A boolean that indicates if the variance of the data have to
#' be forced to unit (variance reduction) or not.
#' 
#' @param span A float between 0 and 1 indicating the span of loess smoothing window.
#' 
#' @author Samuel Wieczorek, Thomas Burger, Helene Borges, Anais Courtier,
#' Enora Fremy
#'
#' @examples
#' 
#' ## Get the list of methods
#' normalizeMethods()
#' 
#' 
#' data(subR25pept)
#' qData <- SummarizedExperiment::assay(subR25pept[[1]])
#' conds <- design.qf(subR25pept)$Condition
#' 
#' 
#' 
#' #normalized <- GlobalQuantileAlignment(qData)
#' 
#' normalized <- SumByColumns(qData, conds,
#'     type = "within conditions",
#'     subset.norm = seq_len(10)
#' )
#' 
#' normalized <- QuantileCentering(
#' SummarizedExperiment::assay(subR25pept), conds,
#' type = "within conditions", subset.norm = seq_len(10)
#' )
#' 
#' normalized <- MeanCentering(qData, conds, type = "overall")
#' 
#' # normalized <- vsn(qData, conds, type = "overall")
#' 
#' normalized <- LOESS(qData, conds, type = "overall")
#'
#' @name normalization_methods
#'
NULL


#' @rdname normalization_methods
#' @export
#' @return xxx
#'
normalizeMethods <- function(target = 'all'){
  switch(target,
    all =  c("GlobalQuantileAlignment", "SumByColumns", "QuantileCentering",
    "MeanCentering", "LOESS",  "vsn"),
  withTracking =  c("SumByColumns", "QuantileCentering", "MeanCentering"),
  withoutTracking = c("GlobalQuantileAlignment", "SumByColumns",
      "QuantileCentering", "MeanCentering", "LOESS", "vsn")
  )

}


#' @return A normalized numeric matrix
#' @export
#' @rdname normalization_methods
#' @importFrom preprocessCore normalize.quantiles
#'
GlobalQuantileAlignment <- function(qData) {

  e <- preprocessCore::normalize.quantiles(as.matrix(qData),
    keep.names = TRUE)
  return(e)
}



#' @return A normalized numeric matrix
#' @export
#' @rdname normalization_methods
#' @importFrom stats median
#'
SumByColumns <- function(qData,
  conds = NULL,
  type = NULL,
  subset.norm = NULL) {
  
  
  if (missing(conds)) {
    stop("'conds' is required")
  }
  if (missing(type)) {
    stop("'type' is required")
  }
  
  
  if (!(type %in% c("overall", "within conditions"))) {
    stop("'type' must have one of the following values: 
            'overall', 'within conditions'")
  }
  
  
  qData <- as.matrix(qData)
  
  e <- 2^qData
  
  if (is.null(subset.norm) || length(subset.norm) < 1) {
    subset.norm <- seq_len(nrow(qData))
  }
  
  if (type == "overall") {
    if (length(subset.norm) == 1) {
      sum_cols <- e[subset.norm, ]
    } else {
      sum_cols <- colSums(e[subset.norm, ], na.rm = TRUE)
    }
    
    for (i in seq_len(nrow(e))) {
      e[i, ] <- (e[i, ] / sum_cols) * (stats::median(sum_cols))
    }
  } else if (type == "within conditions") {
    for (l in unique(conds)) {
      indices <- which(conds == l)
      
      if (length(subset.norm) == 1) {
        sum_cols <- e[subset.norm, indices]
      } else {
        sum_cols <- colSums(e[subset.norm, indices], na.rm = TRUE)
      }
      
      for (i in seq_len(nrow(e))) {
        e[i, indices] <- (e[i, indices] / sum_cols) * 
          stats::median(sum_cols)
      }
    }
  }
  e <- log2(e)
  return(e)
}




#' @return A normalized numeric matrix
#' @export
#' @rdname normalization_methods
#' @importFrom stats quantile
#'
QuantileCentering <- function(qData,
  conds = NULL,
  type = "overall",
  subset.norm = NULL,
  quantile = 0.15) {

  if (missing(conds)) {
    stop("'conds' is required")
  }
  if (missing(type)) {
    stop("'type' is required")
  }
  
  
  if (!(type %in% c("overall", "within conditions"))) {
    stop("'type' must have one of the following values: 'overall', 
            'within conditions'")
  }
  
  
  qData <- as.matrix(qData)
  
  if (is.null(subset.norm) || length(subset.norm) < 1) {
    subset.norm <- seq_len(nrow(qData))
  }
  
  q <- function(x) {
    stats::quantile(x, probs = quantile, na.rm = TRUE)
  }
  
  
  if (length(subset.norm) == 1) {
    quantileOverSamples <- qData[subset.norm, ]
  } else {
    quantileOverSamples <- apply(qData[subset.norm, ], 2, q)
  }
  
  
  if (type == "overall") {
    cOverall <- q(quantileOverSamples)
    qData <- sweep(qData, 2, quantileOverSamples)
    qData <- qData + cOverall
  } else if (type == "within conditions") {
    qData <- sweep(qData, 2, quantileOverSamples)
    cCond <- NULL
    for (l in unique(conds)) {
      indices <- which(conds == l)
      cCond[l] <- q(quantileOverSamples[indices])
      qData[, indices] <- qData[, indices] + cCond[l]
    }
  }
  return(qData)
}


#' @return A normalized numeric matrix
#' @export
#' @rdname normalization_methods
#'
MeanCentering <- function(qData,
  conds,
  type = "overall",
  subset.norm = NULL,
  scaling = FALSE) {
  if (missing(conds)) {
    stop("'conds' is required")
  }
  
  qData <- as.matrix(qData)
  
  if (is.null(subset.norm) || length(subset.norm) < 1) {
    subset.norm <- seq_len(nrow(qData))
  }
  
  if (length(subset.norm) == 1) {
    meanOverSamples <- qData[subset.norm, ]
  } else {
    meanOverSamples <- apply(qData[subset.norm, ], 2, mean, na.rm = TRUE)
  }
  
  if (type == "overall") {
    cOverall <- mean(meanOverSamples)
    qData <- sweep(qData, 2, meanOverSamples)
    if (scaling) {
      qData <- scale(qData, center = FALSE, scale = TRUE)
      attr(qData, "scaled:scale") <- NULL
    }
    qData <- qData + cOverall
  } else if (type == "within conditions") {
    .temp <- sweep(qData, 2, meanOverSamples)
    if (scaling) {
      qData <- scale(qData, center = FALSE, scale = TRUE)
      attr(qData, "scaled:scale") <- NULL
    }
    cCond <- NULL
    for (l in unique(conds)) {
      indices <- which(conds == l)
      cCond[l] <- mean(meanOverSamples[indices])
      qData[, indices] <- .temp[, indices] + cCond[l]
    }
  }
  return(qData)
}

#' @return A normalized numeric matrix
#' @export
#' @rdname normalization_methods
#' 
#'
vsn <- function(qData, 
  conds, 
  type = NULL) {
  pkgs.require('vsn')
  
  
  if (missing(conds)) {
    stop("'conds' is required")
  }
  
  if (type == "overall") {
    vsn.fit <- vsn::vsnMatrix(2^(qData))
    qData <- vsn::predict(vsn.fit, 2^(qData))
  } else if (type == "within conditions") {
    for (l in unique(conds)) {
      indices <- which(conds == l)
      vsn.fit <- vsn::vsnMatrix(2^(qData[, indices]))
      qData[, indices] <- vsn::predict(vsn.fit, 2^(qData[, indices]))
    }
  }
  return(qData)
}


#' @return A normalized numeric matrix
#' @export
#' 
#' @rdname normalization_methods
#'
LOESS <- function(qData, 
  conds, 
  type = "overall", 
  span = 0.7) {
  
  pkgs.require('limma')
  
  
  if (missing(conds)) {
    stop("'conds' is required")
  }
  
  if (type == "overall") {
    qData <- limma::normalizeCyclicLoess(
      x = qData, 
      method = "fast", 
      span = span
    )
  } else if (type == "within conditions") {
    for (l in unique(conds)) {
      indices <- which(conds == l)
      qData[, indices] <- limma::normalizeCyclicLoess(
        x = qData[, indices],
        method = "fast",
        span = span
      )
    }
  }
  return(qData)
}




#' @title Normalisation for QFeatures
#' 
#' @description 
#' This method is a wrapper that provides several methods to normalize quantitative 
#' data from objects of class \code{QFeatures} or \code{SummarizedExperiment}.
#' 
#' They are organized in six main families : GlobalQuantileAlignement,
#' sumByColumns, QuantileCentering, MeanCentering, LOESS, vsn
#' For the first family, there is no type.
#' For the five other families, two type categories are available :
#' "Overall" which means that the value for each protein
#' (ie line in the expression data tab) is computed over all the samples ;
#' "within conditions" which means that the value for each protein
#' (ie line in the \code{SummarizedExperiment::assay()} data tab) is computed condition
#' by condition.
#' The available methods are described in [DaparToolshed::normalizeMethods()].
#'
#' @param obj An object of class \code{QFeatures} or \code{SummarizedExperiment}.
#'            If data is of class \code{QFeatures}, the last assay will be normalized.
#' @param method Define the normalization method used : `"GlobalQuantileAlignment"``, 
#'               `"QuantileCentering"`, `"MeanCentering"`, `"SumByColumns"`, `"LOESS"` or `"vsn"`.
#' @param conditions A vector of conditions in the dataset. 
#'                   If not provided, the vector `"Condition"` from the column metadata will be used.
#' @param type "overall" (shift all the sample distributions at once) or
#'             "within conditions" (shift the sample distributions within each condition at a time).
#' @param subset.norm A vector of index indicating rows to be used for normalization
#' @param quantile A float that corresponds to the quantile used to align the data.
#' @param scaling A boolean that indicates if the variance of the data have to
#'                be forced to unit (variance reduction) or not.
#' @param span xxx
#'
#' @return \code{QFeatures} including a new assay with normalized data or 
#'         \code{SummarizedExperiment} with normalized data.
#'
#' @author Manon Gaudin
#'
#' @examples
#' data(subR25pept)
#' normalized <- normalizeFunction(subR25pept, method = 'GlobalQuantileAlignment')
#' 
#' @export
#'

normalizeFunction <- function(obj,
                              method,
                              conditions = NULL,
                              type = "overall",
                              subset.norm = NULL,
                              quantile = 0.15,
                              scaling = FALSE,
                              span = 0.7
                              ){
  if (missing(obj))
    stop("'obj' is required.")
  if (missing(method))
    stop("'method' is required.")
  
  if (inherits(obj, 'QFeatures')){
    i <- length(obj)
    obj.se <- obj[[i]]
    SummarizedExperiment::colData(obj.se) <- SummarizedExperiment::colData(obj)
  } else if (inherits(obj, 'SummarizedExperiment')){
    obj.se <- obj
  } else { stop("'obj' must be of class 'SummarizedExperiment' or 'QFeatures'.") }
  
  if (is.null(conditions)){
    conditions <- SummarizedExperiment::colData(obj.se)$Condition
  }
  qdata <- SummarizedExperiment::assay(obj.se)
  
  normalized_assay <- switch(method,
                             GlobalQuantileAlignment = {
                               DaparToolshed::GlobalQuantileAlignment(qdata)
                             },
                             QuantileCentering ={
                               DaparToolshed::QuantileCentering(
                                 qData = qdata, 
                                 conds = conditions, 
                                 type = type, 
                                 subset.norm = subset.norm, 
                                 quantile = quantile
                               )
                             },
                             MeanCentering = {
                               DaparToolshed::MeanCentering(
                                 qData = qdata, 
                                 conds = conditions,
                                 type = type,
                                 scaling = scaling,
                                 subset.norm = subset.norm
                               )
                             },
                             SumByColumns = {
                               DaparToolshed::SumByColumns(
                                 qData = qdata,
                                 conds = conditions,
                                 type = type,
                                 subset.norm = subset.norm
                               )
                             }, 
                             LOESS = {
                               DaparToolshed::LOESS(
                                 qData = qdata,
                                 conds = conditions,
                                 type = type,
                                 span = span
                               )
                             },
                             vsn = {
                               DaparToolshed::vsn(
                                 qData = qdata,
                                 conds = conditions,
                                 type = type
                               )
                             })
  
  
  if (inherits(obj, 'QFeatures')){
    new.assay <- obj[[i]]
    SummarizedExperiment::assay(new.assay) <- normalized_assay
    
    new.dataset <- QFeatures::addAssay(obj, new.assay, 'Normalization')
  } else if (inherits(obj, 'SummarizedExperiment')){
    new.dataset <- obj
    SummarizedExperiment::assay(new.dataset) <- normalized_assay
  }
  
  return(new.dataset)
}
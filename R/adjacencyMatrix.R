

#' @title splits an adjacency matrix into specific and shared
#' 
#' @description 
#' Method to split an adjacency matrix into specific and shared
#'
#' @param X An adjacency matrix
#'
#' @return A list of two adjacency matrices
#'
#' @author Samuel Wieczorek
#' 
#' @examples 
#' data(subR25pept)
#' X <- BuildAdjacencyMatrix(subR25pept[[1]])
#' ll <- splitAdjacencyMat(X)
#'
#' @export
#'
splitAdjacencyMat <- function(X) {
  X <- as.matrix(X)
  hasShared <- length(which(rowSums(X) > 1)) > 0
  hasSpec <- length(which(rowSums(X) == 1)) > 0
  
  
  if (hasShared && !hasSpec) {
    tmpShared <- X
    tmpSpec <- X
    tmpSpec[which(rowSums(tmpSpec) > 1), ] <- 0
  } else if (!hasShared && hasSpec) {
    tmpSpec <- X
    tmpShared <- X
    tmpShared[which(rowSums(tmpShared) == 1), ] <- 0
  } else if (hasShared && hasSpec) {
    tmpSpec <- X
    tmpShared <- X
    tmpShared[which(rowSums(tmpShared) == 1), ] <- 0
    tmpSpec[which(rowSums(tmpSpec) > 1), ] <- 0
  } else {
    tmpSpec <- X
    tmpShared <- X
  }
  
  return(
    list(
      Xshared = 
        Matrix::Matrix(tmpShared,
          sparse = TRUE,
          dimnames = list(rownames(tmpShared), colnames(tmpShared))
        ), 
      Xspec = Matrix::Matrix(tmpSpec,
        sparse = TRUE,
        dimnames = list(rownames(tmpSpec), colnames(tmpSpec))
      )))
}






#' @title Function matrix of appartenance group
#' 
#' @description
#' Method to create a binary matrix with proteins in columns and peptides
#' in lines on a `SummarizedExperiment` object (peptides)
#'
#' @param obj.pep An object (peptides) of class `SummarizedExperiment`.
#'
#' @return A binary matrix
#'
#' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
#'
#' @examples
#' data(subR25pept)
#' BuildAdjacencyMatrix(subR25pept[[1]])
#'
#' @export
#' @importFrom Matrix Matrix
#' @importFrom stringr str_trim
#' @import QFeatures
#'

BuildAdjacencyMatrix <- function(obj.pep) {
  
  stopifnot(inherits(obj.pep, 'SummarizedExperiment'))
  protID <- parentProtId(obj.pep)
  
  data <- SummarizedExperiment::assay(obj.pep)
  PG <- SummarizedExperiment::rowData(obj.pep)[, protID]
  
  PG.l <- lapply(
    strsplit(as.character(PG), "[,;]+"),
    function(x) stringr::str_trim(x)
  )
  
  t <- table(data.frame(
    A = rep(seq_along(PG.l), lengths(PG.l)),
    B = unlist(PG.l)
  ))
  
  # if (unique == TRUE) {
  #   ll <- which(rowSums(t) > 1)
  #   if (length(ll) > 0) {
  #     t[ll, ] <- 0
  #   }
  # }
  
  X <- Matrix::Matrix(t,
    sparse = TRUE,
    dimnames = list(rownames(obj.pep), colnames(t))
  )
  
  return(X)
}
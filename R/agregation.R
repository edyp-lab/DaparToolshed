#' This function computes the number of proteins that are only defined by
#' specific peptides, shared peptides or a mixture of two.
#'
#' @title Computes the number of proteins that are only defined by
#' specific peptides, shared peptides or a mixture of two.
#'
#' @param matShared The adjacency matrix with both specific and
#' shared peptides.
#'
#' @return A list
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' obj <- Exp1_R25_pept[seq_len(20)]
#' MShared <- BuildAdjacencyMatrix(obj, protID, FALSE)
#' getProteinsStats(matShared = MShared)
#'
#' @export
#'
getProteinsStats <- function(matShared) {
  if (missing(matShared)) {
    stop("'matShared' is needed.")
  }
  if (is.null(matShared)) {
    stop("'matShared' is NULL")
  }
  
  nbPeptide <- 0
  
  ind.shared.Pep <- which(rowSums(as.matrix(matShared)) > 1)
  M.shared.Pep <- matShared[ind.shared.Pep, ]
  if (length(ind.shared.Pep) == 1) {
    j <- which(as.matrix(M.shared.Pep) == 0)
    M.shared.Pep <- M.shared.Pep[-j]
    pep.names.shared <- names(M.shared.Pep)
  } else {
    j <- which(colSums(as.matrix(M.shared.Pep)) == 0)
    M.shared.Pep <- M.shared.Pep[, -j]
    pep.names.shared <- colnames(M.shared.Pep)
  }
  
  
  ind.unique.Pep <- which(rowSums(as.matrix(matShared)) == 1)
  M.unique.Pep <- matShared[ind.unique.Pep, ]
  if (length(ind.unique.Pep) == 1) {
    j <- which(as.matrix(M.unique.Pep) == 0)
    M.unique.Pep <- M.unique.Pep[-j]
    pep.names.unique <- names(M.unique.Pep)
  } else {
    j <- which(colSums(as.matrix(M.unique.Pep)) == 0)
    M.unique.Pep <- M.unique.Pep[, -j]
    pep.names.unique <- colnames(M.unique.Pep)
  }
  
  
  
  protOnlyShared <- setdiff(
    pep.names.shared,
    intersect(pep.names.shared, pep.names.unique)
  )
  protOnlyUnique <- setdiff(
    pep.names.unique,
    intersect(pep.names.shared, pep.names.unique)
  )
  protMix <- intersect(pep.names.shared, pep.names.unique)
  
  
  
  return(
    list(
      nbPeptides = length(ind.unique.Pep) + length(ind.shared.Pep),
      nbSpecificPeptides = length(ind.unique.Pep),
      nbSharedPeptides = length(ind.shared.Pep),
      nbProt = length(protOnlyShared) + length(protOnlyUnique) + 
        length(protMix),
      protOnlyUniquePep = protOnlyUnique,
      protOnlySharedPep = protOnlyShared,
      protMixPep = protMix
    )
  )
}




#' This function creates a column for the protein dataset after aggregation
#' by using the previous peptide dataset.
#'
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#'
#' @param peptideData A data.frame of meta data of peptides. It is the fData
#' of the MSnset object.
#'
#' @param matAdj The adjacency matrix used to agregate the peptides data.
#'
#' @param columnName The name of the column in Biobase::fData(peptides_MSnset)
#' that the user wants to keep in the new protein data.frame.
#'
#' @param proteinNames The names of the protein in the new dataset
#' (i.e. rownames)
#'
#' @return A vector
#'
#' @author Samuel Wieczorek
#'
#' @example inst/extdata/examples/ex_BuildColumnToProteinDataset.R
#' @export
#'
BuildColumnToProteinDataset <- function(peptideData,
                                        matAdj,
                                        columnName,
                                        proteinNames) {
  nbProt <- ncol(matAdj)
  newCol <- rep("", nbProt)
  i <- 1
  for (p in proteinNames) {
    listeIndicePeptides <- names(which(matAdj[, p] == 1))
    listeData <- unique(
      as.character(
        peptideData[listeIndicePeptides, columnName], ";"
      )
    )
    newCol[i] <- paste0(listeData, collapse = ", ")
    i <- i + 1
  }
  return(newCol)
}




#' This function computes the number of peptides used to aggregate proteins.
#'
#' @title Compute the number of peptides used to aggregate proteins
#'
#' @param M A "valued" adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @return A vector of boolean which is the adjacency matrix
#' but with NA values if they exist in the intensity matrix.
#'
#' @author Alexia Dorffer
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' M <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)], protID, FALSE)
#' CountPep(M)
#'
#' @export
#'
CountPep <- function(M) {
  z <- M
  z[z != 0] <- 1
  return(z)
}


#' Method to compute the number of quantified peptides used for aggregating
#' each protein
#'
#' @title Computes the number of peptides used for aggregating each protein
#'
#' @param X An adjacency matrix
#'
#' @param pepData A data.frame of quantitative data
#'
#' @return A data.frame
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' library(QFeatures)
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' pepData <- assay(obj.pep)
#' GetNbPeptidesUsed(X, pepData)
#' 
GetNbPeptidesUsed <- function(X, pepData) {
  X <- as.matrix(X)
  pepData[!is.na(pepData)] <- 1
  pepData[is.na(pepData)] <- 0
  pep <- t(X) %*% pepData
  
  return(pep)
}




#' @title Computes the detailed number of peptides used for aggregating
#' each protein
#' 
#' @description 
#' Method to compute the detailed number of quantified peptides used for
#' aggregating each protein
#'
#' @param X An adjacency matrix
#'
#' @param qdata.pep A data.frame of quantitative data
#'
#' @return A list of two items
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' library(DaparToolshedData)
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)], protID, FALSE)
#' ll.n <- GetDetailedNbPeptidesUsed(X, 
#' assay(Exp1_R25_pept[seq_len(10)]))
#'
GetDetailedNbPeptidesUsed <- function(X, qdata.pep) {
  X <- as.matrix(X)
  qdata.pep[!is.na(qdata.pep)] <- 1
  qdata.pep[is.na(qdata.pep)] <- 0
  
  mat <- splitAdjacencyMat(X)
  return(list(
    nShared = t(mat$Xshared) %*% qdata.pep,
    nSpec = t(mat$Xspec) %*% qdata.pep
  ))
}



#'
#' @title Computes the detailed number of peptides for each protein
#' 
#' @description
#' Method to compute the detailed number of quantified peptides for each
#' protein
#'
#' @param X An adjacency matrix
#'
#' @return A data.frame
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' n <- GetDetailedNbPeptides(X)
#'
#' @export
#'
GetDetailedNbPeptides <- function(X) {
  mat <- splitAdjacencyMat(as.matrix(X))
  
  
  return(list(
    nTotal = rowSums(t(as.matrix(X))),
    nShared = rowSums(t(mat$Xshared)),
    nSpec = rowSums(t(mat$Xspec))
  ))
}




#' @title Function to create a histogram that shows the repartition of
#' peptides w.r.t. the proteins
#' 
#' @description
#' Method to create a plot with proteins and peptides on
#' a MSnSet object (peptides)
#'
#' @param mat An adjacency matrix.
#'
#' @return A histogram
#'
#' @author Alexia Dorffer, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' mat <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)], "Protein_group_IDs")
#' GraphPepProt(mat)
#'
#' @export
#'
GraphPepProt <- function(mat) {
  
  pkgs.require('graphics')
  
  
  if (is.null(mat)) {
    return(NULL)
  }
  
  mat <- as.matrix(mat)
  t <- t(mat)
  t <- apply(mat, 2, sum, na.rm = TRUE)
  tab <- table(t)
  position <- seq(1, length(tab), by = 3)
  conds <- names(tab)
  
  graphics::barplot(tab,
                    xlim = c(1, length(tab)),
                    xlab = "Nb of peptides",
                    ylab = "Nb of proteins",
                    names.arg = conds,
                    xaxp = c(1, length(tab), 3),
                    las = 1,
                    col = "orange"
  )
}






#' @title Function matrix of appartenance group
#' 
#' @description
#' Method to create a binary matrix with proteins in columns and peptides
#' in lines on a `MSnSet` object (peptides)
#'
#' @param obj.pep An object (peptides) of class `MSnSet`.
#'
#' @param protID The name of proteins ID column
#'
#' @param unique A boolean to indicate whether only the unique peptides must
#' be considered (TRUE) or if the shared peptides have to be integrated (FALSE).
#'
#' @return A binary matrix
#'
#' @author Florence Combes, Samuel Wieczorek, Alexia Dorffer
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protId <- "Protein_group_IDs"
#' obj <- Exp1_R25_pept[[1]]
#' BuildAdjacencyMatrix(obj[seq_len(10)], protId, TRUE)
#'
#' @export
#'

BuildAdjacencyMatrix <- function(obj.pep, protID, unique = TRUE) {
  
  pkgs.require(c("Biobase", "stringr", "Matrix", 'QFeatures'))
  
  data <- assay(obj.pep)
  PG <- rowData(obj.pep)[, protID]
  
  PG.l <- lapply(
    strsplit(as.character(PG), "[,;]+"),
    function(x) stringr::str_trim(x)
  )
  
  t <- table(data.frame(
    A = rep(seq_along(PG.l), lengths(PG.l)),
    B = unlist(PG.l)
  ))
  
  if (unique == TRUE) {
    ll <- which(rowSums(t) > 1)
    if (length(ll) > 0) {
      t[ll, ] <- 0
    }
  }
  
  X <- Matrix::Matrix(t,
                      sparse = TRUE,
                      dimnames = list(rownames(obj.pep), colnames(t))
  )
  
  return(X)
}


#' @title Test
#' @param X xxx
#' @export
#' @examples
#' NULL
#' 
ExtractUniquePeptides <- function(X){

  
  ll <- which(rowSums(X) > 1)
  if (length(ll) > 0) {
    X[ll, ] <- 0
  }

  return(X)
}


#' @title Compute the intensity of proteins with the sum of the intensities
#' of their peptides.
#' 
#' @description 
#' This function computes the intensity of proteins based on the sum of the
#' intensities of their peptides.
#'
#' @param obj.pep A matrix of intensities of peptides
#'
#' @param X An adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @return A matrix of intensities of proteins
#'
#' @author Alexia Dorffer
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(20)]
#' obj.pep.imp <- wrapper.impute.detQuant(obj.pep, na.type = c("Missing POV", "Missing MEC"))
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' ll.agg <- aggregateSum(obj.pep.imp, X)
#'
#' @export
#'
#'
aggregateSum <- function(obj.pep, X) {
  pkgs.require(c("QFeatures", "Biobase"))
  
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else {
    # Agregation of quanti data
    pepData <- 2^(assay(obj.pep))
    protData <- inner.sum(pepData, X)
    # Build protein dataset
    obj.prot <- finalizeAggregation(obj.pep, 
                                    pepData, 
                                    protData, 
                                    metacell$metacell, X
    )
    return(list(
      obj.prot = obj.prot,
      issues = NULL
    ))
  }
}



#' @title xxxx
#'
#' @param obj.pep xxxxx
#'
#' @param X xxxx
#'
#' @param init.method xxxxx
#'
#' @param method xxxxx
#'
#' @param n xxxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' obj.agg <- aggregateIterParallel(obj.pep, X)
#' }
#' 
#' @export
#' 
#' @import foreach
#'
aggregateIterParallel <- function(obj.pep,
                                  X,
                                  init.method = "Sum",
                                  method = "Mean",
                                  n = NULL
) {
  
  
  pkgs.require(c("Msnbase", "parallel", "doParallel", "foreach", "Biobase", 'QFeatures'))
  
  doParallel::registerDoParallel()
  obj.prot <- NULL
  
  # Step 1: Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else { # Step 2 : Agregation of quantitative data
    qData.pep <- 2^(assay(obj.pep))
    protData <- matrix(rep(0, ncol(X) * ncol(obj.pep)),
                       nrow = ncol(X),
                       dimnames = list(colnames(X), rep("cond", ncol(obj.pep)))
    )
    
    protData <- foreach::foreach(
      cond = seq_len(length(unique(Biobase::pData(obj.pep)$Condition))),
      .combine = cbind,
      .packages = "QFeatures"
    ) %dopar% {
      .conds <- Biobase::pData(obj.pep)$Condition
      condsIndices <- which(.conds == unique(.conds)[cond])
      qData <- qData.pep[, condsIndices]
      DAPAR::inner.aggregate.iter(qData, X, init.method, method, n)
    }
    
    protData <- protData[, colnames(assay(obj.pep))]
    colnames(protData) <- colnames(assay(obj.pep))
    
    # Step 3 : Build the protein dataset
    obj.prot <- finalizeAggregation(obj.pep, 
                                    qData.pep, 
                                    protData, 
                                    metacell$metacell, 
                                    X
    )
    
    return(list(
      obj.prot = obj.prot,
      issues = NULL
    ))
  }
}


#' Method to xxxxx
#'
#' @title xxxx
#'
#' @param pepData xxxxx
#'
#' @param X xxxx
#'
#' @param init.method xxx
#'
#' @param method xxx
#'
#' @param n xxxx
#'
#' @return xxxxx
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj[seq_len(10)], protID, FALSE)
#' qdata.agg <- inner.aggregate.iter(assay(obj[seq_len(10)]), X)
#'
#' @export
#'
inner.aggregate.iter <- function(
    pepData,
    X,
    init.method = "Sum",
    method = "Mean",
    n = NULL
) {
  
  if (!(init.method %in% c("Sum", "Mean"))) {
    warning("Wrong parameter init.method")
    return(NULL)
  }
  
  if (!(method %in% c("onlyN", "Mean"))) {
    warning("Wrong parameter method")
    return(NULL)
  }
  
  
  if (method == "onlyN" && is.null(n)) {
    warning("Parameter n is null")
    return(NULL)
  }
  
  yprot <- NULL
  switch(init.method,
         Sum = yprot <- inner.sum(pepData, X),
         Mean = yprot <- inner.mean(pepData, X)
  )
  conv <- 1
  
  while (conv > 10**(-10)) {
    mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot[is.na(mean.prot)] <- 0
    
    X.tmp <- mean.prot * X
    X.new <- X.tmp / rowSums(as.matrix(X.tmp), na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    # l'appel a la fonction ci-dessous depend des parametres choisis par
    # l'utilisateur
    switch(method,
           Mean = yprot <- inner.mean(pepData, X.new),
           onlyN = yprot <- inner.aggregate.topn(pepData, X.new, "Mean", n)
    )
    
    mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
    mean.prot.new[is.na(mean.prot.new)] <- 0
    
    conv <- mean(abs(mean.prot.new - mean.prot))
  }
  return(as.matrix(yprot))
}



#' @title xxxx
#'
#' @param obj.pep xxxxx
#'
#' @param X xxxx
#'
#' @param init.method xxxxx
#'
#' @param method xxxxx
#'
#' @param n xxxx
#'
#' @return A protein object of class `MSnSet`
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)], protID, FALSE)
#' ll.agg <- aggregateIter(Exp1_R25_pept[seq_len(10)], X = X)
#'
#' @export
#'
#'
aggregateIter <- function(
    obj.pep,
    X,
    init.method = "Sum",
    method = "Mean",
    n = NULL
) {
  pkgs.require(c("Biobase", 'QFeatures'))
  
  obj.prot <- NULL
  
  # Step 1 : Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else {
    # Step 2: Agregation of quantitative data
    # For each condition, reproduce iteratively
    # Initialisation: At initialization step, take "sum overall" and  
    # matAdj = X
    # for simplicity.
    # Note : X <- as.matrix(X)
    qData.pep <- 2^(assay(obj.pep))
    
    protData <- matrix(rep(0, ncol(X) * ncol(obj.pep)),
                       nrow = ncol(X),
                       dimnames = list(colnames(X), rep("cond", ncol(obj.pep)))
    )
    
    for (cond in unique(Biobase::pData(obj.pep)$Condition)) {
      condsIndices <- which(Biobase::pData(obj.pep)$Condition == cond)
      qData <- qData.pep[, condsIndices]
      protData[, condsIndices] <- DAPAR::inner.aggregate.iter(
        qData,
        X,
        init.method,
        method,
        n
      )
      .tmp <- assay(obj.pep)
      colnames(protData)[condsIndices] <- colnames(.tmp)[condsIndices]
    }
    
    # Step 3: Build the protein dataset
    obj.prot <- finalizeAggregation(
      obj.pep, 
      qData.pep, 
      protData, 
      metacell$metacell, 
      X
    )
    return(list(
      obj.prot = obj.prot,
      issues = NULL
    ))
  }
}


#' @title Compute the intensity of proteins as the mean of the intensities
#' of their peptides.
#' 
#' @description 
#' #' This function computes the intensity of proteins as the mean of the
#' intensities of their peptides.
#'
#' @param obj.pep A peptide object of class `MSnSet`
#'
#' @param X An adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @return A matrix of intensities of proteins
#'
#' @author Alexia Dorffer
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' obj.pep.imp <- wrapper.impute.detQuant(obj.pep, na.type = c("Missing POV", "Missing MEC"))
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep.imp, protID, FALSE)
#' ll.agg <- aggregateMean(obj.pep.imp, X)
#'
#' @export
#'
#'
aggregateMean <- function(obj.pep, X) {
  pkgs.require(c("QFeatures", "Biobase"))
  
  
  obj.prot <- NULL
  # Agregation of metacell data
  #cat("Aggregate metacell data...\n")
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else {
    # Step 2: Agregation of quantitative data
    #cat("Computing quantitative data for proteins ...\n")
    pepData <- 2^(assay(obj.pep))
    protData <- inner.mean(pepData, as.matrix(X))
    
    # Step 3: Build protein dataset
    #cat("Building the protein dataset...\n")
    obj.prot <- finalizeAggregation(obj.pep, 
                                    pepData, 
                                    protData, 
                                    metacell$metacell, 
                                    X)
    
    return(list(
      obj.prot = obj.prot,
      issues = NULL
    ))
  }
}



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
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
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
  
  return(list(Xshared = tmpShared, Xspec = tmpSpec))
}



#' @title xxxx
#'
#' @param pepData A data.frame of quantitative data
#'
#' @param X An adjacency matrix
#'
#' @return A matrix
#'
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @examples 
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj, protID, FALSE)
#' inner.sum(assay(obj), X)

inner.sum <- function(pepData, X) {
  X <- as.matrix(X)
  pepData[is.na(pepData)] <- 0
  Mp <- t(as.matrix(X)) %*% pepData
  return(Mp)
}


#' @title xxxx
#'
#' @param pepData A data.frame of quantitative data
#'
#' @param X An adjacency matrix
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#' 
#' @export
#' 
#' @examples 
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj, protID, FALSE)
#' inner.mean(assay(obj), X)
#' 
inner.mean <- function(pepData, X) {
  X <- as.matrix(X)
  Mp <- inner.sum(pepData, X)
  Mp <- Mp / GetNbPeptidesUsed(X, pepData)
  
  return(Mp)
}




#' @title xxxx
#'
#' @param pepData A data.frame of quantitative data
#'
#' @param X An adjacency matrix
#'
#' @param method xxxxx
#'
#' @param n xxxxx
#'
#' @return xxxxx
#'
#' @author Samuel Wieczorek
#' 
#' @export
#'
#' @examples 
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj, protID, FALSE)
#' inner.aggregate.topn(assay(obj), X)
#' 
inner.aggregate.topn <- function(pepData, X, method = "Mean", n = 10) {
  
  
  pkgs.require("stats")
  
  
  X <- as.matrix(X)
  med <- apply(pepData, 1, stats::median)
  xmed <- as(X * med, "dgCMatrix")
  for (c in seq_len(ncol(X))) {
    v <- order(xmed[, c], decreasing = TRUE)[seq_len(n)]
    l <- v[which((xmed[, c])[v] != 0)]
    
    if (length(l) > 0) {
      diff <- setdiff(which(X[, c] == 1), l)
      if (length(diff)) {
        X[diff, c] <- 0
      }
    }
  }
  
  Mp <- NULL
  switch(method,
         Mean = Mp <- inner.mean(pepData, X),
         Sum = Mp <- inner.sum(pepData, X)
  )
  
  return(Mp)
}

#' This function computes the intensity of proteins as the sum of the
#' intensities of their n best peptides.
#'
#' @title Compute the intensity of proteins as the sum of the
#' intensities of their n best peptides.
#'
#' @param obj.pep A matrix of intensities of peptides
#'
#' @param X An adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @param method xxx
#'
#' @param n The maximum number of peptides used to aggregate a protein.
#'
#' @return A matrix of intensities of proteins
#'
#' @author Alexia Dorffer, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' protID <- "Protein_group_IDs"
#' X <- BuildAdjacencyMatrix(obj.pep, protID, FALSE)
#' ll.agg <- aggregateTopn(obj.pep, X, n = 3)
#'
#' @export
#'
#'
aggregateTopn <- function(obj.pep,
                          X,
                          method = "Mean",
                          n = 10) {
  pkgs.require(c("QFeatures", "Biobase"))
  
  
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep)
  
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else {
    # Step 2 : Agregation of quantitative data
    pepData <- 2^(assay(obj.pep))
    protData <- inner.aggregate.topn(pepData, X, method = method, n)
    
    # Step 3: Build the protein dataset
    obj.prot <- finalizeAggregation(
      obj.pep, 
      pepData, 
      protData, 
      metacell$metacell, 
      X)
    
    return(list(
      obj.prot = obj.prot,
      issues = NULL
    ))
  }
}




#' Method to finalize the aggregation process
#'
#' @title Finalizes the aggregation process
#'
#' @param obj.pep A peptide object of class `MSnSet`
#'
#' @param pepData xxxx
#'
#' @param X An adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @param protData xxxxx
#'
#' @param protMetacell xxx
#'
#' @return A protein object of class `MSnSet`
#'
#' @author Samuel Wieczorek
#'
#' @export
#'
#' @examples 
#' NULL
#'
#'
finalizeAggregation <- function(obj.pep, pepData, protData, protMetacell, X) {
  
  pkgs.require("utils")
  
  
  
  if (missing(obj.pep)) {
    stop("'obj.pep' is missing")
  }
  if (missing(pepData)) {
    stop("'pepData' is missing")
  }
  if (missing(protData)) {
    stop("'protData' is missing")
  }
  if (missing(protMetacell)) {
    stop("'protMetacell' is missing")
  }
  if (missing(X)) {
    stop("'X' is missing")
  }
  
  # (obj.pep, pepData, protData, metacell, X)
  
  protData <- as.matrix(protData)
  X <- as.matrix(X)
  protData[protData == 0] <- NA
  protData[is.nan(protData)] <- NA
  protData[is.infinite(protData)] <- NA
  
  temp <- GetDetailedNbPeptidesUsed(X, pepData)
  
  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", 
                                   colnames(pepData), 
                                   sep = "")
  rownames(pepSharedUsed) <- colnames(X)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", 
                                 colnames(pepData), 
                                 sep = "")
  rownames(pepSpecUsed) <- colnames(X)
  
  pepTotalUsed <- as.matrix(GetNbPeptidesUsed(X, pepData))
  colnames(pepTotalUsed) <- paste("pepTotal.used.", 
                                  colnames(pepData),
                                  sep = "")
  rownames(pepTotalUsed) <- colnames(X)
  
  n <- GetDetailedNbPeptides(X)
  
  fd <- data.frame(
    proteinId = rownames(protData),
    nPepTotal = n$nTotal,
    nPepShared = n$nShared,
    nPepSpec = n$nSpec
  )
  fd <- cbind(fd,
              pepSpecUsed,
              pepSharedUsed,
              pepTotalUsed,
              protMetacell,
              stringsAsFactors = FALSE
  )
  
  rownames(fd) <- colnames(X)
  obj.prot <- MSnSet(
    exprs = log2(protData),
    fData = fd,
    pData = Biobase::pData(obj.pep)
  )
  
  obj.prot@experimentData@other <- obj.pep@experimentData@other
  obj.prot@experimentData@other$typeOfData <- "protein"
  obj.prot@experimentData@other$keyId <- "proteinId"
  obj.prot@experimentData@other$proteinId <- "proteinId"
  
  
  obj.prot@experimentData@other$Prostar_Version <- NA
  obj.prot@experimentData@other$DAPAR_Version <- NA
  
  tryCatch(
    {
      find.package("Prostar")
      find.package("DAPAR")
      prostar.ver <- Biobase::package.version("Prostar")
      dapar.ver <- Biobase::package.version("DAPAR")
      
      
      obj.prot@experimentData@other$Prostar_Version <- prostar.ver
      obj.prot@experimentData@other$DAPAR_Version <- dapar.ver
    },
    error = function(e) {
      obj.prot@experimentData@other$Prostar_Version <- NA
      obj.prot@experimentData@other$DAPAR_Version <- NA
    }
  )
  
  return(obj.prot)
}




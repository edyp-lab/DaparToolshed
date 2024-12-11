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
#' obj <- Exp1_R25_pept[seq_len(20)]
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' getProteinsStats(X)
#'
#' @export
#'
getProteinsStats <- function(X) {
  if (missing(matShared)) {
    stop("'matShared' is needed.")
  }
  
  stopifnot(!is.null(matShared))
  
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
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' CountPep(X)
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
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' GetNbPeptidesUsed(obj.pep)
#' 
GetNbPeptidesUsed <- function(pepData) {
  
  stopifnot(inherits(pepData, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(pepData)))
  
  assay(pepData)[!is.na(assay(pepData))] <- 1
  assay(pepData)[is.na(assay(pepData))] <- 0
  
  X <- as.matrix(adjacencyMatrix(pepData))
  pep <- t(X) %*% assay(pepData)
  
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
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' ll.n <- GetDetailedNbPeptidesUsed(obj.pep[[length(obj.pep)]])
#'
GetDetailedNbPeptidesUsed <- function(obj.pep) {
  stopifnot(inherits(obj.pep, 'SummarizedExperiment'))
  stopifnot('adjacencyMatrix' %in% names(rowData(obj.pep)))
  
  X <- as.matrix(adjacencyMatrix(obj.pep))
  assay(obj.pep)[!is.na(assay(obj.pep))] <- 1
  assay(obj.pep)[is.na(assay(obj.pep))] <- 0
  
  mat <- splitAdjacencyMat(X)
  return(list(
    nShared = t(mat$Xshared) %*% assay(obj.pep),
    nSpec = t(mat$Xspec) %*% assay(obj.pep)
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
#' X <- BuildAdjacencyMatrix(obj.pep)
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
#' X <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)])
#' GraphPepProt(X)
#'
#' @export
#' 
#' @import graphics
#'
GraphPepProt <- function(mat) {
  
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







#' @title Test
#' @param X xxx
#' @export
#' @examples
#' #' data(Exp1_R25_pept, package="DaparToolshedData")
#' X <- BuildAdjacencyMatrix(Exp1_R25_pept[seq_len(10)])
#' ExtractUniquePeptides(X)
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
#' #obj.pep.imp <- wrapper.impute.detQuant(obj.pep[[length(obj.pep)]], na.type = c("Missing POV", "Missing MEC"))
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' ll.agg <- aggregateProstar2(obj.pep, length(obj.pep))
#'
#' @export
#'
#' @import QFeatures
#' 
aggregateProstar2 <- function(obj, i, FUN = 'Sum') {
  stopifnot(inherits(obj, "QFeatures"))
  stopifnot('adjacencyMatrix' %in% names(rowData(obj.pep[[i]])))
   
    X <- adjacencyMatrix(obj[[i]])
    qMeta <- qMetacell(obj[[i]])
    level <- typeDataset(obj[[i]])
    conds <- colData(obj)$Condition
    
    obj.pep <- obj[[i]]
    
    obj.prot <- NULL
    # Aggregation of metacell data
    metacell <- aggQmetacell(qMeta, X, level, conds)
    if (!is.null(metacell$issues)) {
      return(list(
        obj.prot = NULL,
        issues = metacell$issues
      ))
    } else {
      # Step 2: Agregation of quantitative data
      #cat("Computing quantitative data for proteins ...\n")
      protData <- switch(FUN,
        Sum = inner.sum(obj.pep),
        Mean = inner.mean(obj.pep)
      )
      
      # Step 3: Build protein dataset
      obj.prot <- finalizeAggregation(obj, 
        protData, 
        metacell$metacell, 
        X)
      
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
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep)
#' level <- 'protein'
#' obj.agg <- aggregateIterParallel(obj.pep, X, level)
#' }
#' 
#' @export
#' 
#' @import foreach
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import QFeatures
#'
aggregateIterParallel <- function(obj.pep,
  X,
  init.method = "Sum",
  method = "Mean",
  n = NULL,
  level
) {
  
  
  #pkgs.require(c("Msnbase", "parallel", "doParallel", "foreach", "Biobase", 'QFeatures'))
  
  doParallel::registerDoParallel()
  obj.prot <- NULL
  
  # Step 1: Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep, level)
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
      inner.aggregate.iter(qData, X, init.method, method, n)
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
#' X <- BuildAdjacencyMatrix(obj[seq_len(10)])
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
#' obj <- Exp1_R25_pept[[2]]
#' obj <- obj[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj)
#' ll.agg <- aggregateIter(obj.pep = obj)
#'
#' @export
#'
#'
aggregateIter <- function(
    obj.pep,
    X,
    init.method = "Sum",
    method = "Mean",
    n = NULL,
  level
) {
  
  obj.prot <- NULL
  
  # Step 1 : Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep, level)
  
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
      protData[, condsIndices] <- inner.aggregate.iter(
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


#' @title xxxx
#'
#' @param pepData A data.frame of quantitative data
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
#' adjacencyMatrix(obj[[length(obj)]]) <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' inner.sum(obj[[length(obj)]])

inner.sum <- function(pepData) {
  stopifnot(inherits(pepData, 'SummarizedExperiment'))
  
  X <- as.matrix(adjacencyMatrix(pepData))
  assay(pepData)[is.na(assay(pepData))] <- 0
  
  Mp <- t(as.matrix(X)) %*% (2^assay(pepData))
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
#' X <- BuildAdjacencyMatrix(obj)
#' inner.mean(assay(obj), X)
#' 
inner.mean <- function(pepData) {
  stopifnot(inherits(pepData, 'SummarizedExperiment'))
  
  X <- as.matrix(adjacencyMatrix(pepData))
  
  Mp <- inner.sum(pepData)
  Mp <- Mp / GetNbPeptidesUsed(pepData)
  
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
#' X <- BuildAdjacencyMatrix(obj)
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
#' X <- BuildAdjacencyMatrix(obj.pep)
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' ll.agg <- aggregateTopn(obj.pep, X, n = 3)
#'
#' @export
#'
#'
aggregateTopn <- function(obj.pep,
  X,
  method = "Mean",
  n = 10,
  level) {
  pkgs.require(c("QFeatures", "Biobase"))
  
  
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- AggregateMetacell(X = X, obj.pep = obj.pep, level)
  
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
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#'
#' @import utils
#' 
finalizeAggregation <- function(obj.pep, protData, protMetacell, X) {

  if (missing(obj.pep)) {
    stop("'obj.pep' is missing")
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
  
  temp <- GetDetailedNbPeptidesUsed(obj.pep[[length(obj.pep)]])

  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", 
                                   colnames(assay(obj.pep[[length(obj.pep)]])), 
                                   sep = "")
  rownames(pepSharedUsed) <- colnames(X)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", 
                                 colnames(assay(obj.pep[[length(obj.pep)]])), 
                                 sep = "")
  rownames(pepSpecUsed) <- colnames(X)
  
  pepTotalUsed <- as.matrix(GetNbPeptidesUsed(obj.pep[[length(obj.pep)]]))
  colnames(pepTotalUsed) <- paste("pepTotal.used.", 
                                  colnames(assay(obj.pep[[length(obj.pep)]])),
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
  
  
  
  prot.se <- SummarizedExperiment(
    assays = protData,
    rowData = fd,
    colData = colData(obj.pep)
  )
  ## Add the qMetacell to the new assay
  qMetacell(prot.se) <- protMetacell
  
  
  
  ## Remove the qMetacell that should be dropped anyway and not be aggregated
  ## within QFeatures::aggregateFeatures
  
  
  
  
  X.spec <- X.shared <- X
  
  X.spec[which(rowSums(as.matrix(X.spec)) > 1), ] <- 0
  X.shared[which(rowSums(as.matrix(X.shared)) == 1), ] <- 0
  
  .allPep <- t(as.matrix(X)) %*% !is.na(assay(obj.pep[[length(obj.pep)]]))
  .specPep <- t(as.matrix(X.spec)) %*% !is.na(assay(obj.pep[[length(obj.pep)]]))
  .sharedPep <- t(as.matrix(X.shared)) %*% !is.na(assay(obj.pep[[length(obj.pep)]]))
  rowData(prot.se)[["allPeptidesUsed"]] <-.allPep
  rowData(prot.se)[["specPeptidesUsed"]] <- .specPep
  rowData(prot.se)[["sharedPeptidesUsed"]] <- .sharedPep
  
  
  ## Enrich the new assay
  typeDataset(prot.se) <- "protein"
  idcol(prot.se) <- "proteinId"
  metadata(prot.se)[['Prostar2_Version']] <- NA
  metadata(prot.se)[['DaparToolshed_Version']] <- NA
  
  
  
  tryCatch(
    {
      find.package("Prostar")
      find.package("DaparToolshed")
      prostar.ver <- Biobase::package.version("Prostar2")
      dapar.ver <- Biobase::package.version("DaparToolshed")
      metadata(prot.se)[['Prostar_Version']] <- prostar.ver
      metadata(prot.se)[['DaparToolshed_Version']] <- dapar.ver
    },
    error = function(e) {
      metadata(prot.se)[['Prostar2_Version']] <- NA
      metadata(prot.se)[['DaparToolshed_Version']] <- NA
    }
  )
  
  ## Add the assay to the QFeatures object
  obj.prot <- addAssay(obj.pep,
    prot.se,
    name = 'Aggregated'
  )
  
  ## Link the input assay to the aggregated assay
  addAssayLink(obj.prot,
    from = length(obj.prot) - 1,
    to = 'Aggregated',
    varFrom = 'adjacencyMatrix',
    varTo = 'adjacencyMatrix'
  )
  

  return(obj.prot)
}




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
#' obj.last <- obj[[length(obj)]]
#' X <- BuildAdjacencyMatrix(obj.last)
#' getProteinsStats(X)
#'
#' @export
#'
getProteinsStats <- function(X) {
  if (missing(X)) {
    stop("'X' is needed.")
  }
  
  stopifnot(!is.null(X))
  
  nbPeptide <- 0
  
  ind.shared.Pep <- which(rowSums(as.matrix(X)) > 1)
  M.shared.Pep <- X[ind.shared.Pep, ]
  if (length(ind.shared.Pep) == 1) {
    j <- which(as.matrix(M.shared.Pep) == 0)
    M.shared.Pep <- M.shared.Pep[-j]
    pep.names.shared <- names(M.shared.Pep)
  } else {
    j <- which(colSums(as.matrix(M.shared.Pep)) == 0)
    M.shared.Pep <- M.shared.Pep[, -j]
    pep.names.shared <- colnames(M.shared.Pep)
  }
  
  
  ind.unique.Pep <- which(rowSums(as.matrix(X)) == 1)
  M.unique.Pep <- X[ind.unique.Pep, ]
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
      nbProt = length(protOnlyShared) + length(protOnlyUnique) + length(protMix),
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
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' CountPep(X)
#'
#' @export
#'
CountPep <- function(X) {
  #z <- M
  X[X != 0] <- 1
  return(X)
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
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' GetNbPeptidesUsed(assay(last.obj), X)
#' 
GetNbPeptidesUsed <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  
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
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' ll.n <- GetDetailedNbPeptidesUsed(assay(last.obj), X)
#'
GetDetailedNbPeptidesUsed <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  
  pepData[!is.na(pepData)] <- 1
  pepData[is.na(pepData)] <- 0
  
  mat <- splitAdjacencyMat(X)
  return(list(
    nShared = t(mat$Xshared) %*% pepData,
    nSpec = t(mat$Xspec) %*% pepData
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
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
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
#' obj <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- Exp1_R25_pept[[length(Exp1_R25_pept)]]
#' X <- BuildAdjacencyMatrix(last.obj)
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
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- Exp1_R25_pept[[length(Exp1_R25_pept)]]
#' X <- BuildAdjacencyMatrix(last.obj)
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
#' last.obj <- obj.pep[[length(obj.pep)]]
#' #obj.pep.imp <- wrapper.impute.detQuant(obj.pep[[length(obj.pep)]], na.type = c("Missing POV", "Missing MEC"))
#' X <- BuildAdjacencyMatrix(last.obj)
#' ll.agg <- aggregateProstar2(obj.pep, length(obj.pep), X = X)
#'
#' @export
#'
#' @import QFeatures
#' 
aggregateProstar2 <- function(obj, i, FUN = 'Sum', X) {
  stopifnot(inherits(obj, "QFeatures"))
   
  

    qMeta <- qMetacell(obj[[i]])
    level <- typeDataset(obj[[i]])
    conds <- colData(obj)$Condition
    pepData <- assay(obj[[i]])
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
      pepData <- 2^pepData
      protData <- switch(FUN,
        Sum = inner.sum(pepData, X),
        Mean = inner.mean(pepData, X)
      )
      
      # Step 3: Build protein dataset
      obj.prot <- finalizeAggregation(obj, 
        i,
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
#' last.se <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.se)
#' conds <- colData(obj.pep)$Condition
#' obj.agg <- aggregateIterParallel(obj.pep = obj.pep, i = length(obj.pep), X = X, n = 3)
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
  i,
  X,
  init.method = "Sum",
  method = "Mean",
  n = NULL) {
  
  stopifnot(inherits(obj.pep, "QFeatures"))
  #pkgs.require(c("Msnbase", "parallel", "doParallel", "foreach", "Biobase", 'QFeatures'))
  
  doParallel::registerDoParallel()
  obj.prot <- NULL
  level <- typeDataset(obj.pep[[i]])
  
  # Step 1: Agregation of metacell data
  qMeta <- qMetacell(obj.pep[[i]])
  level <- typeDataset(obj.pep[[i]])
  pepData <- assay(obj.pep[[i]])
  conds <- colData(obj.pep)$Condition
  
  obj.prot <- NULL
  # Aggregation of metacell data
  metacell <- aggQmetacell(qMeta, X, level, conds)
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else { # Step 2 : Agregation of quantitative data
    qData.pep <- 2^(assay(pepData))
    
    protData <- matrix(rep(0, ncol(X) * ncol(qData.pep)),
                       nrow = ncol(X),
                       dimnames = list(colnames(X), rep("cond", ncol(qData.pep)))
    )
    
    protData <- foreach::foreach(
      cond = seq_len(length(unique(conds))),
      .combine = cbind,
      .packages = "QFeatures"
    ) %dopar% {
      condsIndices <- which(conds == unique(conds)[cond])
      qData <- qData.pep[, condsIndices]

      DaparToolshed::inner.aggregate.iter(qData, X, init.method, method, n)
    }
    
    
    protData <- protData[, colnames(pepData)]
    colnames(protData) <- colnames(pepData)
    
    
    
    # Step 3 : Build the protein dataset
    obj.prot <- finalizeAggregation(
      obj.pep, 
      i,
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
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' qdata.agg <- inner.aggregate.iter(assay(obj[[length(obj)]]), X)
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
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' i.sum <- inner.sum(assay(obj[[length(obj)]]), X)

inner.sum <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  
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
#' X <- BuildAdjacencyMatrix(obj)
#' i.mean <- inner.mean(assay(obj), X)
#' 
inner.mean <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  
  pepData[is.na(pepData)] <- 0
  Mp <- inner.sum(pepData, X)
  Mp <- Mp / GetNbPeptidesUsed(pepData, X)
  
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
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' inner.aggregate.topn(assay(obj[[length(obj)]]), X, n=3)
#' 
#' @import stats
#' 
inner.aggregate.topn <- function(pepData, X, method = "Mean", n = 10) {
  #stopifnot(inherits(pepData, 'SummarizedExperiment'))
  
 # pkgs.require("stats")
  
  
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
#' X <- BuildAdjacencyMatrix(obj.pep[[length(obj.pep)]])
#' adjacencyMatrix(obj.pep[[length(obj.pep)]]) <- X
#' conds <- colData(obj.pep)$Condition
#' ll.agg <- aggregateTopn(obj.pep[[length(obj.pep)]], n = 3, conds = conds)
#'
#' @export
#'
#' @import QFeatures
#' 
aggregateTopn <- function(obj.pep,
  i,
  method = "Mean",
  n = 10,
  X) {
  #pkgs.require(c("QFeatures", "Biobase"))
  stopifnot(inherits(obj.pep, 'QFeatures'))
  
  #obj.pep <- obj.pep[[i]]
  obj.prot <- NULL
  
  # Agregation of metacell data
  metacell <- aggQmetacell(
    qMeta = qMetacell(obj.pep[[i]]),
    X = X,
    level = typeDataset(obj.pep[[i]]),
    conds = colData(obj.pep)$Condition
  )
  
  if (!is.null(metacell$issues)) {
    return(list(
      obj.prot = NULL,
      issues = metacell$issues
    ))
  } else {
    # Step 2 : Agregation of quantitative data
    protData <- inner.aggregate.topn(2^assay(obj.pep[[i]]), X, method = method, n)

    # Step 3: Build the protein dataset
    obj.prot <- finalizeAggregation(
      obj.pep, 
      i,
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
finalizeAggregation <- function(obj.pep, i, protData, protMetacell, X) {
stopifnot(inherits(obj.pep, "QFeatures"))
  
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

  protData <- as.matrix(protData)
  X <- as.matrix(X)
  protData[protData == 0] <- NA
  protData[is.nan(protData)] <- NA
  protData[is.infinite(protData)] <- NA
  
  temp <- GetDetailedNbPeptidesUsed(assay(obj.pep[[i]]), X)

  pepSharedUsed <- as.matrix(temp$nShared)
  colnames(pepSharedUsed) <- paste("pepShared.used.", 
                                   colnames(assay(obj.pep[[i]])), 
                                   sep = "")
  rownames(pepSharedUsed) <- colnames(X)
  
  pepSpecUsed <- as.matrix(temp$nSpec)
  colnames(pepSpecUsed) <- paste("pepSpec.used.", 
                                 colnames(assay(obj.pep[[i]])), 
                                 sep = "")
  rownames(pepSpecUsed) <- colnames(X)
  
  pepTotalUsed <- as.matrix(GetNbPeptidesUsed(assay(obj.pep[[i]]), X))
  colnames(pepTotalUsed) <- paste("pepTotal.used.", 
                                  colnames(assay(obj.pep[[i]])),
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
  
  .allPep <- t(as.matrix(X)) %*% !is.na(assay(obj.pep[[i]]))
  .specPep <- t(as.matrix(X.spec)) %*% !is.na(assay(obj.pep[[i]]))
  .sharedPep <- t(as.matrix(X.shared)) %*% !is.na(assay(obj.pep[[i]]))
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
      find.package("Prostar2")
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
  

  rowData(obj.prot[[length(obj.prot) - 1]])['adjacencyMatrix'] <- NULL
  adjacencyMatrix(obj.prot[[length(obj.prot) - 1]]) <- BuildAdjacencyMatrix(obj.prot[[length(obj.prot) - 1]])
  
  
  ## Link the input assay to the aggregated assay
  addAssayLink(obj.prot,
    from = length(obj.prot) - 1,
    to = 'Aggregated',
    varFrom = 'adjacencyMatrix',
    varTo = 'adjacencyMatrix'
  )
  

  return(obj.prot)
}




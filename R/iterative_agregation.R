#'
#'
#'
#'
#'
#' #' Method to aggregate peptides into proteins with the iterative approach.
#' #'
#' #' @title Method to aggregate peptides into proteins with the iterative approach.
#' #'
#' #' @param qPepData A data.frame of quantitative data
#' #'
#' #' @param X An adjacency matrix
#' #'
#' #' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' #'
#' #' @param method The method used for xxx. Default value is 'Mean'.
#' #'
#' #' @param n An integer which is xxx
#' #'
#' #' @return A  matrix containing the quantitative aggregated data for proteins.
#' #'
#' #' @author Samuel Wieczorek, Thomas Burger
#' #'
#' #' @examples
#' #' library(QFeatures)
#' #' datat(ft)
#' #' obj <- Exp1_R25_pept[seq_len(1000),]
#' #' obj <- addListAdjacencyMatrices(obj, 2)
#' #' X <- GetAdjMat(obj[[2]])$all
#' #' qPepData <- assay(obj[[2]])
#' #' inner.aggregate.iter(qPepData, X)
#' #'
#' #' @rdname iterative-aggregation
#' #'
#' inner.aggregate.iter <- function(qData,
#'                                  X,
#'                                  init.method = 'sum',
#'                                  iter.method = 'mean',
#'                                  ...){
#'   #browser()
#'
#'   stopifnot(init.method %in% c('sum', 'mean'))
#'   stopifnot(iter.method %in% c("mean", "topn"))
#'
#'   yprot <- switch(init.method,
#'                   mean = colMeansMat(qData, X),
#'                   sum = colSumsMat(qData, X)
#'                   )
#'   conv <- 1
#'
#'   while(conv > 10**(-10)){
#'     mean.prot <- rowMeans(as.matrix(yprot), na.rm = TRUE)
#'     mean.prot[is.na(mean.prot)] <- 0
#'
#'     X.tmp <- mean.prot*X
#'     X.new <- X.tmp/rowSums(as.matrix(X.tmp), na.rm = TRUE)
#'     X.new[is.na(X.new)] <- 0
#'
#'     # l'appel a la fonction ci-dessous depend des parametres choisis par l'utilisateur
#'    yprot <- switch(iter.method,
#'                    mean = colMeansMat(qData, X.new),
#'                    topn = {
#'                      X.topn <- subAdjMat_topnPeptides(X.new, qData, 'rowMeans', ...)
#'                    }
#'                    )
#'     #yprot <- do.call(iter.method, list(qData, X, ...))
#'
#'     mean.prot.new <- rowMeans(as.matrix(yprot), na.rm = TRUE)
#'     mean.prot.new[is.na(mean.prot.new)] <- 0
#'
#'     conv <- mean(abs(mean.prot.new - mean.prot))
#'   }
#'   return(as.matrix(yprot))
#' }
#'
#'
#'
#' #' Method to aggregate peptides into proteins with the iterative approach with use of parallelism.
#' #'
#' #' @title Method to aggregate peptides into proteins with the iterative approach with use of parallelism.
#' #'
#' #' @param qPepData A matrix of intensities of peptides.
#' #'
#' #' @param X An adjacency matrix in which lines and columns correspond
#' #' respectively to peptides and proteins.
#' #'
#' #' @param conditions xxxx
#' #'
#' #' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' #'
#' #' @param method The method used for xxx. Default value is 'Mean'.
#' #'
#' #' @param n xxxx
#' #'
#' #' @return A matrix of intensities of proteins
#' #'
#' #' @author Samuel Wieczorek
#' #'
#' #' @examples
#' #' library(doParallel)
#' #' library(QFeatures)
#' #' data(Exp1_R25_pept, package="DaparToolshedData")
#' #' obj <- Exp1_R25_pept[seq_len(1000),]
#' #' obj <- addListAdjacencyMatrices(obj, 2)
#' #' X <- as.matrix(GetAdjMat(obj[[2]])$all)
#' #' conditions <- colData(obj)$Condition
#' #' aggIterParallel(assay(obj,2), X, conditions)
#' #'
#' #' @export
#' #'
#' #' @import doParallel
#' #' @import foreach
#' #'
#' #' @rdname iterative-aggregation
#' #'
#' aggIterParallel <- function(qPepData, X, conditions=NULL, init.method='Sum', method='Mean', n=NULL){
#'   if (is.null(conditions)){
#'     warning('The parameter \'conditions\' is NULL: the aggregation cannot be process.')
#'     return(NULL)
#'   }
#'
#'
#'   pkgs.require("Matrix")
#'
#'
#'   doParallel::registerDoParallel()
#'
#'   #qPepData <- 2^(qPepData)
#'   protData <- matrix(rep(0,ncol(X)*nrow(X)), nrow=ncol(X))
#'   cond <- NULL
#'   protData <- foreach(cond = seq_len(length(unique(conditions))),
#'                       .combine=cbind,
#'                       .export=c("inner.aggregate.iter", "inner.sum", "inner.mean","GetNbPeptidesUsed"),
#'                       .packages = "Matrix") %dopar% {
#'                         condsIndices <- which(conditions == unique(conditions)[cond])
#'                         qData <- qPepData[,condsIndices]
#'                         inner.aggregate.iter(qData, X, init.method, method, n)
#'                       }
#'
#'   protData <- protData[,colnames(qPepData)]
#'   return(protData)
#' }
#'
#'
#'
#'
#' #' Method to aggregate peptides into proteins with the iterative approach
#' #'
#' #' @title Method to aggregate peptides into proteins with the iterative approach.
#' #'
#' #' @param qPepData A matrix of intensities of peptides.
#' #'
#' #' @param X An adjacency matrix
#' #'
#' #' @param conditions xxx
#' #'
#' #' @param init.method The method used to initialize the iterative algotirhm. Default is 'Sum'.
#' #'
#' #' @param method The method used for xxx. Default value is 'Mean'.
#' #'
#' #' @param n xxxx
#' #'
#' #' @return A matrix of protein intensities.
#' #'
#' #' @author Samuel Wieczorek
#' #'
#' #' @examples
#' #' ft <- create_ft_example(with.na = TRUE)
#' #' X <- QFeatures::adjacencyMatrix(ft[[1]])
#' #' conds <- design.qf(ft)$Condition
#' #' aggIterative(assay(ft,1), X, conditions=conds)
#' #'
#' #' @export
#' #'
#' #' @rdname iterative-aggregation
#' #'
#' aggIterative <- function(qData,
#'                          X,
#'                          conditions = NULL,
#'                          init.method = 'sum',
#'                          iter.method = 'mean',
#'                          n = NULL){
#'   if (is.null(conditions)){
#'     warning("The parameter 'conditions' is NULL: the aggregation cannot be process.")
#'     return(NULL)
#'   }
#'   protData <- matrix(rep(0,ncol(qData)*ncol(X)),
#'                      nrow=ncol(X),
#'                      dimnames=list(colnames(X), rep("cond", ncol(qData))))
#'
#'   for (cond in unique(conditions)){
#'     ind <- which(conditions == cond)
#'     protData[,ind]  <- inner.aggregate.iter(qData[,ind],
#'                                             X,
#'                                             init.method,
#'                                             iter.method,
#'                                             n)
#'     colnames(protData)[ind] <- colnames(qData)[ind]
#'   }
#'
#'   return(protData)
#'
#' }
#'

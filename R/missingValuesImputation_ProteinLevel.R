#' @title Finds the LAPALA 
#' 
#' @description
#' 
#' Methods available are:
#' 
#' * wrapper.impute.detQuant(): 
#' This method is a wrapper of the function `impute.detQuant()` for objects
#' of class \code{MSnSet}
#' 
#' * wrapper.impute.KNN(): Can impute only POV missing values. This method is 
#' a wrapper for objects of class `QFeatures` and imputes missing values with 
#' a fixed value. This function imputes the missing values condition by 
#' condition.
#' 
#' * wrapper.impute.slsa(): 
#' Imputation of peptides having no values in a biological condition. This 
#' method is a wrapper to the function \code{impute.slsa()} of the package
#' \code{imp4p} adapted to an object of class \code{MSnSet}.
#' 
#' * wrapper.impute.fixedValue():
#' This method is a wrapper to objects of class \code{MSnSet} and imputes 
#' missing values with a fixed value.
#' 
#' * wrapper.impute.pa(): 
#' Imputation of peptides having no values in a biological condition.
#' This method is a wrapper to the function \code{impute.pa} of the package
#' \code{imp4p} adapted to an object of class \code{MSnSet}.
#' 
#' # Utilities functions
#' 
#' * findMECBlock(): xxx
#' * reIntroduceMEC(): xxx
#' * getQuantile4Imp(): 
#' Quantile imputation value definition. This method returns the q-th quantile 
#' of each column of an expression set, up to a scaling factor
#' 
#' @param obj An object of class `QFeatures`.
#' @param grp xxx
#' @param MECIndex A data.frame that contains index of MEC (see findMECBlock) 
#' @param q.min Same as the function `impute.pa()` in the package `imp4p`
#' @param K the number of neighbors.
#' @param fixVal A float.
#' @param qval An expression set containing quantitative values of various
#' replicates
#' @param factor A scaling factor to multiply the imputation value with
#' @param na.type A string which indicates the type of missing values to impute.
#' Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.
#'
#' @param na.type A string which indicates the type of missing values to impute.
#' Available values are: `NA` (for both POV and MEC), `POV`, `MEC`.
#' @param qdata xxx
#' @param coldata xxx
#' 
#' 
#' @return A data.frame containing the indexes of LAPALA
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[[1]][seq_len(100)]
#' grp <- get_group(Exp1_R25_pept)
#' lapala <- findMECBlock(obj, grp)
#' na.type = c("Missing POV", "Missing MEC")
#' obj.imp.pov <- wrapper.impute.detQuant(obj, na.type = na.type)
#' obj.imp.pov <- reIntroduceMEC(obj, grp, lapala)
#' 
#' obj.imp.pov <- wrapper.impute.KNN(obj, grp, 3)
#' 
#' obj.imp.pov <- wrapper.impute.fixedValue(obj, grp, 0.001, na.type = "Missing POV")
#' obj.imp.mec <- wrapper.impute.fixedValue(obj, grp, 0.001, na.type = "Missing MEC")
#' obj.imp.na <- wrapper.impute.fixedValue(obj, grp, 0.001, na.type = c("Missing MEC", "Missing POV"))
#'
#' obj.imp.pov <- wrapper.impute.pa(obj, grp)
#' 
#' qdata <- SummarizedExperiment::assay(obj)
#' quant <- getQuantile4Imp(qdata)
#' 
#' coldata <- colData(Exp1_R25_pept)
#' obj.slsa.pov <- wrapper.impute.slsa(obj, grp, coldata)
#'
#'
#'
#'
#' @name mv_imputation_protein
#' 
NULL



#' @export
#' @rdname mv_imputation_protein
#' @import SummarizedExperiment
#'
findMECBlock <- function(obj, grp) {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  conditions <- unique(grp)
  nbCond <- length(conditions)

    s <- data.frame()
    qdata <- SummarizedExperiment::assay(obj)
    
    for (cond in seq_len(nbCond)) {
        ind <- which(grp == conditions[cond])
        lNA <- which(
            apply(is.na(qdata[, ind]), 1, sum) == length(ind))
        if (length(lNA) > 0) {
            tmp <- data.frame(
                cond, 
                which(
                    apply(
                        is.na(qdata[, ind]), 1, sum) == length(ind)))
            names(tmp) <- c("Condition", "Line")
            s <- rbind(s, tmp)
        }
    }
    return(s)
}


#'
#' @export
#' @rdname mv_imputation_protein
#'
reIntroduceMEC <- function(obj, grp, MECIndex) {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  conditions <- unique(grp)
  
  for (.row in seq(nrow(MECIndex)))
    {
    replicates <- which(grp == conditions[MECIndex[.row, "Condition"]])
    SummarizedExperiment::assay(obj)[MECIndex[.row, "Line"], as.vector(replicates)] <- NA
    }
    return(obj)
}




#'
#' @export
#' @rdname mv_imputation_protein
#' @import omXplore
#' @import SummarizedExperiment
#'
wrapper.impute.KNN <- function(obj = NULL, grp, K) {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
    pkgs.require('impute')
    
    if (missing(obj)) {
        stop("'obj' is required.")
    } else if (is.null(obj)) {
        stop("'obj' is NULL")
    }
    
    data <- SummarizedExperiment::assay(obj)

    conditions <- unique(grp)
    nbCond <- length(conditions)

    for (cond in seq_len(nbCond)) {
        ind <- which(grp== conditions[cond])
        resKNN <- impute::impute.knn(
          SummarizedExperiment::assay(obj)[, ind],
          k = K, 
          rowmax = 0.99, 
          colmax = 0.99, 
          maxp = 1500,
          rng.seed = sample(seq_len(1000), 1)
          )

        SummarizedExperiment::assay(obj)[, ind] <- resKNN$data
    }

    SummarizedExperiment::assay(obj)[data == 0] <- NA
    .metacell <- omXplore::get_metacell(obj)
    SummarizedExperiment::assay(obj)[.metacell == 'Missing MEC'] <- NA
    
    # Transform all previously tagged 'na.type' as 'Imputed'
    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}




#' @rdname mv_imputation_protein
#' @export
#' @import omXplore
#' @import SummarizedExperiment
#'
wrapper.impute.fixedValue <- function(obj, 
  grp,
  fixVal = 0, 
  na.type) {
    if (missing(obj))
        stop("'obj' is required.")
  stopifnot(inherits(obj, 'SummarizedExperiment'))
    if (fixVal == 0)
        warning("Be aware that fixVal = 0. No imputation will be realize.")

    level <- omXplore::get_type(obj)

    if (missing(na.type)) {
        stop(paste0("'na.type' is required. Available values are: ", 
            paste0(metacell.def(level)$node, collapse = " ")))
    } else if (!is.subset(na.type, metacell.def(level)$node)) {
        stop(paste0("Available values for na.type are: ", 
            paste0(metacell.def(level)$node, collapse = " ")))
    }

    ind.na.type <- DaparToolshed::match.metacell(
      omXplore::get_metacell(obj),
      na.type,
      level = level
      )
    
    .ind <- is.na(SummarizedExperiment::assay(obj)) & ind.na.type
    SummarizedExperiment::assay(obj)[.ind] <- fixVal
    obj <- UpdateMetacellAfterImputation(obj)
    return(obj)
}


#'
#' @rdname mv_imputation_protein
#' @export
#' @import SummarizedExperiment
#'
wrapper.impute.pa <- function(
    obj = NULL, 
    grp,
    q.min = 0.025) {
  stopifnot(inherits(obj, 'SummarizedExperiment'))
    pkgs.require('imp4p')

    if (is.null(obj)) 
        stop("'obj' is required.")

    res <- imp4p::impute.pa(
      SummarizedExperiment::assay(obj), 
      conditions = as.factor(grp), 
      q.min = q.min,
      q.norm = 3,
      eps = 0)
    
    SummarizedExperiment::assay(obj) <- res[["tab.imp"]]

    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}




#' @export
#' @rdname mv_imputation_protein
#' @import omXplore
#' @import SummarizedExperiment
#'
wrapper.impute.detQuant <- function(
    obj, 
    qval = 0.025, 
    factor = 1, 
    na.type) {
  
    if (missing(obj))
        stop("'obj' is required.")
  stopifnot(inherits(obj, 'SummarizedExperiment'))
    if (missing(na.type)){
      na.type <- c('Missing POV', 'Missing MEC')
    } else {
      if (!is.subset(na.type, c('Missing POV', 'Missing MEC'))) {
      stop("'na.type' is required. Available values are: 'Missing POV', 'Missing MEC'")
      }
    }

    qdata <- SummarizedExperiment::assay(obj)
    values <- getQuantile4Imp(qdata, qval, factor)
    for (iter in seq_len(ncol(qdata))) {
        col <- qdata[, iter]
        ind.na.type <- DaparToolshed::match.metacell(
          omXplore::get_metacell(obj)[,iter],
                                      pattern = na.type,
                                      level = omXplore::get_type(obj)
                                      )

        col[ind.na.type] <- values$shiftedImpVal[iter]
        qdata[, iter] <- col
    }

    SummarizedExperiment::assay(obj) <- qdata
    #msg <- "Missing values imputation using deterministic quantile"
    #metadata(obj)$processing <- c(metadata(obj)$processing, msg)
    #metadata(obj)$processing$imputation.method <- "detQuantile"
    obj <- UpdateMetacellAfterImputation(obj)

    return(obj)
}



#'
#' @return A list of two vectors, respectively containing the imputation values
#' and the rescaled imputation values
#'
#' @rdname mv_imputation_protein
#' 
#' @export
#'
getQuantile4Imp <- function(qdata, qval = 0.025, factor = 1) {
    
    pkgs.require('stats')
    r1 <- apply(qdata, 2, stats::quantile, qval, na.rm = TRUE)
    r2 <- r1 * factor
    return(list(ImpVal = r1, shiftedImpVal = r2))
}





#' @export
#' 
#' @rdname mv_imputation_protein
#' @import SummarizedExperiment
#'
wrapper.impute.slsa <- function(
    obj = NULL,
    grp,
    coldata) {
    
    pkgs.require('imp4p')
    
    if (is.null(obj))
        stop("'obj' is required.")
  stopifnot(inherits(obj, 'SummarizedExperiment'))
  
    MECIndex <- findMECBlock(obj, grp)

    # sort conditions to be compliant with impute.slsa
    conds <- factor(grp, levels = unique(grp))
    sample.names.old <- coldata[, 'quantCols']
    sTab <- as.data.frame(coldata)
    new.order <- unlist(
      lapply(split(sTab, conds), 
        function(x) {x["quantCols"]})
      )
    qdata <- SummarizedExperiment::assay(obj)[, new.order]
    res <- imp4p::impute.slsa(qdata,
      conditions = conds,
      nknn = 10,
      selec = "all",
      weight = 1,
      ind.comp = 1
    )

    # restore old order
    res <- res[, sample.names.old]

    SummarizedExperiment::assay(obj) <- res
    obj <- reIntroduceMEC(obj, grp, MECIndex)

    obj <- UpdateMetacellAfterImputation(obj)
    return(obj)
}

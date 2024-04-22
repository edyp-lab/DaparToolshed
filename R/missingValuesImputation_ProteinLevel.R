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
#' @param i xxx
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
#'
#' @return A data.frame containing the indexes of LAPALA
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(100)]
#' lapala <- findMECBlock(obj, 1)
#' na.type = c("Missing POV", "Missing MEC")
#' obj <- wrapper.impute.detQuant(obj, 1, na.type = na.type)
#' obj <- reIntroduceMEC(obj, 1, lapala)
#' 
#' obj <- Exp1_R25_pept[seq_len(10), ]
#' obj.imp.pov <- wrapper.impute.KNN(obj, 1, 3)
#' 
#' obj.imp.pov <- wrapper.impute.fixedValue(obj, 1, 0.001, na.type = "Missing POV")
#' obj.imp.mec <- wrapper.impute.fixedValue(obj, 1, 0.001, na.type = "Missing MEC")
#' obj.imp.na <- wrapper.impute.fixedValue(obj, 1, 0.001, na.type = c("Missing MEC", "Missing POV"))
#'
#' obj.imp.pov <- wrapper.impute.pa(obj, 1)
#' 
#' qdata <- SummarizedExperiment::assay(obj[[1]])
#' quant <- getQuantile4Imp(qdata)
#' 
#' obj <- Exp1_R25_pept[seq_len(100)]
#' obj.slsa.pov <- wrapper.impute.slsa(obj, 1)
#'
#'
#'
#'
#' @name mv_imputation_protein
#' 
NULL



#' @export
#' @rdname mv_imputation_protein
#'
findMECBlock <- function(obj, i) {
    groups <- unique(get_group(obj))
    nbCond <- length(groups)

    s <- data.frame()
    qdata <- SummarizedExperiment::assay(x = obj, i = i)
    for (cond in seq_len(nbCond)) {
        ind <- which(get_group(obj) == groups[cond])
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
reIntroduceMEC <- function(obj, i, MECIndex) {
  .cond <- get_group(obj)
  conditions <- unique(.cond)
  
  for (.row in seq(nrow(MECIndex)))
    {
    replicates <- which(.cond == conditions[MECIndex[.row, "Condition"]])
    SummarizedExperiment::assay(obj[[i]])[MECIndex[.row, "Line"], as.vector(replicates)] <- NA
    }
    return(obj)
}




#'
#' @export
#' @rdname mv_imputation_protein
#'
wrapper.impute.KNN <- function(obj = NULL, i, K) {

    pkgs.require('impute')
    
    stopifnot(inherits(obj, 'QFeatures'))
    if (missing(obj)) {
        stop("'obj' is required.")
    } else if (is.null(obj)) {
        stop("'obj' is NULL")
    }


    data <- SummarizedExperiment::assay(obj[[i]])

    conditions <- unique(get_group(obj))
    nbCond <- length(conditions)

    for (cond in seq_len(nbCond)) {
        ind <- which(get_group(obj) == conditions[cond])
        resKNN <- impute::impute.knn(
          SummarizedExperiment::assay(obj[[i]])[, ind],
          k = K, 
          rowmax = 0.99, 
          colmax = 0.99, 
          maxp = 1500,
          rng.seed = sample(seq_len(1000), 1)
          )
        SummarizedExperiment::assay(obj[[i]])[, ind] <- resKNN[[1]]
    }

    SummarizedExperiment::assay(obj[[i]])[data == 0] <- NA
    .metacell <- get_metacell(obj[[i]])
    SummarizedExperiment::assay(obj[[i]])[.metacell == 'Missing MEC'] <- NA
    
    # Transform all previously tagged 'na.type' as 'Imputed'
    obj[[i]] <- UpdateMetacellAfterImputation(obj[[i]])

    return(obj)
}




#' @rdname mv_imputation_protein
#' @export
#'
wrapper.impute.fixedValue <- function(obj, 
  i,
  fixVal = 0, 
  na.type) {
    if (missing(obj))
        stop("'obj' is required.")

    if (fixVal == 0)
        warning("Be aware that fixVal = 0. No imputation will be realize.")

    level <- get_type(obj[[i]])

    if (missing(na.type)) {
        stop(paste0("'na.type' is required. Available values are: ", 
            paste0(metacell.def(level)$node, collapse = " ")))
    } else if (!is.subset(na.type, metacell.def(level)$node)) {
        stop(paste0("Available values for na.type are: ", 
            paste0(metacell.def(level)$node, collapse = " ")))
    }

    ind.na.type <- match.metacell(
      get_metacell(obj[[i]]),
      na.type,
      level = level
      )
    
    .ind <- is.na(SummarizedExperiment::assay(obj[[i]])) & ind.na.type
    SummarizedExperiment::assay(obj[[i]])[.ind] <- fixVal
    obj[[i]] <- UpdateMetacellAfterImputation(obj[[i]])
    return(obj)
}


#'
#' @rdname mv_imputation_protein
#' @export
#'
wrapper.impute.pa <- function(
    obj = NULL, 
    i,
    q.min = 0.025) {
    
    pkgs.require('imp4p')

    if (is.null(obj)) 
        stop("'obj' is required.")

    res <- imp4p::impute.pa(
      SummarizedExperiment::assay(obj[[i]]), 
      conditions = as.factor(get_group(obj)), 
      q.min = q.min,
      q.norm = 3,
      eps = 0)
    
    SummarizedExperiment::assay(obj[[i]]) <- res[["tab.imp"]]

    obj[[i]] <- UpdateMetacellAfterImputation(obj[[i]])

    return(obj)
}




#' @export
#' @rdname mv_imputation_protein
#'
wrapper.impute.detQuant <- function(
    obj, 
    i,
    qval = 0.025, 
    factor = 1, 
    na.type) {
  
    if (missing(obj))
        stop("'obj' is required.")
    
    if (missing(na.type)){
      na.type <- c('Missing POV', 'Missing MEC')
    } else {
      if (!is.subset(na.type, c('Missing POV', 'Missing MEC'))) {
      stop("'na.type' is required. Available values are: 'Missing POV', 'Missing MEC'")
      }
    }

    qdata <- SummarizedExperiment::assay(obj[[i]])
    values <- getQuantile4Imp(qdata, qval, factor)
    for (iter in seq_len(ncol(qdata))) {
        col <- qdata[, iter]
        ind.na.type <- match.metacell(get_metacell(obj[[i]])[,iter],
                                      pattern = na.type,
                                      level = get_type(obj[[i]])
                                      )

        col[ind.na.type] <- values$shiftedImpVal[iter]
        qdata[, iter] <- col
    }

    SummarizedExperiment::assay(obj[[i]]) <- qdata
    msg <- "Missing values imputation using deterministic quantile"
    metadata(obj[[i]])$processing <- c(metadata(obj[[i]])$processing, msg)

    metadata(obj[[i]])$processing$imputation.method <- "detQuantile"
    obj[[i]] <- UpdateMetacellAfterImputation(obj[[i]])

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
#'
wrapper.impute.slsa <- function(
    obj = NULL, 
  i = 1) {
    
    pkgs.require('imp4p')
    
    if (is.null(obj))
        stop("'obj' is required.")

    MECIndex <- findMECBlock(obj, i)

    # sort conditions to be compliant with impute.slsa
    conds <- factor(get_group(obj), levels = unique(get_group(obj)))
    sample.names.old <- colData(obj)$Sample.name
    sTab <- colData(obj)
    new.order <- unname(unlist(
      lapply(split(as.data.frame(sTab), conds), function(x) {x["Sample.name"]})
      ))
    qdata <- SummarizedExperiment::assay(obj[[i]])[, new.order]

    res <- imp4p::impute.slsa(qdata,
                              conditions = conds,
                              nknn = 10,
                              selec = "all",
                              weight = 1,
                              ind.comp = 1
                              )

    # restore old order
    res <- res[, sample.names.old]

    SummarizedExperiment::assay(obj[[i]]) <- res
    obj <- reIntroduceMEC(obj, i, MECIndex)

    obj[[i]] <- UpdateMetacellAfterImputation(obj[[i]])
    return(obj)
}

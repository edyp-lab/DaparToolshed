#' @title Delete the lines in the matrix of intensities and the metadata table
#' given their indices.
#'
#' @param obj An object of class `SummarizedExperiment` containing
#' quantitative data.
#' 
#' @param conds xxx
#'
#' @param level A vector of integers which are the indices of lines to
#' delete.
#'
#' @param pattern A string to be included in the `SummarizedExperiment`
#' object for log.
#'
#' @param type xxx
#'
#' @param percent xxx
#'
#' @param op xxx
#'
#' @param th xxx
#'
#' @return An instance of class `SummarizedExperiment` that have been filtered.
#'
#' @author Samuel Wieczorek
#'
#' @examplesIf interactive()
#' NA
#'
#' @export
#'
GetIndices_FunFiltering <- function(obj,
  conds, 
  level, 
  pattern = NULL,
  type = NULL,
  percent, 
  op, 
  th) {
    if (missing(obj))
        stop("'obj' is required.")
    
    if (missing(level))
        stop("'level' is required.")

    if (missing(pattern))
        stop("'pattern' is required.")

    if (missing(type))
        stop("'type' is required.")

    if (missing(percent))
        stop("'percent' is required.")

    if (missing(op))
        stop("'op' is required.")

    if (missing(th))
        stop("'th' is required.")



    indices <- NULL
    is.subset <- sum(pattern == intersect(pattern,  metacell.def(level)$node))==length(pattern)
    if (!is.subset) {
        warning("Available values for pattern are: ", paste0(metacell.def(level)$node, collapse=', ' ))
        return(NULL)
    }

    mask <- match.metacell(metadata = get_metacell(obj),
                           pattern = pattern,
                           level = level
                           )

    indices <- switch(type,
        WholeLine = GetIndices_WholeLine(metacell.mask = mask),
        WholeMatrix = GetIndices_WholeMatrix(
            metacell.mask = mask,
            op = op,
            percent = percent,
            th = th
        ),
        AllCond = GetIndices_BasedOnConditions(
            metacell.mask = mask,
            type = type,
            conds = conds,
            percent = percent,
            op = op,
            th = th
        ),
        AtLeastOneCond = GetIndices_BasedOnConditions(
            metacell.mask = mask,
            type = type,
            conds = conds,
            percent = percent,
            op = op,
            th = th
        )
    )

    return(indices)
}



#' @title
#' Lists the metacell scopes for filtering
#'
#' @export
#' 
#' @return xxx
#' 
#' @examples 
#' MetacellFilteringScope()
#'
MetacellFilteringScope <- function() {
    c("None", "WholeLine", "WholeMatrix", "AllCond", "AtLeastOneCond")
}



#' @title xxx
#'
#' @export
#' 
#' @return A `character()`
#' 
#' @examples 
#' SymFilteringOperators()
#'
SymFilteringOperators <- function() {
    c("<=", "<", ">=", ">", "==", "!=")
}


#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either all 
#' conditions or at least one condition.
#'
#' @param metacell.mask xxx
#'
#' @param op  String for operator to use. List of operators is available with 
#' SymFilteringOperators().
#'
#' @param percent A boolean to indicate whether the threshold represent an 
#' absolute value (percent = FALSE) or
#' a percentage (percent=TRUE).
#'
#' @param th A floating number which is in the interval [0, 1]
#'
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- 'peptide'
#' pattern <- "Missing"
#' metacell.mask <- match.metacell(metadata = GetMetacell(obj), 
#' pattern = pattern, level = level)
#' percent <- FALSE
#' th <- 3
#' op <- ">="
#' ind <- GetIndices_WholeMatrix(metacell.mask, op, percent, th)
#'
#' @export
#' 
#' @return xxx
#'
GetIndices_WholeMatrix <- function(metacell.mask,
    op = "==",
    percent = FALSE,
    th = 0) {

    # Check parameters
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is required.")
    }
    if (isTRUE(percent)) {
        if (th < 0 || th > 1) {
            warning("With percent=TRUE, the threshold 'th' must be in the 
                interval [0, 1].")
            return(NULL)
        }
    } else {
        th.upbound <- ncol(metacell.mask)
        if (th > th.upbound) {
            warn.txt <- paste0(
                "Param `th` is not correct. It must be an integer greater 
                than or equal to 0 and less or equal than ",
                th.upbound
            )
            warning(warn.txt)
            return(NULL)
        }
    }

    if (!(op %in% SymFilteringOperators())) {
        warning(paste0(
            "'op' must be one of the following values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
        return(NULL)
    }

    indices <- NULL
    if (isTRUE(percent)) {
        inter <- rowSums(metacell.mask) / ncol(metacell.mask)
        indices <- which(eval(parse(text = paste0("inter", op, th))))
    } else {
        inter <- apply(metacell.mask, 1, sum)
        indices <- which(eval(parse(text = paste0("inter", op, th))))
    }


    if (length(indices) == 0) indices <- NULL
    return(indices)
}


#' @title
#' Search lines which respects query on all their elements.
#'
#' @description
#' This function looks for the lines where each element respect the query.
#'
#' @param metacell.mask xxx
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq.int(from=20, to=30)]
#' level <- 'peptide'
#' pattern <- "Missing POV"
#' metacell.mask <- match.metacell(metadata = GetMetacell(obj), 
#' pattern = pattern, level = level)
#' ind <- GetIndices_WholeLine(metacell.mask)
#'
#' @export
#' 
#' @return xxx
#'
GetIndices_WholeLine <- function(metacell.mask) {
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is missing.")
    }

    indices <- which(rowSums(metacell.mask) == ncol(metacell.mask))
    if (length(indices) == 0) indices <- NULL
    return(indices)
}


#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either 
#' all conditions
#' or at least one condition.
#'
#' @param metacell.mask xxx
#'
#' @param type Available values are:
#' * 'AllCond' (the query is valid in all the conditions),
#' * 'AtLeatOneCond' (the query is valid in at leat one condition.
#'
#' @param conds xxx
#'
#' @param percent xxx
#'
#' @param op  String for operator to use. List of operators is available 
#' with SymFilteringOperators().
#'
#' @param th The theshold to apply
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' level <- GetTypeofData(obj)
#' pattern <- 'Missing'
#' metacell.mask <- match.metacell(metadata=GetMetacell(obj), 
#' pattern=pattern, level=level)
#' type <- 'AllCond'
#' conds <- Biobase::pData(obj)$Condition
#' op <- '>='
#' th <- 0.5
#' percent <- TRUE
#' ind <- GetIndices_BasedOnConditions(metacell.mask, type, conds, 
#' percent, op, th)
#'
#' @return xxx
#'
#' @export
#'
GetIndices_BasedOnConditions <- function(metacell.mask,
    type,
    conds,
    percent,
    op,
    th) {

    # Check parameters
    if (missing(metacell.mask)) {
        stop("'metacell.mask' is missing.")
    }
    if (missing(conds)) {
        stop("'conds' is missing.")
    }
    if (missing(type)) {
        stop("'type' is missing.")
    } else if (!(type %in% c("AllCond", "AtLeastOneCond"))) {
        stop("'type' must be one of the following: AllCond, AtLeastOneCond.")
    }
    if (missing(percent)) {
        stop("'percent' is missing.")
    }
    if (missing(op)) {
        stop("'op' is missing.")
    }
    if (missing(th)) {
        stop("'th' is missing.")
    } else if (!(op %in% SymFilteringOperators())) {
        stop(paste0(
            "'op' must be one of the following values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
    }

    u_conds <- unique(conds)
    nbCond <- length(u_conds)

    if (isTRUE(percent)) {
        if (th < 0 || th > 1) {
            warning("With percent=TRUE, the threshold 'th' must be in the 
                interval [0, 1].")
            return(NULL)
        }
    } else {
        th.upbound <- min(unlist(lapply(u_conds, 
            function(x) length(which(conds == x)))))
        if (th > th.upbound) {
            warn.txt <- paste0(
                "Param `th` is not correct. It must be an integer greater 
                than or equal to 0 and less or equal than ",
                th.upbound
            )
            warning(warn.txt)
            return(NULL)
        }
    }

    indices <- NULL
    s <- matrix(rep(0, nrow(metacell.mask) * nbCond),
        nrow = nrow(metacell.mask),
        ncol = nbCond
    )

    for (c in seq_len(nbCond)) {
        ind.cond <- which(conds == u_conds[c])
        inter <- rowSums(metacell.mask[, ind.cond])
        if (isTRUE(percent)) {
            inter <- inter / length(ind.cond)
        }

        s[, c] <- eval(parse(text = paste0("inter", op, th)))
    }

    indices <- switch(type,
        AllCond = which(rowSums(s) == nbCond),
        AtLeastOneCond = which(rowSums(s) >= 1)
    )

    return(indices)
}

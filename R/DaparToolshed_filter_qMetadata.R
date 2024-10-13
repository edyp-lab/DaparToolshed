#' @title
#' Search lines which respects request on one or more conditions.
#'
#' @description
#' This function looks for the lines that respect the request in either all 
#' conditions or at least one condition.
#'
#' @name qMetacell-filter
#' 
#' @return NA
#'
#' @examples
#' data(ft, package='DaparToolshed')
#' obj <- ft[[1]]
#' level <- typeOfData(ft, 1)
#' pattern <- "missing"
#' mask <- match.metacell(
#'     metadata = qMetacell(obj),
#'     pattern = pattern,
#'     level = level
#' )
#' percent <- FALSE
#' th <- 3
#' op <- ">="
#' ind <- qMetacellWholeMatrix(mask, op, percent, th)
#'
#' data(ft, package='DaparToolshed')
#' obj <- ft[[1]]
#' mask <- match.metacell(
#'     metadata = qMetacell(obj),
#'     pattern = "missing POV",
#'     level = typeOfData(obj)
#' )
#' ind <- qMetacellWholeLine(mask)
#'
#' #' data(ft, package='DaparToolshed')
#' obj <- ft[[1]]
#' mask <- match.metacell(
#'     metadata = qMetacell(obj),
#'     pattern = "missing",
#'     level = typeOfData(obj)
#' )
#' type <- "AllCond"
#' conds <- design.qf(ft)$Condition
#' op <- ">="
#' th <- 2
#' percent <- FALSE
#' ind <- qMetacellOnConditions(mask, type, conds, percent, op, th)
#'
NULL


#' @export
#'
#' @rdname qMetacell-filter
#'
#' @return A vector of filtering scopes
#' 
#' @examples
#' qMetacellFilteringScope()
#'
qMetacellFilteringScope <- function() {
    c("None", "WholeLine", "WholeMatrix", "AllCond",  "AtLeastOneCond")
}



#' @export
#'
#' @rdname qMetacell-filter
#'
#' @return A vector of operators
#' 
#' @examples 
#' SymFilteringOperators()
#'
SymFilteringOperators <- function() {
    c( "<=", "<", ">=", ">", "==", "!=")
}




#' @param object An instance of the class `SummarizedExperiment`
#' @param cmd A `character(1)` xxx
#' @param pattern A `character(1)` xxx
#' @param percent A boolean to indicate whether the threshold represent an 
#' absolute value (percent = FALSE) or a percentage (percent=TRUE).
#' @param th A floating number which is in the interval [0, 1]
#' @param operator String for operator to use. List of operators is available 
#' with `SymFilteringOperators()`.
#'
#'
#' @return NA
#'
#' @export
#'
#' @rdname qMetacell-filter
#' 
#' @examples 
#' data(ft_na)
#' obj <- ft_na[[1]]
#'
qMetacellWholeMatrix <- function(object, 
  cmd, 
  pattern, 
  percent = "Percentage", 
  th, 
  operator) {
    if (missing(object)) {
        return(NULL)
    } else
      stopifnot(inherits(object, 'SummarizedExperiment'))
  
    if (missing(cmd) || missing(pattern) || 
        missing(percent) || missing(th) ||
        missing(operator)) {
        return(object)
    }


    stopifnot(inherits(object, "SummarizedExperiment"))
    stopifnot("qMetacell" %in% names(rowData(object)))

    stopifnot(!is.null(cmd) || (cmd %in% c("keep", "delete")))

    indices <- NULL
    level <- typeDataset(object)


    # Check parameters
    mask <- match.metacell(
      metadata = qMetacell(object),
        pattern = pattern,
        level = level
    )
    if (percent == "Percentage") {
        if (th < 0 || th > 1) {
            stop("With percent=TRUE, the threshold 'th' must be in the 
                interval [0, 1].")
        }
    } else {
        if (th > ncol(mask)) {
            stop(paste0(
                "Param `th` is not correct. It must be an integer greater 
                than or equal to 0 and less or equal than ",
                ncol(mask)
            ))
        }
    }

    if (!(operator %in% SymFilteringOperators())) {
        stop(paste0(
            "'op' must be one of the followinf values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
    }



    indices <- NULL
    if (percent == "Percentage") {
        inter <- rowSums(mask) / ncol(mask)
    } else {
        inter <- apply(mask, 1, sum)
    }

    indices <- which(eval(parse(text = paste0("inter", operator, th))))
    if (length(indices) > 0) {
        object <- switch(cmd,
            keep = object[indices, ],
            delete = object[-indices, ]
        )
    }

    return(object)
}



#' @param object xxx
#' @param cmd xxx
#' @param pattern xxx
#'
#' @return NA
#' @export
#'
#' @rdname qMetacell-filter
#' 
#' @import omXplore
#'
qMetacellWholeLine <- function(object, cmd, pattern) {
    if (missing(object)) {
        return(NULL)
    }
    if (missing(cmd)) {
        return(object)
    }
    if (missing(pattern)) {
        return(object)
    }


    stopifnot(inherits(object, "SummarizedExperiment"))
    stopifnot("qMetacell" %in% names(rowData(object)))
    stopifnot(!is.null(cmd) || (cmd %in% c("keep", "delete")))

    indices <- NULL
    level <- typeDataset(object)

    if (!all(pattern %in% omXplore::metacell.def(level)$node)) {
        warning("Either 'pattern' nor 'type' are equal to 'None'")
        return(NULL)
    }

    mask <- match.metacell(
        metadata = qMetacell(object),
        pattern = pattern,
        level = level
    )


    indices <- which(rowSums(mask) == ncol(mask))
    if (length(indices) > 0) {
        object <- switch(cmd,
            keep = object[indices, ],
            delete = object[-indices, ]
        )
    }

    return(object)
}



#' @param object xxx
#'
#' @param cmd Available values are:
#'
#' @param mode xxx
#' @param pattern ss
#' @param conds xxx
#'
#' @param percent xxx
#'
#' @param operator  String for operator to use. List of operators is available 
#' with `SymFilteringOperators()`.
#'
#' @param th The threshold to apply
#'
#' @return NA
#' @export
#'
#' @rdname qMetacell-filter
#'
qMetacellOnConditions <- function(object,
    cmd,
    mode,
    pattern,
    conds,
    percent = "Percentage",
    operator,
    th) {

    # Check parameters
    if (missing(pattern)) {
        stop("'pattern' is missing.")
    }
    if (missing(conds)) {
        stop("'conds' is missing.")
    }
    if (missing(mode)) {
        stop("'type' is missing.")
    } else if (!(mode %in% c("AllCond", "AtLeastOneCond"))) {
        stop("'mode' must be one of the following: AllCond, AtLeastOneCond.")
    }
    if (missing(percent)) {
        stop("'percent' is missing.")
    }
    if (missing(th)) {
        stop("'th' is missing.")
    } else if (percent == "Percentage"){
        stopifnot(th >= 0 && th <= 1)
    }

    if (missing(operator)) {
        stop("'operator' is missing.")
    } else if (!(operator %in% SymFilteringOperators())) {
        stop(paste0(
            "'operator' must be one of the followinf values: ",
            paste0(SymFilteringOperators(), collapse = " ")
        ))
    }
    u_conds <- unique(conds)
    nbCond <- length(u_conds)

    if (percent == "Count") {
        th.upbound <- min(unlist(lapply(
            u_conds,
            function(x) length(which(conds == x))
        )))
        if (th > th.upbound) {
            stop(paste0(
                "Param `th` is not correct. It must be an integer greater than 
                or equal to 0 and less or equal than ",
                th.upbound
            ))
        }
    }


    mask <- match.metacell(
      metadata = qMetacell(object),
        pattern = pattern,
        level = typeDataset(object)
    )

    indices <- NULL
    temp <- matrix(rep(NA, nrow(mask) * nbCond),
        nrow = nrow(mask),
        ncol = nbCond
    )

    for (c in seq_len(nbCond)) {
        ind.cond <- which(conds == u_conds[c])
        inter <- rowSums(mask[, ind.cond])
        if (percent == "Percentage") {
            inter <- inter / length(ind.cond)
        }
        temp[, c] <- eval(parse(text = paste0("inter", operator, th)))
    }

    indices <- switch(mode,
        AllCond = which(rowSums(temp) == nbCond),
        AtLeastOneCond = which(rowSums(temp) >= 1)
    )

    if (length(indices) > 0) {
        object <- switch(cmd,
            keep = object[indices, ],
            delete = object[-indices, ]
        )
    }

    return(object)
}

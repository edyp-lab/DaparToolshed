#' @title Check if xxxxxx
#'
#' @param tab A data.frame which correspond to xxxxxx
#'
#' @return A list of two items
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' test.design(design.qf(Exp1_R25_pept)[, seq_len(3)])
#'
#' @export
#'
test.design <- function(tab) {
  valid <- TRUE
  txt <- NULL
  level <- NULL
  
  level.a <- factor(tab[, 1], ordered = TRUE)
  level.b <- factor(tab[, 2], ordered = TRUE)
  name.level.a <- colnames(tab)[1]
  name.level.b <- colnames(tab)[2]
  
  level.c <- NULL
  if (ncol(tab) == 3) {
    level.c <- factor(tab[, 3], ordered = TRUE)
    name.level.c <- colnames(tab)[3]
  }
  
  # verification intersection sur B
  # verification de la non redondance'intersection
  # vide entre les groupes
  uniqueA <- unique(level.a)
  ll <- lapply(
    uniqueA,
    function(x) {
      as.character(level.b)[which(level.a == x)]
    }
  )
  n <- NULL
  for (i in seq_len(length(uniqueA) - 1)) {
    for (j in seq.int(from=(i + 1), to = length(uniqueA))) {
      n <- c(n, intersect(ll[[i]], ll[[j]]))
    }
  }
  
  if (length(n) > 0) {
    valid <- FALSE
    txt <- c(txt, paste0(
      "The value ",
      n,
      " in column '",
      colnames(tab)[2],
      "' is not correctly set.\n"
    ))
  }
  
  
  # verification si niveau hierarchique inf
  if (length(levels(level.a)) == length(levels(level.b))) {
    ## c'est un design de niveau n-1 en fait
    valid <- FALSE
    txt <- c(
      txt,
      paste0(
        "The column ",
        name.level.b,
        " is not informative. ",
        "Thus, the design is not of level (n-1).\n"
      )
    )
  } else if (!is.null(level.c)) {
    if (length(levels(level.b)) == length(levels(level.c))) {
      ## c'est un design de niveau n-1 en fait
      valid <- FALSE
      txt <- c(
        txt,
        paste0(
          "The column ",
          name.level.c,
          " is not informative. ",
          "Thus, the design is of level (n-1).\n"
        )
      )
    }
  }
  
  # verification si niveau non informatif
  return(list(
    valid = valid,
    warn = txt
  ))
}


#' @title Check if the design is valid
#'
#' @param conds A vector
#'
#' @return A list
#'
#' @author Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' check.conditions(design.qf(Exp1_R25_pept)$Condition)
#'
#' @export
#'
check.conditions <- function(conds) {
  res <- list(valid = TRUE, warn = NULL)
  
  if (("" %in% conds) || (NA %in% conds)) {
    res <- list(valid = FALSE, 
                warn = "The conditions are note full filled.")
    return(res)
  }
  
  # Check if there is at least two conditions
  if (length(unique(conds)) < 2) {
    res <- list(valid = FALSE, 
                warn = "The design must contain at least two conditions.")
    return(res)
  }
  
  
  # check if each condition has at least two values
  nValPerCond <- unlist(lapply(unique(conds), function(x) {
    length(conds[which(conds == x)])
  }))
  if (all(nValPerCond < 2)) {
    res <- list(valid = FALSE, 
                warn = "The design must contain at least two values per condition.")
    return(res)
  }
  
  return(res)
}

#' @title Check if the design is valid
#'
#' @param sTab The data.frame which correspond to the `pData()` function 
#' of package `MSnbase`.
#'
#' @return A boolean
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' check.design(design.qf(Exp1_R25_pept)[, seq_len(3)])
#'
#' @export
#'
check.design <- function(sTab) {
  res <- list(valid = FALSE, warn = NULL)
  
  names <- colnames(sTab)
  level.design <- ncol(sTab) - 2
  
  
  res <- check.conditions(sTab$Condition)
  if (!res$valid) {
    return(res)
  }
  # Check if all the column are fullfilled
  
  if (level.design == 1) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    }
  } else if (level.design == 2) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Tech.Rep colmumn are not full filled.")
      return(res)
    }
  } else if (level.design == 3) {
    if (("" %in% sTab$Bio.Rep) || (NA %in% sTab$Bio.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Bio.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Tech.Rep) || (NA %in% sTab$Tech.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Tech.Rep colmumn are not full filled.")
      return(res)
    } else if (("" %in% sTab$Analyt.Rep) || (NA %in% sTab$Analyt.Rep)) {
      res <- list(valid = FALSE, 
                  warn = "The Analyt.Rep colmumn are not full filled.")
      return(res)
    }
  }
  
  # Check if the hierarchy of the design is correct
  if (level.design == 1) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep")])
  } else if (level.design == 2) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep", "Tech.Rep")])
  } else if (level.design == 3) {
    res <- test.design(sTab[, c("Condition", "Bio.Rep", "Tech.Rep")])
    if (res$valid) {
      res <- test.design(sTab[, c("Bio.Rep", "Tech.Rep", "Analyt.Rep")])
    }
  }
  
  return(res)
}

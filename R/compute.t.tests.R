
#' @title xxxxxx
#'
#' @param obj A matrix of quantitative data, without any missing values.
#' @param i xxx
#' @param contrast Indicates if the test consists of the comparison of each
#' biological condition versus
#' each of the other ones (contrast=1;
#' for example H0:"C1=C2" vs H1:"C1!=C2", etc.)
#' or each condition versus all others (contrast=2; e.g.  H0:"C1=(C2+C3)/2" vs
#' H1:"C1!=(C2+C3)/2", etc. if there are three conditions).
#'
#' @param type xxxxx
#'
#' @return A list of two items : logFC and P_Value; both are dataframe. The
#' first one contains the logFC values of all the comparisons (one column for
#' one comparison), the second one contains the pvalue of all the comparisons
#' (one column for one comparison). The names of the columns for those two
#' dataframes are identical and correspond to the description of the comparison.
#'
#' @author Florence Combes, Samuel Wieczorek
#'
#' @examples
#' data(Exp1_R25_prot, package = "DaparToolshedData")
#' obj <- Exp1_R25_prot[seq_len(1000)]
#' level <- 'protein'
#' metacell.mask <- match.metacell(omXplore::get_metacell(obj[[1]]), 
#' c("Missing POV", "Missing MEC"), level)
#' # Simulate imputation
#' obj <- NAIsZero(obj, 1)
#' obj <- NAIsZero(obj, 2)
#' ttest <- compute_t_tests(obj, 2)
#' 
#' ttest <- compute_t_tests(obj, 1)
#'
#' @import omXplore
#'
#' @export
#'
compute_t_tests <- function(obj, 
  i = 1, 
  contrast = c("OnevsOne", "OnevsAll"),
  type = c("Student", "Welch")
  ){
  
  pkgs.require('stats')
  
  if(missing(contrast))
    contrast <- match.arg(contrast)[1]
  else
    contrast <- match.arg(contrast)
  
  
  if(missing(type))
    type <- match.arg(type)[1]
  else
    type <- match.arg(type)
  
  
  switch(type,
         Student = .type <- TRUE,
         Welch = .type <- FALSE
  )
  
  
  qData <- as.matrix(SummarizedExperiment::assay(obj, i))
  sTab <- as.data.frame(MultiAssayExperiment::colData(obj))
  res <- list()
  logFC <- list()
  P_Value <- list()
  nbComp <- NULL
  Conditions <- omXplore::get_group(obj)
  
  Conditionsf <- factor(Conditions, levels = unique(Conditions))
  Cond.Nb <- length(levels(Conditionsf))
  
  
  if (contrast == "OnevsOne") {
    res.tmp <- NULL
    nbComp <- Cond.Nb * (Cond.Nb - 1) / 2
    
    for (i in seq_len(Cond.Nb - 1)) {
      for (j in seq.int(from = (i + 1), to = Cond.Nb)) {
        
        c1Indice <- which(Conditions == levels(Conditionsf)[i])
        c2Indice <- which(Conditions == levels(Conditionsf)[j])
        res.tmp <- apply(
          qData[, c(c1Indice, c2Indice)], 1,
          function(x) {
            stats::t.test(x ~ Conditionsf[c(c1Indice, c2Indice)], var.equal = .type)
          }
        )
        p.tmp <- unlist(lapply(res.tmp, function(x) x$p.value))
        m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[1])))
        m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[2])))
        m1.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[1])))[1]
        m2.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[2])))[1]
        logFC.tmp <- m1.tmp - m2.tmp
        

        txt <- paste(levels(Conditionsf)[i], "_vs_", 
                     levels(Conditionsf)[j], sep = "")
        
        logFC[[paste(txt, "logFC", sep = "_")]] <- logFC.tmp
        P_Value[[paste(txt, "pval", sep = "_")]] <- p.tmp
      }
    }
  } ## end Contrast==1
  
  if (contrast == "OnevsAll") {
    nbComp <- Cond.Nb
    
    for (i in seq(nbComp)) {
      c1 <- which(Conditions == levels(Conditionsf)[i])
      
      Cond.t.all <- seq_len(length(Conditions))
      Cond.t.all[c1] <- levels(Conditionsf)[i]
      Cond.t.all[-c1] <- "all"
      
      res.tmp <- apply(qData, 1, function(x) {
        stats::t.test(x ~ Cond.t.all, var.equal = .type)
      }
      )
      
      p.tmp <- unlist(lapply(res.tmp, function(x) x$p.value))
      m1.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[1])))
      m2.tmp <- unlist(lapply(res.tmp, function(x) as.numeric(x$estimate[2])))
      m1.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[1])))[1]
      m2.name <- names(unlist(lapply(res.tmp, function(x) x$estimate[2])))[1]
      logFC.tmp <- m1.tmp - m2.tmp

      txt <- paste(levels(Conditionsf)[i], "_vs_(all-", 
                   levels(Conditionsf)[i], ")", sep = "")
      
      logFC[[paste(txt, "logFC", sep = "_")]] <- logFC.tmp
      P_Value[[paste(txt, "pval", sep = "_")]] <- p.tmp
    }
  } # End Contrast=2
  
  
  res.l <- list(
    logFC = as.data.frame(logFC),
    P_Value = as.data.frame(P_Value)
  )
  colnames(res.l$logFC) <- names(logFC)
  colnames(res.l$P_Value) <- names(P_Value)
  
  return(res.l)
}

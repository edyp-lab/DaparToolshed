#' @title Aggregate an assay's quantitative features
#'
#' @description
#' This function aggregates the quantitative features of an assay,
#' applying a summarization function (`fun`) to sets of features.
#' The `fcol` variable name points to a rowData column that defines
#' how to group the features during aggregate. This variable has to 
#' be an adjacency matrix. This function uses [QFeatures::aggregateFeatures()]
#' to aggregate quantitative data.
#'
#' The list of agregation methods can be obtained with the function
#' [aggregateMethods()]. This function compiles both methods from the
#' packages `DaparToolshed` and `QFeatures`.
#'
#' @param object An instance of class `QFeatures` or `SummarizedExperiment`
#' @param i The index or name of the assay which features will be aggregated the create the new assay.
#' @param fcol A `character(1)` naming a rowdata variable (of assay `i` in case of a `QFeatures`) 
#'             defining how to aggregate the features of the assay. 
#'             This variable is a (possibly sparse) matrix. See below for details.
#' @param name A `character(1)` naming the new assay. Default is `newAssay`. 
#'             Note that the function will fail if there's already an assay with `name`.
#' @param fun A function used for quantitative feature aggregation. 
#'            See details for examples.
#' @param shared A `boolean` indication if shared peptides should be considered. If `TRUE`, shared peptides 
#' @param n A `numeric(1)` specifying the number of peptides to use for each protein. If `NULL`, all peptides are considered. 
#' @param ... Additional parameters passed the `fun`.
#'
#' @return A `QFeatures` object with an additional assay or a `SummarizedExperiment` object (or subclass thereof).
#'
#' @details 
#' This function uses [QFeatures::aggregateFeatures()] to aggregate quantitative data.
#'
#' @section Iterative aggregation function:
#' xxxxxx
#' xxxxx
#'
#' @section Quantitative metadata aggregation:
#' The function to aggregate the quantitative metadata is `aggQmetadat()`.
#'
#' @seealso The *QFeatures* vignette provides an extended example and
#'     the *Aggregation* vignette, for a complete quantitative
#'     proteomics data processing pipeline.
#'
#' @name DaparToolshed-aggregate
#'
#' @return NULL
#'
#' @examples
#'
#' ## ---------------------------------------
#' ## An example QFeatures with PSM-level data
#' ## ---------------------------------------
#' \dontrun{
#' data(ft, package='DaparToolshed')
#' ft
#'
#' ## Aggregate peptides into proteins
#' ## using the adjacency matrix
#' feat1 <- aggregateFeatures4Prostar(object = ft,
#' i = 1,
#' name = 'aggregated',
#' fcol = 'adjacencyMatrix',
#' fun = 'colSumsMat')
#' feat1
#'
#' assay(feat1[[1]])
#' assay(feat1[[2]])
#' aggcounts(feat1[[2]])
#' assay(feat1[[3]])
#' aggcounts(feat1[[3]])
#' rowData(feat1[[2]])
#' }
NULL

#' @exportMethod aggregateFeatures4Prostar
#' @rdname DaparToolshed-aggregate
#' @importFrom MsCoreUtils robustSummary
setMethod(
  "aggregateFeatures4Prostar", "QFeatures",
  function(object, i, fcol, name = "newAssay",
           fun = MsCoreUtils::robustSummary, shared = TRUE, n = NULL, ...) {
    if (length(object) == 0) {
      return(object)
    }
    if (name %in% names(object)) {
      stop("There's already an assay named '", name, "'.")
    }
    if (missing(i)) {
      i <- length(object)
    }
    
    # Add stats on agregation
    aggAssay <- aggregateFeatures4Prostar(
      object = object[[i]],
      fcol = fcol,
      fun = fun,
      conds = design.qf(object)$Condition,
      shared = shared,
      n = n,
      ...
    )
    colData(aggAssay) <- colData(object)
    
    ## Add the assay to the QFeatures object
    object <- addAssay(object,
                       aggAssay,
                       name = name
    )
    
    ## Link the input assay to the aggregated assay
    addAssayLink(object,
                 from = i,
                 to = name,
                 varFrom = fcol,
                 varTo = fcol
    )
  }
)

#' @exportMethod aggregateFeatures4Prostar
#' @rdname DaparToolshed-aggregate
#' @importFrom MsCoreUtils robustSummary
setMethod(
  "aggregateFeatures4Prostar", "SummarizedExperiment",
  function(object, fcol, fun = MsCoreUtils::robustSummary, conds, shared = TRUE, n = NULL, ...) {
    .aggregateFeatures4Prostar(object, fcol, fun, conds, shared, n, ...)
  }
)

.aggregateFeatures4Prostar <- function(object, fcol, fun, conds, shared = TRUE, n = NULL, ...) {
  
  if (!is(rowData(object)[[fcol]], "Matrix")){stop("'fcol' must refer to a matrix. 
                                                    You can create one from a vector using the function PSMatch::makeAdjacencyMatrix.")}
  
  if (fun != "robustSummary"){
    assay(object) <-  2^(assay(object))
  }
  
  matadj <- rowData(object)[[fcol]]
  peptdata <- assay(object)
  
  if (!shared){
    if (length(which(rowSums(matadj) > 1)) != 0){
      matadj[which(rowSums(matadj) > 1), ] <- 0
    }
  }
  
  if (!is.null(n)){
    if (!is.numeric(n)){stop("'n' must be numeric.")}
    matadj <- select_topn(peptdata, matadj, n, funpept = "Mean") 
  }
  
  repeat_counts <- rowSums(matadj)
  
  assay_SE <- peptdata[rep(1:nrow(peptdata), repeat_counts), ]
  rownames(assay_SE) <- NULL
  
  vect_prot <- apply(matadj, 1, function(row) paste(colnames(matadj)[row == 1], collapse = ";"))
  rowdata_SE <- data.frame(unlist(strsplit(vect_prot, ";")))
  colnames(rowdata_SE) <- fcol
  
  copy_object <- SummarizedExperiment(assays = SimpleList(assay = assay_SE),
                                      colData = colData(object),
                                      rowData = rowdata_SE)
  
  ###QUANTITATIVE DATA
  # Create the aggregated assay
  aggAssay <- QFeatures::aggregateFeatures(copy_object, fcol, fun, na.rm = TRUE, ...)
  assays(aggAssay)<- assays(aggAssay)[1]
  assay(aggAssay)[assay(aggAssay) == 0 | is.nan(assay(aggAssay))] <- NA
  
  # If some proteins do not contain any peptide but are included in the adjacency matrix
  if (nrow(assay(aggAssay)) != ncol(matadj)){
    '%ni%' <- Negate('%in%')
    missprot <- colnames(matadj)[which(colnames(matadj) %ni% rownames(assay(aggAssay)))]
    
    # Add missing proteins to assay
    missprot_mat_a <- matrix(0, nrow = length(missprot), ncol = ncol(assay(aggAssay))) 
    rownames(missprot_mat_a) <- missprot
    colnames(missprot_mat_a) <- colnames(assay(aggAssay))
    assay_missprot <- rbind(assay(aggAssay), missprot_mat_a)
    
    # Add missing proteins to  rowData
    missprot_mat_rd <- matrix(NA, nrow = length(missprot), ncol = ncol(rowData(aggAssay))) 
    rownames(missprot_mat_rd) <- missprot
    colnames(missprot_mat_rd) <- colnames(rowData(aggAssay))
    rowdata_missprot <- rbind(rowData(aggAssay), missprot_mat_rd)
    
    # Make new SummarizedExperiment, including missing proteins
    aggAssay <- SummarizedExperiment(assays = SimpleList(assay = assay_missprot),
                                     colData = colData(aggAssay),
                                     rowData = rowdata_missprot)
  }
  
  ###METACELL DATA
  # Removing aggregated rowdata
  rowData(aggAssay) <- NULL
  # Add rowdata
  protname_order <- rownames(assay(aggAssay))
  aggAssay <- metacell_agg(aggAssay, object, matadj, conds, protname_order)
  
  ## Enrich the new assay
  typeDataset(aggAssay) <- "protein"
  idcol(aggAssay) <- NULL
  
  return(aggAssay)
}


#' @title Aggregate an assay's quantitative features with shared peptide redistribution
#'
#' @description
#' This function aggregates the quantitative features of an assay,
#' applying a summarization function (`fun`) to sets of features.
#' The `fcol` variable name points to a rowData column that defines
#' how to group the features during aggregate. This variable has to 
#' be an adjacency matrix. This function uses [DaparToolshed::inner.aggregate.iter()]
#' to aggregate quantitative data.
#'
#' The list of agregation methods can be obtained with the function
#' [aggregateMethods()]. This function compiles both methods from the
#' packages `DaparToolshed` and `QFeatures`.
#'
#' @param object An instance of class `QFeatures` or `SummarizedExperiment`
#' @param i The index or name of the assay which features will be aggregated the create the new assay.
#' @param fcol A `character(1)` naming a rowdata variable (of assay `i` in case of a `QFeatures`) 
#'             defining how to aggregate the features of the assay. 
#'             This variable is a (possibly sparse) matrix. See below for details.
#' @param name A `character(1)` naming the new assay. Default is `newAssay`. 
#'             Note that the function will fail if there's already an assay with `name`.
#' @param init.method A function used for initializing the aggregation. 
#'                    Available functions are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. 
#'                    See [DaparToolshed::inner.aggregate.iter()] for details.
#' @param method A function used for the aggregation. 
#'               Available functions are `Sum`, `Mean`, `Median` or `medianPolish`. 
#'               See [DaparToolshed::inner.aggregate.iter()] for details.
#' @param ponderation A `character(1)` defining what to consider to create the coefficient for redistribution of shared peptides. 
#'                    Available values are `Global` (default), `Condition` or `Sample`. 
#' @param n A `numeric(1)` specifying the number of peptides to use for each protein. If `NULL`, all peptides are considered. 
#' @param max_iter A `numeric(1)` setting the maximum number of iteration.
#'
#' @return A `QFeatures` object with an additional assay or a `SummarizedExperiment` object (or subclass thereof).
#'
#' @details 
#' This function uses [DaparToolshed::inner.aggregate.iter()] to aggregate quantitative data.
#'
#' @section Iterative aggregation function:
#' xxxxxx
#' xxxxx
#'
#' @section Quantitative metadata aggregation:
#' The function to aggregate the quantitative metadata is `aggQmetadat()` 
#'
#' @seealso The *QFeatures* vignette provides an extended example and the *Aggregation* vignette, 
#'          for a complete quantitative proteomics data processing pipeline.
#'
#' @name DaparToolshed-aggregateRedistribution
#'
#' @return NULL
#'
#' @examples
#'
#' ## ---------------------------------------
#' ## An example QFeatures with PSM-level data
#' ## ---------------------------------------
#' \dontrun{
#' data(ft, package='DaparToolshed')
#' ft
#'
#' ## Aggregate peptides into proteins
#' ## using the adjacency matrix
#' feat1 <- aggregateRedistribution(object = ft,
#' i = 1,
#' name = 'aggregated',
#' fcol = 'adjacencyMatrix',
#' init.method = 'Mean',
#' method = 'Mean')
#' feat1
#'
#' assay(feat1[[1]])
#' assay(feat1[[2]])
#' aggcounts(feat1[[2]])
#' assay(feat1[[3]])
#' aggcounts(feat1[[3]])
#' rowData(feat1[[2]])
#' }
NULL

#' @exportMethod aggregateRedistribution
#' @rdname DaparToolshed-aggregateRedistribution
#' @importFrom MsCoreUtils robustSummary
setMethod(
  "aggregateRedistribution", "QFeatures",
  function(object, i, name = "newAssay", fcol, init.method = "Mean", method = "Mean", ponderation = "Global", n = NULL, uniqueiter = FALSE, max_iter = 500) {
    if (length(object) == 0) {
      return(object)
    }
    if (name %in% names(object)) {
      stop("There's already an assay named '", name, "'.")
    }
    if (missing(i)) {
      i <- length(object)
    }
    
    # Add stats on agregation
    aggAssay <- aggregateRedistribution(
      object = object[[i]],
      fcol = fcol,
      init.method = init.method,
      method = method,
      ponderation = ponderation,
      n = n,
      uniqueiter = uniqueiter,
      conds = design.qf(object)$Condition, 
      max_iter = max_iter
    )
    colData(aggAssay) <- colData(object)
    
    ## Add the assay to the QFeatures object
    object <- addAssay(object,
                       aggAssay,
                       name = name
    )
    
    ## Link the input assay to the aggregated assay
    addAssayLink(object,
                 from = i,
                 to = name,
                 varFrom = fcol,
                 varTo = fcol
    )
  }
)

#' @exportMethod aggregateRedistribution
#' @rdname DaparToolshed-aggregateRedistribution
#' @importFrom MsCoreUtils robustSummary
setMethod(
  "aggregateRedistribution", "SummarizedExperiment",
  function(object, fcol, init.method = "Mean", method = "Mean", ponderation = "Global", n = NULL, uniqueiter = FALSE, conds, max_iter = 500) {
    .aggregateRedistribution(object, fcol, init.method, method, ponderation, n, uniqueiter, conds, max_iter)
  }
)


.aggregateRedistribution <- function(object,
                                     fcol,
                                     init.method = "Mean",
                                     method = "Mean",
                                     ponderation = "Global",
                                     n = NULL, 
                                     uniqueiter = FALSE,
                                     conds,
                                     max_iter = 500) {
  if (!is(rowData(object)[[fcol]], "Matrix")){stop("'fcol' must refer to a sparse matrix.")}
  if (!(init.method %in% c("Sum", "Mean", "Median", "medianPolish", "robustSummary"))) {
    stop("Wrong parameter init.method")
  }
  if (!(method %in% c("Sum", "Mean", "Median", "medianPolish", "robustSummary"))) {
    stop("Wrong parameter method")
  }
  if (!(ponderation %in% c("Global", "Condition", "Sample"))) {
    stop("Wrong parameter ponderation")
  }
  if (!is.null(n) & !is.numeric(n)){
    stop("'n' must be NULL or numeric.")
  }
  
  ###QUANTITATIVE DATA
  ## Create the aggregated assay
  quanti_data <- 2^assay(object)
  switch(ponderation,
         Global = {aggAssay <- inner.aggregate.iter(quanti_data, rowData(object)[[fcol]], init.method, method, n, uniqueiter, max_iter = max_iter)},
         Condition = {
           aggAssay <- NULL
           for (condition in unique(conds)){
             pepDatasample <- quanti_data[, which(conds==condition), drop = FALSE]
             aggregsample <- inner.aggregate.iter(as.matrix(pepDatasample),  rowData(object)[[fcol]], init.method, method, n, uniqueiter, max_iter = max_iter)
             aggAssay <- cbind(aggAssay, aggregsample)
           }
           colnames(aggAssay) <- colnames(quanti_data)},
         Sample = {
           aggAssay <- NULL
           for (sample in colnames(quanti_data)){
             pepDatasample <- quanti_data[, sample, drop = FALSE]
             aggregsample <- inner.aggregate.iter(as.matrix(pepDatasample), rowData(object)[[fcol]], init.method, method, n, uniqueiter, max_iter = max_iter)
             aggAssay <- cbind(aggAssay, aggregsample)
           }
           colnames(aggAssay) <- colnames(quanti_data)}
  )
  assay(aggAssay)[assay(aggAssay) == 0 | is.nan(assay(aggAssay))] <- NA
  # Create aggregated SummarizedExperiment
  aggSE <- SummarizedExperiment(assays=SimpleList(assay = aggAssay), colData = colData(object))
  
  ###METACELL DATA
  # Add rowdata
  protname_order <- rownames(aggAssay)
  aggSE <- metacell_agg(aggSE, object, rowData(object)[[fcol]], conds, protname_order)
  
  ## Enrich the new assay
  typeDataset(aggSE) <- "protein"
  idcol(aggSE) <- NULL
  
  return(aggSE)
}


#' @title Get the type of dataset
#' @description xxx
#'
#' @param qMeta  An object of class 'SummarizedExperiment'
#' @param X xxxx
#' @param level A `character(1)` which is the type of dataset
#' @param conds A `character()` vector which is the names of conditions
#' for each sample in the dataset.
#'
#' @return xxxxx
#'
#' @examples
#' data(ft, package='DaparToolshed')
#' qMeta <- qMetacell(ft, 1)
#' X <- adjacencyMatrix(ft[[1]])
#' level <- typeDataset(ft[[1]])
#' conds <- colData(ft)$Condition
#' aggQmeta <- aggQmetacell(qMeta, X, level, conds)
#'
#' @rdname DaparToolshed-aggregate
#' 
#' @export
#'
aggQmetacell <- function(qMeta, X, level, conds) {
  # stopifnot(inherits(object, "SummarizedExperiment"))
  
  res <- list(
    metacell = data.frame(stringsAsFactors = TRUE),
    issues = NULL
  )
  rowcol <- function(meta.col, X.col)  meta.col[X.col > 0]
  
  res$metacell <- as.data.frame(matrix(rep('', ncol(qMeta)*ncol(X)), nrow = ncol(X)), stringsAsFactors = FALSE)
  
  for (j in seq(ncol(qMeta))) {
    for (i in seq(ncol(X))) {
      res$metacell[i, j] <- metacombine(rowcol(qMeta[, j], X[, i]), level)
    }
  }
  
  dimnames(res$metacell) <- list(colnames(X), colnames(qMeta))
  res$metacell[res$metacell == "NA"] <- NA
  # Delete protein with only NA
  
  # Post processing of metacell to discover 'imputed POV', 'imputed MEC'
  res$metacell <- Set_POV_MEC_tags(conds, res$metacell, level)
  
  # Search for issues
  prot.ind <- unique(rownames(which(res$metacell == "STOP", arr.ind = TRUE)))
  if (!is.null(prot.ind)) {
    res$issues <- stats::setNames(
      lapply(
        prot.ind,
        function(x) {
          rownames(X)[which(X[, which(colnames(X) == x)] == 1)]
        }
      ),
      prot.ind
    )
  }
  
  return(res)
}


#' @export
#' @rdname DaparToolshed-aggregate
aggregateMethods <- function() {
  stats::setNames(
    c("medianPolish",
      "robustSummary",
      "colMeansMat",
      "colSumsMat"
    ),
    nm = c(
      "median Polish",
      "robust Summary",
      "col Means Mat",
      "col Sums Mat"
    )
  )
}



#' @title Aggregation of peptide-level assay QFeatures
#' 
#' @description This function aggregate both quantitative and rowdata from the last assay contained in a `QFeatures`. 
#' Note that the function assumes that the intensities in the QFeatures are already log-transformed.
#' 
#' @param qf An instance of class [QFeatures]. The last assay contained in `qf` will be aggregated. Intensities are assumed to already be log-transformed.
#' @param includeSharedPeptides How shared peptides are handled. Either `Yes_As_Specific` (default), `Yes_Iterative_Redistribution`, `Yes_Simple_Redistribution` or `No`.
#'                              See below for details.
#' @param operator A function used for quantitative feature aggregation. 
#'                 Available functions are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. 
#'                 See below for details.
#' @param considerPeptides A `character(1)` defining what peptide to consider. Available values are `allPeptides` (default) and `topN`.
#' @param adjMatrix A `character(1)` naming a rowdata variable from the last assay of `qf` containing an adjacency matrix. 
#' @param ponderation A `character(1)` defining what to consider to create the coefficient for redistribution of shared peptides. 
#'                    Available values are `Global` (default), `Condition` or `Sample`. 
#' @param n If `topN`, specify the number of peptides to use for each protein.
#' @param aggregated_col A `character()` of column names from rowdata to be aggregated. 
#' @param max_iter A `numeric(1)` setting the maximum number of iteration.
#' 
#' @return A [QFeatures] with an aggregated assay added.
#' 
#' @details
#' Aggregation of quantitative data is performed using [aggregateFeatures], or [inner.aggregate.iter] if `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution` is selected.
#' 
#' The handling of shared peptide is as follow : 
#' - `Yes_As_Specific` : Shared peptides are used multiple times. 
#'   Each peptide is duplicated as many times as the number of proteins in which they are present, and thus are considered as if they are specific to each protein. 
#' - `Yes_Simple_Redistribution` : Intensity of shared peptides are redistributed proportionally to each protein. See [inner.aggregate.iter] for more information.
#' - `Yes_Iterative_Redistribution` : Intensity of shared peptides are redistributed proportionally to each protein. See [inner.aggregate.iter] for more information.
#' - `No` : No shared peptides are used. If a peptide contained only shared peptides, its intensity is set as 0 for every sample. 
#' 
#' Available functions are : 
#' - `Sum` : [base::colSums()] or [base::rowSums()] if `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.
#' - `Mean` : [base::colMeans()] or [base::rowMeans()] if `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.
#' - `Median` : [matrixStats::colMedians()] or [matrixStats::rowMedians()] if `Yes_Iterative_Redistribution` or `Yes_Simple_Redistribution`.
#' - `medianPolish` : [MsCoreUtils::medianPolish()].
#' - `robustSummary` : [MsCoreUtils::robustSummary()].
#' 
#' @author Samuel Wieczorek, Manon Gaudin
#' 
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' ft <- Exp1_R25_pept[1:100]
#' obj.agg <- RunAggregation(ft, "Yes_As_Specific", "Sum", "allPeptides", aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' obj.agg <- RunAggregation(ft, "Yes_As_Specific", "Mean", "allPeptides", aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' obj.agg <- RunAggregation(ft, "Yes_As_Specific", "Sum", "topN", n = 4, aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' obj.agg <- RunAggregation(ft, "Yes_As_Specific", "Mean", "topN", n = 4, aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' 
#' obj.agg <- RunAggregation(ft, "No", "Sum", "allPeptides")
#' obj.agg <- RunAggregation(ft, "No", "Sum", "topN", n = 4)
#' 
#' obj.agg <- RunAggregation(ft, "Yes_Redistribution", "Sum", "allPeptides", aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' obj.agg <- RunAggregation(ft, "Yes_Redistribution", "Sum", "topN", n = 4, aggregated_col = colnames(rowData(ft[[length(ft)]])))
#' }
#' 
#' @export
#'
RunAggregation <- function(qf,
                           includeSharedPeptides = 'Yes_As_Specific',
                           operator = 'Mean',
                           considerPeptides = 'allPeptides',
                           adjMatrix = 'adjacencyMatrix',
                           ponderation = 'Global',
                           n = NULL,
                           aggregated_col = NULL,
                           max_iter = 500){
  
  stopifnot(inherits(qf, "QFeatures"))
  if (!(includeSharedPeptides %in% c("Yes_As_Specific", "Yes_Iterative_Redistribution", "Yes_Simple_Redistribution", "No"))) {
    stop(includeSharedPeptides, " is not an option for 'includeSharedPeptides'.")
  }
  if (!(operator %in% c("Sum", "Mean", "Median","medianPolish", "robustSummary"))) {
    stop(operator, " is not an option for 'operator'.")
  }
  if (!(considerPeptides %in% c("allPeptides", "topN"))) {
    stop(considerPeptides, " is not an option for 'considerPeptides'.")
  }
  if (!is.null(n) & !is.numeric(n)) {
    stop("n is not numeric or NULL.")
  }
  
  # Each possible situation
  caseA <- includeSharedPeptides == "Yes_Iterative_Redistribution" && considerPeptides == "allPeptides"
  caseB <- includeSharedPeptides == "Yes_Iterative_Redistribution" && considerPeptides == "topN"
  caseC <- includeSharedPeptides == "Yes_Simple_Redistribution" && considerPeptides == "allPeptides"
  caseD <- includeSharedPeptides == "Yes_Simple_Redistribution" && considerPeptides == "topN"
  caseE <- includeSharedPeptides == "Yes_As_Specific" && considerPeptides == "allPeptides"
  caseF <- includeSharedPeptides == "Yes_As_Specific" && considerPeptides == "topN"
  caseG <- includeSharedPeptides == "No" && considerPeptides == "allPeptides"
  caseH <- includeSharedPeptides == "No" && considerPeptides == "topN"
  
  case <- setNames(c(caseA, caseB, caseC, caseD, caseE, caseF, caseG, caseH),
                   nm = c('caseA', 'caseB', 'caseC', 'caseD', 'caseE', 'caseF', 'caseG', 'caseH'))
  
  message('Aggregating data')
  switch (as.character(names(which(case))),
          caseA = {
            # Redistribution of shared peptides and all peptides
            ll.agg <- aggregateRedistribution(
              object = qf, 
              i = length(qf), 
              name = 'aggregated', 
              fcol = adjMatrix, 
              init.method = operator, 
              method = operator,
              ponderation = ponderation,
              max_iter = max_iter
            )
          },
          caseB = {
            # Redistribution of shared peptides and top n peptides
            ll.agg <- aggregateRedistribution(
              object = qf, 
              i = length(qf), 
              name = 'aggregated', 
              fcol = adjMatrix, 
              init.method = operator, 
              method = operator, 
              ponderation = ponderation,
              n = n,
              max_iter = max_iter
            )
          },
          caseC = {
            # Redistribution of shared peptides and all peptides
            ll.agg <- aggregateRedistribution(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              init.method = operator,
              method = operator,
              ponderation = ponderation,
              uniqueiter = TRUE
            )
          },
          caseD = {
            # Redistribution of shared peptides and top n peptides
            ll.agg <- aggregateRedistribution(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              init.method = operator,
              method = operator,
              ponderation = ponderation,
              n = n,
              uniqueiter = TRUE
            )
          },
          caseE = {
            # Shared peptides as specific and all peptides
            switch(operator,
                   Sum = operator <- "colSums",
                   Mean = operator <- "colMeans", 
                   Median = operator <- "colMedians",
                   medianPolish = operator <- "medianPolish",
                   robustSummary = operator <- "robustSummary"
            )
            
            ll.agg <- aggregateFeatures4Prostar(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              fun = operator
            )
          },
          caseF = {
            # Shared peptides as specific and top n peptides
            switch(operator,
                   Sum = operator <- "colSums",
                   Mean = operator <- "colMeans", 
                   Median = operator <- "colMedians",
                   medianPolish = operator <- "medianPolish",
                   robustSummary = operator <- "robustSummary"
            )
            
            ll.agg <- aggregateFeatures4Prostar(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              fun = operator,
              n = n
            )
          },
          caseG = {
            # Only unique peptides and all peptides
            switch(operator,
                   Sum = operator <- "colSums",
                   Mean = operator <- "colMeans", 
                   Median = operator <- "colMedians",
                   medianPolish = operator <- "medianPolish",
                   robustSummary = operator <- "robustSummary"
            )
            
            ll.agg <- aggregateFeatures4Prostar(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              fun = operator,
              shared = FALSE
            )
          },
          caseH = {
            # Only unique peptides and top n peptides
            switch(operator,
                   Sum = operator <- "colSums",
                   Mean = operator <- "colMeans", 
                   Median = operator <- "colMedians",
                   medianPolish = operator <- "medianPolish",
                   robustSummary = operator <- "robustSummary"
            )
            
            ll.agg <- aggregateFeatures4Prostar(
              object = qf,
              i = length(qf),
              name = 'aggregated',
              fcol = adjMatrix,
              fun = operator,
              shared = FALSE,
              n = n
            )
          }
  )
  
  if(!is.null(aggregated_col)){
    message('Adding aggregated metadata')
    ll.agg <- Add_Aggregated_rowData(ll.agg, aggregated_col, length(ll.agg))
  }
  
  return(ll.agg)
}



#' This function creates a column for the protein dataset after aggregation
#' by using the previous peptide dataset.
#'
#' @title creates a column for the protein dataset after agregation by
#'  using the previous peptide dataset.
#'
#' @param peptideData A data.frame of meta data of peptides. It is the rowData
#' of the SummarizedExperiment object.
#'
#' @param matAdj The adjacency matrix used to agregate the peptides data.
#'
#' @param columnName The name(s) of the column in Biobase::rowData(peptides_MSnset)
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
BuildColumnToProteinDataset <- function(
    peptideData,
    matAdj,
    columnName,
    proteinNames) {
  
  stopifnot(length(columnName) == 1)
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





#' @title Add aggregated rowData
#' 
#' @description Aggregation of rowData of a `QFeatures` assay.
#'
#' @param obj An instance of class [QFeatures].
#' @param col A `character()` of column names from rowdata to be aggregated. 
#' @param i.agg A `numeric(1)` indicating the index of the assay to which add the aggregated rowData, using the previous assay's rowData.  
#'
#' @return An instance of `QFeatures` class with aggregated rowData in specified assay.
#'
#' @author Samuel Wieczorek, Manon Gaudin
#'
#' @export
#'
Add_Aggregated_rowData <- function(obj, col, i.agg){
  stopifnot(inherits(obj, "QFeatures"))
  stopifnot(is.numeric(i.agg))
  stopifnot(i.agg > 1)
  stopifnot(is.character(col))
  
  col <- col[-match(c('qMetacell', 'adjacencyMatrix'), col)]
  if (length(col) == 0) {return(obj)}
  
  matadj <- adjacencyMatrix(obj[[i.agg - 1]])
  protnames <- rownames(rowData((obj[[i.agg]])))
  peptrowdata <- rowData((obj[[i.agg - 1]]))
  
  # For each protein, create the aggregated column
  for (col.name in col) {
    newCol <- BuildColumnToProteinDataset(
      peptideData = peptrowdata,
      matAdj = matadj,
      columnName = col.name,
      proteinNames = protnames
    )
    
    rowData(obj[[i.agg]])[[paste0("agg_", col.name)]] <- newCol
  }
  
  return(obj)
}


#' @title Metacell aggregation
#' 
#' @description  xxx
#' 
#' @param aggregatedSE An instance of class [SummarizedExperiment] containing the aggregated data.
#' @param originalSE An instance of class [SummarizedExperiment] containing the non-aggregated data. 
#' @param adj_mat An adjacency matrix.
#' @param conds 
#' @param protname_order 
#' 
#' @return A `SummarizedExperiment` containing the aggregated data.
#' 
#' @author Samuel Wieczorek, Manon Gaudin 
#' 
#' @export 
#'
metacell_agg <- function(aggregatedSE, originalSE, adj_mat, conds, protname_order){
  # Get metacell from original data
  aggQ <- aggQmetacell(
    qMeta = qMetacell(originalSE),
    X = adj_mat,
    level = typeDataset(originalSE),
    conds = conds
  )
  aggQ$metacell <- aggQ$metacell[protname_order, ]
  
  # Get matrix with number of peptide
  X.spec <- X.shared <- X <- adj_mat
  if (!isEmpty(which(rowSums(as.matrix(X.spec)) > 1))){
    X.spec[which(rowSums(as.matrix(X.spec)) > 1), ] <- 0
  }
  if (!isEmpty(which(rowSums(as.matrix(X.shared)) == 1))){
    X.shared[which(rowSums(as.matrix(X.shared)) == 1), ] <- 0
  }
  .allPep <- t(as.matrix(X)) %*% !is.na(assay(originalSE))
  .allPep <- .allPep[match(protname_order, rownames(.allPep)), , drop = FALSE]
  .specPep <- t(as.matrix(X.spec)) %*% !is.na(assay(originalSE))
  .specPep <- .specPep[match(protname_order, rownames(.specPep)), , drop = FALSE]
  .sharedPep <- t(as.matrix(X.shared)) %*% !is.na(assay(originalSE))
  .sharedPep <- .sharedPep[match(protname_order, rownames(.sharedPep)), , drop = FALSE]
  
  # Add number of peptide
  rowData(aggregatedSE)[["proteinId"]] <- protname_order
  rowData(aggregatedSE)[["nPepTotal"]] <- .allPep[, 1]
  rowData(aggregatedSE)[["nPepShared"]] <-  .sharedPep[, 1]
  rowData(aggregatedSE)[["nPepSpec"]] <- .specPep[, 1]
  
  # Add number of peptide for each sample
  cols <- colnames(aggregatedSE)
  rowData(aggregatedSE)[paste0("pepSpec.used.", cols)] <- as.data.frame(.specPep)
  rowData(aggregatedSE)[paste0("pepShared.used.", cols)] <- as.data.frame(.sharedPep)
  rowData(aggregatedSE)[paste0("pepTotal.used.", cols)] <- as.data.frame(.allPep)
  
  # Add metacell tag
  rowData(aggregatedSE) <- cbind(rowData(aggregatedSE), aggQ$metacell)
  qMetacell(aggregatedSE) <- aggQ$metacell
  metadata(aggregatedSE)[['aggQmetacell_issues']] <- aggQ$issues
  
  # Add complete matrix with number of peptide
  rowData(aggregatedSE)[["allPeptidesUsed"]] <- .allPep
  rowData(aggregatedSE)[["specPeptidesUsed"]] <- .specPep
  rowData(aggregatedSE)[["sharedPeptidesUsed"]] <- .sharedPep
  
  return(aggregatedSE)
}


#' @title Selection of top n peptides
#' 
#' @description  xxx
#' 
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix. 
#' @param n  A `numeric(1)` specifying the number of peptides to use for each protein.
#' @param funpept A function used for determining a peptide's value. 
#'                Available functions are `Sum`, `Mean` or `Median`. 
#' 
#' @return An adjacency matrix with only the top n peptides selected. 
#' 
#' @author Manon Gaudin
#' 
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(50)]
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' X.topn <- select_topn(assay(obj[[length(obj)]]), X, n = 3)
#' }
#'  
#' @export
#'
select_topn <- function(pepData, X, n = 10, funpept = "Mean") {
  
  switch(funpept,
         Sum =  valpept <- rowSums(pepData, na.rm = TRUE),
         Mean =  valpept <- rowMeans(pepData, na.rm = TRUE),
         Median =  valpept <- rowMedians(pepData, na.rm = TRUE)
  )
  names(valpept) <- rownames(pepData)
  
  for (prot in colnames(X)){
    selectpeptname <- names(X[which(X[, prot] == 1), prot])
    selectpeptval <- valpept[selectpeptname]
    
    if (length(selectpeptval) > n){
      selectpeptval <- selectpeptval[order(selectpeptval, decreasing = TRUE)]
      selectpeptname <- names(selectpeptval)[(n+1):length(selectpeptval)]
      X[selectpeptname, prot] <- 0
    }
  }
  
  return(X)
}


#' This function computes the number of proteins that are only defined by
#' specific peptides, shared peptides or a mixture of two.
#'
#' @title Computes the number of proteins that are only defined by
#' specific peptides, shared peptides or a mixture of two.
#'
#' @param X The adjacency matrix with both specific and
#' shared peptides.
#'
#' @return A list
#'
#' @author Samuel Wieczorek, Manon Gaudin
#'
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(20)]
#' obj.last <- obj[[length(obj)]]
#' X <- BuildAdjacencyMatrix(obj.last)
#' getProteinsStats(X)
#' }
#'
#' @export
#'
getProteinsStats <- function(X) {
  if (missing(X)) {
    stop("'X' is needed.")
  }
print("TODO : verifier pourquoi on passe plusieurs fois dedans !!!!!!!")
  stopifnot(!is.null(X))
  
  nbPeptide <- 0
  
  ind.shared.Pep <- which(rowSums(as.matrix(X)) > 1)
  M.shared.Pep <- X[ind.shared.Pep, ]
  if (length(ind.shared.Pep) == 1) {
    j <- which(as.matrix(M.shared.Pep) == 0)
    if (length(j) != 0) {M.shared.Pep <- M.shared.Pep[-j]}
    pep.names.shared <- names(M.shared.Pep)
  } else {
    j <- which(colSums(as.matrix(M.shared.Pep)) == 0)
    if (length(j) != 0) {M.shared.Pep <- M.shared.Pep[, -j]}
    pep.names.shared <- colnames(M.shared.Pep)
  }
  
  
  ind.unique.Pep <- which(rowSums(as.matrix(X)) == 1)
  M.unique.Pep <- X[ind.unique.Pep, ]
  if (length(ind.unique.Pep) == 1) {
    j <- which(as.matrix(M.unique.Pep) == 0)
    if (length(j) != 0) {M.unique.Pep <- M.unique.Pep[-j]}
    pep.names.unique <- names(M.unique.Pep)
  } else {
    j <- which(colSums(as.matrix(M.unique.Pep)) == 0)
    if (length(j) != 0) {M.unique.Pep <- M.unique.Pep[, -j]}
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
#' @param X A "valued" adjacency matrix in which lines and columns correspond
#' respectively to peptides and proteins.
#'
#' @return A vector of boolean which is the adjacency matrix
#' but with NA values if they exist in the intensity matrix.
#'
#' @author Alexia Dorffer
#'
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' CountPep(X)
#' }
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
#' \dontrun{
#' library(QFeatures)
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' GetNbPeptidesUsed(assay(last.obj), X)
#' }
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
#' @param pepData A data.frame of quantitative data
#'
#' @return A list of two items
#'
#' @author Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' \dontrun{
#' library(DaparToolshedData)
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' ll.n <- GetDetailedNbPeptidesUsed(assay(last.obj), X)
#' }
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
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj.pep <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- obj.pep[[length(obj.pep)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' n <- GetDetailedNbPeptides(X)
#' }
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
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- Exp1_R25_pept[[length(Exp1_R25_pept)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' GraphPepProt(X)
#' }
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
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' last.obj <- Exp1_R25_pept[[length(Exp1_R25_pept)]]
#' X <- BuildAdjacencyMatrix(last.obj)
#' ExtractUniquePeptides(X)
#' }
#' 
ExtractUniquePeptides <- function(X){
  ll <- which(rowSums(X) > 1)
  if (length(ll) > 0) {
    X[ll, ] <- 0
  }
  
  return(X)
}


#' @title Aggregation of quantitative data with redistribution of shared peptides
#' 
#' @description This function aggregate quantitative data using a method of redistribution of shared peptides. Intensity of shared peptides are redistributed proportionally to each protein. 
#' Note that the function assumes that the intensities are not log-transformed.
#'
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix. 
#' @param init.method A function used for initializing the aggregation. 
#'                    Available functions are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. 
#'                    See below for details.
#' @param method A function used for the aggregation. 
#'               Available functions are `Sum`, `Mean`, `Median`, `medianPolish` or `robustSummary`. 
#'               See below for details.
#' @param n A `numeric(1)` specifying the number of peptides to use for each protein. If `NULL`, all peptides are considered. 
#' @param uniqueiter A bole
#' @param topn_fun A function used to determine how to choose the top n peptides. 
#'                 Available functions are `Sum`, `Mean` or `Median`. 
#'                 See below for details.
#' @param max_iter A `numeric(1)` setting the maximum number of iteration.
#'
#' @return A `matrix` containing the aggregated values. 
#' 
#' @details
#' Available functions are : 
#' - `Sum` : [base::rowSums()]
#' - `Mean` : [base::rowMeans()] 
#' - `Median` : [matrixStats::rowMedians()]
#' - `medianPolish` : [MsCoreUtils::medianPolish()], not available for `topn_fun`. 
#'   Note that this method takes significantly more time than the others, and is parallelized to be more efficient. 
#' - `robustSummary` : [MsCoreUtils::robustSummary()], not available for `topn_fun`. 
#'   Note that this method takes significantly more time than the others, and is parallelized to be more efficient. 
#' 
#' @author Samuel Wieczorek, Manon Gaudin
#'
#' @examples
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' qdata.agg <- inner.aggregate.iter(assay(obj[[length(obj)]]), X)
#' }
#' 
#' @export
#'
inner.aggregate.iter <- function(
    pepData,
    X,
    init.method = "Mean",
    method = "Mean",
    n = NULL,
    uniqueiter = FALSE,
    topn_fun = "Mean",
    max_iter = 500
) {
  if (!is(X, "Matrix") & !is(X, "matrix")){stop("'X' must refer to a matrix.")}
  if (!(init.method %in% c("Sum", "Mean", "Median", "medianPolish", "robustSummary"))) {
    stop("Wrong parameter init.method")
  }
  if (!(method %in% c("Sum", "Mean", "Median", "medianPolish", "robustSummary"))) {
    stop("Wrong parameter method")
  }
  if (!is.null(n) && !is.numeric(n)) {
    stop("Parameter n is not numeric")
  }
  if (!is.numeric(max_iter)) {
    stop("Parameter max_iter must be a strictly positive integer.")
  }
  if ((max_iter < 1) | (max_iter != as.integer(max_iter))) {
    stop("Parameter max_iter must be a strictly positive integer.")
  }
  
  X.split <- DaparToolshed::splitAdjacencyMat(X)
  ProtSharedPept <- names(which(colSums(X.split$Xshared) !=0 & colSums(X.split$Xspec) !=0))
  
  yprot <- NULL
  # Initialisation
  switch(init.method,
         Sum = yprot <- inner.sum(pepData, X.split$Xspec),
         Mean = yprot <- inner.mean(pepData, X.split$Xspec), 
         Median = yprot <- inner.median(pepData, X.split$Xspec),
         medianPolish = yprot <- inner.medianpolish(pepData, X.split$Xspec),
         robustSummary = yprot <- 2^inner.robustsummary(log2(pepData), X.split$Xspec)
  )
  
  if (!is.null(n)){
    stopifnot(is.numeric(n))
    if (!(topn_fun %in% c("Sum", "Mean", "Median"))) {
      stop("Can only use 'Sum', 'Mean' or 'Median' with top n peptides.")
    }
    # Select the method used
    switch(topn_fun,
           Sum = topnfun <- match.fun("rowSums"),
           Mean = topnfun <- match.fun("rowMeans"),
           Median = topnfun <- match.fun("rowMedians")
    )
    
    # Select get value for each peptide
    valpept <- topnfun(pepData, na.rm = TRUE)
  }
  
  conv <- 1
  n_iter <- 1
  while (conv > 10**(-10)) {
    # print("____________iter___________")
    # print(n_iter)
    
    # Get the value of each protein
    val.prot <- rowMeans(yprot, na.rm = TRUE)
    val.prot[is.na(val.prot)] <- 0
    # Calculate each peptide coeficient
    X.tmp <- t(t(X) * val.prot) #val.prot * X
    X.new <- X.tmp / rowSums(X.tmp, na.rm = TRUE)
    X.new[is.na(X.new)] <- 0
    
    # If only consider top n peptides
    if (!is.null(n)){
      # Select top n peptide for each protein
      X.pept <- X.new
      # @x vector for values, @i for row index
      X.pept@x <- X.pept@x * valpept[X.pept@i + 1]  
      
      for (j in seq_len(ncol(X.pept))) {
        col_start <- X.pept@p[j] + 1
        col_end <- X.pept@p[j + 1]
        
        if (col_end > col_start) {
          vals <- X.pept@x[col_start:col_end]
          rows <- X.pept@i[col_start:col_end] + 1
          if (length(vals) > n) {
            keep_idx <- order(vals, decreasing = TRUE)[1:n]
            drop_rows <- rows[-keep_idx]
            X.new[drop_rows, j] <- 0
          }
        }
      }
    }
    
    # Get aggregated value with peptide redistribution
    switch(method,
           Mean = yprot <- inner.mean(pepData, X.new),
           Sum = yprot <- inner.sum(pepData, X.new),
           Median = {X.new <- X.new[, ProtSharedPept]
           protval <- inner.median(pepData, X.new)
           yprot[rownames(protval),] <- protval},
           medianPolish = {X.new <- X.new[, ProtSharedPept]
           protval <- inner.medianpolish(pepData, X.new)
           yprot[rownames(protval),] <- protval},
           robustSummary = {X.new <- X.new[, ProtSharedPept]
           pepDatalog <- log2(pepData)
           pepDatalog[which(pepDatalog == "-Inf")] <- NA
           protval <- inner.robustsummary(pepDatalog, X.new)
           protval <- 2^protval
           protval[which(protval == 1)] <- NA
           yprot[rownames(protval),] <- protval}
    )
    
    # Get the value of each protein with the new coefficient 
    val.prot.new <- rowMeans(yprot, na.rm = TRUE)
    val.prot.new[is.na(val.prot.new)] <- 0
    if (uniqueiter | (n_iter == max_iter)){
      conv <- 0
    } else {
      convN <- conv
      conv <- mean(abs(val.prot-val.prot.new)/val.prot.new, na.rm = T) 
      if ((convN - conv) == 0){
        conv <- 0
      }
    }
    n_iter <- n_iter + 1
  }
  
  return(as.matrix(yprot))
}



#' @title Sum-based aggregation
#' 
#' @description Aggregation using sum method. 
#'
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix. 
#'
#' @return  A `matrix` containing the aggregated values. 
#'
#' @author Samuel Wieczorek
#' 
#' @examples 
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj[[length(obj)]])
#' i.sum <- inner.sum(assay(obj[[length(obj)]]), X)
#' }
#' 
#' @export
#' 
inner.sum <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  #stopifnot(inherits(X, 'matrix'))
  
  pepData[is.na(pepData)] <- 0
  
  Mp <- t(as.matrix(X)) %*% pepData
  return(Mp)
}



#' @title Mean-based aggregation
#' 
#' @description Aggregation using mean method. 
#'
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix. 
#'
#' @return A `matrix` containing the aggregated values.
#'
#' @author Samuel Wieczorek
#' 
#' @examples 
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj)
#' i.mean <- inner.mean(assay(obj), X)
#' }
#' 
#' @export
#' 
inner.mean <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  #stopifnot(inherits(X, 'matrix'))
  
  pepData[is.na(pepData)] <- 0
  Mp <- inner.sum(pepData, X)
  Mp <- Mp / GetNbPeptidesUsed(pepData, X)
  
  return(Mp)
}


#' @title Median-based aggregation
#' 
#' @description Aggregation using median method. 
#' 
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix.  
#' 
#' @return A `matrix` containing the aggregated values.
#' 
#' @author Manon Gaudin
#' 
#' @examples 
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj)
#' i.mean <- inner.median(assay(obj), X)
#' }
#' 
#' @export
#' 
inner.median <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  #stopifnot(inherits(X, 'matrix'))
  
  protval <- NULL
  for (prot in colnames(X)){
    coef <- X[, prot]
    peptval <- pepData*coef
    peptval[peptval == 0] <- NA
    medprot <- colMedians(peptval, na.rm = TRUE)
    protval <- rbind(protval, medprot)
  }
  rownames(protval) <- colnames(X)
  
  return(protval)
}


#' @title medianPolish-based aggregation
#' 
#' @description Aggregation using medianPolish method. 
#' Note that this method is parallelized to be more efficient. 
#' 
#' @param pepData A `matrix` containing the peptide intensities. 
#' @param X A `matrix` acting as an adjacency matrix. 
#' 
#' @return A `matrix` containing the aggregated values. 
#' 
#' @author Manon Gaudin
#' 
#' @examples 
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj)
#' i.mean <- inner.medianpolish(assay(obj), X)
#' }
#' 
#' @export
#' 
inner.medianpolish <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  #stopifnot(inherits(X, 'matrix'))
  
  protval <- NULL
  for (prot in colnames(X)){
    coef <- X[, prot]
    peptval <- pepData*coef
    peptval[peptval == 0] <- NA
    peptvalnoNA <- peptval[which(rowSums(peptval, na.rm = TRUE) != 0), , drop = FALSE]
    if (nrow(peptvalnoNA) > 1) { # If more than 1 peptide
      if (ncol(peptvalnoNA) == 1) { # If only 1 column
        # If only 1 column, medianPolish act the same as colMedians but colMedians is faster
        medprot <- colMedians(peptvalnoNA, na.rm = TRUE)
      } else { 
        medprot <- MsCoreUtils::medianPolish(peptvalnoNA, na.rm = TRUE)
      }
    } else if (nrow(peptvalnoNA) == 1) { # If only 1 peptide
      medprot <- peptvalnoNA[,]
    } else { # If no peptide
      medprot <- rep(0, ncol(peptval))
    }
    protval <- rbind(protval, medprot)
  }
  rownames(protval) <- colnames(X)
  
  # Ncpus <- parallel::detectCores() - 1
  # cl <- parallel::makeCluster(Ncpus)
  # doParallel::registerDoParallel(cl)
  # 
  # protval <- NULL
  # X <- as.matrix(X)
  # protval <- foreach::foreach(prot = colnames(X), .combine = rbind, .packages = "matrixStats") %dopar% {
  #   coef <- X[, prot]
  #   peptval <- pepData * coef
  #   peptval[peptval == 0] <- NA
  #   peptvalnoNA <- peptval[which(rowSums(peptval, na.rm = TRUE) != 0),]
  #   medprot <- MsCoreUtils::medianPolish(peptvalnoNA, na.rm = TRUE)
  #   medprot
  # }
  # parallel::stopCluster(cl)
  # rownames(protval) <- colnames(X)
  
  return(protval)
}



#' @title robustSummary-based aggregation
#' 
#' @description Aggregation using robustSummary method. 
#' 
#' @param pepData A `matrix` containing the peptide intensities. 
#'                Note that the function assume that data is already log-transformed.
#' @param X A `matrix` acting as an adjacency matrix.  
#' 
#' @return A `matrix` containing the aggregated values.
#' 
#' @author Manon Gaudin
#' 
#' @examples 
#' \dontrun{
#' data(Exp1_R25_pept, package="DaparToolshedData")
#' obj <- Exp1_R25_pept[seq_len(10)]
#' X <- BuildAdjacencyMatrix(obj)
#' i.mean <- inner.robustSummary(assay(obj), X)
#' }
#' 
#' @export
#' 
inner.robustsummary <- function(pepData, X) {
  stopifnot(inherits(pepData, 'matrix'))
  #stopifnot(inherits(X, 'matrix'))
  
  protval <- NULL
  for (prot in colnames(X)){
    coef <- X[, prot]
    peptval <- pepData*coef
    peptval[peptval == 0] <- NA
    peptvalnoNA <- peptval[which(rowSums(peptval, na.rm = TRUE) != 0), , drop = FALSE]
    if (nrow(peptvalnoNA) > 1) { # If more than 1 peptide
      if (ncol(peptvalnoNA) == 1) { # If only 1 column
        # Duplicate the column and give them different names as robust summary can't work if there is only 1 column 
        peptvaldup <- cbind(peptvalnoNA, peptvalnoNA)
        colnames(peptvaldup) <- c(colnames(peptvalnoNA), "dup")
        medprot <- MsCoreUtils::robustSummary(peptvaldup)
        medprot <- medprot[1]
      } else { 
        medprot <- MsCoreUtils::robustSummary(peptvalnoNA)
      }
    } else if (nrow(peptvalnoNA) == 1) { # If only 1 peptide
      medprot <- peptvalnoNA[,]
    } else { # If no peptide
      medprot <- rep(0, ncol(peptval))
    }
    protval <- rbind(protval, medprot)
  }
  rownames(protval) <- colnames(X)
  
  return(protval)
}

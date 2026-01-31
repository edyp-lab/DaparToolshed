## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setups, include=FALSE----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(123)

## ----setup, message = FALSE---------------------------------------------------
library(DaparToolshed)

## ----dataset------------------------------------------------------------------
data.file <- system.file("extdata/data", "Exp1_R25_pept_500.txt", package="DaparToolshed")
data <- read.table(data.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
  
sample.file <- system.file("extdata/data", "samples_Exp1_R25.txt", package="DaparToolshed")
sample <- read.table(sample.file, header=TRUE, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)

## ----importdata, message = FALSE----------------------------------------------
obj <- createQFeatures(data = data,
                       sample = sample,
                       indQData = 56:61,
                       keyId = "Sequence",
                       indexForMetacell = 43:48,
                       logData = TRUE,
                       force.na = TRUE,
                       typeDataset = "peptide",
                       parentProtId = "Protein_group_IDs",
                       analysis = "Pept_Data",
                       software = "maxquant",
                       name = "original")

obj

## ----metadataaccess, message = FALSE------------------------------------------
 #quantitative data tags
head(qMetacell(obj[[2]]))

 #metadata
head(SummarizedExperiment::rowData(obj[[2]]), n = 3)

## ----filters, message = FALSE-------------------------------------------------
 #create filter to remove empty lines
filter_emptyline <- FunctionFilter("qMetacellWholeLine",
                                  cmd = 'delete',
                                  pattern = 'Missing MEC')
  
 #create filter to remove contaminant
filter_contaminant <- QFeatures::VariableFilter(field = "Potential_contaminant",
                                                value = "+",
                                                condition = "==",
                                                not = TRUE)

## ----filtering, message = FALSE-----------------------------------------------
 #apply all filters and create new assay
obj <- filterFeaturesOneSE(object = obj,
                           i = length(obj),
                           name = "Filtered",
                           filters = list(filter_emptyline, filter_contaminant))

## ----removeprot, message = FALSE----------------------------------------------
 #remove proteins with no associated peptides
X <- SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix
SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix <- X[, -which(Matrix::colSums(X) == 0)]

obj

## ----normalizationmethods-----------------------------------------------------
 #list of available methods
normalizeMethods()

## ----normalization------------------------------------------------------------
obj <- normalizeFunction(obj,
                         method = "MeanCentering",
                         scaling = TRUE,
                         type = "overall")

obj

## ----normalizationplot, message = FALSE, warning=FALSE------------------------
obj1 <- obj[[length(obj)]]
obj2 <- obj[[length(obj)-1]]
protId <- DaparToolshed::idcol(obj1)

.n <- floor(0.02 * nrow(obj1))
.subset <- seq(nrow(obj1))

compareNormalizationD_HC(
  qDataBefore = SummarizedExperiment::assay(obj1),
  qDataAfter = SummarizedExperiment::assay(obj2),
  keyId = SummarizedExperiment::rowData(obj1)[, protId],
  conds = DaparToolshed::design.qf(obj)$Condition,
  n = .n,
  subset.view = .subset
)

## ----missingvalueplots, warning=FALSE-----------------------------------------
pal <- unique(GetColorsForConditions(design.qf(obj)$Condition))
pattern <- c("Missing MEC", "Missing POV")
grp <- design.qf(obj)$Condition

 #number of line with different amount of NA
metacellPerLinesHisto_HC(obj[[length(obj)]], group = grp, pattern = pattern)
 #number of line with different amount of NA per condition
metacellPerLinesHistoPerCondition_HC(obj[[length(obj)]], group = grp, pattern = pattern)

## ----missingvaluedistribplots-------------------------------------------------
hc_mvTypePlot2(obj[[length(obj)]], group = grp, pattern = pattern, pal = pal)

## ----imputation, message = FALSE----------------------------------------------
obj <- wrapper.pirat(data = obj,
                     adjmat = SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix,
                     extension = "base")

obj

## ----aggregationadjmatrix-----------------------------------------------------
X <- SummarizedExperiment::rowData(obj[[length(obj)]])$adjacencyMatrix
info_pept <- getProteinsStats(X)

cat("Total number of peptides : ", info_pept$nbPeptides, "\n",
    "Number of specific peptides : ", info_pept$nbSpecificPeptides, "\n",
    "Number of shared peptides : ", info_pept$nbSharedPeptides, "\n")

cat("Total number of proteins : ", info_pept$nbProt, "\n",
    "Number of protein with only specific peptides : ", length(info_pept$protOnlyUniquePep), "\n",
    "Number of protein with only shared peptides : ", length(info_pept$protOnlySharedPep), "\n",
    "Number of protein with both : ", length(info_pept$protMixPep), "\n")

## ----aggregation, message = FALSE---------------------------------------------
obj <- RunAggregation(obj,
                      adjMatrix = 'adjacencyMatrix',
                      includeSharedPeptides = 'Yes_Iterative_Redistribution',
                      operator = 'Mean',
                      considerPeptides = 'allPeptides',
                      ponderation = "Global",
                      max_iter = 500
                      )

obj

## ----aggregationsummary-------------------------------------------------------
summary(SummarizedExperiment::assay(obj[["aggregated"]]))

## ----filteringprot------------------------------------------------------------
 #apply the filter to remove empty lines created
obj <- filterFeaturesOneSE(object = obj,
                           i = length(obj),
                           name = "FilterProt",
                           filters = list(filter_emptyline))

obj

## ----limma--------------------------------------------------------------------
res_pval_FC <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), 
                                 design.qf(obj), 
                                 comp.type = "OnevsOne")

str(res_pval_FC)

## ----foldchange, message = FALSE, warning=FALSE-------------------------------
 #logFC threshold
Foldchange_thlogFC <- 0.04
        
hc_logFC_DensityPlot(
  df_logFC = as.data.frame(res_pval_FC$logFC),
  th_logFC = Foldchange_thlogFC
  )

## ----pushpvalue---------------------------------------------------------------
pval <- unlist(res_pval_FC$P_Value)

 #push p-values for proteins with more than 60% of values derived from imputation
pval <- pushpvalue(obj,
                   pval,
                   scope = "WholeMatrix",
                   pattern = c("Imputed POV", "Imputed MEC"),
                   percent = FALSE,
                   threshold = 5,
                   operator = ">=")

cat("Number of pushed p-values : ", length(which(pval == 1)), "\n")

## ----pvalcalibration----------------------------------------------------------
 #remove protein with logFC under threshold
pval_logfc_ind <- which(abs(res_pval_FC$logFC) >= Foldchange_thlogFC)
 #remove pushed p-values
pushedpval_ind <- which(pval == 1)
pval_logfc_ind <- setdiff(pval_logfc_ind, pushedpval_ind)
pval_logfc <- pval[pval_logfc_ind]

 #calibration plot with all methods 
calibration_all <- wrapperCalibrationPlot(vPVal = pval_logfc, 
                                          pi0Method  = "ALL")
calibration_all$pi0

## ----pvalcalibrationmethod----------------------------------------------------
 #chosen pi0
pi0 <- 1

 #calibration plot with chosen method
wrapperCalibrationPlot(vPVal = pval_logfc, 
                       pi0Method  = pi0)

 #histogram of p-values density
histPValue_HC(pval_logfc,
              bins = 50,
              pi0 = pi0
)

 #calculate adjusted p-values
adjusted_pvalues <- diffAnaComputeAdjustedPValues(
  pval_logfc,
  pi0)

adjusted_pvalues_comp <- unlist(res_pval_FC$P_Value)
adjusted_pvalues_comp[pval_logfc_ind] <- adjusted_pvalues

## ----FDRcontrol---------------------------------------------------------------
 #define p-value threshold
pval_threshold <- 0.01
FDRcontrol_thpval <- -log10(pval_threshold)

 #get FDR
logpval <- -log10(pval)
logpval_thpval_ind <- which(logpval >= FDRcontrol_thpval)
fdr <- diffAnaComputeFDR(adjusted_pvalues_comp[logpval_thpval_ind])

 #determine differentially abundant proteins
isDifferential <- is.differential(pvalue = pval,
                                 logFC = res_pval_FC$logFC,
                                 thpvalue = FDRcontrol_thpval,
                                 thlogFC = Foldchange_thlogFC)
NbSignificant <- length(which(isDifferential == 1))

cat("pvalue threshold : ", pval_threshold, "\n",
    "logpvalue threshold : ", FDRcontrol_thpval, "\n",
    "FDR : ", round(100*fdr, 2), "%\n",
    "Number significant proteins : ", NbSignificant, "\n",
    "Number of expected false discovery : ", fdr*NbSignificant, "\n")

## ----volcanoplot--------------------------------------------------------------
 #create dataframe for volcano plot
df <- data.frame(
  x = unlist(res_pval_FC$logFC),
  y = logpval,
  index = as.character(rownames(obj[[length(obj)]]))
  )
tooltip <- "proteinId"
df <- cbind(df, SummarizedExperiment::rowData(obj[[length(obj)]])[, tooltip, drop = F])
colnames(df) <- gsub(".", "_", colnames(df), fixed = TRUE)
if (ncol(df) > 3) {
  colnames(df)[4:ncol(df)] <- paste0("tooltip_",
                                      colnames(df)[4:ncol(df)])
}
cond <- unique(design.qf(obj)$Condition)

 #volcano plot
diffAnaVolcanoplot_rCharts(
  df,
  th_pval = FDRcontrol_thpval,
  th_logfc = Foldchange_thlogFC,
  conditions = cond
)

## ----DAassay------------------------------------------------------------------
 #create new assay with all values of interest
obj <- QFeatures::addAssay(obj, obj[[length(obj)]], 'DifferentialAnalysis')
SummarizedExperiment::rowData(obj[[length(obj)]])$pval <- unlist(res_pval_FC$P_Value)
SummarizedExperiment::rowData(obj[[length(obj)]])$logpval <- logpval
SummarizedExperiment::rowData(obj[[length(obj)]])$logFC <- unlist(res_pval_FC$logFC)

SummarizedExperiment::rowData(obj[[length(obj)]])$adjusted_pval <- adjusted_pvalues_comp

SummarizedExperiment::rowData(obj[[length(obj)]])$isDifferential <- isDifferential

obj

## ----sessioninfo, echo=FALSE--------------------------------------------------
sessionInfo()


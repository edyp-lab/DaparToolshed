library(testthat)
library(DaparToolshed)

# prepare data
## peptide
peptname <- paste0("Pept", seq_len(10))
protname <- paste0("Prot", seq_len(6))

numdatapept <- data.frame(S1 = c(2, 5, 9, 10, 7, 8, 6, 8, 6, 7),
                          S2 = c(7, 8, 5, 5, 7, 9, 9, 4, 9, 1),
                          S3 = c(4, 1, 4, 9, 6, 8, 5, 3, 4, 3),
                          S4 = c(9, 8, 9, 8, 3, 6, 4, 1, 4, 5),
                          names = peptname)
datapept <- QFeatures::readQFeatures(numdatapept, quantCols = seq_len(4), fnames = "names", name = "datapept")
S4Vectors::metadata(datapept[[1]])[["typeDataset"]] <- "peptide"

coldata <- data.frame(quantCols = c("S1", "S2", "S3", "S4"),
                      Condition = c("C1", "C1", "C2", "C2"),
                      Bio.Rep   = as.character(seq_len(4)))
rownames(coldata) <- coldata$quantCols
SummarizedExperiment::colData(datapept) <- coldata

adjmat <- matrix(
  c(1, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 1, 0,
    0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1),
  byrow = TRUE,
  nrow = 10,
  ncol = 6,
  dimnames = list(peptname, protname))
adjmat <- as(adjmat, "dgCMatrix")
SummarizedExperiment::rowData(datapept[[1]])$adjacencyMatrix <- adjmat

metacellpept <- data.frame(metacell_S1 = rep("Quant. by direct id", 10),
                           metacell_S2 = rep("Quant. by direct id", 10),
                           metacell_S3 = rep("Quant. by direct id", 10),
                           metacell_S4 = rep("Quant. by direct id", 10))
rownames(metacellpept) <- peptname
SummarizedExperiment::rowData(datapept[[1]])$qMetacell <- metacellpept

S4Vectors::metadata(datapept[[1]])$idcol <- "names"
SummarizedExperiment::rowData(datapept[[1]])$testcol <- paste0("V", seq_len(10))

### with NA
datapeptNA <- datapept
SummarizedExperiment::assay(datapeptNA[[1]])[which(SummarizedExperiment::assay(datapeptNA[[1]]) %in% c(1, 8))] <- NA
metacellpeptna <- data.frame(metacell_S1 = c(rep("Quant. by direct id", 5), "Missing POV", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 2)),
                             metacell_S2 = c("Quant. by direct id", "Missing POV", rep("Quant. by direct id", 7), "Missing POV"),
                             metacell_S3 = c("Quant. by direct id", "Missing MEC", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 4)),
                             metacell_S4 = c("Quant. by direct id", "Missing MEC", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 2)))
rownames(metacellpeptna) <- peptname
SummarizedExperiment::rowData(datapeptNA[[1]])$qMetacell <- metacellpeptna



## protein
protname <- paste0("Prot", seq_len(10))

numdataprot <- data.frame(S1 = c(2, 5, 9, 10, 7, 8, 6, 8, 6, 7),
                          S2 = c(7, 8, 5, 5, 7, 9, 9, 4, 9, 1),
                          S3 = c(4, 1, 4, 9, 6, 8, 5, 3, 4, 3),
                          S4 = c(9, 8, 9, 8, 3, 6, 4, 1, 4, 5),
                          names = protname)
dataprot <- QFeatures::readQFeatures(numdataprot, quantCols = seq_len(4), fnames = "names", name = "dataprot")
S4Vectors::metadata(dataprot[[1]])[["typeDataset"]] <- "protein"

coldata <- data.frame(quantCols = c("S1", "S2", "S3", "S4"),
                      Condition = c("C1", "C1", "C2", "C2"),
                      Bio.Rep   = as.character(seq_len(4)))
rownames(coldata) <- coldata$quantCols
SummarizedExperiment::colData(dataprot) <- coldata

metacellprot <- data.frame(metacell_S1 = rep("Quant. by direct id", 10),
                           metacell_S2 = rep("Quant. by direct id", 10),
                           metacell_S3 = rep("Quant. by direct id", 10),
                           metacell_S4 = rep("Quant. by direct id", 10))
rownames(metacellprot) <- protname
SummarizedExperiment::rowData(dataprot[[1]])$qMetacell <- metacellprot

S4Vectors::metadata(dataprot[[1]])$idcol <- "names"

### with NA
dataprotNA <- dataprot
SummarizedExperiment::assay(dataprotNA[[1]])[which(SummarizedExperiment::assay(dataprotNA[[1]]) %in% c(1, 8))] <- NA
metacellprotna <- data.frame(metacell_S1 = c(rep("Quant. by direct id", 5), "Missing POV", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 2)),
                             metacell_S2 = c("Quant. by direct id", "Missing POV", rep("Quant. by direct id", 7), "Missing POV"),
                             metacell_S3 = c("Quant. by direct id", "Missing MEC", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 4)),
                             metacell_S4 = c("Quant. by direct id", "Missing MEC", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 2)))
rownames(metacellprotna) <- protname
SummarizedExperiment::rowData(dataprotNA[[1]])$qMetacell <- metacellprotna



test_that("check history", {
  # Initialization
  resHini <- InitializeHistory()
  
  expect_true(is.data.frame(resHini))
  
  # Add
  resHadd <- Add2History(resHini, "Process1", "Step1", "Parameter1", "Value1")
  
  expect_true(is.data.frame(resHadd))
  
  # Add to SE
  resHaddSE <- SetHistory(dataprot[[1]], resHadd)
  
  expect_true(is(resHaddSE, "SummarizedExperiment"))
  
  # Get
  resHget <- GetHistory(resHaddSE)
  
  expect_true(is.data.frame(resHget))
})

test_that("check filtering", {
  # Filters
  resFilter1 <- FunctionFilter('qMetacellWholeLine', cmd = 'delete', pattern = 'Missing POV')
  resFilter2 <- FunctionFilter('specPeptides', list())
  resFilter3 <- FunctionFilter('topnPeptides', fun = 'rowSums', top = 2)
  resFilter4 <- FunctionFilter('qMetacellOnConditions',
                               cmd = 'delete',
                               mode = 'AtLeastOneCond',
                               pattern = 'Missing POV',
                               conds = SummarizedExperiment::colData(dataprotNA)$Condition,
                               percent = TRUE,
                               th = 0.4,
                               operator = '<')
    
  expect_true(is(resFilter1, "FunctionFilter"))
  expect_true(is(resFilter2, "FunctionFilter"))
  expect_true(is(resFilter3, "FunctionFilter"))
  expect_true(is(resFilter4, "FunctionFilter"))
  
  
  # Functions
  resFilterFun1 <- filterFeaturesOneSE(dataprotNA, 1, filters = list(resFilter4))
  resFilterFun2 <- filterFeaturesOneSE(datapeptNA, 1, filters = list(resFilter2))

  expect_true(is(resFilterFun1, "QFeatures"))
  expect_true(is(resFilterFun2, "QFeatures"))
  expect_equal(length(resFilterFun1), 2)
  expect_equal(length(resFilterFun2), 2)
})


test_that("check normalization", {
  # Function
  ## GlobalQuantileAlignment
  resNormGQA <- normalizeFunction(dataprot, "GlobalQuantileAlignment")
  
  expect_true(is(resNormGQA, "QFeatures"))
  expect_equal(length(resNormGQA), 2)
  
  ## QuantileCentering
  resNormQC <- normalizeFunction(dataprot, "QuantileCentering")
  
  expect_true(is(resNormQC, "QFeatures"))
  expect_equal(length(resNormQC), 2)
  
  ## MeanCentering
  resNormMC <- normalizeFunction(dataprot, "MeanCentering")
  
  expect_true(is(resNormMC, "QFeatures"))
  expect_equal(length(resNormMC), 2)
  
  ## SumByColumns
  resNormSBC <- normalizeFunction(dataprot, "SumByColumns")
  
  expect_true(is(resNormSBC, "QFeatures"))
  expect_equal(length(resNormSBC), 2)
  
  ## LOESS
  resNormLOESS <- normalizeFunction(dataprot, "LOESS")
  
  expect_true(is(resNormLOESS, "QFeatures"))
  expect_equal(length(resNormLOESS), 2)
})


test_that("check imputation", {
  grp <- design_qf(dataprotNA)$Condition
  
  # Protein level
  ## Fixed value
  resprotFV <- wrapperImputeFixedValue(dataprotNA[[1]], grp, 5, na.type = "Missing MEC")
  
  expect_true(is(resprotFV, "SummarizedExperiment"))
  expect_true("Imputed MEC" %in% unique(unlist(qMetacell(resprotFV))))
  expect_true("Missing POV" %in% unique(unlist(qMetacell(resprotFV))))
  expect_equal(SummarizedExperiment::assay(resprotFV)[which(qMetacell(resprotFV) == "Imputed MEC")], c(5, 5))
  
  ## KNN
  dataprotNAImpMEC <- QFeatures::addAssay(dataprotNA, resprotFV, name = "impMEC")
  
  resprotKNN <- wrapperImputeKNN(dataprotNAImpMEC[[2]], grp, 3)
  
  expect_true(is(resprotKNN, "SummarizedExperiment"))
  expect_true("Imputed POV" %in% unique(unlist(qMetacell(resprotKNN))))
  
  ## pa
  resprotPA <- wrapperImputePA(dataprotNA[[1]], grp)
  
  expect_true(is(resprotPA, "SummarizedExperiment"))
  expect_true("Imputed MEC" %in% unique(unlist(qMetacell(resprotPA))))
  expect_true("Imputed POV" %in% unique(unlist(qMetacell(resprotPA))))
  
  ## Det Quantile
  resprotDQ <- wrapperImputeDetQuant(dataprotNA[[1]], na.type = "Missing POV")
  
  expect_true(is(resprotDQ, "SummarizedExperiment"))
  expect_true("Imputed POV" %in% unique(unlist(qMetacell(resprotDQ))))
  expect_true("Missing MEC" %in% unique(unlist(qMetacell(resprotDQ))))
  expect_equal(SummarizedExperiment::assay(resprotDQ)[which(qMetacell(resprotDQ) == "Imputed POV")], c(2.525, 2.525, 4.175, 4.175, 3, 3.15, 3.15), tolerance = 0.001)
  
  
  # Peptide level
  ## mle
  respeptMLE <- wrapperImputeMLE(datapeptNA[[1]], grp)
  
  expect_true(is(respeptMLE, "SummarizedExperiment"))
  expect_true("Imputed POV" %in% unique(unlist(qMetacell(respeptMLE))))
  expect_true("Missing MEC" %in% unique(unlist(qMetacell(respeptMLE))))
  
  ## pa
  respeptPA <- wrapperImputePA2(datapeptNA[[1]], design_qf(datapeptNA))

  expect_true(is(respeptPA, "SummarizedExperiment"))
  expect_true("Imputed POV" %in% unique(unlist(qMetacell(respeptPA))))
  expect_true("Imputed MEC" %in% unique(unlist(qMetacell(respeptPA))))
})


test_that("check aggregation", {
  # complete dataset
  ## basic
  resAggSpe <- RunAggregation(datapept,
                        includeSharedPeptides = 'Yes_As_Specific',
                        operator = 'Mean',
                        considerPeptides = 'allPeptides',
                        adjMatrix = 'adjacencyMatrix',
                        aggregated_col = "testcol")

  expect_true(is(resAggSpe, "QFeatures"))
  expect_equal(length(resAggSpe), 2)
  expect_equal(SummarizedExperiment::assay(resAggSpe[[2]])[2, ], log2(colMeans(2^numdatapept[2:4, -5])))
  
  ## redistribution
  ### mean
  resAggREDmean <- RunAggregation(datapept,
                           includeSharedPeptides = 'Yes_Simple_Redistribution',
                           operator = 'Mean',
                           considerPeptides = 'allPeptides',
                           adjMatrix = 'adjacencyMatrix',
                           ponderation = 'Global')
  
  expect_true(is(resAggREDmean, "QFeatures"))
  expect_equal(length(resAggREDmean), 2)
  expect_equal(SummarizedExperiment::assay(resAggREDmean[[2]])[2, ], log2(colMeans(2^numdatapept[2:4, -5])))
  
  ### sum
  resAggREDsum <- RunAggregation(datapept,
                              includeSharedPeptides = 'Yes_Simple_Redistribution',
                              operator = 'Sum',
                              considerPeptides = 'allPeptides',
                              adjMatrix = 'adjacencyMatrix',
                              ponderation = 'Condition')
  
  expect_true(is(resAggREDsum, "QFeatures"))
  expect_equal(length(resAggREDsum), 2)
  
  ### median
  resAggREDmed <- RunAggregation(datapept,
                               includeSharedPeptides = 'Yes_Simple_Redistribution',
                               operator = 'Median',
                               considerPeptides = 'allPeptides',
                               adjMatrix = 'adjacencyMatrix',
                               ponderation = 'Sample')
  
  expect_true(is(resAggREDmed, "QFeatures"))
  expect_equal(length(resAggREDmed), 2)
  
  ### medianPolish
  resAggREDmedPol <- RunAggregation(datapept,
                                 includeSharedPeptides = 'Yes_Simple_Redistribution',
                                 operator = 'medianPolish',
                                 considerPeptides = 'allPeptides',
                                 adjMatrix = 'adjacencyMatrix',
                                 ponderation = 'Global')
  
  expect_true(is(resAggREDmedPol, "QFeatures"))
  expect_equal(length(resAggREDmedPol), 2)
  
  # protein with no associated peptide
  datapeptNO <- datapept
  SummarizedExperiment::rowData(datapeptNO[[1]])$adjacencyMatrix[10,6] <- 0
  
  resAggNO <- RunAggregation(datapeptNO,
                          includeSharedPeptides = 'Yes_As_Specific',
                          operator = 'Mean',
                          considerPeptides = 'allPeptides',
                          adjMatrix = 'adjacencyMatrix')
  
  expect_true(is(resAggNO, "QFeatures"))
  expect_equal(length(resAggNO), 2)
  expect_equal(nrow(resAggNO[[2]]), 6)
  expect_true(all(is.na(SummarizedExperiment::assay(resAggNO[[2]])[6, ])))
  
  # misc
  adjmat <- SummarizedExperiment::rowData(datapept[[1]])$adjacencyMatrix
  
  resProtStat <- getProteinsStats(adjmat)
  resAggMet <- aggregateMethods()
  resAggGraph <- GraphPepProt(adjmat)
  resAggExtract <- ExtractUniquePeptides(adjmat)
  resAggDetPept <- GetDetailedNbPeptides(adjmat)
  resAggDetNbPept <- GetDetailedNbPeptidesUsed(SummarizedExperiment::assay(datapept[[1]]), adjmat)
  resAggCountPept <- CountPep(adjmat)
  resAggTopn <- select_topn(SummarizedExperiment::assay(datapept[[1]]), adjmat, n = 1)
  
  expect_true(is.list(resProtStat))
  expect_equal(resProtStat$nbSpecificPeptides, 9)
  expect_equal(length(resProtStat$protMixPep), 2)
  expect_true(is.character(resAggMet))
  expect_true(is.matrix(resAggGraph))
  expect_true(is(resAggExtract, "Matrix"))
  expect_true(is.list(resAggDetPept))
  expect_true(is.list(resAggDetNbPept))
  expect_true(is(resAggCountPept, "Matrix"))
  expect_true(is(resAggTopn, "Matrix"))
})


test_that("check pept prot cc", {
  # get pept prot cc
  resCC <- getPepProtCC(SummarizedExperiment::rowData(datapept[[1]])$adjacencyMatrix)
  
  expect_true(is.list(resCC))
  
  # plot cc
  resCCplot <- plotJitter(resCC)
  
  expect_true(is.null(resCCplot))
})


test_that("check metacell", {
  resMetaNb <- GetNbTags(datapept[[1]])
  resMetaParent <- Parent('protein', c('Missing', 'Missing POV', 'Missing MEC'))
  resMetaChildren <- Children('protein', 'Missing')
  resMetaGet <- GetMetacellTags(dataprot, 1, level = 1)
  resMetaSet <- Set_POV_MEC_tags(dataprot[[1]], design_qf(dataprot)$Condition)
  resMetaSoft <- GetSoftAvailables()
  
  expect_true(is.integer(resMetaNb))
  expect_true(is.character(resMetaParent))
  expect_true(is.character(resMetaChildren))
  expect_true(is.character(resMetaGet))
  expect_true(is.data.frame(resMetaSet))
  expect_true(is.character(resMetaSoft))
})


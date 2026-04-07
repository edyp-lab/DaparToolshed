library(testthat)
library(DaparToolshed)

# prepare data
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



test_that("check limma", {
  tmpcoldata <- SummarizedExperiment::colData(dataprot)
  
  # Check and test
  resTestDesign <- testDesign(tmpcoldata[, -1])
  resCheckDesign <- checkDesign(tmpcoldata[, -1])
  resCheckCondition <- checkConditions(design_qf(dataprot)$Condition)
  
  expect_true(resTestDesign$valid)
  expect_true(resCheckDesign$valid)
  expect_true(resCheckCondition$valid)
  
  # Design 
  resMakeD <- makeDesign(tmpcoldata)
  resMakeD1 <- makeDesign1(tmpcoldata)
  resMakeD2 <- makeDesign2(tmpcoldata)
  resMakeD3 <- makeDesign3(cbind(tmpcoldata, Tech.Rep = seq_len(4)))
  resGetDlvl <- getDesignLevel(tmpcoldata)
  
  expect_equal(ncol(resMakeD), 2)
  expect_equal(ncol(resMakeD1), 2)
  expect_equal(ncol(resMakeD2), 4)
  expect_equal(ncol(resMakeD3), 4)
  expect_equal(resGetDlvl, 1)
  
  # Contrast
  resContrast <- makeContrast(resMakeD, design_qf(dataprot)$Condition)
  
  expect_equal(resContrast, "( ConditionA )/ 1 -( ConditionB )/ 1")
  
  # limma
  resLimma <- limmaCompleteTest(as.matrix(SummarizedExperiment::assay(dataprot[[1]])),
                                tmpcoldata,
                                comp.type = "anova1way")
  
  expect_true(is.list(resLimma))
  expect_equal(names(resLimma), c("logFC", "P_Value"))
})


test_that("check t-test", {
  # OnevsOne
  resTTestOvO <- compute_t_tests(dataprot, 1, "OnevsOne")
  
  expect_true(is.list(resTTestOvO))
  expect_equal(names(resTTestOvO), c("logFC", "P_Value"))
  expect_equal(resTTestOvO$logFC[, , TRUE], c(-2, 2, 0.5, -1, 2.5, 1.5, 3, 4, 3.5, 0))
  
  # OnevsAll
  # here only 2 conditions
  resTTestOvA <- compute_t_tests(dataprot, 1, "OnevsAll")
  
  expect_true(is.list(resTTestOvA))
  expect_equal(names(resTTestOvA), c("logFC", "P_Value"))
  expect_equal(resTTestOvA$logFC[, 1], -resTTestOvA$logFC[, 2])
  expect_equal(resTTestOvA$P_Value[, 1], resTTestOvA$P_Value[, 2])
})


test_that("check anova", {
  # Apply anova
  resAppAnova <- applyAnovasOnProteins(dataprot, 1)
  
  expect_true(is.matrix(resAppAnova))
  expect_true(is.list(resAppAnova[1]))
  expect_equal(length(resAppAnova), 10)
  
  # Test
  ## Omnibus
  resTestAnovaOmni <- testAnovaModels(resAppAnova, test = "Omnibus")
  
  expect_true(is.list(resTestAnovaOmni))
  expect_equal(names(resTestAnovaOmni), c("logFC", "P_Value"))
  expect_true(all(is.na(resTestAnovaOmni$logFC)))
  
  ## TukeyHSD
  resTestAnovaTukeyHSD <- testAnovaModels(resAppAnova, test = "TukeyHSD")
  
  expect_true(is.list(resTestAnovaTukeyHSD))
  expect_equal(names(resTestAnovaTukeyHSD), c("logFC", "P_Value"))
  expect_equal(resTestAnovaTukeyHSD$logFC[, , TRUE], c(2, -2, -0.5, 1, -2.5, -1.5, -3, -4, -3.5, 0))
  
  ## TukeyNoMTC
  resTestAnovaTukeyNmtc <- testAnovaModels(resAppAnova, test = "TukeyNoMTC")
  
  expect_true(is.list(resTestAnovaTukeyNmtc))
  expect_equal(names(resTestAnovaTukeyNmtc), c("logFC", "P_Value"))
  expect_equal(resTestAnovaTukeyNmtc$logFC[, , TRUE], c(2, -2, -0.5, 1, -2.5, -1.5, -3, -4, -3.5, 0))
})


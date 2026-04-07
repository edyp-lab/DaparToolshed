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

### with NA
dataprotNA <- dataprot
SummarizedExperiment::assay(dataprotNA[[1]])[which(SummarizedExperiment::assay(dataprotNA[[1]]) %in% c(1, 8))] <- NA
metacellprotna <- data.frame(metacell_S1 = c(rep("Quant. by direct id", 5), "Missing POV", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 2)),
                             metacell_S2 = c("Quant. by direct id", "Missing POV", rep("Quant. by direct id", 7), "Missing POV"),
                             metacell_S3 = c("Quant. by direct id", "Missing MEC", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 4)),
                             metacell_S4 = c("Quant. by direct id", "Missing MEC", "Quant. by direct id", "Missing POV", rep("Quant. by direct id", 3), "Missing POV", rep("Quant. by direct id", 2)))
rownames(metacellprotna) <- protname
SummarizedExperiment::rowData(dataprotNA[[1]])$qMetacell <- metacellprotna




test_that("check plot metacell", {
  grp <- design_qf(dataprotNA)$Condition
  pal <- unique(GetColorsForConditions(grp, ExtendPalette(2)))
  
  resMetaPlot1 <- metacellPerLinesHisto_HC(dataprotNA[[1]], group = grp, pattern = "Missing POV")
  resMetaPlot2 <- metacellPerLinesHistoPerCondition_HC(dataprotNA[[1]], group = grp, pattern = "Missing POV", pal = pal)
  resMetaPlot3 <- metacellHisto_HC(dataprotNA[[1]], group = grp, pattern = "Missing POV", pal = pal)
  
  expect_true(is(resMetaPlot1, "htmlwidget"))
  expect_true(is(resMetaPlot2, "htmlwidget"))
  expect_true(is(resMetaPlot3, "htmlwidget"))
})


test_that("check normalization plot", {
  resNormLOESS <- normalizeFunction(dataprot, "LOESS")
  
  # 1
  pal <- unique(GetColorsForConditions(SummarizedExperiment::colData(resNormLOESS)$Condition, 
                                       ExtendPalette(2)))
  resNormPlot1 <- compareNormalizationD_HC(SummarizedExperiment::assay(resNormLOESS[[1]]),
                                           SummarizedExperiment::assay(resNormLOESS[[2]]),
                                           conds = SummarizedExperiment::colData(resNormLOESS)$Condition,
                                           n = 10,
                                           pal = pal)
  
  expect_true(is(resNormPlot1, "htmlwidget"))
  
  # 2 
  resNormPlot2 <- plotCompareAssays(resNormLOESS, 1, 2,n = 10)
  
  expect_true(is(resNormPlot2, "htmlwidget"))
})


test_that("check imputation plot", {
  grp <- design_qf(dataprotNA)$Condition
  pal <- ExtendPalette(2, "Dark2")
  
  # 1
  resImpPlot1heat <- heatmapForMissingValues(SummarizedExperiment::assay(dataprotNA[[1]]))
  resImpPlot1mv <- wrapperMVImage(dataprotNA[[1]], group = grp, pattern = "Missing POV")
  
  expect_true(is.null(resImpPlot1heat))
  expect_true(is.null(resImpPlot1mv))
  
  
  # 2
  resImpPlot2 <- hc_mvTypePlot2(dataprotNA[[1]], group = grp, pattern = "Missing POV", pal = pal)
  
  expect_true(is(resImpPlot2, "htmlwidget"))
})


test_that("check hypothesistest plot", {
  pal <- ExtendPalette(1, "Paired")
  
  resLimma <- limmaCompleteTest(as.matrix(SummarizedExperiment::assay(dataprot[[1]])),
                                SummarizedExperiment::colData(dataprot),
                                comp.type = "OnevsOne")
  
  resHTPlot <- hc_logFC_DensityPlot(resLimma$logFC, th_logFC = 1, pal = pal)
  
  expect_true(is(resHTPlot, "htmlwidget"))
})


test_that("check volcano plot", {
  cond <- unique(design_qf(dataprot)$Condition)
  
  resLimma <- limmaCompleteTest(as.matrix(SummarizedExperiment::assay(dataprot[[1]])),
                                SummarizedExperiment::colData(dataprot),
                                comp.type = "OnevsOne")
  
  df <- data.frame(
  x = resLimma$logFC[[1]],
  y = -log10(resLimma$P_Value[[1]]),
  index = as.character(rownames(dataprot[[1]])))
  
  resVolcanoPlot <- diffAnaVolcanoplot_rCharts(df,
                                               th_pval = 0.2,
                                               th_logfc = 0.7,
                                               conditions = cond)
  
  expect_true(is(resVolcanoPlot, "htmlwidget"))
})

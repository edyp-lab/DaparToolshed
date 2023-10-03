
library(QFeatures)
library(DaparViz)

#vList <- BuildExampleDataset('QFeatures')
vList <- BuildExampleDataset('MSnbase')
#vList <- BuildExampleDataset('list')

vData <- GetVizData(vList, 1)
pal <- ExtendPalette(length(unique(vData@conds)), "Dark2")


metacellPerLinesHisto_HC(vData, pattern = "Missing POV")
metacellPerLinesHisto_HC(vData)
metacellPerLinesHisto_HC(vData, pattern = "Quantified")
metacellPerLinesHisto_HC(vData, pattern = "Quant. by direct id")
metacellPerLinesHisto_HC(vData, pattern = "Quant. by recovery")
pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
metacellPerLinesHisto_HC(vData, pattern = pattern)


metacellPerLinesHistoPerCondition_HC(vData, pattern = "Missing POV")
metacellPerLinesHistoPerCondition_HC(vData)
metacellPerLinesHistoPerCondition_HC(vData, pattern = "Quantified")
metacellPerLinesHistoPerCondition_HC(vData, pattern = "Quant. by direct id")
metacellPerLinesHistoPerCondition_HC(vData, pattern = "Quant. by recovery")
pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
metacellPerLinesHistoPerCondition_HC(vData, pattern = pattern)


metacellHisto_HC(vData, pattern = "Missing POV")
metacellHisto_HC(vData)
metacellHisto_HC(vData, pattern = "Quantified")
metacellHisto_HC(vData, pattern = "Quant. by direct id")
metacellHisto_HC(vData, pattern = "Quant. by recovery")
pattern <- c("Quantified", "Quant. by direct id", "Quant. by recovery")
metacellHisto_HC(vData, pattern = pattern)
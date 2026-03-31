
#' @title Builds sub_Exp1_R25_pept dataset
#' 
#' @description This script allows to build the sub_Exp1_R25_pept dataset.
#' This dataset is a subset of the final outcome of a quantitative mass
#' spectrometry-based proteomic analysis of two samples containing different
#' concentrations of 48 human proteins (UPS1 standard from Sigma-Aldrich)
#' within a constant yeast background (see Giai Gianetto et al. (2016) for
#' details). It contains the abundance values of the different human and
#' yeast proteins identified and quantified in these two conditions. The two
#' conditions represent the measured abundances of peptides when respectively
#' 5 fmol and 10 fmol of UPS1 human proteins were mixed with the yeast extract
#' before mass spectrometry analyses. This results in a concentration ratio of 2.
#' Three technical replicates were acquired for each condition.
#'
#' The original dataset is available as a CSV file
#' (see inst/extdata/Exp1_R25_pept_100.txt). In the latter case, 
#' the quantitative data are those of the raw intensities.
#' 
#' This dataset is a subset containing the first 100 peptides from the original 
#' dataset, which comes from: https://doi.org/10.1002/pmic.201500189
#' 
#' @rdname sub_Exp1_R25_pept
#' @examples 
#' builds_sub_Exp1_R25_pept()
#' 
#' @importFrom utils read.table
#' 
#' @export

builds_sub_Exp1_R25_pept <- function(){
  
  if (! requireNamespace("DaparToolshed", quietly = TRUE)) {
    stop("Please install DaparToolshed: BiocManager::install('DaparToolshed')")
  }
  
  exprsFile <- system.file(
    "extdata/data", 
    "Exp1_R25_pept_100.txt", 
    package = "DaparToolshed"
  )
  
  qdata = read.table(
    exprsFile, 
    header = TRUE, 
    sep = "\t", 
    as.is = TRUE
  )
  
  metadataFile <- system.file(
    "extdata/data", 
    "samples_Exp1_R25.txt", 
    package = "DaparToolshed"
  )
  
  metadata = read.table(
    metadataFile, 
    header = TRUE, 
    sep = "\t", 
    as.is = TRUE
  )
  
  
  
  subR25pept <- createQFeatures(
    data = qdata,
    file = 'subExp1R25_100_pept',
    sample = metadata,
    indQData = seq.int(from = 56, to = 61),
    keyId = "Sequence",
    indexForMetacell = seq.int(from = 43, to = 48),
    logData = TRUE,
    force.na = TRUE,
    typeDataset = 'peptide',
    parentProtId = "Protein_group_IDs",
    analysis = "Pept_data",
    description = NULL,
    processes = NULL,
    typePipeline = NULL,
    name.pipeline = NULL,
    software = 'maxquant') 
  
  save(subR25pept, file = 'data/subR25pept.RData', compress='xz')

}
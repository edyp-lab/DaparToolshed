
#' @title Builds sub_Exp1_R25_pept dataset
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
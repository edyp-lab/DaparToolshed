
#' @title Builds sub_Exp1_R25_prot dataset
#' 
#' @rdname sub_Exp1_R25_prot
#' @examples 
#' builds_sub_Exp1_R25_prot()
#' 
#' @importFrom utils read.table
#' 
#' @export

builds_sub_Exp1_R25_prot <- function(){
  
  if (! requireNamespace("DaparToolshed", quietly = TRUE)) {
    stop("Please install DaparToolshed: BiocManager::install('DaparToolshed')")
  }
  
  exprsFile <- system.file(
    "extdata/data", 
    "Exp1_R25_prot_100.txt", 
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
  
  
  
  subR25prot <- createQFeatures(
    data = qdata,
    file = 'subExp1R25_100_prot',
    sample = metadata,
    indQData = seq.int(from = 49, to = 54),
    keyId = "Majority_protein_IDs",
    indexForMetacell = seq.int(from = 36, to = 41),
    logData = TRUE,
    force.na = TRUE,
    typeDataset = 'protein',
    parentProtId = NULL,
    analysis = "Protein_data",
    description = NULL,
    processes = NULL,
    typePipeline = NULL,
    name.pipeline = NULL,
    software = 'maxquant') 
  
  
  save(subR25prot, file = 'data/subR25prot.RData', compress='xz')
  
}
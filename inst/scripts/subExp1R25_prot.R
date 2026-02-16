
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
  library(DaparToolshedData)
  data(Exp1_R25_prot)
  subR25prot <- Exp1_R25_prot[1:100]
  
  
  save(subR25prot, file = 'data/subR25prot.RData', compress='xz')
  
}
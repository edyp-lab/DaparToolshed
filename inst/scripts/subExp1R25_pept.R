
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
  library(DaparToolshedData)
  data(Exp1_R25_pept)
  subR25pept <- Exp1_R25_pept[1:100]
  
  save(subR25pept, file = 'data/subR25pept.RData', compress='xz')

}
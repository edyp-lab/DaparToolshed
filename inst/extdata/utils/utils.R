#' @export
load_dataset <- function(path_to_file){
  demo <- readRDS(path_to_file)
  
  return(demo)
}
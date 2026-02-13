

#' @title Get the last validated step before current position.
#'
#' @description This function returns the indice of the last validated step before
#' the current step.
#'
#' @param history xxx
#' @param process xxx
#' @param step.name xxx
#' @param param.name xxx
#' @param value xxx
#' @return A `integer(1)`
#'
#' @export
#' @examples
#' NULL
Add2History <- function(history, process, step.name, param.name, value){
  if (inherits(value, 'list'))
    value <- unlist(value)
  
  if (is.null(value))
    value <- NA
  
  history[nrow(history)+1, ] <- c(process, step.name, param.name, value)
  
  return(history)
}




#' @title Get the last validated step before current position.
#'
#' @description This function returns the indice of the last validated step before
#' the current step.
#'
#' @return A `integer(1)`
#'
#' @export
#' @examples
#' InitializeHistory()
#' 
InitializeHistory <- function(){
  
  history <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
    c('Process', 'Step', 'Parameter', 'Value'))
  
  return(history)
}

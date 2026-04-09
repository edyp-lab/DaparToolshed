

#' @title Add row to history
#'
#' @description This function adds a row to the history.
#'
#' @param history A `data.frame` corresponding to the current history.
#' @param process A `character(1)` corresponding to the process name.
#' @param step.name A `character(1)` corresponding to the step name.
#' @param param.name A `character(1)` corresponding to the parameter name.
#' @param value The value of the corresponding parameter.
#' 
#' @return A `data.frame` with one added row
#'
#' @export
#' @examples
#' history <- InitializeHistory()
#' Add2History(history, "Process1", "Step1", "Parameter1", "Value1")
#' 
Add2History <- function(history, process, step.name, param.name, value){
  if (inherits(value, 'list'))
    value <- unlist(value)
  
  if (is.null(value))
    value <- NA
  
  history[nrow(history)+1, ] <- c(process, step.name, param.name, value)
  
  return(history)
}




#' @title Initialize the history
#'
#' @description This function initializes the history.
#'
#' @return An empty `data.frame` with 4 columns ('Process', 'Step', 'Parameter' and 'Value')
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

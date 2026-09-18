#' @title Add row to history
#'
#' @description This function adds a row to the history.
#'
#' @param history A `data.frame` corresponding to the current history.
#' @param step A `character(1)` corresponding to the step name.
#' @param substep A `character(1)` corresponding to the substep name.
#' @param param.name A `character(1)` corresponding to the parameter name.
#' @param value The value of the corresponding parameter.
#' 
#' @return A `data.frame` with one added row
#'
#' @examples
#' history <- InitializeHistory()
#' Add2History(history, "Example step", "First sub-step", "my param", "THE value")
#' 
#' @export
#' 
Add2History <- function(history, step, substep, param.name, value){
  if (inherits(value, "list")) {
    value <- paste(names(value), unlist(value), collapse = ", ", sep = "=")
  }
  
  if (is.null(value)) {
    value <- NA
  }
  
  history[nrow(history) + 1, ] <- c(step, substep, param.name, value)
  
  return(history)
}



#' @title Initialize the history
#'
#' @description This function initializes the history.
#'
#' @return An empty `data.frame` with 4 columns ('Step', 'Substep', 'Parameter' and 'Value')
#'
#' @examples
#' InitializeHistory()
#' 
#' @export
#' 
InitializeHistory <- function() {
  history <- NULL
  history <- setNames(
    data.frame(matrix(ncol = 4, nrow = 0)),
    c("Step", "Substep", "Parameter", "Value")
  )
  
  return(history)
}

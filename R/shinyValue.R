
shinyOutput <- function(FUN, id, num, ...) {
  inputs <- character(num)
  for (i in seq_len(num)) {
    inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
  }
  inputs
}


# function for dynamic inputs in DT
shinyInput <- function(FUN, id, num, ...) {
  inputs <- character(num)
  for (i in seq_len(num)) {
    inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
  }
  inputs
}


# function to read DT inputs
shinyValue <- function(input, id, num) {
  unlist(lapply(seq_len(num), function(i) {
    value <- input[[paste0(id, i)]]
    if (is.null(value)) NA else value
  }))
}
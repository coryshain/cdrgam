#' Report message to the console
#'
#' Write a string to standard error. Useful alternative to ``message()``
#' that does not append a newline, allowing it to be used for progress
#' reporting.
#'
#' @param x A string
#' @export
report <- function(x) {
    cat(x, file = stderr())
}

#' Get the names of the smooth terms in a GAM object
#'
#' Get the names of the smooth terms in a GAM object
#' @param object A GAM object
#' @return A character vector
#' @export
get_smooth_names <- function(object) {
    return(vapply(object[["smooth"]], FUN = `[[`, FUN.VALUE = character(1), "label"))
}

#' Check whether a variable is a single string
#'
#' Check whether a variable is a single string (i.e., a character vector
#' of length 1)
#' @param x A variable
#' @return `TRUE` if `x` is a single string, `FALSE` otherwise
#' @export
is.string <- function(x){
    return(is.character(x) && length(x) == 1)
}
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
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

#' Save a CDR-GAM object to a file
#' @param obj A CDR-GAM object, with the mgcv model object
#'   stored in the `m` field and additional metadata stored
#'   in other fields.
#' @param file A file path
#' @param clean A logical value indicating whether to clean
#'   the data from the GAM object before saving, which can
#'   dramatically reduce disk and memory consumption at the
#'   loss of some functionality (see `clean_data_from_gam()`
#'   for details)
#' @param keep_model A logical value indicating whether to keep
#'   the main model matrix in the object. If `FALSE`, the
#'   `model` field of the GAM object will be removed, resulting
#'   in significant additional savings but breaking some
#'   native `mgcv` functionality like plotting or prediction
#'   without supplying a dataset. Ignored if `clean=FALSE`. See
#'   `clean_data_from_gam()` for details.
#' @param ... Additional arguments to `saveRDS()`
#' @export
save.cdrgam <- function(obj, file, clean=TRUE, keep_model=FALSE, ...) {
    if (clean) {
        obj$m <- clean_data_from_gam(obj$m, keep_model=keep_model)
    }
    saveRDS(obj, file=file, ...)
}

#' Load a CDR-GAM object from a file
#' @param file A file path
#' @return A CDR-GAM object, as saved by `save.cdrgam()`,
#'   with the mgcv model object stored in the `m` field
#'   and additional metadata stored in other fields.
#' @param ... Additional arguments to `readRDS()`
#' @export
load.cdrgam <- function(file, ...) {
    obj <- readRDS(file, ...)
    return(obj)
}

#' Clean data from a GAM object
#'
#' Clean data from a GAM object by removing model matrices
#' (training data) that `mgcv` automatically stores in the
#' object alongside the model itself.
#'
#' GAMs fitted by `mgcv` store all data involved in fitting
#' the model, including copies of transformed data associated
#' with smooth terms. Since CDR-GAMs are fitted using
#' (potentially very large and redundant) matrices of data
#' lagged over time, this can result in huge filesizes when
#' saving in the standard way. Removing this data can yield
#' large (orders-of-magnitude) reductions in object size
#' both in memory and on disk.
#'
#' This function removes the stored training data
#' while keeping the model object itself and its fitted
#' values intact. How much savings this provides depends on
#' the `keep_model` feature, which governs whether the main
#' model matrix is retained (`keep_model=TRUE`). Keeping the
#' model preserves most (all?) native `mgcv` functionality,
#' including prediction using `predict(m, newdata=NULL)` and
#' native `mgcv` plotting functions. However, the main model
#' matrix is still potentially quite large. Using
#' `keep_model=FALSE` further reduces the object size at the
#' cost of breaking all native `mgcv` functionality that
#' depends on access to the model matrix (including the
#' prediction and plotting functions mentioned above). This
#' is usually not a problem because `cdrgam` provides its
#' own prediction and plotting functions, which do not
#' rely on the model matrix.
#'
#' @param m A GAM object
#' @param keep_model A logical value indicating whether to keep
#'   the main model matrix in the object. If `FALSE`, the
#'   `model` field of the GAM object will be removed, resulting
#'   in significant additional savings but breaking some
#'   native `mgcv` functionality like plotting or prediction
#'   without supplying a dataset.
#' @export
clean_data_from_gam <- function(m, keep_model=FALSE) {
    # Remove main model matrix copies
    m$call$data <- NULL
    if (!keep_model) {
        m$model <- NULL
    }
    # Remove smooth model matrices
    for (s_ix in seq_along(m$smooth)) {
        m$smooth[[s_ix]]$ind <- NULL
        for (m_ix in seq_along(m$smooth[[s_ix]]$margin)) {
            m$smooth[[s_ix]]$margin[[m_ix]]$X <- NULL
        }
    }
    return(m)
}

#' Print the memory size of all child elements of an object
#'
#' Print the memory size of all child elements of an object.
#' Useful for debugging memory usage.
#' @param obj An object
#' @export
size_by_element <- function(obj) {
    if (is.list(obj)) {
        if (is.null(names(obj))) {
            names_ <- 1:length(obj)
        } else {
            names_ <- names(obj)
        }
        for (i in seq_along(obj)) {
            name_ <- names_[i]
            child <- obj[[i]]
            cat(paste0(name_, ': ', capture.output(print(pryr::object_size(child))), '\n'))
        }
    } else if (!is.null(names(obj))) {
        for (name_ in names(obj)) {
            child <- obj[[name_]]
            cat(paste0(name_, ': ', capture.output(print(pryr::object_size(child))), '\n'))
        }
    }
}
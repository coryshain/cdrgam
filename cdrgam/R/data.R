#' Get time windows for each observation in Y
#'
#' Compute local time windows (row start and end indices in impulse
#' matrix `X`) for each row in response matrix `Y`, given a specific window
#' size and direction. Used to constructed the time-lagged impulse matrix
#' need to fit a CDR-GAM.
#'
#' @param X A data.frame with M rows. Must contain a column matching `time_col`
#'   representing timestamps.
#' @param Y A data.frame with N rows. Must contain a column matching `time_col`
#'   representing timestamps.
#' @param series_ids A character vector of column names in `X` and `Y` that
#'   identify series groups (a series group is a sequence of observations
#'   that forms a contiguous time series). If NULL, all rows are treated
#'   as a single series.
#' @param forward A logical indicating whether to compute the window
#'   forward in time (TRUE) or backward in time (FALSE).
#' @param window_length An integer specifying the maximum number of rows
#'   in X to include in each window.
#' @param t_delta_cutoff An integer specifying the maximum time difference
#'   between the start (if `forward=TRUE`) or end (if `forward=FALSE`) of
#'   a window and the timestamp of the corresponding row in `Y`. If NULL,
#'   no cutoff is applied.
#' @param time_col A string specifying the column name in `X` and `Y` that
#'   contains timestamps.
#' @param verbose A logical indicating whether to print progress messages.
#' @return A list with two elements: 'first_obs' and 'last_obs'. Each
#'   is a vector with N elements, where N is the number of rows in `Y`.
#'   'first_obs' contains the row index in `X` of the first observation in
#'   the window for each row in `Y`. 'last_obs' contains the row index in `X`
#'   of the last observation in the window for each row in `Y`.
#' @export
get_time_windows <- function(
        X,
        Y,
        series_ids=NULL,
        forward=FALSE,
        window_length=32,
        t_delta_cutoff=NULL,
        time_col='time',
        verbose=TRUE
) {
    if (is.null(window_length)) {
        window_length <- 0
    }

    m <- nrow(X)
    n <- nrow(Y)

    X_time <- X[[time_col]]
    Y_time <- Y[[time_col]]
    if (forward) {
        X_time <- -X_time
        Y_time <- -Y_time
    }
    if (is.null(series_ids)) {
        X_ix_sort <- do.call(order, X_time)
        Y_ix_sort <- do.call(order, Y_time)
    } else {
        X_ix_sort <- do.call(order, cbind(data.frame(X[,series_ids]), X_time))
        Y_ix_sort <- do.call(order, cbind(data.frame(Y[,series_ids]), Y_time))
    }
    X_ix_sort_to_orig <- c((1:m)[X_ix_sort], m + 1)  # Add extra count for endpoint
    Y_ix_sort_to_orig <- c((1:n)[Y_ix_sort], n + 1)  # Add extra count for endpoint
    X <- X[X_ix_sort,]
    Y <- Y[Y_ix_sort,]
    X_time <- X_time[X_ix_sort]
    Y_time <- Y_time[Y_ix_sort]

    X_id_vectors <- NULL
    Y_id_vectors <- NULL

    if (!is.null(series_ids)) {
        for (col in series_ids) {
            X_id_vectors <- cbind(X_id_vectors, as.character(X[[col]]))
            Y_id_vectors <- cbind(Y_id_vectors, as.character(Y[[col]]))
        }
    } else {
        X_id_vectors <- matrix(1, nrow = nrow(X), ncol = 1)
        Y_id_vectors <- matrix(1, nrow = nrow(Y), ncol = 1)
    }

    Y_id <- Y_id_vectors[1, ]

    first_obs <- rep(1, n)
    last_obs <- rep(1, n)

    i <- 1
    j <- 1
    start <- 1
    end <- 1
    epsilon <- 1.1920929e-07

    while (i <= n) {
        if (verbose && (i == 1 || i %% 1000 == 999 || i == n)) {
            report(paste("\r", i, "/", n))
        }

        if (any(Y_id_vectors[i, ] != Y_id)) {
            start <- j
            end <- j
            Y_id <- Y_id_vectors[i, ]
        }

        if (j == 1 || any(X_id_vectors[j - 1, ] != Y_id)) {
            while (j <= m && any(X_id_vectors[j, ] != Y_id)) {
                j <- j + 1
                start <- j
                end <- j
            }
        }

        while ((j <= m) && (X_time[j] <= (Y_time[i] + epsilon)) && all(X_id_vectors[j, ] == Y_id)) {
            j <- j + 1
            while (!is.null(t_delta_cutoff) && abs(X_time[start] - Y_time[i]) > t_delta_cutoff && start < j) {
                start <- start + 1
            }
            end <- j
        }

        if (forward) {
            first_obs[i] <- max(end - 1, start)
            last_obs[i] <- start
        } else {
            first_obs[i] <- start
            last_obs[i] <- end
        }

        i <- i + 1
    }

    # Unsort X indices
    first_obs <- X_ix_sort_to_orig[first_obs]
    last_obs <- X_ix_sort_to_orig[last_obs]

    # Unsort Y indices
    first_obs <- first_obs[Y_ix_sort_to_orig]
    last_obs <- last_obs[Y_ix_sort_to_orig]

    if (forward) {
        last_obs <- last_obs + 1
        if (!is.null(window_length)) {
            last_obs <- pmin(last_obs, first_obs + window_length)
        }
    } else if (!is.null(window_length)) {
        first_obs <- pmax(first_obs, last_obs - window_length)
    }

    report('\n')

    return(list(first_obs = first_obs, last_obs = last_obs))
}

#' Get CDR data for GAM fitting
#'
#' Transform data frames of impulses/predictors `X` and responses `Y`
#' into the input format needed by a CDR-GAM model. This function should
#' be applied to both training data and held-out evaluation data before
#' passing either to the model. An important nuance is that held-out data
#' must be constructed using the same vector of factor levels as that used
#' for training, for all categorical variables in the model. Otherwise,
#' predictions will be incorrect unless all factors in both datasets
#' contain the same unique levels in the same sequence. To support this
#' constraint, the function allows the user to provide a vector of
#' factor levels for use in constructing factor-valued variables
#' (random grouping factors). These factor levels can be obtained
#' from a fitted model by calling `levels(model$var.summary[[factor]])`
#' and then supplied to this function to construct a new dataset for
#' prediction from the fitted model.
#'
#' @param X A data.frame with M rows. Must contain a column matching `time_col`
#'   representing timestamps.
#' @param Y A data.frame with N rows. Must contain a column matching `time_col`
#'   representing timestamps.
#' @param response_name A string specifying the column name of `Y` that
#'   contains the response (dependent) variable.
#' @param predictor_names A character vector of column names in `X` that
#'   contain the predictor (independent) variables.
#' @param series_ids A character vector of column names in `X` and `Y`
#'   that identify series groups (a series group is a sequence of
#'   observations that forms a contiguous time series). If NULL, all rows
#'   are treated as a single series.
#' @param ranef_names A character vector of column names in `X` and `Y` that
#'   contain random grouping factors. If NULL, no random grouping factors
#'   are included in the dataset.
#' @param ranef_levels A character vector of factor levels to use for
#'   constructing random grouping factors. If NULL, factor levels are
#'   inferred automatically. IMPORTANT: This parameter is strongly
#'   encouraged for prediction datasets to ensure that factor levels are
#'   matched to the model. Otherwise, predictions may be incorrect.
#' @param history_length An integer specifying the number of rows in `X`
#'   to include in the history (backward) window for each row in `Y`.
#' @param future_length An integer specifying the number of rows in `X`
#'   to include in the future (forward) window for each row in `Y`.
#' @param t_delta_cutoff An integer specifying the maximum time difference
#'   between the start/end of a window and the timestamp of the corresponding
#'   row in `Y`. If NULL, no cutoff is applied.
#' @param time_col A string specifying the column name in `X` and `Y` that
#'   contains timestamps.
#' @param verbose A logical indicating whether to print progress messages.
#' @return A list containing matrix-valued (NxT, where T represents
#'   `history_length`+`future_length` and N represents the number of rows in Y),
#'   windowed representations of each predictor in the input data, along with
#'   three special matrix-valued variables of the same dimensionality:
#'   - `time_col`: The timestamp of each observation in the window.
#'   - 't_delta': The time difference between the timestamp of each
#'     observation in the window and the timestamp of the corresponding
#'     row in `Y`.
#'   - 'mask': A binary mask indicating whether each observation in the
#'     window is valid (1) or invalid (0). An observation is invalid if
#'     the time difference between the timestamp of the observation and
#'     the timestamp of the corresponding row in `Y` exceeds `t_delta_cutoff`,
#'     or if the observation is outside the bounds of the input data.
#'   This list can be supplied directly to the `data` argument of the
#'   `mgcv::gam()` (and related) functions.
#' @export
get_cdr_data <- function(
        X,
        Y,
        response_name,
        predictor_names,
        series_ids=NULL,
        ranef_names=NULL,
        ranef_levels=NULL,
        history_length=16,
        future_length=0,
        t_delta_cutoff=NULL,
        time_col='time',
        verbose=TRUE
) {
    # Distinguish series groups
    if (is.null(series_ids)) {
        series_id <- '!!!SERIES_ID!!!'
        X[[series_id]] <- as.factor(1)
        Y[[series_id]] <- as.factor(1)
    } else {
        series_id <- paste(series_ids, collapse='_')
        X[[series_id]] <- as.factor(apply(X[, series_ids], 1, paste, collapse='_'))
        Y[[series_id]] <- factor(
            apply(Y[, series_ids], 1, paste, collapse='_'),
            levels=levels(X[[series_id]])
        )
    }

    # Ensure sort order
    if (is.null(series_id)) {
        X_ix_sort <- do.call(order, X[[time_col]])
        Y_ix_sort <- do.call(order, Y[[time_col]])
    } else {
        X_ix_sort <- do.call(order, cbind(data.frame(X[,series_ids]), X[[time_col]]))
        Y_ix_sort <- do.call(order, cbind(data.frame(Y[,series_ids]), Y[[time_col]]))
    }
    X <- X[X_ix_sort,]
    Y <- Y[Y_ix_sort,]

    n <- nrow(Y)
    w <- history_length + future_length

    # Initialize outputs
    out <- list()
    out[[series_id]] <- Y[[series_id]]
    out[[response_name]] <- Y[[response_name]]
    if (!is.null(ranef_names)) {
        for (ranef_name in ranef_names) {
            # This uses Simon Wood's trick for creating a factor matrix:
            # https://rdrr.io/cran/mgcv/man/linear.functional.terms.html
            ranef <- NULL
            for (i in 1:w) {
                 ranef_ <- Y[[ranef_name]]
                 ranef <- c(ranef, ranef_)
            }
            if (is.null(ranef_levels)) {
                ranef <- as.factor(ranef)
                if (get_cdr_population_code() %in% levels(ranef)) {
                    stop(paste0('Reserved keyword ', get_cdr_population_code(),
                                ' found as a level of random effect ', ranef_name, '. Please rename.'))
                }
            } else {
                ranef <- factor(ranef, levels=ranef_levels)
            }
            dim(ranef) <- c(n, w)
            if (!(get_cdr_population_code() %in% levels(ranef))) {
                levels(ranef) <- c(levels(ranef), get_cdr_population_code())
            }
            out[[ranef_name]] <- ranef
        }
    }
    out[[time_col]] <- matrix(0, n, w)
    out[['t_delta']] <- matrix(0, n, w)
    out[['mask']] <- matrix(0, n, w)
    for (predictor_name in predictor_names) {
        pred <- matrix(0, n, w)
        class(pred) <- class(X[[predictor_name]])
        out[[predictor_name]] <- pred
    }

    # Compute windows
    message('Computing and applying local time windows')
    first_obs <- NULL
    last_obs <- NULL
    if (history_length > 0) {
      windows <- get_time_windows(
        X,
        Y,
        series_ids = series_id,
        window_length = history_length,
        t_delta_cutoff = t_delta_cutoff,
        time_col=time_col,
        verbose=verbose
      )
      first_obs <- windows[['first_obs']]
      last_obs <- windows[['last_obs']]
    }
    if (future_length > 0) {
      windows <- get_time_windows(
        X,
        Y,
        series_ids = series_id,
        window_length = future_length,
        forward = TRUE,
        t_delta_cutoff = t_delta_cutoff,
        time_col=time_col,
        verbose=verbose
      )
      last_obs <- windows[['last_obs']]
      if (is.null(first_obs)) {
        first_obs <- windows[['first_obs']]
      }
    }

    # Assign windowed data to outputs
    for (i in seq_len(n)) {
        s <- first_obs[i]
        e <- last_obs[i]
        t <- e - s
        ix <- seq_len(t)
        in_ix <- s + ix - 1
        out_ix <- w - length(ix) + ix
        X_time <- X[[time_col]][in_ix]
        out[[time_col]][i, out_ix] <- X_time
        out[['t_delta']][i, out_ix] <- Y[[time_col]][i] - X_time
        out[['mask']][i, out_ix] <- 1
        for (predictor_name in predictor_names) {
            out[[predictor_name]][i, out_ix] <- X[[predictor_name]][in_ix]
        }
    }

    return(out)
}

#' Get means of numeric columns in a CDR dataset
#'
#' Compute the means of numeric predictors in a CDR dataset, excluding
#' values marked as invalid by the 'mask' column.
#'
#' @param x A list containing matrix-valued variables, including a 'mask'
#'   variable indicating whether each observation is valid (1) or invalid (0).
#'   Each other variable should be a matrix of the same dimensionality as
#'   'mask'. This list is typically the output of `get_cdr_data()`.
#' @return A list containing the means of each numeric predictor in `x`.
#' @export
get_cdr_means <- function(x) {
    out <- list()
    sel <- !x[['mask']]
    for (col in names(x)) {
        if (col != 'mask' & is.numeric(x[[col]])) {
            val <- x[[col]]
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- NA
            }
	        out[[col]] <- mean(val, na.rm=TRUE)
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- 0
            }
        }
    }
    return(out)
}

#' Get standard deviations of numeric columns in a CDR dataset
#'
#' Compute the SDs of numeric predictors in a CDR dataset, excluding
#' values marked as invalid by the 'mask' column.
#'
#' @param x A list containing matrix-valued variables, including a 'mask'
#'   variable indicating whether each observation is valid (1) or invalid (0).
#'   Each other variable should be a matrix of the same dimensionality as
#'   'mask'. This list is typically the output of `get_cdr_data()`.
#' @return A list containing the SDs of each numeric predictor in `x`.
#' @export
get_cdr_sds <- function(x) {
    out <- list()
    sel <- !x[['mask']]
    for (col in names(x)) {
        if (col != 'mask' & is.numeric(x[[col]])) {
            val <- x[[col]]
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- NA
            }
	        out[[col]] <- sd(val, na.rm=TRUE)
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- 0
            }
        }
    }
    return(out)
}

#' Get vectors of quantiles for numeric columns in a CDR dataset
#'
#' Compute quantile vectors of numeric predictors in a CDR dataset, excluding
#' values marked as invalid by the 'mask' column.
#'
#' @param x A list containing matrix-valued variables, including a 'mask'
#'   variable indicating whether each observation is valid (1) or invalid (0).
#'   Each other variable should be a matrix of the same dimensionality as
#'   'mask'. This list is typically the output of `get_cdr_data()`.
#' @param n An integer specifying the number of quantiles to compute.
#' @return A list containing the quantile vectors of each numeric predictor
#'   in `x`.
#' @export
get_cdr_quantiles <- function(x, n=101) {
    out <- list()
    sel <- !x[['mask']]
    for (col in names(x)) {
        if (col != 'mask' & is.numeric(x[[col]])) {
            val <- x[[col]]
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- NA
            }
	        out[[col]] <- quantile(val, seq(0, 1, length=n), na.rm=TRUE)
            if (NCOL(val) == NCOL(sel)) {
                val[sel] <- 0
            }
        }
    }
    return(out)
}
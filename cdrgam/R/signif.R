#' Permutation test a difference between models
#'
#' This function performs a permutation test to determine if the difference
#' between two models is significant.
#'
#' @param a0 A numeric vector of the base model's values.
#' @param a1 A numeric vector of the alternative model's values.
#' @param n_iter The number of iterations to perform.
#' @param n_tails The number of tails to use in the test.
#' @param agg The aggregation function to use. Either 'sum' or 'mean'.
#' @param statistic The statistic to use. Either 'logLik' or 'error'.
#' @param nested A logical value indicating if the test is nested. If TRUE, the
#'   test is only performed if the observed difference is positive (i.e., if
#'   the alternative model is better than the base model). Otherwise, the test
#'   is always performed.
#' @param verbose A logical value indicating if progress should be printed.
#' @return A list with the p-value, the observed difference, and the samples.
#' @export
permutation_test <- function(
        a0,
        a1,
        n_iter=10000,
        n_tails=2,
        agg='sum',
        statistic='logLik',
        nested=FALSE,
        verbose=TRUE
) {
    # Check if the input is a vector
    if (!is.vector(a0) || !is.vector(a1)) {
        stop('Input must be a vector')
    }
    # Check if the input is numeric
    if (!is.numeric(a0) || !is.numeric(a1)) {
        stop('Input must be numeric')
    }
    # Check if the input is not empty
    if (length(a0) == 0 || length(a1) == 0) {
        stop('Input must not be empty')
    }
    # Check if input lengths differ
    if (length(a0) != length(a1)) {
        stop('Input lengths must be equal')
    }
    # Check if the number of tails is 1 or 2
    if (n_tails != 1 && n_tails != 2) {
        stop('Number of tails must be 1 or 2')
    }
    # Check if the aggregation function is valid
    if (agg != 'sum' && agg != 'mean') {
        stop('Aggregation function must be sum or mean')
    }
    # Check if the nested parameter is logical
    if (!is.logical(nested)) {
        stop('nested parameter must be logical')
    }
    # Check if the mode parameter is logLik or error
    if (statistic != 'logLik' && statistic != 'error') {
        stop('mode parameter must be logLik or error')
    }
    # Check if the verbose parameter is logical
    if (!is.logical(verbose)) {
        stop('verbose parameter must be logical')
    }
    # Check if the number of iterations is positive
    if (n_iter <= 0) {
        stop('Number of iterations must be positive')
    }
    # Compute the observed difference
    agg_fn <- get(agg)
    if (statistic == 'logLik') {
        obs_diff <- agg_fn(a1) - agg_fn(a0)
    } else {
        obs_diff <- agg_fn(a0) - agg_fn(a1)
    }
    if (nested && obs_diff < 0) {
        return(list(
            p_value=1,
            observed_difference=obs_diff,
            samples=rep(0, n_iter)
        ))
    }
    if (n_tails == 2) {
        obs_diff <- abs(obs_diff)
    }

    if (verbose) {
        message(sprintf('Difference in test statistic: %s', obs_diff))
        message('Permutation testing...')
    }
    src <- cbind(a0, a1)
    n <- nrow(src)
    samples <- numeric(n_iter)
    hits <- 0
    for (i in 1:n_iter) {
        if (verbose && (i == 1 || i %% 100 == 0 || i == n_iter)) {
            report(sprintf('\r%d/%d', i, n_iter))
        }
        ix <- as.integer(runif(n) > 0.5)
        ix0 <- 1 + ix
        ix1 <- 2 - ix
        a0_ <- src[cbind(1:n, ix0)]
        a1_ <- src[cbind(1:n, ix1)]
        if (statistic == 'logLik') {
            sample_diff <- agg_fn(a1_) - agg_fn(a0_)
        } else {
            sample_diff <- agg_fn(a0_) - agg_fn(a1_)
        }
        if (n_tails == 1) {
            if (obs_diff < 0) {
                if (sample_diff <= obs_diff) {
                    hits <- hits + 1
                }
            } else {
                if (sample_diff >= obs_diff) {
                    hits <- hits + 1
                }
            }
        } else {
            sample_diff <- abs(sample_diff)
        }
        samples[i] <- sample_diff
        if (sample_diff >= obs_diff) {
            hits <- hits + 1
        }
    }
    if (verbose) {
        message('')
    }
    p_value <- (hits + 1) / (n_iter + 1)

    return(list(
        p_value=p_value,
        observed_difference=obs_diff,
        samples=samples
    ))
}

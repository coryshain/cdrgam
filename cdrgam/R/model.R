#' Compute a standard CDR-GAM formula string
#'
#' Compute a string representation of the right-hand side (RHS) of a standard
#' CDR-GAM formula from a set of preditors and optional parameters. The
#' model will contain a "rate" term (deconvolutional intercepts), additive
#' impulse response functions (IRFs) for each predictor individually, and
#' intercept, rate, and IRF terms for each random grouping factor, if
#' applicable. Model formulas can of course be written by the modeler,
#' so this is simply a convenience function for typical cases, which
#' guarantees that the resulting model is a valid CDR-GAM. More complex
#' terms (e.g., terms involving IRFs of predictor interactions) currently
#' must be added to the output by hand using string concatenation
#' (`paste()` function).
#'
#' @param preds A character vector of predictor names
#' @param k A numeric or list of numerics, the degree of the IRF splines for
#'   each predictor. If a single numeric, the same value is used for all
#'   predictors.
#' @param k_t A numeric or list of numerics, the degree of the IRF splines for
#'   the time delta variable of each IRF. If a single numeric, the same value
#'   is used for all predictors.
#' @param bs A character or list of characters, the basis function for the
#'   IRF splines for each predictor. If a single character, the same value is
#'   used for all IRFs. See `?mgcv::smooth.terms` for details.
#' @param bs_t A character or list of characters, the basis function for the
#'   IRF splines for the time delta variable of each IRF. If a single character,
#'   the same value is used for all IRFs. See `?mgcv::smooth.terms` for details.
#' @param ran_gfs A character vector of random grouping factor names. If
#'   provided, the formula will include terms for the random grouping factors.
#'   If NULL, no random effects terms are included.
#' @param random_intercept A logical, whether to include random intercept terms
#'   for the random grouping factors. Ignored if `ran_gfs` is NULL.
#' @param random_rate A logical, whether to include random rate (deconvolutional
#'   intercept) terms for the random grouping factors. Ignored if `ran_gfs` is
#'   NULL.
#' @param random_IRFs A logical, whether to include random IRF terms for the
#'   random grouping factors. Ignored if `ran_gfs` is NULL.
#' @return A string representing the RHS of the model formula
#' @export
get_formula_string <- function(
        preds,
        k=5,
        k_t=10,
        bs='cr',
        bs_t='cr',
        ran_gfs=NULL,
        random_intercept=TRUE,
        random_rate=TRUE,
        random_IRFs=TRUE
) {
    if (is.numeric(k)) {
        k_ <- list()
        for (pred in preds) {
            k_[[pred]] <- k
        }
        k <- k_
    }
    if (is.numeric(k_t)) {
        k_t_ <- list()
        for (pred in c('t_delta', preds)) {
            k_t_[[pred]] <- k_t
        }
        k_t <- k_t_
    }
    if (is.character(bs)) {
        bs_ <- list()
        for (pred in preds) {
            bs_[[pred]] <- bs
        }
        bs <- bs_
    }
    if (is.character(bs_t)) {
        bs_t_ <- list()
        for (pred in c('t_delta', preds)) {
            bs_t_[[pred]] <- bs_t
        }
        bs_t <- bs_t_
    }
    if (is.null(ran_gfs)) {
        ran_gfs <- character()
    }
    f <- paste0('~ te(t_delta, k=', k_t[['t_delta']], ', bs="', bs_t[['t_delta']], '", by=mask)')
    for (pred in preds) {
        f <- paste0(f, ' + te(t_delta, I(', pred, '), ', pred, ', k=c(', k_t[[pred]], ', ', k[[pred]],
                      '), bs=c("', bs_t[[pred]], '", "', bs[[pred]], '", "re"), by=mask)')
    }
    for (ran_gf in ran_gfs) {
        if (random_intercept) {
            f <- paste0(f, ' + te(', ran_gf, ', bs="re", by=mask)')
        }
        if (random_rate) {
            f <- paste0(f, ' + te(t_delta, ', ran_gf, ', bs=c("', bs_t[['t_delta']], '", "re"), by=mask)')
        }
        if (random_IRFs) {
            for (pred in preds) {
                f <- paste0(f, ' + te(t_delta, I(', pred, '), ', pred, ', ', ran_gf,
                             ', bs=c("', bs_t[[pred]], '", "', bs[[pred]], '", "re", "re"), by=mask)')
            }
        }
    }
    return(f)
}

#' Evaluate a CDR-GAM model on a dataset
#'
#' Evaluate a CDR-GAM model on a dataset by computing the log-likelihood of the
#' data given the model. The model must be a GAM object from the `mgcv` package
#' fitted using one of the following families: "gaussian", "gaulss", or "shash".
#'
#' @param model A CDR-GAM object (i.e., an `mgcv` model fitted with a CDR-GAM
#'   model structure on CDR-GAM-structured data).
#' @param X A list containing matrix-valued variables, including a 'mask'
#'   variable indicating whether each observation is valid (1) or invalid (0).
#'   Each other variable should be a matrix of the same dimensionality as
#'   'mask'. This list is typically the output of `get_cdr_data()`.
#' @param y A numeric vector of response values. Must be the same length as
#'   the number of rows in the matrices in `X`.
#' @return A number representing the log-likelihood of the data given the
#'   model
#' @export
evaluate <- function(
        model,
        X,
        y
) {
    # TODO: Expand to other families supported by mgcv
    family <- model$family$family
    response <- predict(model, newdata=X, type='response')
    if (family == 'gaussian') {
        mu <- response
        sigma <- sqrt(mean(residuals(model)^2))
        ll <- sum(dnorm(y, mean=mu, sd=sigma, log=TRUE))
    } else if (family == 'gaulss') {
        mu <- response[, 1]
        sigma <- 1 / response[, 2]
        ll <- sum(dnorm(y, mean=mu, sd=sigma, log=TRUE))
    } else if (family == 'shash') {
        mu <- response[, 1]
        sigma <- exp(response[, 2])  # mgcv represents on a log scale
        nu <- response[, 3]
        tau <- exp(response[, 4])  # mgcv represents on a log scale
        ll <- sum(gamlss.dist::dSHASHo2(y, mean=mu, sd=sigma, nu=nu, tau=tau, log=TRUE))
    } else {
        stop(paste0('Unknown family: ', family))
    }

    return(list(
        logLik = ll
    ))
}
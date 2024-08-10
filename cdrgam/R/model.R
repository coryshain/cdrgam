#' Compute a standard CDR-GAM formula string
#'
#' Compute a string representation of the right-hand side (RHS) of a standard
#' CDR-GAM formula from a set of preditors (or predictor sets in the case of
#' interactions) and optional parameters. The model will contain a "rate"
#' term (deconvolutional intercepts), additive impulse response functions
#' (IRFs) for each predictor set individually, and corresponding
#' intercept, rate, and IRF terms for each random grouping factor, if
#' applicable. Model formulas can of course be written by the modeler,
#' so this is simply a convenience function for typical cases, which
#' guarantees that the resulting model is a valid CDR-GAM. However, since
#' the function returns a string, the output can be edited as needed by
#' the modeler in special cases.
#'
#' @param preds A character vector of predictor names, or a list of character
#'   vectors of predictor names. If a list, each element of the list--i.e.,
#'   each (possibly singleton) set of predictors--defines a separate smooth
#'   that will be added to the model. A smooth involving multiple predictors
#'   will represent an IRF of the interaction of all predictors in the set.
#'   For example, `c('a', 'b')` will create independent IRFs for `a` and `b`,
#'   whereas `list('a', 'b', c('a', 'b'))` will additionally create an IRF
#'   of the interaction of `a` and `b`.
#' @param k A numeric or list of numerics, the degree of the IRF splines for
#'   each predictor. If a single numeric, the same value is used for all
#'   predictors.
#' @param k_t A numeric or list of numerics, the degree of the IRF splines for
#'   the time delta variable of each IRF. If a single numeric, the same value
#'   is used for all predictors.
#' @param bs A character or list of characters, the basis function for the
#'   IRF splines for each predictor. If a single character, the same value is
#'   used for all IRFs. See `?mgcv::smooth.terms` for details. If NULL,
#'   the IRF will be constrained to be linear in the predictor.
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
#' @param random_irf A logical, whether to include random IRF terms for the
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
  random_irf=TRUE
) {
    defaults <- list(
        k=5,
        k_t=10,
        bs='cr',
        bs_t='cr',
        random_intercept=TRUE,
        random_rate=TRUE,
        random_irf=TRUE
    )

    # Expand function arguments as needed into lists with full coverage of relevant input values
    expand_arg <- function(x, preds, argname, type='numeric', add_t_delta=FALSE) {
        if (is(x, type)) { # Single value of the correct type
            default <- x
            x <- list()
        } else { # By assumption, named list of values
            default <- defaults[[argname]]
        }
        if (is.null(preds)) { # No predictors to expand to
            preds <- list()
        }
        x_ <- list()
        for (pred in preds) {
            for (pred_ in pred) {  # Allows nesting
                if (pred_ %in% names(x)) {  # If a specific value is provided, use it
                    x_[[pred_]] <- x[[pred_]]
                } else {  # Otherwise, use the default
                    x_[[pred_]] <- default
                }
            }
        }
        if (add_t_delta && !('t_delta' %in% names(x_))) {  # Add a value for t_delta if required
            x_[['t_delta']] <- default
        }
        x <- x_
        return(x)
    }

    k <- expand_arg(k, preds, 'k', type='numeric')
    k_t <- expand_arg(k_t, preds, 'k_t', type='numeric', add_t_delta=TRUE)
    bs <- expand_arg(bs, preds, 'bs', type='character')
    bs_t <- expand_arg(bs_t, preds, 'bs_t', type='character', add_t_delta=TRUE)
    if (is.null(ran_gfs)) {
        ran_gfs <- character()
    }
    random_intercept <- expand_arg(random_intercept, ran_gfs, 'random_intercept', type='logical')
    random_rate <- expand_arg(random_rate, ran_gfs, 'random_rate', type='logical')
    random_irf <- expand_arg(random_irf, ran_gfs, 'random_irf', type='logical')

    # Helper function to simplify per-predictor code
    get_pred_formula <- function(
        pred,
        k,
        k_t,
        bs,
        bs_t,
        ran_gf=NULL
    ) {
        smooth_in <- 't_delta'
        linear_in <- character()
        bs_arg <- paste0('"', bs_t[['t_delta']], '"')
        k_arg <- as.character(k_t[['t_delta']])
        for (pred_ in pred) {  # Allows nesting, permitting multiple preds to interact
            if (!is.null(bs[[pred_]])) {
                smooth_in <- c(smooth_in, paste0('I(', pred_, ')'))
                bs_arg <- c(bs_arg, paste0('"', bs_t[[pred_]], '"'))
                k_arg <- c(k_arg, as.character(k[[pred_]]))
            }
            linear_in <- c(linear_in, pred_)
        }
        bs_arg <- c(bs_arg, rep('"re"', length(pred))) # Add re bases for linear terms
        if (!(is.null(ran_gf))) {
            # Add random grouping factor
            linear_in <- c(linear_in, ran_gf)
            bs_arg <- c(bs_arg, '"re"')
        }
        smooth_in <- paste(smooth_in, collapse=', ')
        linear_in <- paste(linear_in, collapse=', ')
        bs_arg <- paste(bs_arg, collapse=', ')
        k_arg <- paste(k_arg, collapse=', ')
        pred_f <- paste0(' + te(', smooth_in, ', ', linear_in, ', k=c(', k_arg, '), bs=c(', bs_arg, '), by=mask)')
        return(pred_f)
    }

    f <- paste0('~ te(t_delta, k=', k_t[['t_delta']], ', bs="', bs_t[['t_delta']], '", by=mask)')
    for (pred in preds) {
        pred_f <- get_pred_formula(pred, k, k_t, bs, bs_t)
        f <- paste0(f, pred_f)
    }
    for (ran_gf in ran_gfs) {
        if (random_intercept[[ran_gf]]) {
            f <- paste0(f, ' + te(', ran_gf, ', bs="re", by=mask)')
        }
        if (random_rate[[ran_gf]]) {
            f <- paste0(f, ' + te(t_delta, ', ran_gf, ', bs=c("', bs_t[['t_delta']], '", "re"), by=mask)')
        }
        if (random_irf[[ran_gf]]) {
            for (pred in preds) {
                pred_f <- get_pred_formula(pred, k, k_t, bs, bs_t, ran_gf=ran_gf)
                f <- paste0(f, pred_f)
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
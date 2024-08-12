fit_cdrnn <- function(
        cfg,
        model_name
) {
    if (!('data' %in% names(cfg))) {
        stop('Required field "data" missing from config')
    }
    if (!(model_name %in% names(cfg))) {
        stop(paste0('Required field "', model_name, '" missing from config'))
    }
    data_cfg <- cfg$data
    model_cfg <- cfg[[model_name]]

    # Load data
    X_train <- data_cfg$X_train
    if (is.string(X_train)) {  # Provided as a filepath
        X_train <- get(load(X_train))
    }
    Y_train <- data_cfg$Y_train
    if (is.string(Y_train)) {  # Provided as a filepath
        Y_train <- get(load(Y_train))
    }
    train <- get_cdr_data(
        data_cfg$X_train,
        data_cfg$Y_train,
        response_name=model_cfg$response_name,
        predictor_names=cfg$predictor_names,
        ranef_names=cfg$predictor_names,
        history_length=model_cfg$history_length,
        future_length=model_cfg$future_length,
        t_delta_cutoff=model_cfg$t_delta_cutoff
    )
    means <- get_cdr_means(train)
    sds <- get_cdr_sds(train)
    quantiles <- get_cdr_quantiles(train)

    # Construct formula
    f_fixed <- get_formula_string(cfg$PRED_COLS, k=cfg$K, k_t=cfg$K_T, bs=cfg$BS, bs_t=cfg$BS_T)
    f_subj <- get_formula_string(cfg$PRED_COLS, k=cfg$K, k_t=cfg$K_T, bs=cfg$BS, bs_t=cfg$BS_T, ran_gf='subject')
    f_item <- get_formula_string(NULL, k=cfg$K, k_t=cfg$K_T, bs=cfg$BS, bs_t=cfg$BS_T, ran_gf='item', use_rate=FALSE)
    f <- paste(f_fixed, f_subj, f_item, sep=' + ')

    # Fit model
    if (cfg$FIT) {
        message('Fitting')
        message(paste0('  Model: ', f))
        m <- mgcv::gam(
            as.formula(paste(cfg$DV, f, sep=' ~ ')),
            data=train,
            gamma=0.1,
            optimizer=c('outer', 'optim'),
            drop.unused.levels=FALSE
        )
        out <- list(
            m=m,
            means=means,
            sds=sds,
            quantiles=quantiles
        )
        save(out, file=paste0(cfg$MODEL_NAME, '.RData'))
    } else {
        message('Loading')
        out <- get(load(paste0(cfg$MODEL_NAME, '.RData')))
        m <- out$m
        means <- out$means
        sds <- out$sds
        quantiles <- out$quantiles
    }

    print(summary(m))
    sink(paste0(cfg$MODEL_NAME, '_eval.txt'))
    print(summary(m))
    sink()

    if (cfg$EVAL) {
        ll_train <- evaluate(m, train, train[[cfg$DV]])
    }
}

get_columns_from_cfg <- function(
        cfg
) {
    #TODO: Implement
}

#' Compute a standard CDR-GAM formula string
#'
#' Compute a string representation of the right-hand side (RHS) of a standard
#' CDR-GAM formula from a set of preditors (or predictor sets in the case of
#' interactions) and optional parameters. Model formulas can of course be
#' hand written, so this is simply a convenience function for typical cases,
#' which guarantees that the resulting model is a valid CDR-GAM. However,
#' since the function returns a string, the output can be edited as needed by
#' the modeler in special cases. The function returns either a set of fixed
#' effects terms or a set of random effects terms (if `ran_gf` is provided).
#' Thus, to construct mixed models, simply call this function multiple times,
#' once for the fixed effects specification and once for each random effects
#' specification, and concatenate the results with a `+` separator.
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
#' @param use_intercept A logical, whether to include an intercept term
#' @param use_rate A logical, whether to a rate (deconvolutional
#'   intercept) term
#' @param ran_gf A string containing the name of a random grouping factor. If
#'   provided, the formula will interact all terms with the random grouping factor.
#'   If NULL, the string will represent fixed effects terms.
#' @param t_delta_col A string specifying the name of the column
#'   containing the difference in time between impulses and response.
#' @param mask_col A string specifying the name of the column
#'   containing the mask over valid timepoints.
#' @return A string representing the RHS of the model formula
#' @export
get_formula_string <- function(
        preds,
        k=5,
        k_t=10,
        bs='cr',
        bs_t='cr',
        use_intercept=TRUE,
        use_rate=TRUE,
        ran_gf=NULL,
        t_delta_col='t_delta',
        mask_col='mask'
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
        if (is(x, type) | is.null(x)) { # Single value of the correct type
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
        if (add_t_delta && !(t_delta_col %in% names(x_))) {  # Add a value for t_delta if required
            x_[[t_delta_col]] <- default
        }
        x <- x_
        return(x)
    }

    k <- expand_arg(k, preds, 'k', type='numeric')
    k_t <- expand_arg(k_t, preds, 'k_t', type='numeric', add_t_delta=TRUE)
    bs <- expand_arg(bs, preds, 'bs', type='character')
    bs_t <- expand_arg(bs_t, preds, 'bs_t', type='character', add_t_delta=TRUE)

    # Helper function to simplify per-predictor code
    get_pred_formula <- function(
        pred,
        k,
        k_t,
        bs,
        bs_t,
        ran_gf=NULL
    ) {
        smooth_in <- t_delta_col
        linear_in <- character()
        bs_arg <- paste0('"', bs_t[[t_delta_col]], '"')
        k_arg <- as.character(k_t[[t_delta_col]])
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
        pred_f <- paste0(
            'te(', smooth_in, ', ', linear_in, ', k=c(', k_arg, '), bs=c(', bs_arg, '), by=', mask_col, ')'
        )
        return(pred_f)
    }

    f <- character()
    if (use_intercept)  {
        if (is.null(ran_gf)) {
            f <- c(f, '1')
        } else {
            f <- c(f, paste0('te(', ran_gf, ', bs="re", by=', mask_col, ')'))
        }
    }
    if (use_rate) {
        if (is.null(ran_gf)) {
            f <- c(f, paste0(
                'te(',t_delta_col, ', k=', k_t[[t_delta_col]], ', bs="', bs_t[[t_delta_col]], '", by=', mask_col, ')'
            ))
        } else {
            f <- c(f, paste0(
                'te(', t_delta_col, ', ', ran_gf, ', k=', k_t[[t_delta_col]],
                ', bs=c("', bs_t[[t_delta_col]], '", "re"), by=', mask_col, ')'
            ))
        }
    }
    for (pred in preds) {
        f <- c(f, get_pred_formula(pred, k, k_t, bs, bs_t, ran_gf=ran_gf))
    }
    f <- paste(f, collapse=' + ')
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
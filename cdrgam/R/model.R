#' High-level CDR-GAM wrapper for fitting, evaluating, and plotting
#'
#' A convenience function for fitting, evaluating, and plotting a CDR-GAM model
#' using a configuration file. The function will fit the model if `fit` is
#' `TRUE`, evaluate the model if `eval_partition` is provided, and visualize
#' the model if `plot` is `TRUE`.
#'
#' @param cfg A configuration object or a string with the path to a YAML file
#' @param model_name A string with the name of the model to fit
#' @param fit A logical, whether to fit the model
#' @param eval_partition A string or a list of strings with the names of the
#'    partition(s) (e.g., train, val, test) on which to evaluate
#' @param plot A logical, whether to visualize the model estimates
#' @param plot_cfg A configuration object or a string with the path to a YAML
#'    file with additional plot configuration settings (useful for overriding
#'    defaults used in the main config)
#' @param overwrite A logical, whether to refit and overwrite an existing
#'    model file
#' @export
main <- function(
        cfg,
        model_name,
        fit=TRUE,
        eval_partition=c('train', 'val'),
        plot=TRUE,
        plot_cfg=NULL,
        overwrite=FALSE
) {
    if (is.string(cfg)) {  # Provided as a filepath
        cfg <- get_cfg(cfg)
    }
    validate_cfg(cfg, model_name)
    if (fit) {
        fit_cdrnn(cfg, model_name=model_name, overwrite=overwrite)
    }
    if (!is.null(eval_partition)) {
        for (eval_partition_ in eval_partition) {
            X_part <- paste0('X_', eval_partition_)
            Y_part <- paste0('Y_', eval_partition_)
            if (!(is.null(cfg$data[[X_part]]) || is.null(cfg$data[[Y_part]]))) {
                evaluate_cdrnn(cfg, model_name=model_name, eval_partition=eval_partition_)
            } else {
                message(paste0('No data provided for partition "', eval_partition_, '". Skipping evaluation.'))
            }
        }
    }
    if (plot) {
        plot_cdrnn(cfg, model_name=model_name, plot_cfg=plot_cfg)
    }
}

#' Fit a CDR-GAM model
#'
#' Fit a CDR-GAM model based on a configuration file. The function will load
#' the data, fit the model, and save the model to disk. If the model already
#' exists, the function will skip fitting unless the `overwrite` flag is set.
#'
#' @param cfg A configuration object or a string with the path to a YAML file
#' @param model_name A string with the name of the model to fit
#' @param overwrite A logical, whether to overwrite an existing model
#' @export
fit_cdrnn <- function(
        cfg,
        model_name,
        overwrite=FALSE
) {
    message(rep('=', 80))
    message('FITTING CDR-GAM MODEL')

    # Get cfg
    if (is.string(cfg)) {  # Provided as a filepath
        cfg <- get_cfg(cfg)
    }
    validate_cfg(cfg, model_name)
    data_cfg <- cfg$data
    model_cfg <- cfg$models[[model_name]]
    name_map <- cfg$name_map

    # Get paths
    output_dir <- file.path(cfg$output_dir, model_name)
    model_path <- file.path(output_dir, 'model.RData')
    summary_path <- file.path(output_dir, 'summary.txt')
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }

    # Load data
    X <- data_cfg$X_train
    sep <- data_cfg$sep
    if (is.string(X)) {  # Provided as a filepath
        X <- read.csv(X, sep=sep, header=TRUE)
    }
    Y <- data_cfg$Y_train
    if (is.string(Y)) {  # Provided as a filepath
        Y <- read.csv(Y, sep=sep, header=TRUE)
    }
    response_name <- all.vars(as.formula(paste('~', model_cfg$response)))
    predictor_names <- get_columns_from_cfg(model_cfg$formula)
    ranef_names <- get_ranefs_from_cfg(model_cfg$formula)
    other_names <- get_others_from_cfg(model_cfg$formula)
    cdrgam_data <- get_cdr_data(
        X,
        Y,
        response_name=response_name,
        predictor_names=predictor_names,
        series_ids=data_cfg$series_ids,
        ranef_names=ranef_names,
        other_names=other_names,
        filters=data_cfg$filters,
        history_length=data_cfg$history_length,
        future_length=data_cfg$future_length,
        t_delta_cutoff=data_cfg$t_delta_cutoff
    )
    means <- get_cdr_means(cdrgam_data)
    sds <- get_cdr_sds(cdrgam_data)
    quantiles <- get_cdr_quantiles(cdrgam_data)

    f <- get_formula_from_config(model_cfg$formula, model_cfg$response)
    response_params <- names(f)

    # Fit model
    fit <- overwrite || !file.exists(model_path)
    if (fit) {
        message('  Fitting')
        message(paste0('    n: ', length(cdrgam_data[[response_name]])))
        message('    Model:')
        for (response_param in response_params) {
            f_str <- deparse(f[[response_param]])
            f_str <- gsub("^\\s+|\\s+$", "", f_str)
            f_str <- paste(f_str, collapse=' ')
            message(paste0('      ', response_param, ': ', f_str, sep=''))
        }
        if (length(f) == 1) {
            f <- f[[1]]
        }
        fit_kwargs <- list(f, data=cdrgam_data, drop.unused.levels=FALSE)
        keys <- names(model_cfg)
        keys <- keys[keys %in% names(GLOBAL.CDRGAM$model)]
        fit_kwargs <- c(fit_kwargs, model_cfg[keys])
        m <- do.call(mgcv::gam, fit_kwargs)
        model <- list(
            m=m,
            means=means,
            sds=sds,
            quantiles=quantiles,
            response_params=response_params,
            name_map=name_map,
            cfg=cfg
        )
    } else {
        message('  Model already exists, skipping fitting and reloading')
        model <- get(load(model_path))
        m <- model$m
        model <- list(
            m=m,
            means=means,
            sds=sds,
            quantiles=quantiles,
            response_params=response_params,
            name_map=name_map,
            cfg=cfg
        )
    }
    message('  Saving')
    save(model, file=model_path)

    print(summary(m))
    sink(summary_path)
    print(summary(m))
    sink()
}

#' Evaluate a CDR-GAM model
#'
#' Evaluate a CDR-GAM model based on a configuration file. The function will
#' load the model, compute the log-likelihood of the data given the model, and
#' save the results to disk.
#'
#' @param cfg A configuration object or a string with the path to a YAML file
#' @param model_name A string with the name of the model to evaluate
#' @param eval_partition A string with the name of the partition on which to
#'   evaluate the model (e.g., 'val', 'test')
#' @param extra_cols A logical, whether to include all columns from `Y` in
#'   the output
#' @export
evaluate_cdrnn <- function(
        cfg,
        model_name,
        eval_partition="val",
        extra_cols=FALSE
) {
    message(rep('=', 80))
    message('EVALUATING CDR-GAM MODEL')
    if (is.string(cfg)) {  # Provided as a filepath
        cfg <- get_cfg(cfg)
    }
    validate_cfg(cfg, model_name)
    data_cfg <- cfg$data
    model_cfg <- cfg$models[[model_name]]

    # Load data
    sep <- data_cfg$sep
    X_part <- paste0('X_', eval_partition)
    Y_part <- paste0('Y_', eval_partition)
    X <- read.csv(data_cfg[[X_part]], sep=sep, header=TRUE)
    Y <- read.csv(data_cfg[[Y_part]], sep=sep, header=TRUE)
    response_name <- all.vars(as.formula(paste('~', model_cfg$response)))
    predictor_names <- get_columns_from_cfg(model_cfg$formula)
    ranef_names <- get_ranefs_from_cfg(model_cfg$formula)
    other_names <- get_others_from_cfg(model_cfg$formula)
    cdrgam_data <- get_cdr_data(
        X,
        Y,
        response_name=response_name,
        predictor_names=predictor_names,
        series_ids=data_cfg$series_ids,
        ranef_names=ranef_names,
        other_names=other_names,
        filters=data_cfg$filters,
        history_length=data_cfg$history_length,
        future_length=data_cfg$future_length,
        t_delta_cutoff=data_cfg$t_delta_cutoff
    )
    n <- length(cdrgam_data[[response_name]])

    output_dir <- file.path(cfg$output_dir, model_name)
    model_path <- file.path(output_dir, 'model.RData')
    output_path <- file.path(output_dir, paste0('output_', eval_partition, '.csv'))
    eval_path <- file.path(output_dir, paste0('eval_', eval_partition, '.txt'))

    message('  Loading model')
    model <- get(load(model_path))
    m <- model$m
    # Use model.matrix to handle responses containing transformations
    obs <- model.matrix(as.formula(paste('~', response_name)), data=cdrgam_data)[, response_name]
    message('  Computing likelihoods')
    lls <- evaluate(m, cdrgam_data, obs)$logLik
    ll <- sum(lls)

    message(paste('  Model:', model_name))
    message(paste('  partition:', eval_partition))
    message(paste('  n:', n))
    message(paste('  logLik:', ll))
    sink(eval_path)
    cat(paste('Model:', model_name, '\n'))
    cat(paste('partition:', eval_partition, '\n'))
    cat(paste('n:', n, '\n'))
    cat(paste('logLik:', ll, '\n'))
    sink()

    if (extra_cols) {
        out <- Y
        out[['cdrgam.logLik']] <- lls
    } else {
        out <- data.frame(cdrgam.logLik=lls)
    }
    write.table(out, output_path, row.names=FALSE, col.names=TRUE, sep=sep)
    return(NULL)
}

#' Plot a CDR-GAM model
#'
#' Visualize the estimates of a CDR-GAM model based on a configuration file.
#' The function will load the model, produce plots, and save the plots to disk.
#'
#' @param cfg A configuration object or a string with the path to a YAML file
#' @param model_name A string with the name of the model to plot
#' @param plot_cfg A configuration object or a string with the path to a YAML
#'   file with additional plot configuration settings (useful for overriding
#'   defaults used in the main config)
#' @export
plot_cdrnn <- function(
        cfg,
        model_name,
        plot_cfg=NULL
) {
    message(rep('=', 80))
    message('PLOTTING CDR-GAM MODEL')
    if (is.string(cfg)) {  # Provided as a filepath
        cfg <- get_cfg(cfg)
    }
    validate_cfg(cfg, model_name)
    model_cfg <- cfg$models[[model_name]]
    response_name <- model_cfg$response
    f <- get_formula_from_config(model_cfg$formula, model_cfg$response)
    response_params <- names(f)
    if (is.null(plot_cfg)) {
        plot_cfg <- cfg$plot
    } else {
        if (is.string(plot_cfg)) {
            plot_cfg <- get_plot_cfg(plot_cfg)
        }
        plot_defaults <- cfg$plot
        keys <- names(plot_defaults)
        keys <- keys[!(keys %in% names(plot_cfg))]
        plot_cfg[keys] <- plot_defaults[keys]
    }
    t_delta_ref <- plot_cfg$t_delta_ref
    if (is.null(t_delta_ref)) {
        t_delta_ref <- 0
    }

    output_dir <- file.path(cfg$output_dir, model_name)
    message('  Loading model')
    model <- get(load(file.path(output_dir, 'model.RData')))
    m <- model$m
    means <- model$means
    sds <- model$sds
    quantiles <- model$quantiles
    if ('name_map' %in% names(plot_cfg)) {
        name_map <- plot_cfg$name_map
    } else {
        name_map <- model$name_map
    }

    message('  Generating plots')
    for (i in seq_along(response_params)) {
        plot_response_name <- response_name
        for (key in names(name_map)) {
            if (grepl(key, plot_response_name, fixed=TRUE)) {
                plot_response_name <- name_map[[key]]
                break
            }
        }
        response_param <- response_params[[i]]
        plot_response_name <- paste(plot_response_name, response_param, sep=', ')
        plot_path <- file.path(output_dir, paste0('irf_univariate_', response_name, '_', response_param, '.png'))
        if (is.null(plot_cfg$t_delta_xlim)) {
            xlim <- c(0, 1)
        } else {
            xlim <- plot_cfg$t_delta_xlim
        }
        p <- plot_irfs(
            m,
            means=means,
            sds=sds,
            xlim=xlim,
            irf_name_map=name_map,
            response_param=i,
            ylabel=plot_response_name,
            legend=plot_cfg$legend
        )
        ggplot2::ggsave(plot_path, plot=p, width=plot_cfg$width, height=plot_cfg$height, scale=plot_cfg$scale)
        plot_path <- file.path(
            output_dir,
            paste0('curvature_', response_name, '_', response_param, '_at_delay', t_delta_ref, '.png')
        )
        p <- plot_curvature(
            m,
            means=means,
            sds=sds,
            quantiles=quantiles,
            range=0.9,
            irf_name_map=name_map,
            response_param=i,
            t_delta_ref=t_delta_ref,
            ylabel=plot_response_name,
            legend=plot_cfg$legend
        )
        ggplot2::ggsave(plot_path, plot=p, width=plot_cfg$width, height=plot_cfg$height, scale=plot_cfg$scale)
    }
}

#' Validate a CDR-GAM configuration
#'
#' Validate a CDR-GAM configuration object to ensure that it contains the
#' required fields.
#'
#' @param cfg A configuration object or a string with the path to a YAML file
#' @param model_name A string with the name of the model to validate
#' @export
validate_cfg <- function(cfg, model_name) {
    if (is.string(cfg)) {  # Provided as a filepath
        cfg <- get_cfg(cfg)
    }
    if (!('data' %in% names(cfg))) {
        stop('Required field "data" missing from config')
    }
    if (!('models' %in% names(cfg))) {
        stop(paste0('Required field "models" missing from config'))
    }
    if (!(model_name %in% names(cfg$models))) {
        stop(paste0('Required field "', model_name, '" missing from config'))
    }
}

#' Get a CDR-GAM formula from a configuration
#'
#' Get a CDR-GAM formula from a formula configuration object. The
#' function will return a list of formulas, one for each response parameter.
#'
#' @param formula_cfg A formula configuration object
#' @param response_name A string with the name of the response variable
#' @return A list of formulas
#' @export
get_formula_from_config <- function(
        formula_cfg,
        response_name
) {
    f <- list()
    for (formula_name in names(formula_cfg)) {
        formula <- formula_cfg[[formula_name]]
        f_ <- list()
        for (formula_ in formula) {
            f_ <- c(f_, do.call(get_formula_string, formula_))
        }
        f_ <- paste(f_, collapse=' + ')
        f_ <- paste('~', f_)
        f[[formula_name]] <- f_
    }
    f[[1]] <- paste(response_name, f[[1]])
    response_params <- names(f)
    for (key in response_params) {
        f[[key]] <- as.formula(f[[key]], env=baseenv())  # Setting env avoids exploding memory footprint for closures
    }
    return(f)
}

get_columns_from_cfg <- function(
        formula_cfg
) {
    columns <- NULL
    for (formula in formula_cfg) {  # Iterate over response parameters
        for (formula_ in formula) {  # Iterate over formula components
            for (irf in formula_$irfs) {  # Iterate over IRFs within formula component
                if (!is.null(names(irf)) && names(irf) == c('inputs', 'impulses')) {
                    # Flatten IRFs that distinguish inputs and impulses
                    irf <- c(irf[['inputs']], irf[['impulses']])
                }
                for (pred in irf) {  # Iterate over predictors within IRF
                    new_columns <- all.vars(as.formula(paste('~', pred)))
                    columns <- c(columns, new_columns)
                }
            }
        }
    }
    columns <- unique(columns)
    return(columns)
}


get_ranefs_from_cfg <- function(
        formula_cfg
) {
    columns <- NULL
    for (formula in formula_cfg) {  # Iterate over response parameters
        for (formula_ in formula) {  # Iterate over formula components
            columns <- c(columns, formula_$ran_gf)
        }
    }
    columns <- unique(columns)
    return(columns)
}

get_others_from_cfg <- function(
        formula_cfg
) {
    columns <- NULL
    for (formula in formula_cfg) {  # Iterate over response parameters
        for (formula_ in formula) {  # Iterate over formula components
            others <- c(columns, formula_$others)
            for (other in others) {  # Get all variables in the formula component
                new_columns <- all.vars(as.formula(paste('~', other)))
                columns <- c(columns, new_columns)
            }
        }
    }
    columns <- unique(columns)
    return(columns)
}

#' Compute a standard CDR-GAM formula string
#'
#' Compute a string representation of the right-hand side (RHS) of a standard
#' CDR-GAM formula from a set of preditors (or predictor sets in the case of
#' interactions) and optional parameters.
#'
#' Model formulas can of course be hand written, so this is simply a
#' convenience function for typical cases, which guarantees that the
#' resulting model is a valid CDR-GAM. However, since the function returns
#' a string, the output can be edited as needed by the modeler in special
#' cases. The function returns either a set of fixed effects terms or a set
#' of random effects terms (if `ran_gf` is provided). Thus, to construct
#' mixed models, simply call this function multiple times, once for the
#' fixed effects specification and once for each random effects specification,
#' and concatenate the results with a `+` separator.
#'
#' This leverages the `linear.functional.terms` feature of `mgcv`, which
#' *sums* the results of applying the same smooth to multiple values of a
#' predictor for a given datapoint, weighting the elements in the
#' summation by the value of the variable supplied to the `by` term of the
#' smooth. When such a smooth is parameterized by the time offset *t_delta*
#' and `by` is the impulse (predictor), a `linear.functional.term` becomes
#' a discrete convolution over time, where the smooth represents the
#' impulse response. If the smooth is also parameterized by the impulse
#' itself, then the IRF can be non-linear in the impulse.
#'
#' For example, the following `mgcv` term defines an IRF that is linear
#' in predictor `a`, assuming that `t_delta` and `a` are both supplied
#' as matrices as required by `linear.functional.terms` (see
#' `get_cdr_data()` for details):
#'     `te(t_delta, k=10, bs='cr', by=a)`
#' By contrast, the following defines an IRF that is nonlinear in `a`
#' (because `a` also parameterizes the smooth itself):
#'     `te(t_delta, a, k=c(10, 5), bs=c('cr', 'cr'), by=a)`
#'
#' @param irfs A character vector of predictor names, or a list of character
#'   vectors of predictor names, or a list of lists of character vectors of
#'   predictor names. If a list, each element of the list--i.e., each
#'   (possibly singleton) set of predictors--defines a separate smooth
#'   that will be added to the model. If a list of lists, each element must
#'   be a pair (list of length 2) in which the first subelement (named
#'   "inputs") is a character vector of predictors to include in the
#'   inputs to the IRF, and the second subelement (named "impulses") is
#'   a character vector of impulses to convolve using the IRF.
#'   A smooth involving multiple predictors will represent an IRF of the
#'   interaction of all predictors in the set. For example, `c('a', 'b')` will
#'   create independent IRFs for `a` and `b`, whereas
#'   `list('a', 'b', c('a', 'b'))` will additionally create an IRF of the
#'  interaction of `a` and `b`.
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
#' @param others A character vector of additional `mgcv` GAM terms to include
#'   in the model without convolving them. All predictors in these terms must
#'   be found in the response data (`Y`), since only the response is guaranteed
#'   to be conformable without convolution. To distinguish these terms (encoded
#'   as vectors) from potentially identically-named convolved predictors (encoded
#'   as matrices), the column must be suffixed with '_Y' in the term. Thus, to
#'   add a term that fits a non-convolutional smooth to a predictor called `time`,
#'   the term should use `time_Y`, e.g.:
#'       `te(time_Y, k=10, bs='cr')`
#' @param t_delta_col A string specifying the name of the column
#'   containing the difference in time between impulses and response.
#' @param mask_col A string specifying the name of the column
#'   containing the mask over valid timepoints.
#' @return A string representing the RHS of the model formula
#' @export
get_formula_string <- function(
        irfs=NULL,
        k=5,
        k_t=10,
        bs='cr',
        bs_t='cr',
        use_intercept=TRUE,
        use_rate=TRUE,
        ran_gf=NULL,
        others=NULL,
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
    expand_arg <- function(x, variables, argname, type='numeric', add_t_delta=FALSE) {
        if (is(x, type) | is.null(x)) { # Single value of the correct type
            default <- x
            x <- list()
        } else { # By assumption, named list of values
            default <- defaults[[argname]]
        }
        if (is.null(variables)) { # No predictors to expand to
            variables <- list()
        }
        x_ <- list()
        for (variable in variables) {
            for (variable_ in variable) {  # Allows nesting
                sel <- variable_ %in% names(x)
                keys_in <- variable_[sel]
                keys_out <- variable_[!sel]
                x_[keys_in] <- x[keys_in]  # If a specific value is provided, use it
                x_[keys_out] <- list(default)  # Otherwise, use the default
            }
        }
        if (add_t_delta && !(t_delta_col %in% names(x_))) {  # Add a value for t_delta if required
            x_[t_delta_col] <- list(default)
        }
        x <- x_
        return(x)
    }

    k <- expand_arg(x=k, variables=irfs, argname='k', type='numeric')
    k_t <- expand_arg(x=k_t, variables=irfs, argname='k_t', type='numeric', add_t_delta=TRUE)
    bs <- expand_arg(x=bs, variables=irfs, argname='bs', type='character')
    bs_t <- expand_arg(x=bs_t, variables=irfs, argname='bs_t', type='character', add_t_delta=TRUE)

    # Helper function to simplify per-predictor code
    get_irf_formula <- function(
            impulses,
            k,
            k_t,
            bs,
            bs_t,
            inputs=NULL,
            ran_gf=NULL
    ) {
        smooth_in <- t_delta_col
        ranef_in <- character()
        by_in <- character()
        bs_arg <- paste0('"', bs_t[[t_delta_col]], '"')
        k_arg <- as.character(k_t[[t_delta_col]])
        if (is.null(inputs)) {  # Assume identity between IRF inputs and impulses
            inputs <- impulses
        }
        for (input in inputs) {  # Allows nesting, permitting multiple preds to interact
            if (!is.null(bs[[input]])) {
                # smooth_in <- c(smooth_in, paste0('I(', input, ')'))
                smooth_in <- c(smooth_in, input)
                bs_arg <- c(bs_arg, paste0('"', bs[[input]], '"'))
                k_arg <- c(k_arg, as.character(k[[input]]))
            }
        }
        for (impulse in impulses) {  # Allows nesting, permitting multiple preds to interact
            by_in <- c(by_in, impulse)
        }
        # bs_arg <- c(bs_arg, rep('"re"', length(impulses))) # Add re bases for linear terms
        if (!(is.null(ran_gf))) {
            # Add random grouping factor
            for (ran_gf_ in ran_gf) {  # Allows nesting, for crossed random grouping factors
                ranef_in <- c(ranef_in, ran_gf)
                bs_arg <- c(bs_arg, '"re"')
            }
        }
        smooth_in <- paste(smooth_in, collapse=', ')
        ranef_in <- paste(ranef_in, collapse=', ')
        if (is.null(ran_gf)) {
            preds <- smooth_in
        } else {
            preds <- paste(smooth_in, ranef_in, sep=', ')
        }
        by_in <- paste(by_in, collapse='*')
        bs_arg <- paste(bs_arg, collapse=', ')
        k_arg <- paste(k_arg, collapse=', ')
        pred_f <- paste0(
            'te(', preds, ', k=c(', k_arg, '), bs=c(', bs_arg, '), by=', by_in, ')'
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
    for (irf in irfs) {
        if (!is.null(names(irf)) && names(irf) == c('inputs', 'impulses')) {
            impulses <- irf[['impulses']]
            inputs <- irf[['inputs']]
        } else {
            impulses <- irf
            inputs <- NULL
        }
        f <- c(f, get_irf_formula(impulses, k, k_t, bs, bs_t, inputs=inputs, ran_gf=ran_gf))
    }
    for (other in others) {
        f <- c(f, other)
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
#' @return A numeric vector representing the log-likelihood of each datapoint
#'   given the model. Dataset likelihood can be computed by summing this vector.
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
        ll <- dnorm(y, mean=mu, sd=sigma, log=TRUE)
    } else if (family == 'gaulss') {
        mu <- response[, 1]
        sigma <- 1 / response[, 2]
        ll <- dnorm(y, mean=mu, sd=sigma, log=TRUE)
    } else if (family == 'shash') {
        mu <- response[, 1]
        sigma <- exp(response[, 2])  # mgcv represents on a log scale
        nu <- response[, 3]
        tau <- exp(response[, 4])  # mgcv represents on a log scale
        ll <- gamlss.dist::dSHASHo2(y, mean=mu, sd=sigma, nu=nu, tau=tau, log=TRUE)
    } else {
        stop(paste0('Unknown family: ', family))
    }

    return(list(
        logLik = ll
    ))
}
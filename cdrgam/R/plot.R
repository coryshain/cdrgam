#' Compute CDR-GAM plot data
#'
#' Compute data (numeric arrays representing lines or surfaces) from a CDR-GAM
#' model by predicting different terms under systematic manipulations of the
#' input. This function is very flexible at the expense of user-friendliness.
#' It can be called directly to create arbitrary plots, but it is typically
#' called by higher-level functions for generating commonly-needed plot types
#' (e.g., `plot_irfs()`).
#'
#' @param model A CDR-GAM object (i.e., an `mgcv` model fitted with a CDR-GAM
#'   model structure on CDR-GAM-structured data).
#' @param smooth_name The name of the smooth term to plot.
#' @param xvar The name of the variable to manipulate along the x-axis.
#' @param yvar The name of the variable to manipulate along the y-axis. If NULL,
#'   the function will return a line plot (1D manipulation). If provided,
#'   the function will return a surface plot (2D manipulation).
#' @param X_ref A list of reference values for some subset of variables in the
#'   model. If NULL, or if a variable is missing, the function will default to
#'   a value of 0 for numeric variables, the first level for factors that define
#'   fixed effects, and the value of `cdr_population_code()` for factors that
#'   define random grouping factors.
#' @param xaxis A numeric vector of values to use for the x-axis. If NULL, the
#'   function will construct the axis from the value of `xlim`.
#' @param xlim A numeric vector of length 2 specifying the minimum and maximum
#'   values for the x-axis. If NULL, the function will default to 0 and 1.
#' @param xres The number of points to use for the x-axis. If NULL, the function
#'   will default to 1024 for line plots and 32 for surface plots.
#' @param yaxis A numeric vector of values to use for the y-axis. If NULL, the
#'   function will construct the axis from the value of `ylim`.
#' @param ylim A numeric vector of length 2 specifying the minimum and maximum
#'   values for the y-axis. If NULL, the function will default to 0 and 1.
#' @param yres The number of points to use for the y-axis. If NULL, the function
#'   will default to 32.
#' @param reshape_3d A logical value indicating whether to reshape the output
#'   data into the matrix format expected by `plotly`'s 3D plotting
#'   functions. If FALSE, all values will be returned as flattened vectors,
#'   which can be used to e.g., construct heatmaps using `ggplot2`.
#' @return A list containing the following elements:
#'   - smooth_name: The name of the smooth term.
#'   - xvar: The name of the variable manipulated along the x-axis.
#'   - yvar: The name of the variable manipulated along the y-axis if appicable
#'       (otherwise NULL)
#'   - x: A numeric vector of x-axis values.
#'   - y: A numeric vector of y-axis values, if appicable (otherwise NULL)
#'   - resp: A numeric vector of response values.
#'   - resp.se: A numeric vector of standard errors for the response values.
#' @export
get_plot_data <- function(
        model,
        smooth_name,
        xvar,
        yvar=NULL,
        X_ref=NULL,
        xaxis=NULL,
        xlim=NULL,
        xres=NULL,
        yaxis=NULL,
        ylim=NULL,
        yres=NULL,
        reshape_3d=FALSE
) {
    if (is.null(model)) {stop('Value must be provided for model')}
    if (is.null(smooth_name)) {stop('Value must be provided for smooth_name')}
    if (is.null(xvar)) {stop('Value must be provided for xvar')}
    if ((!is.null(yvar)) && (xvar == yvar)) {stop('Cannot vary two axes along the same variable')}

    is_3d <- !is.null(yvar)

    # Initialize plotting resolutions
    if (is.null(xaxis)) {
        if (is_3d) {
            if (is.null(xres)) {
                xres <- 32
            }
        } else {
            if (is.null(xres)) {
                xres <- 1024
            }
        }
        xvar_base <- seq(0., 1., length.out=xres)
    } else {
        xres <- length(xaxis)
        xvar_base <- NULL
    }

    if (is_3d) {
        if (is.null(yaxis)) {
            if (is.null(yres)) {
                yres <- 32
            }
            yvar_base <- seq(0., 1., length.out=yres)
        } else {
            yres <- length(yaxis)
            yvar_base <- NULL
        }
        T <- xres * yres
    } else {
        T <- xres
    }

    # Initialize predictor reference
    term_names <- names(model$var.summary)
    if (is.null(X_ref)) {
        X_ref <- list()
    }
    for (term_name in term_names) {
        var.summary <- model$var.summary[[term_name]]
        if (!(term_name %in% names(X_ref))) {
            if (is.factor(var.summary)) {
                if (get_cdr_population_code() %in% levels(var.summary)) {  # Random grouping factor
                    X_ref[[term_name]] <- get_cdr_population_code()
                } else {  # Other factor
                    X_ref[[term_name]] <- levels(var.summary)[1]
                }
            } else if (term_name == 'mask') {  # Mask, unused for plotting
                X_ref[[term_name]] <- 1
            } else {  # Numeric
                X_ref[[term_name]] <- 0
            }
        }
    }

    if (is.null(xlim)) {
        xmin <- NULL
        xmax <- NULL
    } else {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }

    if (is.null(ylim)) {
        ymin <- NULL
        ymax <- NULL
    } else {
        ymin <- ylim[1]
        ymax <- ylim[2]
    }

    # Construct x-axis manipulation
    xdict <- list(
        axis_var=xvar,
        axis=xaxis,
        axis_ix=1,
        ax_min=xmin,
        ax_max=xmax,
        base=xvar_base,
        tile_3d=c(TRUE, yres)
    )
    params <- list(xdict)

    if (is_3d) {
        ydict <- list(
            axis_var=yvar,
            axis=yaxis,
            axis_ix=2,
            ax_min=ymin,
            ax_max=ymax,
            base=yvar_base,
            tile_3d=c(FALSE, xres)
        )
        params[[2]] <- ydict
    }

    plot_axes <- list(NULL, NULL)

    X <- NULL

    for (par in params) {
        axis_var <- par[['axis_var']]
        axis <- par[['axis']]
        axis_ix <- par[['axis_ix']]
        ax_min <- par[['ax_min']]
        ax_max <- par[['ax_max']]
        base <- par[['base']]
        tile_3d <- par[['tile_3d']]

        if (is.null(X)) {
            X <- list()
            for (term_name in term_names) {
                X[[term_name]] <- rep(X_ref[[term_name]], T)
            }
        }

        if (axis_var %in% term_names) {
            if (is.null(axis)) {
                if (is.null(ax_min)) {
                    ax_min <- 0
                }
                if (is.null(ax_max)) {
                    ax_max <- 1
                }
                axis <- (base * (ax_max - ax_min) + ax_min)
            }
            if (is_3d) {
                each <- tile_3d[1]
                n_tile <- tile_3d[2]
                if (each) {
                    axis <- rep(axis, each=n_tile)
                } else {
                    axis <- rep(axis, times=n_tile)
                }
            }
            X[[axis_var]] <- axis
            plot_axes[[axis_ix]] <- axis
        } else {
            stop(paste0('Unrecognized value for axis_var: ', axis_var))
        }
    }

    # Ensure correct shape and type for variables

    for (term_name in term_names) {
        dim(X[[term_name]]) <- c(length(X[[term_name]]), 1)
    }

    preds <- predict(model, newdata=X, type='terms', se.fit=TRUE, terms=smooth_name)

    resp <- preds$fit[,smooth_name]
    resp.se <- (preds$se.fit[,smooth_name]) * 2  # Approximate 95% CI

    dim(resp) <- T
    dim(resp.se) <- T

    x <- plot_axes[[1]]
    y <- plot_axes[[2]]

    if (reshape_3d) {
        df_ <- reshape_3d_plot(plot_axes[[1]], plot_axes[[2]], resp, resp.se)
        x <- df_$x
        y <- df_$y
        resp <- df_$z
        resp.se <- df_$z.se
    }

    return(list(
      smooth_name=smooth_name,
      xvar=xvar,
      yvar=yvar,
      x=x,
      y=y,
      resp=resp,
      resp.se=resp.se
    ))
}

#' Reshape 3D plot data
#'
#' Reshape flat vector-valued 3D plot data into a list of 2D matrices for
#' plotting with `plotly`. Vector lengths must either be perfect squares
#' (if `nx` is not provided) or must be divisible by `nx`.
#'
#' @param x A numeric vector of x-axis values
#' @param y A numeric vector of y-axis values
#' @param z A numeric vector of response values
#' @param z.se A numeric vector of standard errors for the response values
#' @param nx An integer specifying the number of x-axis values
#' @return A list containing reshaped 2D matrices for for all inputs
#' @export
reshape_3d_plot <- function(x, y, z, z.se=NULL, nx=NULL) {
    if (is.null(nx)) {
        nx <- as.integer(sqrt(length(x)))
    }
    ny <- length(y) / nx
    x_ <- x
    y_ <- y
    z_ <- z
    if (!is.null(z.se)) {
        z.se_ <- z.se
    }

    dim(x_) <- c(nx, ny)
    dim(y_) <- c(nx, ny)
    dim(z_) <- c(nx, ny)
    if (!is.null(z.se)) {
        dim(z.se_) <- c(nx, ny)
    }

    return(list(
        x=x_,
        y=y_,
        z=z_,
        z.se=z.se_
    ))
}

#' Plot lines
#'
#' Plot lines with confidence intervals. Each line is represented by a
#' numeric vector of x-axis values, a numeric vector of y-axis values
#' (called "resp"), and a numeric vector of standard errors for the
#' y-axis values (called "resp.se").
#'
#' @param data_2d A list of lists, each containing the following
#'   keys: "x", "resp", and (optionally) "resp.se". This is typically
#'   the output of `get_plot_data()`.
#' @param xlabel A string specifying the x-axis label
#' @param ylabel A string specifying the y-axis label
#' @param legend A logical value indicating whether to include a legend
#' @param add_origin A logical value indicating whether to add lines
#'   at x=0 and y=0
#' @return A ggplot2 object
#' @export
plot_lines <- function(
        data_2d,
        xlabel='x',
        ylabel='y',
        legend=TRUE,
        add_origin=TRUE
) {
    X <- data.frame(matrix(ncol=5, nrow=0))
    colnames(X) <- c('x', 'y', 'lq', 'uq', 'color')
    n_colors <- length(data_2d)
    palette <- rainbow(n_colors)

    include_ci <- FALSE
    for (line in data_2d) {
        X_ <- data.frame(
          x=line$x,
          y=line$resp
        )
        if (!(is.null(line$resp.se))) {
            include_ci <- TRUE
            X_$lq <- line$resp - line$resp.se
            X_$uq <- line$resp + line$resp.se
        }
        if ('name' %in% names(line)) {
            X_$color <- line$name
        } else {
            X_$color <- rep(line$smooth_name, length(line$x))
        }
        X <- rbind(X, X_)
    }

    if (length(X) > 0) {
        color_levels <- unique(X$color)
        X$color <- factor(X$color, levels=color_levels)
    }

    p <- ggplot2::ggplot(X, ggplot2::aes(x=x, y=y, lq=lq, uq=uq, color=color)) + ggplot2::theme_classic()
    if (add_origin) {
        p <- p + ggplot2::geom_vline(xintercept=0, color='black')
        p <- p + ggplot2::geom_hline(yintercept=0, color='black')
    }
    p <- p + ggplot2::geom_line()
    if (include_ci) {
        p <- p + ggplot2::geom_line(ggplot2::aes(x=x, y=lq), linetype='dashed')
        p <- p + ggplot2::geom_line(ggplot2::aes(x=x, y=uq), linetype='dashed')
    }

    p <- p + ggplot2::theme(legend.title=ggplot2::element_blank())
    p <- p + ggplot2::scale_colour_manual(values=palette)
    p <- p + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)

    if (!legend) {
        p <- p + ggplot2::theme(legend.position='none')
    }

    return(p)
}

#' Get impulse response function metadata for plots
#'
#' Compute metadata about each impulse response function (IRF) in a CDR-GAM
#' model. This function is useful because IRFs are implicit in the
#' specification of `mgcv` smooths and are not directly accessible from the
#' model itself, although they can be inferred from standard `mgcv` model
#' metadata. This function performs this inference and is used by
#' downstream functions (e.g., high-level plotting functions like
#' `plot_irfs()`).
#'
#' @param model A CDR-GAM object (i.e., an `mgcv` model fitted with a CDR-GAM
#'   model structure on CDR-GAM-structured data).
#' @param response_param An integer specifying the index of the response
#'   parameter to use for the IRF. E.g., if the response is Gaussian,
#'   `response_param=1` refers to the mean parameter, while `response_param=2`
#'   refers to the standard deviation parameter.
#' @param means A list of means for each variable in the model. This is
#'   typically the output of `get_cdr_means()`. If NULL, the function will
#'   default to 0 for all numeric variables and the first level for all
#'   factor variables.
#' @param sds A list of standard deviations for each variable in the model.
#'   This is typically the output of `get_cdr_sds()`. If NULL, the function
#'   will default to 1.
#' @param xlim A numeric vector of length 2 specifying the minimum and maximum
#'   values for the x-axis. If NULL, the function will default to 0 and 1.
#' @param irf_name_map A list of mappings from IRF names to human-readable
#'   names. This is useful for renaming IRFs that are not easily interpretable
#'   from the default names generated by `mgcv` (i.e. most of them).
#' @param add_rate A logical value indicating whether to add a rate term to the
#'   IRF list. This is useful for models that include a rate term.
#' @param exclude A character vector of term names to exclude from the IRF list.
#' @return A list of IRF metadata, each containing the following elements:
#'   - irf_name: The name of the IRF
#'   - term_name: The name of the term (predictor) in the model
#'   - smooth_name: The name of the smooth term in the model
#'   - ref: A list of reference values for the IRF
#'   - xlim: A numeric vector of length 2 specifying the minimum and maximum
#'     values for the x-axis
#' @export
get_irf_metadata <- function(
        model,
        response_param=1,
        means=NULL,
        sds=NULL,
        xlim=c(0, 1),
        irf_name_map=NULL,
        add_rate=TRUE,
        exclude=c('t_delta', 'mask')
) {
    if (is.null(sds)) {
        sds <- list()
    }
    if (is.null(irf_name_map)) {
        irf_name_map <- list()
    }
    irf_keys <- names(irf_name_map)
    irf_keys <- irf_keys[order(nchar(irf_keys), decreasing=TRUE)]

    term_dtypes <- attr(model$terms, 'dataClasses')
    term_labels <- names(term_dtypes)
    for (term_label in term_labels) {
        if (grepl('I\\([^)]+\\)', term_label)) {
            reduced_label <- gsub('I\\(([^)]+)\\)', '\\1', term_label)
            if (reduced_label %in% term_labels) {
                exclude <- c(exclude, term_label)
            }
        }
    }
    factors <- term_labels[term_dtypes == 'factor']
    gf <- character()
    for (factor in factors) {
        if (get_cdr_population_code() %in% levels(model$var.summary[[factor]])) {
            gf <- c(gf, factor)
        }
    }
    exclude <- c(exclude, factors)

    irfs <- list()
    for (smooth in model$smooth) {
        smooth_name <- smooth$label
        response_param_ix <- as.integer(gsub('[a-zA-Z]+(\\.(\\d*))?\\(.+', '\\2', smooth_name))
        if (is.na(response_param_ix)) {
            response_param_ix <- 0
        }
        response_param_ix <- response_param_ix + 1
        if (response_param_ix != response_param) {
            next
        }
        sel <- !(smooth$term %in% exclude)
        term_names <- smooth$term
        if (sum(term_names %in% gf) > 0) { # Skip random effects
            next
        }
        if (add_rate && length(term_names) == 1 && term_names == 't_delta') {  # Rate term
            irf <- list(
                irf_name='Rate',
                term_name='t_delta',
                smooth_name=smooth_name,
                ref=means,
                xlim=xlim
            )
            irfs <- c(irfs, list(irf))
        } else if ('t_delta' %in% term_names) {
            term_names <- term_names[sel]
            for (term_name in term_names) {
                irf_name <- paste0(term_name, ' | ', smooth_name)
                for (irf_key in irf_keys) {
                    if (grepl(irf_key, irf_name, fixed=TRUE)) {
                        irf_name <- irf_name_map[[irf_key]]
                        break
                    }
                }
                if (term_name %in% names(sds)) {
                    delta <- sds[[term_name]]
                } else {
                    delta <- 1
                }
                ref <- means
                ref[[term_name]] <- ref[[term_name]] + delta

                irf <- list(
                    irf_name=irf_name,
                    term_name=term_name,
                    smooth_name=smooth_name,
                    ref=ref,
                    xlim=xlim
                )
                irfs <- c(irfs, list(irf))
            }
        }
    }
    return(irfs)
}

#' Generate impulse response function line plots
#'
#' Generate line plots of impulse response functions (IRFs) from a CDR-GAM
#' model. This function is a high-level wrapper around `get_irf_metadata()`
#' and `get_plot_data()`.
#'
#' @param model A CDR-GAM object (i.e., an `mgcv` model fitted with a CDR-GAM
#'   model structure on CDR-GAM-structured data).
#' @param response_param An integer specifying the index of the response
#'   parameter to use for the IRF. E.g., if the response is Gaussian,
#'   `response_param=1` refers to the mean parameter, while `response_param=2`
#'   refers to the standard deviation parameter.
#' @param means A list of means for each variable in the model. This is
#'   typically the output of `get_cdr_means()`. If NULL, the function will
#'   default to 0 for all numeric variables and the first level for all
#'   factor variables.
#' @param sds A list of standard deviations for each variable in the model.
#'   This is typically the output of `get_cdr_sds()`. If NULL, the function
#'   will default to 1.
#' @param xlim A numeric vector of length 2 specifying the minimum and maximum
#'   values for the x-axis. If NULL, the function will default to 0 and 1.
#' @param irf_name_map A list of mappings from IRF names to human-readable
#'   names. This is useful for renaming IRFs that are not easily interpretable
#'   from the default names generated by `mgcv` (i.e. most of them).
#' @param xlabel A string specifying the x-axis label
#' @param ylabel A string specifying the y-axis label
#' @param legend A logical value indicating whether to include a legend
#' @param exclude A character vector of term names to exclude from the IRF list.
#' @return A ggplot2 object
#' @export
plot_irfs <- function(
        model,
        response_param=1,
        means=NULL,
        sds=NULL,
        xlim=c(0, 1),
        irf_name_map=NULL,
        xlabel='Delay (s)',
        ylabel='y',
        legend=TRUE,
        exclude=c('t_delta', 'mask')
) {
    irfs <- get_irf_metadata(
        model,
        response_param=response_param,
        means=means,
        sds=sds,
        xlim=xlim,
        irf_name_map=irf_name_map,
        exclude=exclude
    )
    data_2d <- list()
    for (irf in irfs) {
        plot_data <- get_plot_data(
            model,
            smooth_name=irf$smooth_name,
            xvar='t_delta',
            X_ref=irf$ref,
            xlim=irf$xlim
        )
        plot_data['name'] <- irf$irf_name
        data_2d <- c(data_2d, list(plot_data))
    }

    p <- plot_lines(
        data_2d,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=legend
    )
    return(p)
}

#' Generate curvature plots
#'
#' Generate plots of the curvature of the response surface of a CDR-GAM model
#' at a specific delay (t_delta) value. This function is a high-level wrapper
#' around `get_irf_metadata()` and `get_plot_data()`.
#'
#' @param model A CDR-GAM object (i.e., an `mgcv` model fitted with a CDR-GAM
#'   model structure on CDR-GAM-structured data).
#' @param response_param An integer specifying the index of the response
#'   parameter to use for the IRF. E.g., if the response is Gaussian,
#'   `response_param=1` refers to the mean parameter, while `response_param=2`
#'   refers to the standard deviation parameter.
#' @param means A list of means for each variable in the model. This is
#'   typically the output of `get_cdr_means()`. If NULL, the function will
#'   default to 0 for all numeric variables and the first level for all
#'   factor variables.
#' @param sds A list of standard deviations for each variable in the model.
#'   This is typically the output of `get_cdr_sds()`. If NULL, the function
#'   will default to 1.
#' @param quantiles A list of quantiles for each variable in the model. This is
#'   typically the output of `get_cdr_quantiles()`.
#' @param range A numeric value specifying the range of quantiles to use for
#'   the x-axis. For example, a value of 0.9 would use the 5th to 95th
#'   quantiles for the x-axis. Overridden by `xlim`.
#' @param xlim A list of numeric vectors of length 2 specifying the minimum and
#'   maximum values for the x-axis for each term in the model. If NULL, the
#'   function will default to 0 and 1.
#' @param t_delta_ref A numeric value specifying the delay (t_delta) value at
#'   which to evaluate the curvature.
#' @param irf_name_map A list of mappings from IRF names to human-readable
#'   names. This is useful for renaming IRFs that are not easily interpretable
#'   from the default names generated by `mgcv` (i.e. most of them).
#' @param xlabel A string specifying the x-axis label
#' @param ylabel A string specifying the y-axis label
#' @param legend A logical value indicating whether to include a legend
#' @param exclude A character vector of term names to exclude from the IRF list.
#' @return A ggplot2 object
#' @export
plot_curvature <- function(
        model,
        response_param=1,
        means=NULL,
        sds=NULL,
        quantiles=NULL,
        range=1,
        xlim=NULL,
        t_delta_ref=0,
        irf_name_map=NULL,
        xlabel='Predictor',
        ylabel='y',
        legend=TRUE,
        exclude=c('t_delta', 'mask')
) {

    irfs <- get_irf_metadata(
        model,
        response_param=response_param,
        means=means,
        sds=sds,
        irf_name_map=irf_name_map,
        exclude=exclude,
        add_rate=FALSE
    )
    data_2d <- list()
    for (irf in irfs) {
        term_name <- irf$term_name
        if (!is.null(xlim) && term_name %in% names(xlim)) {
            irf$xlim <- xlim[[term_name]]
        } else if (!is.null(quantiles) && term_name %in% names(quantiles)) {
            ix_min <- round((1 - range) / 2 * length(quantiles[[term_name]])) + 1
            ix_max <- round((1 - ((1 - range) / 2)) * length(quantiles[[term_name]]))
            irf$xlim <- c(quantiles[[term_name]][ix_min], quantiles[[term_name]][ix_max])
        } else {
            irf$xlim <- c(0, 1)
        }
        irf$X_ref$t_delta <- t_delta_ref
        plot_data <- get_plot_data(
            model,
            smooth_name=irf$smooth_name,
            xvar=term_name,
            X_ref=irf$ref,
            xlim=irf$xlim
        )
        plot_data['name'] <- irf$irf_name
        data_2d <- c(data_2d, list(plot_data))
    }

    p <- plot_lines(
        data_2d,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=legend
    )
    return(p)
}
#' Get CDR-GAM model configuration object from YAML file
#'
#' Get configuration object from YAML file. This function reads a YAML file
#' and returns a list with the configuration. The configuration is populated
#' as needed against a set of defaults.
#'
#' @param path A string with the path to the YAML file
#' @return A list with the configuration
#' @export
get_cfg <- function(path) {
    defaults <- GLOBAL.CDRGAM
    cfg <- yaml::read_yaml(path)
    if (!('data' %in% names(cfg))) {
        stop('Required field "data" missing from config')
    }
    if (!('models' %in% names(cfg))) {
        stop(paste0('Required field "models" missing from config'))
    }

    if (!('output_dir' %in% names(cfg))) {
        output_dir <- file.path('cdrgam_results', gsub('.yml', '', basename(path), fixed=TRUE))
        cfg$output_dir <- output_dir
    }

    keys <- names(defaults$data)
    keys <- keys[!(keys %in% names(cfg$data))]
    cfg$data[keys] <- defaults$data[keys]
    if (!('plot' %in% names(cfg))) {
        cfg$plot <- list()
    }
    keys <- names(defaults$plot)
    keys <- keys[!(keys %in% names(cfg$plot))]
    cfg$plot[keys] <- defaults$plot[keys]
    if (is.null(cfg$plot$t_delta_xlim) && !is.null(cfg$data$t_delta_cutoff)) {
        cfg$plot$t_delta_xlim <- c(0, cfg$data$t_delta_cutoff)
    }
    if (is.null(cfg$plot$output_dir)) {
        cfg$plot$output_dir <- cfg$output_dir
    }
    if ('formula' %in% names(cfg)) {
        keys <- names(cfg$formula)
        defaults$formula[keys] <- cfg$formula[keys]
    }
    if ('model' %in% names(cfg)) {
        keys <- names(cfg$model)
        defaults$model[keys] <- cfg$model[keys]
    }

    for (model in names(cfg$models)) {
        keys <- names(defaults$model)
        keys <- keys[!(keys %in% names(cfg$models[[model]]))]
        cfg$models[[model]][keys] <- defaults$model[keys]
        formula_defaults <- defaults$formula
        keys <- names(cfg$models[[model]])
        keys <- keys[keys %in% names(formula_defaults)]
        formula_defaults[keys] <- cfg$models[[model]][keys]
        for (formula_key in names(cfg$models[[model]]$formula)) {
            formula <- cfg$models[[model]]$formula[[formula_key]]
            for (ix in seq_along(formula)) {
                keys <- names(formula_defaults)
                keys <- keys[!(keys %in% names(cfg$models[[model]]$formula[[formula_key]][[ix]]))]
                cfg$models[[model]]$formula[[formula_key]][[ix]][keys] <- formula_defaults[keys]
            }
        }
    }
    return(cfg)
}

#' Get CDR-GAM plot configuration object from YAML file
#'
#' Get configuration object from YAML file. This function reads a YAML file
#' and returns a list with the configuration. The configuration is populated
#' as needed against a set of defaults.
#'
#' @param path A string with the path to the YAML file
#' @return A list with the configuration
#' @export
get_plot_cfg <- function(path) {
    defaults <- GLOBAL.CDRGAM
    cfg <- yaml::read_yaml(path)
    keys <- names(defaults$plot)
    keys <- keys[!(keys %in% names(cfg))]
    cfg[keys] <- defaults$plot[keys]
    return(cfg)
}
#' Command-line utility for running CDR-GAM models
#'
#' A command-line utility for fitting, plotting, and evaluating CDR-GAM models.
#' Wraps the `main()` function for use in the command line as follows:
#'     `Rscript -e "cdrgam::cli()" <cfg> <model_name> --<option1> <value1> --<option2> <value2> ...`
#' where `<cfg>` is the path to a YAML configuration file, `<model_name>` is the
#' name of the model to fit, and `<option1>`, `<option2>`, etc., are any optional
#' parameters to `main()`. Note that, due to limitations of command line parsing in R,
#' options that accept lists (e.g., `eval_partition`) must be provided as comma-delimited strings.
#' @export
cli <- function() {
    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for running CDR-GAM models.",
            "Two positional arguments are required: <cfg> and <model_name>."
        )
    )
    # Optional arguments to `main()`
    parser <- optparse::add_option(parser, '--fit', type='logical', default=TRUE, help='Whether to fit the model')
    parser <- optparse::add_option(parser, '--plot', type='logical', default=TRUE, help='Whether to plot the model')
    parser <- optparse::add_option(parser, '--evaluate', type='logical', default=TRUE,
                                   help='Whether to evaluate the model')
    parser <- optparse::add_option(parser, '--overwrite', type='logical',
                                   help='Whether to overwrite an existing model', action='store_true')
    parser <- optparse::add_option(parser, '--clean', type='logical', help='Whether to clean data from fitted model')
    parser <- optparse::add_option(parser, '--keep_model',  type='logical',
                                   help='Whether to save model matrix in fitted model')
    parser <- optparse::add_option(parser, '--eval_partition', default='train,val',
                                   help='Data partition(s) on which to evaluate (comma-delimited)')
    parser <- optparse::add_option(parser, '--extra_cols', type='logical',
                                   help='Whether to add extra columns from Y to the prediction output')
    parser <- optparse::add_option(parser, '--plot_cfg', help='Path to additional plot config file')
    parser <- optparse::add_option(parser, '--dump_plot_data', type='logical',
                                   help='Whether to save plot data to CSV table')
    cliargs <- optparse::parse_args(parser, positional_arguments=2)
    args <- list(cfg=cliargs$args[[1]], model_name=cliargs$args[[2]])
    for (x in names(cliargs$options)) {
        if (x != 'help' && !is.null(cliargs$options[[x]])) {
            args[[x]] <- cliargs$options[[x]]
            if (x == 'eval_partition') {
                args[[x]] <- strsplit(cliargs$options[[x]], ',')[[1]]
            }
        }
    }
    do.call(main, args)
}

#' Command-line utility for making slurm jobs from CDR-GAM config files
#'
#' A command-line utility for generating slurm job scripts from CDR-GAM config
#' files. The utility generates one job script per model in each config file.
#' The job scripts are written to the current working directory.  Note that,
#' due to limitations of command line parsing in R, options that accept lists
#' (e.g., `eval_partition`) must be provided as comma-delimited strings.
#'
#' @export
make_jobs <- function() {
    base <- paste(
        '#!/bin/bash',
        '#',
        '#SBATCH --job-name=%s',
        '#SBATCH --output="%s-%%N-%%j.out"',
        '#SBATCH --time=%d:00:00',
        '#SBATCH --mem=%dgb',
        '#SBATCH --ntasks=%d',
        sep='\n'
    )
    slurm_names <- c('time', 'mem', 'ntasks', 'exclude', 'partition')

    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for making slurm jobs from CDR-GAM config files.",
            "One or more positional arguments (each corresponding to a config file)",
            "are required."
        )
    )
    # SLURM options
    parser <- optparse::add_option(parser, c('-t', '--time'), default=24, help='Max runtime (in hours)')
    parser <- optparse::add_option(parser, c('-m', '--mem'), default=8, help='Memory allocation (in GB)')
    parser <- optparse::add_option(parser, c('-n', '--ntasks'), default=2, help='Number of cores ("tasks")')
    parser <- optparse::add_option(parser, c('-e', '--exclude'), help='Comma-delimited list of nodes to exclude')
    parser <- optparse::add_option(parser, c('-P', '--partition'), help='SLURM partition to use')
    # Optional arguments to `main()`
    parser <- optparse::add_option(parser, '--fit', type='logical', default=TRUE, help='Whether to fit the model')
    parser <- optparse::add_option(parser, '--plot', type='logical', default=TRUE, help='Whether to plot the model')
    parser <- optparse::add_option(parser, '--evaluate', type='logical', default=TRUE,
                                   help='Whether to evaluate the model')
    parser <- optparse::add_option(parser, '--overwrite', type='logical',
                                   help='Whether to overwrite an existing model', action='store_true')
    parser <- optparse::add_option(parser, '--clean', type='logical', help='Whether to clean data from fitted model')
    parser <- optparse::add_option(parser, '--keep_model',  type='logical',
                                   help='Whether to save model matrix in fitted model')
    parser <- optparse::add_option(parser, '--eval_partition', default='train,val',
                                   help='Data partition(s) on which to evaluate (comma-delimited)')
    parser <- optparse::add_option(parser, '--extra_cols', type='logical',
                                   help='Whether to add extra columns from Y to the prediction output')
    parser <- optparse::add_option(parser, '--plot_cfg', help='Path to additional plot config file')
    parser <- optparse::add_option(parser, '--dump_plot_data', type='logical',
                                   help='Whether to save plot data to CSV table')
    cliargs <- optparse::parse_args(parser, positional_arguments=c(1, Inf))
    cfg_paths <- cliargs$args
    options <- cliargs$options

    kwargs <- options[!(names(options) %in% c(slurm_names, 'help'))]
    job_name <- character()
    if (options$fit) {
        job_name <- c(job_name, 'fit')
    }
    if (options$plot) {
        job_name <- c(job_name, 'plot')
    }
    if (options$evaluate) {
        job_name <- c(job_name, 'evaluate')
    }
    job_name <- paste(job_name, collapse='-')

    for (cfg_path in cfg_paths) {
        cfg <- get_cfg(cfg_path)
        models <- names(cfg$models)
        for (model in models) {
            basename_ <- gsub('\\.yml$', '', basename(cfg_path))
            job_name_ <- paste0(basename_, '-', model, '_', job_name)
            job_path <- paste0(job_name_, '.pbs', sep='')
            script <- sprintf(base, job_name_, job_name_, options$time, options$mem, options$ntasks)
            if (!is.null(options$exclude)) {
                script <- paste(script, sprintf('#SBATCH --exclude=%s', options$exclude), sep='\n')
            }
            if (!is.null(options$partition)) {
                script <- paste(script, sprintf('#SBATCH --partition=%s', options$partition), sep='\n')
            }
            script <- paste0(script, '\n\n')
            script <- paste0(
                script,
                sprintf('Rscript -e "cdrgam::cli()" %s %s %s', cfg_path, model,
                        paste0('--', paste(names(kwargs), kwargs, sep='='), collapse=' '))
            )
            script <- paste0(script, '\n')
            cat(script, file=job_path)
        }
    }
}

#' Command-line utility for statistically comparing two CDR-GAM models
#'
#' A command-line utility for statistically comparing two CDR-GAM models. The
#' utility wraps the `test_cdrgam()` function for use in the command line as
#' follows:
#'    `Rscript -e "cdrgam::test()" <cfg0> <model0> <model1> --<option1> <value1> --<option2> <value2> ...`
#' where `<cfg0>` is the path to the first model configuration file, `<model0>`
#' is the name of the first model, `<model1>` is the name of the second model, and
#' `<option1>`, `<option2>`, etc., are any optional parameters to `test_cdrgam()`.
#' Note that, due to limitations of command line parsing in R, the `--eval_partition`
#' option can accept lists, but only as comma-separated (rather than space-separated) strings.
#' @export
test <- function() {
    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for statistically comparing two CDR-GAM models.",
            "Three positional arguments are required: <cfg0>, <model0>, and <model1>."
        )
    )
    parser <- optparse::add_option(parser, '--cfg1',
                                   help='Path to the second model configuration file (if distinct from cfg1)')
    parser <- optparse::add_option(parser, '--eval_partition', default='val',
                                   help='Data partition(s) on which to evaluate (comma-delimited)')
    parser <- optparse::add_option(parser, '--statistic', default='logLik', help='Statistic to test')
    parser <- optparse::add_option(parser, '--output_dir', help='Path to output directory')
    parser <- optparse::add_option(parser, '--n_iter', help='Number of resampling iterations to run')
    parser <- optparse::add_option(parser, '--n_tails', help='Number of tails in the test')
    parser <- optparse::add_option(parser, '--nested', help='Whether model1 is to be treated as nested inside model0',
                                   action='store_true')
    cliargs <- optparse::parse_args(parser, positional_arguments=3)
    args <- list(cfg0=cliargs$args[[1]], model_name0=cliargs$args[[2]], model_name1=cliargs$args[[3]])
    for (x in names(cliargs$options)) {
        if (x != 'help' && !is.null(cliargs$options[[x]])) {
            args[[x]] <- cliargs$options[[x]]
            if (x == 'eval_partition') {
                args[[x]] <- strsplit(cliargs$options[[x]], ',')[[1]]
            }
        }
    }
    do.call(test_cdrgam, args)
}

#' Command-line utility for comparing configs of saved CDR-GAM models
#'
#' This utility compares the configuration files for two
#' saved CDR-GAM models (fitted via `fit_cdrgam()` and saved via
#' `save_cdrgam()`).
#' @export
diff_cfgs <- function() {
    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for comparing configs of saved CDR-GAM models.",
            "Two positional arguments are required: <model_path1> and <model_path2>."
        )
    )
    cliargs <- optparse::parse_args(parser, positional_arguments=2)
    model_path1 <- cliargs$args[[1]]
    model_path2 <- cliargs$args[[2 ]]
    cfg1 <- readRDS(model_path1)$cfg
    cfg2 <- readRDS(model_path2)$cfg

    cfg1_path <- tempfile()
    cfg2_path <- tempfile()

    yaml::write_yaml(cfg1, cfg1_path)
    yaml::write_yaml(cfg2, cfg2_path)
    out <- capture.output(system(sprintf('diff %s %s', cfg1_path, cfg2_path)))
    cat(paste(out, '\n'), stdout())
}
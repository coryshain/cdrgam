#' Command-line utility for running CDR-GAM models
#'
#' A command-line utility for fitting, evaluating, and plotting CDR-GAM models.
#' Wraps the `main()` function for use in the command line as follows:
#'     `Rscript -e "cdrgam::cli()" <cfg> <model_name> --<option1> <value1> --<option2> <value2> ...`
#' where `<cfg>` is the path to a YAML configuration file, `<model_name>` is the
#' name of the model to fit, and `<option1>`, `<option2>`, etc., are any optional
#' parameters to `main()`. Note that, due to limitations of command line parsing in R,
#' the `--eval_partition` option can accept lists, but only as comma-separated (rather
#' than space-separated) strings.
#' @export
cli <- function() {
    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for running CDR-GAM models.",
            "Two positional arguments are required: <cfg> and <model_name>."
        )
    )
    parser <- optparse::add_option(parser, '--fit', default=TRUE, help='Whether to fit the model')
    parser <- optparse::add_option(parser, '--eval_partition', default='train,val',
                                   help='Partition(s) on which to evaluate')
    parser <- optparse::add_option(parser, '--plot', default=TRUE, help='Whether to plot the model')
    parser <- optparse::add_option(parser, '--plot_cfg', help='Path to the plot configuration file')
    parser <- optparse::add_option(parser, '--overwrite', help='Whether to overwrite an existing model',
                                   action='store_true')
    cliargs <- optparse::parse_args(parser, positional_arguments=2)
    args <- list(cfg=cliargs$args[[1]], model_name=cliargs$args[[2]])
    for (x in names(cliargs$options)) {
        if (x != 'help') {
            args[[x]] <- cliargs$options[[x]]
        } else if (x == 'eval_partition') {
            args[[x]] <- strsplit(cliargs$options[[x]], ',')
        }
    }
    do.call(main, args)
}

#' Command-line utility for making slurm jobs from CDR-GAM config files
#'
#' A command-line utility for generating slurm job scripts from CDR-GAM config
#' files. The utility generates one job script per model in each config file.
#' The job scripts are written to the current working directory.
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

    parser <- optparse::OptionParser(
        description=paste(
            "Command-line utility for making slurm jobs from CDR-GAM config files.",
            "One or more positional arguments (each corresponding to a config file)",
            "are required."
        )
    )
    parser <- optparse::add_option(parser, c('-j', '--jobtype'), default='all',
                                   help='Type of job to run (all, fit, eval, plot), or comma-delimite list of these')
    parser <- optparse::add_option(parser, c('-t', '--time'), default=24, help='Max runtime (in hours)')
    parser <- optparse::add_option(parser, c('-m', '--mem'), default=8, help='Memory allocation (in GB)')
    parser <- optparse::add_option(parser, c('-n', '--ntasks'), default=2, help='Number of cores ("tasks")')
    parser <- optparse::add_option(parser, c('-e', '--exclude'), help='Comma-delimited list of nodes to exclude')
    parser <- optparse::add_option(parser, c('-p', '--eval_partition'), default='val',
                                   help='Data partition to use for evaluation')
    parser <- optparse::add_option(parser, c('-P', '--partition'), help='SLURM partition to use')
    parser <- optparse::add_option(parser, c('-c', '--plot_cfg'), help='Path to additional plot config file')
    cliargs <- optparse::parse_args(parser, positional_arguments=c(1, Inf))
    cfg_paths <- cliargs$args
    options <- cliargs$options

    kwargs <- list()
    if (options$jobtype == 'all' || options$jobtype == 'fit') {
        kwargs$fit <- TRUE
    } else {
        kwargs$fit <- FALSE
    }
    if (options$jobtype == 'all' || options$jobtype == 'eval') {
        kwargs$eval_partition <- options$eval_partition
    } else {
        kwargs$eval_partition <- ''
    }
    if (options$jobtype == 'all' || options$jobtype == 'plot') {
        kwargs$plot <- TRUE
    } else {
        kwargs$plot <- FALSE
    }

    for (cfg_path in cfg_paths) {
        cfg <- get_cfg(cfg_path)
        models <- names(cfg$models)
        for (model in models) {
            basename_ <- gsub('\\.yml$', '', basename(cfg_path))
            job_name <- paste0(basename_, '-', model)
            job_path <- paste0(job_name, '_slurm.sh', sep='')
            script <- sprintf(base, job_name, job_name, options$time, options$mem, options$ntasks)
            if (!is.null(options$exclude)) {
                script <- paste(script, sprintf('#SBATCH --exclude=%s', options$exclude), sep='\n')
            }
            if (!is.null(options$partition)) {
                script <- paste(script, sprintf('#SBATCH --partition=%s', options$partition), sep='\n')
            }
            script <- paste0(script, '\n\n')
            script <- paste0(
                script,
                sprintf('Rscript -e "cdrgam::cli()" %s %s %s\n', cfg_path, model,
                        paste0('--', paste(names(kwargs), kwargs, sep='='), collapse=' '))
            )
            cat(script, file=job_path)
        }
    }
}
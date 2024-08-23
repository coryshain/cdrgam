GLOBAL.CDRGAM <- list(
    population_code='!!!POPULATION!!!',
    data=list(
        sep=',',
        history_length=16,
        future_length=0,
        t_delta_cutoff=NULL
    ),
    formula=list(
        k_t=20,
        k=5,
        bs_t='cr',
        bs=NULL
    ),
    model=list(
        gam=list(
            gamma=1
        ),
        family='gaussian'
    ),
    plot=list(
        t_delta_xlim=NULL,
        width=2,
        height=1.6,
        scale=4,
        legend=TRUE
    )
)

#' Get the CDR-GAM population code
#'
#' Get the string code for the reserved value (representing the population
#' level)  of CDR-GAM random effects ("!!!POPULATION!!!"). Use for
#' prediction and plotting of fixed effects.
#'
#' @returns A string
#' @export
get_cdr_population_code <- function() {
    return(GLOBAL.CDRGAM[['population_code']])
}

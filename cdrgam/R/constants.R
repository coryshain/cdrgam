#' Get the CDR-GAM population code
#'
#' Get the string code for the reserved value (representing the population
#' level)  of CDR-GAM random effects ("!!!POPULATION!!!"). Use for
#' prediction and plotting of fixed effects.
#'
#' @returns A string
#' @export
cdr_population_code <- function() {
    return('!!!POPULATION!!!')
}
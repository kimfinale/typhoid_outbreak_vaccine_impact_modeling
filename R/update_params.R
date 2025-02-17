#' Update params
#' Updata parameters such that they balance. Changing popoulation size will
#' lead to the chanages in the initials values of susceptibles, infecteds, etc.
#' Similarly, changing the proportion of immune in the beginning of the population
#' will lead to chantges in initial values of the state variables
#' Change in the date of introduction change the simulation period
#' Updating occurs by way of the name of the parameter
#'
#' @param pars # parameter values to be varied
#' @param pars_baseline # This represents a complete set of parameters needed to run the model including initial values for the state variables. Some of the parameter values may be updated based on the pars
#' @param model # model (e.g., seirw)
#' @param output_time #
#' @param output_state
#'
#' @return
#' @export
#'
#' @examples
#'
update_params <- function(pars=NULL,
                          pars_baseline=NULL) {

  if (is.null(pars) | is.null(pars_baseline)) {
    stop("Both pars and pars_baseline must not be NULL")
  }

  params = pars_baseline

  nms = names(pars)
  for (nm in nms) {
    params[[nm]] = pars[[nm]]
  }

  # simulation times changes by Day 1 (introduction of the virus)
  # index starts from 1 and therefore one must be added to get the correct
  params[["ndays"]] <- round(params[["day1"]]) +
    (params[["obs_length"]] * params[["output_days"]]) + 1

  pop = params[["population"]]
  pop = pop * params[["n0"]] # adjusted initial population size
  s0 = params[["s0"]] # proportion susceptible
  # initial values for the state variables
  params[["susceptible"]] <- pop*s0
  i0 = params[["i0"]] # proportion infectious

  if (i0 > 1e-12) { # set to i0=0 if you want to use initial number of I
    fa = params[["fA"]] # proportion asymptomatic
    latent_pd = 1 / params[["epsilon"]] # mean latent period = 1/epsilon
    infect_pd = 1 / params[["gamma"]] # mean infectious period = 1/gamma
    fe = latent_pd/(latent_pd + infect_pd)
    fi = infect_pd / (latent_pd + infect_pd) # I / (E + I)

    params[["exposed"]] <- pop*i0*(1-fi)
    params[["infectious"]] <- pop*i0*fi*(1-fa)
    params[["asymptomatic"]] <- pop*i0*fi*fa
    params[["recovered"]] <- pop*(1-s0-i0)
  }
  # params[["water"]] <- 0
  # params[["cumul_exposed"]] <- 0
  # params[["cumul_infectious"]] <- 0

  return (params)
}

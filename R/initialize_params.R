initialize_params <- function(...){

  params <- list() # input parameters for the SEIR model

  # parameters to be estimates based on the outbreak data
  params$day1 <- 0 #
  params$s0 <- 1.0 - 1e-6 #
  params$i0 <- 1e-6 # prop initially infected (E + I + A)
  params$n0 <- 1 # ratio of the underlying population compared to the population size of the admin in which the outbreak reported
  # params$r0 <- 0.0 #
  params$R0 <- 3.0 #
  params$prop_report <- 1 #
  # baseline population size is critical for the outbreak size
  # and the currently is the the population size of the admin in which the
  # outbreak was reported.  Actual population is likely to be a fraction
  # params$prop_eff_pop <- 1 #

  # time window over which the number of cases is tracked (i.e., weekly reported)
  params$output_days <- 7
  params$output_state <- "CI"
  # obs_length refers the duration during which data are available in days
  # this is updated for each data set
  params$obs_length <- 365 # 20 weeks
  # this refers to the total simulation days and has to be larger than or equal
  # to the obs_length
  params$ndays <- 365 # total number of days for output

  params$model <- seiarw_2ag_erlang_vacc #seiarw #

  params$tau <- 0.01 # time step size for numerical integration

  params$epsilon <- 1/1.4 # mean latent period = 1/epsilon
  params$gamma <- 1/2 # mean infectious period = 1/gamma
  params$fA <- 0.5 # fraction of asymptomatic state
  params$bA <- 0.05 # relative infectiousness of asymptomatic state
  params$kappa <- 575 # excretion rate cells per person per day
  params$prop_children <- 0.193 #
# exponent used to model sub-exponential growth (0 < expon < 1) (foi=I^{expon}*S/N)
  params$expon <- 1 # exponential growth when expo = 1

  params$xi <- 1/21 # mean decay rate of Vibrio cholerae
  params$K <- 10000 # half-infective bacteria dose (10,000 cells/ml)

  params$R0W <- 0.0 #
  params$sigma <- 1/(4*365) # 1/sigma = mean duration of natural immunity
  # 1/sigma_v1 = mean duration of OCV-induced immunity (1st dose)
  params$sigma_v1 <- 1/(2*365)
  # 1/sigma_v2 = mean duration of OCV-induced immunity (2nd dose)
  params$sigma_v2 <- 1/(4*365)

  # vaccination-campaign associated
  params$vacc_cov_v1 <- c(0.8, 0.8) # two age groups
  params$vacc_cov_v2 <- c(0.8, 0.8)
  params$campaign_1_start <- 1e6
  params$campaign_2_start <- 1e6
  params$campaign_dur <- 14 # 14 days of vaccination campaign
  # params$delay_until_2nd_campaign <- 14 # 30 days of delay between the 1st and the 2nd campaign
  params$vacc_eff_v1 <- c(0.3, 0.5) # two age groups
  params$vacc_eff_v2 <- c(0.5, 0.7)
  params$immunity_dev_rate <- 1/7 # vaccine-induced immunity develops after 7 d on average
  # case threshold over which intervention will be implemented
  params$case_threshold <- 10

  params$alpha <- 0.0; # proportional reduction in R0
  params$case_track_window <- 7;

  # initial values for the state variables
  params$population <- 1
  params$susceptible <- params$population * params$s0
  params$exposed <- 0
  params$infectious <- params$population * params$i0
  params$asymptomatic <- 0
  params$recovered <- params$population * (1 - params$s0 - params$i0)
  params$water <- 0
  params$cumul_exposed <- 0
  params$cumul_infectious <- 0

  # update params
  params = update_params(pars=list(...), pars_baseline=params)

  return(params)
}

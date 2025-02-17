#' Sets the parameter valurotective factor (e.g.,
#' at least basic sanitation) and those without for a given overall incidence
#' rate and the proportion of population for each category

draw_parameter_samples <- function(nruns=200,
                       parameter_data = NULL){

  # adopted from
  # EDIT 18 Jan 2024
  # vaccination_overlap
  # this parameter represents the proportion of the vaccine recipients of the
  # first round that received the second dose. This proportion was extracted
  # by the fact that the mean coverage during the first round was 92.1% and the
  # mean coverage who received the two doses was 69.9% [range: 27.5% to 95.3%]
  # Therefore, the proportion of the first-dose recipients who received the
  # dose as well = 69.9/92.1 = 0.7589577
  # lower and upper bound were derived by the fact
  # lower bound = 27.5/92.1 = 0.2985885
  # upper bound = 95.3/92.1 = 1.034745
  # since the upper bound could not exceed 1,
  # we arbitrarily set as the 95.3%, which is somewhat high and appears to be
  # consistent with coverages at the high end


  # param_names <- c("run_id",
  #                  "vacc_cov",
  #                  "vaccination_overlap",
  #                  "vacc_eff_direct_U5",
  #                  "vacc_eff_direct_5plus",
  #                  "vacc_eff_id_indirect", # select the curve fitting result
  #                  "dur_illness",
  #                  "prop_mild",
  #                  "prop_moderate",
  #                  "prop_severe",
  #                  "dis_weight_mild",
  #                  "dis_weight_modeerate",
  #                  "dis_weight_severe",
  #                  "discount_rate",
  #                  "reference_year")
  #
  # params <- data.frame(matrix(NA, nrow = nruns, ncol = length(param_names)))
  #
  # names(params) <- param_names
  # params$run_id <- 1:nruns
  #
  # # create local variables for convenience
  # ir <- incidence_rate_data
  # rm(incidence_rate_data)
  # cfr <- case_fatality_ratio_data
  # rm(case_fatality_ratio_data)
  # pars <- parameter_data
  # rm(parameter_data)

  # sobol design
  sobol_seq <-
    pomp::sobol_design(
      lower = c(camp_dur=0, ve_delay=0, dve_U5=0, dve_5p=0, ive_id=0, dw_mild=0,
                dw_mod=0, dw_sev=0, cfr=0, vo=0),
      upper = c(camp_dur=1, ve_delay=1, dve_U5=1, dve_5p=1, ive_id=1, dw_mild=1,
                dw_mod=1, dw_sev=1, cfr=1, vo=1), nseq = nruns)


    # parameter transformation
    # extract central, lower bound, and upper bounds
    # take the central value as the mean
    # determine variance of the beta distribution as the wider
    # interval between central-lower and central-upper to account for
    # uncertainty sufficiently
    # determine the alpha and beta parameter
  if (is.null(parameter_data)){
    parameter_data <- data.table::fread("data/parameters.csv")
  }
  dat <- parameter_data

  dur_vcamp_est <- dat[dat$Parameter == "Duration_Campaign",]$Value
  dur_vcamp_lb <- dat[dat$Parameter == "Duration_Campaign",]$Lower_95
  dur_vcamp_ub <- dat[dat$Parameter == "Duration_Campaign",]$Upper_95
  dur_vcamp_min <- dat[dat$Parameter == "Duration_Campaign",]$Min
  dur_vcamp_max <- dat[dat$Parameter == "Duration_Campaign",]$Max
  dur_vcamp_sd_approx <- (dur_vcamp_ub - dur_vcamp_est) / qnorm(0.975)
  dur_vcamp_sample <- truncnorm::qtruncnorm(sobol_seq$camp_dur, a=dur_vcamp_min,
                                           b=dur_vcamp_max, mean=dur_vcamp_est,
                                           sd=ceiling(dur_vcamp_sd_approx))

  # calculate the standard deviation given mean and p, which is the value you want to have at the lower alpha percentile.
  # assume a truncated normal distribution
  ve_delay_lb <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Lower_95
  ve_delay_min <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Min
  ve_delay_max <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Max
  ve_delay_sd_approx <- (ve_delay_est - ve_delay_lb) / qnorm(0.975)
  ve_delay_sample <- truncnorm::qtruncnorm(sobol_seq$ve_delay, a=ve_delay_min,
                                           b=ve_delay_max, mean=ve_delay_est,
                                           sd=ceiling(ve_delay_sd_approx))


  ## vaccine efficacy
  ss <- cbind(sobol_seq$dve_U5, sobol_seq$dve_5p)

  str <- c("Efficacy_Vaccine_U5", "Efficacy_Vaccine")
  dve_sample <- data.frame(matrix(rep(NA,nrow(ss)*2), ncol=2))

  for (i in 1:2) {
    # vaccine efficacy is
    dve_est <- dat[dat$Parameter == str[i],]$Value
    dve_lb <- dat[dat$Parameter == str[i],]$Lower_95
    dve_ub <- dat[dat$Parameter == str[i],]$Upper_95
    #
    irr_est <- 1 - dve_est
    irr_lb <- 1 - dve_lb
    irr_ub <- 1 - dve_ub
    # make use of that log(IRR) is approximately normal
    logirr_se <- (log(irr_lb) - log(irr_ub))/2/1.96
    # make use of that incidence rate ratio (irr) is approximately normally distributed
    dve_sample[, i] <- 1 - exp(qnorm(ss[,i], mean = log(irr_est), sd = logirr_se))
  }


  p_trans <- list(dur_vcamp = dur_vcamp_sample,
                 ve_delay = ve_delay_sample,
                 dve_U5 = dve_sample[,1],
                 dve_5plus = dve_sample[,2],
                 ive_id = round(sobol_seq$ive_id*nruns))


  params_transformed <- do.call('cbind', p_trans)
  return(params_transformed)
}

draw_parameter_samples <- function(nruns = 200, parameter_data = NULL){

  param_names <- c("run_id",
                 "dur_campaign",
                 "delay_vacc_eff",
                 "vacc_eff_direct_U5",
                 "vacc_eff_direct_5plus",
                 "vacc_eff_indirect_id")

# sobol design
  sobol_seq <-
    pomp::sobol_design(
      lower = c(camp_dur=0, ve_delay=0, dve_U5=0, dve_5p=0, ive_id=0,
                dw_mild=0, dw_mod=0, dw_sev=0, cfr=0, vo=0),
      upper = c(camp_dur=1, ve_delay=1, dve_U5=1, dve_5p=1, ive_id=1,
                dw_mild=1, dw_mod=1, dw_sev=1, cfr=1, vo=1), nseq = nruns)

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
  ve_delay_est <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Value
  ve_delay_lb <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Lower_95
  ve_delay_min <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Min
  ve_delay_max <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Max
  ve_delay_sd_approx <- (ve_delay_est - ve_delay_lb) / qnorm(0.975)
  ve_delay_sample <- truncnorm::qtruncnorm(sobol_seq$ve_delay, a=ve_delay_min,
                                           b=ve_delay_max, mean=ve_delay_est,
                                           sd=ceiling(ve_delay_sd_approx))

  ## samples for vaccine efficacy are drawn making use of that log-transformed incidence rate ratios
  ## are approximately normally distributed
  ## vaccine efficacy
  ss <- cbind(sobol_seq$dve_U5, sobol_seq$dve_5p)

  str <- c("Efficacy_Vaccine_U5", "Efficacy_Vaccine")
  dve_sample <- data.frame(matrix(rep(NA,nrow(ss)*2), ncol=2))

  for (i in 1:2) {
    dve_est <- dat[dat$Parameter == str[i],]$Value
    dve_lb <- dat[dat$Parameter == str[i],]$Lower_95
    dve_ub <- dat[dat$Parameter == str[i],]$Upper_95
    # incidence rate ratio
    irr_est <- 1 - dve_est
    irr_lb <- 1 - dve_lb
    irr_ub <- 1 - dve_ub
    # log(IRR) is approximately normal
    logirr_se <- (log(irr_lb) - log(irr_ub))/2/qnorm(0.975)
    dve_sample[, i] <- 1 - exp(qnorm(ss[,i], mean = log(irr_est), sd = logirr_se))
  }

  # transformed parameters
  p_trans <- list(run_id = 1:nruns,
                  dur_vcamp = dur_vcamp_sample,
                  ve_delay = ve_delay_sample,
                  dve_U5 = dve_sample[,1],
                  dve_5plus = dve_sample[,2],
                  ive_id = round(sobol_seq$ive_id * nruns))

  params_transformed <- as.data.frame(do.call('cbind', p_trans))
  names(params_transformed) <- param_names

  return(params_transformed)
}

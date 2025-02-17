# tstamp
# create the time stamp
tstamp <- function(year=TRUE, month=TRUE, day=TRUE,
                   hour=FALSE, minute=FALSE, second=FALSE) {
  stamp1 <- c()
  stamp2 <- c()
  if (year & !month & !day) {
    stamp <- format(Sys.time(), "%Y")
  } else if (year & month & !day) {
    stamp1 <- format(Sys.time(), "%Y%m")
  } else if (year & month & day) {
    stamp1 <- format(Sys.time(), "%Y%m%d")
  } else if (!year & month & day) {
    stamp1 <- format(Sys.time(), "%m%d")
  } else if (year & !month & day) {
    stamp1 <- format(Sys.time(), "%Y%d")
  } else if (!year & month & !day) {
    stamp1 <- format(Sys.time(), "%m")
  } else if (!year & !month & day) {
    stamp1 <- format(Sys.time(), "%d")
  } else{ stamp1 <- "You'd better select parameters well."}

  if (hour & !minute & !second) {
    stamp2 <- format(Sys.time(), "%H")
  } else if (hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%H%M")
  } else if (hour & minute & second) {
    stamp2 <- format(Sys.time(), "%H%M%S")
  } else if (!hour & minute & !second) {
    stamp2 <- format(Sys.time(), "%M")
  } else if (!hour & !minute & second) {
    stamp2 <- format(Sys.time(), "%S")
  } else if (!hour & minute & second) {
    stamp2 <- format(Sys.time(), "%M%S")
  } else{}

  if (!is.null(stamp2)) {
    stamp1 <- paste0(stamp1, "T", stamp2)
  }
  return (stamp1)
}


# figure size for plotting
figure_size <- data.frame(journal=c("Nature","Elsevier","Lancet"),
                          single=c(89,90,75),
                          double=c(183,190,154),
                          unit2=c("mm","mm","mm"))


# clean_country_names
# Clean and standardize country names according to VIMC report templates
clean_country_names <- function(country){
  for (i in 1:length(country)) {
    if (country[i] %in% c("DR Congo", "Democratic Republic of the Congo",
                          "DRC", "Congo, Dem. Rep.", "Congo, DR",
                          "Congo, the Democratic Republic of the",
                          "Congo - Kinshasa")){
      country[i] <- "Congo, Democratic Republic of the"
    }
    if (country[i] %in% c("Congo, Rep.", "Republic of the Congo", "Congo",
                          "Congo - Brazzaville")){
      country[i] <- "Congo, Republic of the"
    }
    if (country[i] %in% c("São Tomé and Príncipe")){
      country[i] <- "Sao Tome e Principe"
    }
    if (country[i] %in% c("Iran", "Iran, Islamic Rep.",
                          "Iran (Islamic Republic of)")){
      country[i] <- "Iran, Islamic Republic of"
    }
    if (country[i] %in% c("North Korea", "Korea:North", "Korea, DPR", "DPRK",
                          "Democratic People's Republic of Korea",
                          "Korea DPR")){
      country[i] <- "Korea, Democratic People's Republic of"
    }
    if (country[i] %in% c("South Korea", "Korea:South", "Korea, Rep.")){
      country[i] <- "Korea, the Republic of"
    }
    if (country[i] %in% c("Sudan: South")){
      country[i] <- "South Sudan"
    }
    if (country[i] %in% c("Sudan: North")){
      country[i] <- "Sudan"
    }
    if (country[i] %in% c("Venezuela", "Venezuela, RB",
                          "Venezuela (Bolivarian Republic of)")){
      country[i] <- "Venezuela, Bolivarian Republic of"
    }
    if (country[i] %in% c("Tanzania", "United Republic of Tanzania")){
      country[i] <- "Tanzania, United Republic of"
    }
    if (country[i] %in% c("Syria")){
      country[i] <- "Syrian Arab Republic"
    }
    if (country[i] %in% c("Moldova", "Republic of Moldova")){
      country[i] <- "Moldova, Republic of"
    }
    if (country[i] %in% c("CAR")){
      country[i] <- "Central African Republic"
    }
    if (country[i] %in% c("Lao", "Laos", "Lao PDR")){
      country[i] <- "Lao People's Democratic Republic"
    }
    if (country[i] %in% c("US", "USA")){
      country[i] <- "United States of America"
    }
    if (country[i] %in% c("C?te d'Ivoire", "CÃ´te d'Ivoire",
                          "Cì²™te d'Ivoire", "Côte d'Ivoire",
                          "Côte d’Ivoire")){
      country[i] <- "Cote d'Ivoire"
    }
    if (country[i] %in% c("Bolivia", "Bolivia (Plurinational State of)")){
      country[i] <- "Bolivia, Plurinational State of"
    }
    if (country[i] %in% c("Cape Verde")){
      country[i] <- "Cabo Verde"
    }
    if (country[i] %in% c("Micronesia", "Micronesia (Federated States of)")){
      country[i] <- "Micronesia, Federated States of"
    }
    if (country[i] %in% c("Sao Tome e Principe")){
      country[i] <- "Sao Tome and Principe"
    }
    if (country[i] %in% c("Vietnam")){
      country[i] <- "Viet Nam"
    }
    if (country[i] %in% c("Eswatini")){ # to be consistent with other data files
      country[i] <- "Swaziland"
    }
  }
  return (country)
}

# draw_parameter_samples
# draw samples of parameter sets with size of nruns
# First create sobol low discrepancy sequence of random numbers
# Second, adjust the sequence based on the known distribution for the parameter
#
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

# run_vacc_impact
# compute the impact of vaccine - static model
# account for age (< 5 yo vs. 5+ yo)
# This function computes vaccine impact from Week 1 as though the vaccine
# recipients have already acquired vaccine-derived immunity
run_vacc_impact <- function(outbreak_data = NULL, # time series
                            ive_data = NULL,
                            parameters = NULL,
                            age_dist = NULL,
                            vacc_week = NULL,
                            vacc_cov = NULL,
                            runid = NULL) {

  outbreak_ids <- unique(outbreak_data$id)
  # lists to store simulation results
  lst <- vector("list", length(outbreak_ids))
  lst2 <- vector("list", length(vacc_cov))
  lst3 <- vector("list", length(vacc_week))

  p <- parameters[runid, c("dur_campaign", "delay_vacc_eff",
                           "vacc_eff_direct_U5", "vacc_eff_direct_5plus")]
  # indirect effectiveness was pre-calculated and
  # here we only randomly select a column indicated a sample
  ive_col <- parameters$vacc_eff_indirect_id[runid] + 1
  col_name <- paste0("X", ive_col)
  dve <- c(p$vacc_eff_direct_U5, p$vacc_eff_direct_5plus)

  for (k in 1:length(vacc_week)) {
    time_to_vaccination <- (vacc_week[k] - 1)*7 + 3.5 # vaccinations start at mid-day
    week_delay <- round((time_to_vaccination + p$dur_campaign/2 + p$delay_vacc_eff)/7)
    for (j in 1:length(vacc_cov)) {
      for (i in 1:length(outbreak_ids)) {
        d <- vacc_impact(
          data = outbreak_data[outbreak_data$id == outbreak_ids[i],],
          vacc_week = vacc_week[k],
          vacc_cov = vacc_cov[j],
          dve = dve,
          ive_data = ive_data,
          ive_colname = col_name,
          week_delay = week_delay)

        lst[[i]] <- d
      }
      lst2[[j]] <- rbindlist(lst)
    }
    lst3[[k]] <- rbindlist(lst2)
  }

  res <- rbindlist(lst3)
  res$runid <- runid

  return (as.data.frame(res))
}

# run_vacc_impact_weekly
# compute the impact of vaccine - static model
# account for age (< 5 yo vs. 5+ yo)
# This function computes vaccine impact from Week 1 as though the vaccine
# recipients have already acquired vaccine-derived immunity
run_vacc_impact_weekly <- function(outbreak_data = NULL, # time series
                            outbreak_data2 = NULL, # summary data
                            ive_data = NULL,
                            parameters = NULL,
                            age_dist = NULL,
                            vacc_week = NULL,
                            vacc_cov = NULL,
                            runid = NULL) {

  outbreak_ids <- unique(outbreak_data$id)
  # lists to store simulation results
  lst <- vector("list", length(outbreak_ids))
  lst2 <- vector("list", length(vacc_cov))
  lst3 <- vector("list", length(vacc_week))

  p <- parameters[runid, c("dur_campaign", "delay_vacc_eff",
                           "vacc_eff_direct_U5", "vacc_eff_direct_5plus")]
  # indirect effectiveness was pre-calculated and
  # here we only randomly select a column indicated a sample
  ive_col <- parameters$vacc_eff_indirect_id[runid] + 1
  col_name <- paste0("X", ive_col)
  dve <- c(p$vacc_eff_direct_U5, p$vacc_eff_direct_5plus)


  for (k in 1:length(vacc_week)) {
    time_to_vaccination <- (vacc_week[k] - 1)*7 + 3.5 # vaccinations start at mid-day
    week_delay <- round((time_to_vaccination + p$dur_campaign/2 + p$delay_vacc_eff)/7)
    for (j in 1:length(vacc_cov)) {
      for (i in 1:length(outbreak_ids)) {
        d <- vacc_impact_weekly(
          data = outbreak_data[outbreak_data$id == outbreak_ids[i],],
          outbreak_data2 = outbreak_data2,
          outbreak_id = outbreak_ids[i],
          age_dist = age_dist,
          vacc_week = vacc_week[k],
          vacc_cov = vacc_cov[j],
          dve = dve,
          ive_data = ive_data,
          ive_colname = col_name,
          week_delay = week_delay)

        lst[[i]] <- d
      }
      lst2[[j]] <- rbindlist(lst)
    }
    lst3[[k]] <- rbindlist(lst2)
  }

  res <- rbindlist(lst3)
  res$runid <- runid

  return (as.data.frame(res))
}




vacc_impact <- function(data = NULL,
                        vacc_week = NULL,
                        vacc_cov = NULL,
                        dve = NULL,
                        ive_data = NULL,
                        ive_colname = NULL,
                        week_delay = NULL) {

  # yr_ch <- strsplit(as.character(data$date[1]),"-")[[1]][1]
  # # prop under 5 according to year and country
  # prop_U5 <-
  #   age_dist[(age_dist$ISO3 == data$country[1] &
  #               as.character(age_dist$year) == yr_ch), "0-4"]

  vc <- vacc_cov
  # this only accounts for vaccine-induced immunity but not existing in the population
  # this could be another reason for underestimation for the vaccine impact
  # variables that may different across simulations
  vc_eff <- vacc_cov * (prop_U5 * dve[1] + (1-prop_U5) * dve[2])
  ive <- as.numeric(ive_data[ive_data$X1 == round(vacc_cov_eff, digits=2),
                             ive_colname])

  # fraction to be averted via vaccination
  f_avert <- 1 - (vacc_cov*(1-dve)*(1-ive)+(1-vacc_cov)*(1-ive))


  # apply delay to vaccine effectiveness
  frac_averted <- rep(f_avert, nrow(data))


  nr <- nrow(data)
  if (week_delay <= nrow(data)) {
    nr <- week_delay
  }
  frac_averted[1:nr] <- 0


  data$case_wk_averted_U5 <-  data$sCh * frac_averted

  df <- data[1, c("country", "week", "id", "data_id", "date","good_fit",
    "OCV_use")]
  df$week_vaccination = vacc_week
  df$week_delay_to_vacc_eff <- week_delay
  df$sCh_tot <- sum(data$sCh)
  df$confirmed_tot <- outbreak_data2[outbreak_data2$ID_outbreak == outbreak_id,
                                     "total_confirmed_cases"]
  df$death_tot <- data$death_total[1]
  df$vacc_cov <- vacc_cov
  df$vacc_cov_eff <- vacc_cov_eff
  df$ive <- ive
  df$prop_U5 <- prop_U5
  df$sCh_averted <- sum(data$case_wk_averted)

  return(df)
}

vacc_impact_weekly <- function(data = NULL,
                        outbreak_data2 = NULL,
                        outbreak_id = NULL,
                        age_dist = NULL,
                        vacc_week = NULL,
                        vacc_cov = NULL,
                        dve = NULL,
                        ive_data = NULL,
                        ive_colname = NULL,
                        week_delay = NULL) {

  yr_ch <- strsplit(as.character(data$date[1]),"-")[[1]][1]
  # prop under 5 according to year and country
  prop_U5 <-
    age_dist[(age_dist$ISO3 == data$country[1] &
                as.character(age_dist$year) == yr_ch), "0-4"]

  vc <- vacc_cov
  # this only accounts for vaccine-induced immunity but not existing in the population
  # this could be another reason for underestimation for the vaccine impact
  # variables that may different across simulations
  vacc_cov_eff <- vacc_cov * (prop_U5 * dve[1] + (1-prop_U5) * dve[2])
  ive <- as.numeric(ive_data[ive_data$X1 == round(vacc_cov_eff, digits=2),
                             ive_colname])

  # fraction to be averted via vaccination
  fa_U5 <- 1 - (vacc_cov*(1-dve[1])*(1-ive)+(1-vacc_cov)*(1-ive))
  fa_5up <- 1 - (vacc_cov*(1-dve[2])*(1-ive)+(1-vacc_cov)*(1-ive))

  # apply delay to vaccine effectiveness
  frac_averted_U5 <- rep(fa_U5, nrow(data))
  frac_averted_5up <- rep(fa_5up, nrow(data))

  nr <- nrow(data)
  if (week_delay <= nrow(data)) {
    nr <- week_delay
  }
  frac_averted_U5[1:nr] <- 0
  frac_averted_5up[1:nr] <- 0

  data$case_wk_averted_U5 <-  data$sCh * prop_U5 * frac_averted_U5
  data$case_wk_averted_5up <- data$sCh * (1-prop_U5) * frac_averted_5up
  data$case_wk_averted_tot <- data$case_wk_averted_U5 + data$case_wk_averted_5up

  df <- data[, c("country", "week", "id", "data_id", "date","good_fit",
                  "OCV_use")]
  df$week_vaccination <- vacc_week
  df$week_delay_to_vacc_eff <- week_delay
  df$sCh_tot <- data$sCh
  df$confirmed_tot <- outbreak_data2[outbreak_data2$ID_outbreak == outbreak_id,
                                     "total_confirmed_cases"]
  df$death_tot <- data$death_total
  df$vacc_cov <- vacc_cov
  df$vacc_cov_eff <- vacc_cov_eff
  df$ive <- ive
  df$prop_U5 <- prop_U5
  df$sCh_averted_U5 <- data$case_wk_averted_U5
  df$sCh_averted_5up <- data$case_wk_averted_5up
  df$sCh_averted_tot <- data$case_wk_averted_tot

  return(df)
}

# vacc_impact <- function(data = NULL,
#                         outbreak_data2 = NULL,
#                         outbreak_id = NULL,
#                         age_dist = NULL,
#                         vacc_week = NULL,
#                         vacc_cov = NULL,
#                         dve = NULL,
#                         ive_data = NULL,
#                         ive_colname = NULL,
#                         week_delay = NULL) {
#
#   yr_ch <- strsplit(as.character(data$date[1]),"-")[[1]][1]
#   # prop under 5 according to year and country
#   prop_U5 <-
#     age_dist[(age_dist$ISO3 == data$country[1] &
#                 as.character(age_dist$year) == yr_ch), "0-4"]
#
#   vc <- vacc_cov
#   # this only accounts for vaccine-induced immunity but not existing in the population
#   # this could be another reason for underestimation for the vaccine impact
#   # variables that may different across simulations
#   vacc_cov_eff <- vc * (prop_U5 * dve[1] + (1-prop_U5) * dve[2])
#   ive <- as.numeric(ive_data[ive_data$X1 == round(vacc_cov_eff, digits=2),
#                              ive_colname])
#
#   frac_averted_U5 <- 1 - (vc*(1-dve[1])*(1-ive)+(1-vc)*(1-ive))
#   frac_averted_5up <- 1 - (vc*(1-dve[2])*(1-ive)+(1-vc)*(1-ive))
#
#   # categorize by confirmed cases
#   data$week_vaccination <- vacc_week
#   data$week_delay_to_vacc_eff <- week_delay
#
#   data$confirmed <-
#     outbreak_data2[outbreak_data2$ID_outbreak == outbreak_id, "total_confirmed_cases"]
#   data$prop_U5 <- prop_U5
#   # reported cases per each week
#   data$case_wk_tot <- data$sCh
#   data$case_wk_U5 <- data$sCh * prop_U5
#   data$case_wk_5up <- data$sCh * (1-prop_U5)
#   # total cases per outbreak  - stays constant across week
#   data$case_total_U5 <- data$case_total * prop_U5
#   data$case_total_5up <- data$case_total * (1-prop_U5)
#   # cases to occur during the remainder of the outbreak
#   # i.e., decreases over the week
#   data$case_rem_U5 <- data$case_rem * prop_U5
#   data$case_rem_5up <- data$case_rem * (1-prop_U5)
#
#   data$vacc_cov <- vc
#
#   data$vacc_cov_eff <- vacc_cov_eff
#
#   data$ive <- ive
#   # weekly cases
#   data$case_wk_averted_U5 <- data$case_wk_U5 * frac_averted_U5
#   data$case_wk_averted_5up <- data$case_wk_5up * frac_averted_5up
#   data$case_wk_averted_tot <- data$case_wk_averted_U5 + data$case_wk_averted_5up
#
#   # data$case_averted_U5 <- sum(data$case_wk_averted_U5)
#   # data$case_averted_5up <- sum(data$case_wk_averted_5up)
#   # data$case_averted_tot <- sum(data$case_wk_averted_tot)
#   # # cases to be expected during the remainder of the outbreak
#   # data$case_rem_averted_U5 <- data$case_rem_U5 * frac_averted_U5
#   # data$case_rem_averted_5up <- data$case_rem_5up * frac_averted_5up
#   # data$case_rem_averted_tot <- data$case_rem_averted_U5 + data$case_rem_averted_5up
#   # # death only reported as cumulative number during the outbreak
#   # data$death_averted_U5 <- data$death_total * prop_U5 * frac_averted_U5
#   # data$death_averted_5up <- data$death_total * (1-prop_U5) * frac_averted_5up
#   # data$death_averted_tot <- data$death_averted_U5 + data$death_averted_5up
#
#   # apply delay in vaccine effect before calculating percent impact
#   d <- apply_vacc_delay(data=data, week_delay=week_delay)
#
#   return(d)
# }

# apply_vacc_week
# This function generates the impact of the vaccine when it was applied
# at Week 2 or later based on the impact estimates when vaccination occurred at
# Week 1
# The result is to simply remove the impact during the Week 1 and (week - 1)
apply_vacc_week <- function(data, week) {
  if (week <= 1) stop("week must be larger than 1")

  outbreak_ids <- unique(data$id)
  vacc_cov <- unique(data$vacc_cov)

  lst <- vector("list", length(outbreak_ids))
  lst2 <- vector("list", length(vacc_cov))

  for (j in 1:length(vacc_cov)) {
    for (i in 1:length(outbreak_ids)) {
      d <- data[(data$id == outbreak_ids[i] & data$vacc_cov == vacc_cov[j]),]
      wk_diff <- week - d$week_vaccination[1]
      if (wk_diff <= 0) {
        stop("Week different must be positive.")
      }
      d$week_vaccination <- week + wk_diff
      d$week_delay_to_vacc_eff <- d$week_delay_to_vacc_eff + wk_diff

      d <- apply_vacc_delay(d, d$week_delay_to_vacc_eff[1])
      lst[[i]] <- d
    }
    lst2[[j]] <- rbindlist(lst)
  }
  return(as.data.frame(rbindlist(lst2)))
}

# vacc_impact <- function(data = NULL,
#                         outbreak_data2 = NULL,
#                         outbreak_id = NULL,
#                         age_dist = NULL,
#                         vacc_week = NULL,
#                         vacc_cov = NULL,
#                         dve = NULL,
#                         ive_data = NULL,
#                         ive_colname = NULL,
#                         week_delay = NULL) {
#
#   yr_ch <- strsplit(as.character(data$date[1]),"-")[[1]][1]
#   # prop under 5 according to year and country
#   prop_U5 <-
#     age_dist[(age_dist$ISO3 == data$country[1] &
#                 as.character(age_dist$year) == yr_ch), "0-4"]
#
#   vc <- vacc_cov
#   # this only accounts for vaccine-induced immunity but not existing in the population
#   # this could be another reason for underestimation for the vaccine impact
#   # variables that may different across simulations
#   vacc_cov_eff <- vc * (prop_U5 * dve[1] + (1-prop_U5) * dve[2])
#   ive <- as.numeric(ive_data[ive_data$X1 == round(vacc_cov_eff, digits=2),
#                              ive_colname])
#
#   frac_averted_U5 <- 1 - (vc*(1-dve[1])*(1-ive)+(1-vc)*(1-ive))
#   frac_averted_5up <- 1 - (vc*(1-dve[2])*(1-ive)+(1-vc)*(1-ive))
#
#   # categorize by confirmed cases
#   data$week_vaccination <- vacc_week
#   data$week_delay_to_vacc_eff <- week_delay
#
#   data$confirmed <-
#     outbreak_data2[outbreak_data2$ID_outbreak == outbreak_id, "total_confirmed_cases"]
#   data$prop_U5 <- prop_U5
#   # reported cases per each week
#   data$case_wk_tot <- data$sCh
#   data$case_wk_U5 <- data$sCh * prop_U5
#   data$case_wk_5up <- data$sCh * (1-prop_U5)
#   # total cases per outbreak  - stays constant across week
#   data$case_total_U5 <- data$case_total * prop_U5
#   data$case_total_5up <- data$case_total * (1-prop_U5)
#   # cases to occur during the remainder of the outbreak
#   # i.e., decreases over the week
#   data$case_rem_U5 <- data$case_rem * prop_U5
#   data$case_rem_5up <- data$case_rem * (1-prop_U5)
#
#   data$vacc_cov <- vc
#
#   data$vacc_cov_eff <- vacc_cov_eff
#
#   data$ive <- ive
#
#   data$case_wk_averted_U5 <- data$case_wk_U5 * frac_averted_U5
#   data$case_wk_averted_5up <- data$case_wk_5up * frac_averted_5up
#   data$case_wk_averted_tot <- data$case_wk_averted_U5 + data$case_wk_averted_5up
#
#   data$case_rem_averted_U5 <- data$case_rem_U5 * frac_averted_U5
#   data$case_rem_averted_5up <- data$case_rem_5up * frac_averted_5up
#   data$case_rem_averted_tot <- data$case_rem_averted_U5 + data$case_rem_averted_5up
#
#   data$death_averted_U5 <- data$death_total * prop_U5 * frac_averted_U5
#   data$death_averted_5up <- data$death_total * (1-prop_U5) * frac_averted_5up
#   data$death_averted_tot <- data$death_averted_U5 + data$death_averted_5up
#
#   # apply delay in vaccine effect before calculating percent impact
#   d <- apply_vacc_delay(data=data, week_delay=week_delay)
#
#   d$pct_case_averted_U5 = 100 * d$case_rem_averted_U5 / d$case_total_U5
#   d$pct_case_averted_5up = 100 * d$case_rem_averted_5up / d$case_total_5up
#   d$pct_case_averted_tot = 100 * d$case_rem_averted_tot / d$case_total
#
#   d$pct_gain_case_averted_U5 = c(-diff(d$pct_case_averted_U5), NA)
#   d$pct_gain_case_averted_5up = c(-diff(d$pct_case_averted_5up), NA)
#   d$pct_gain_case_averted_tot = c(-diff(d$pct_case_averted_tot), NA)
#
#   return(d)
# }

# create_cumul_vars
# Create three additional variables
# 1. case_total: total number of cases (ie, outbreak size)
# 2. case_cum: cumulative number of cases by week
# 3. case_rem: cases to be expected during the remainder of the outbreak.
# case_rem is used to compute the case averted
# (i.e., case averted  = case_rem * (1-(1-DVE)*(1-IVE))...
#
create_cumul_vars <- function(data1, data2){
  outbreak_ids <- unique(data1$id)
  lst <- vector("list", length(outbreak_ids))
  for (i in 1:length(outbreak_ids)) {
    d <- data1[data1$id == outbreak_ids[i],]
    total_deaths <- data2[data2$ID_outbreak == outbreak_ids[i],]$total_deaths
    d$case_cum <- cumsum(d$sCh)
    d$case_total <- max(d$case_cum)
    d$death_total <- total_deaths
    for (j in 1:nrow(d)) {
      d$case_rem[j] <- ifelse(j==1, d$case_total[1],
                              d$case_total[1] - d$case_cum[j-1])
    }
    lst[[i]] <- d
  }
  data1_cumsum <- as.data.frame(rbindlist(lst))

  return(data1_cumsum)
}

# apply_vacc_delay
# this is to remove the vaccine impact that occurred before the delay
# the overall idea is to compute the vaccine impact as though the vaccine
# takes effect at Week 1 and then follow
#
apply_vacc_delay <- function(data = NULL, week_delay = NULL,
           case_vars = c("case_wk_averted_U5", "case_wk_averted_5up",
                         "case_wk_averted_tot")) {

  nr <- nrow(data)
  if (week_delay < nr) {
    nr <- week_delay
  }

  data[1:nr, case_vars] <- 0

  return(data)
}

# apply_vacc_delay <-
#   function(data=NULL, # vaccine impact estimates
#            week_delay=NULL,
#            case_vars = c("case_wk_averted_U5", "case_wk_averted_5up",
#                          "case_wk_averted_tot",  "case_rem_averted_U5",
#                          "case_rem_averted_5up", "case_rem_averted_tot"),
#            death_vars = c("death_averted_U5",
#                           "death_averted_5up", "death_averted_tot")) {
#
#   nr <- nrow(data)
#   if (week_delay < nr) {
#     cols <- c(case_vars, death_vars)
#     data[1:(nr-week_delay), c(case_vars, death_vars)] <-
#       data[(week_delay+1):nr, cols]
#     data[(nr-week_delay+1):nr, case_vars] <- 0
#
#     for (wk in 1:week_delay) {
#       for (dv in death_vars) {
#         if(!is.na(data[wk+(nr-week_delay), dv])){
#           data[wk, dv] <- 0
#         }
#       }
#     }
#   }
#   else {
#     data[, case_vars] <- 0
#
#     for (wk in 1:nr) {
#       for (dv in death_vars) {
#         if(!is.na(data[wk, dv])){
#           data[wk, dv] <- 0
#         }
#       }
#     }
#   }
#
#   return(data)
# }



# adjust_week_pct_effect

# This function fills up the row of vaccine effect of zero for weeks that follow
# the end of the outbreak up to the `week_final`. # For instance, if the outbreak
# ended around Week 8 and we are exploring the impact of the vaccination
# from Week 0 and Week 16, then Week 9 - Week 16, the impact of the vaccine
# should be zero instead of a missing value.
# That is, setting the effect of the vaccine introduced after the outbreak is zero.
# In practice, the function adds empty rows for the outbreaks that ended before
# the week_final parameter such that summary statistics correctly account for
# weeks after the outbreak

# The following columns are NA's at the las row and there don't have to be set to zero
# c("pct_gain_case_averted_U5",
#   "pct_gain_case_averted_5up",
#   "pct_gain_case_averted_tot")

adjust_week_pct_effect <- function(dat, week_final) {
  # these columns will be zero when new rows are added
  case_vars <-
    c("sCh", "case_wk_tot", "case_wk_U5", "case_wk_5up",
      "case_rem", "case_rem_U5", "case_rem_5up",
      "case_wk_averted_U5" , "case_wk_averted_5up", "case_wk_averted_tot",
      "case_rem_averted_U5", "case_rem_averted_5up", "case_rem_averted_tot",
      "pct_case_averted_U5", "pct_case_averted_5up", "pct_case_averted_tot")

  death_vars <- c("death_averted_U5", "death_averted_5up", "death_averted_tot")

  idx <- unique(dat$data_id)
  vc <- unique(dat$vacc_cov)
  lst2 <- vector("list", length(vc))
  for (j in 1:length(vc)) {
    cat("Filling the rows for the vaccine coverage =", vc[j], "\n")
    lst <- vector("list", length(idx))
    for (i in 1:length(idx)) {
      d <- dat[(dat$data_id == idx[i] & dat$vacc_cov == vc[j]),]
      maxwk <- max(d$week)
      if (maxwk < week_final) {
        nr <- week_final - maxwk
        for (k in 1:nr) {
          last_row <- d[nrow(d),]
          if (k == 1) {
            # last_row[,c(6,12,16:37)] <- 0
            last_row[, case_vars] <- 0
            for (dv in death_vars) {
              if(!is.na(last_row[, dv])){
                last_row[, dv] <- 0
              }
            }
          }
          last_row$week <- last_row$week + 1
          last_row$date <- last_row$date + 7
          d <- rbind(d, last_row)
        }
      }
      lst[[i]] <- d
    }
    lst2[[j]] <- data.table::rbindlist(lst)
  }
  newdata <- data.table::rbindlist(lst2)

  return(as.data.frame(newdata))
}



compute_yld <- function(cases, parms=NULL, option=1){
  if (is.null(parms)) {
    stop("Parameters to compute YLD must be provided")
  }
  pa <- parms[parms$Parameter == "Prop_Asymptomatic", "Value"]
  pmild <- parms[parms$Parameter == "Prop_Mild", "Value"]
  pm <- parms[parms$Parameter == "Prop_Moderate", "Value"]
  ps <- parms[parms$Parameter == "Prop_Severe", "Value"]
  dur <- parms[parms$Parameter == "Duration_Illness", "Value"]
  wa <- parms[parms$Parameter == "Disability_Weight_Asymptomatic", "Value"]
  wmild <- parms[parms$Parameter == "Disability_Weight_Mild", "Value"]
  wm <- parms[parms$Parameter == "Disability_Weight_Moderate", "Value"]
  ws <- parms[parms$Parameter == "Disability_Weight_Severe", "Value"]

  # dis_wt_tot <- sum(c(pa, pmild, pm, ps) * c(wa, wmild, wm, ws))
  # to account for moderate and severe cases
  dis_wt_tot <- pm/(pm+ps)*wm + ps/(pm+ps)*ws #0.2163957
  if (option == 2) {
    dis_wt_tot <- pmild/(pmild+pm+ps)*wmild + pm/(pmild+pm+ps)*wm +
      ps/(pmild+pm+ps)*ws #0.1780156
  }
  #   disability weight for severe diarrheal disease as estimated from the 2016 IHME Global Burden of Disease
  # Study (0.247) (Institute for Health Metrics and Evaluation, Global Burden of Disease Collaborative Network
  # 2017).
  yld_per_capita <- (dur/365) * dis_wt_tot
  yld_tot <- cases * yld_per_capita

  return(yld_tot)
}


get_life_expectancy <- function(life_exp_data, age, data){
  life_exp_data |>
    dplyr::filter(`ISO3 Alpha-code` == data$country,
                  Year == as.character(data.table::year(data$date))) |>
    dplyr::select(as.character(age)) |>
    as.numeric() -> life_exp

  life_exp_data[(life_exp_data$`ISO3 Alpha-code` == cntry & life_exp_data$Year == year),
                as.character(age) ]
  return(life_exp)
}

compute_yll <- function(deaths,
                        life_exps=NULL,
                        parms=NULL) {

  dr <- parms[parms$Parameter == "Rate_Discount", "Value"]
  yll <- deaths * (1/dr) * (1 - exp(-dr * life_exps))

  return(yll)
}



theme_2D_contour <- function(){
  theme(panel.background=element_blank(),
        plot.background=element_blank(),
        # panel.border=element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        # axis.line.x=element_line(colour="black"),
        # axis.line.y=element_line(colour="black"),
        axis.text=element_text(size=11),
        axis.title=element_text(size=11),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vetical",
        legend.key.size = unit(1.2,"mm"),
        legend.text = element_text(size=10),
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.margin = margin(-2, 0, 0, 0, unit="mm"),
        legend.title = element_text(face="plain",size=10))
}

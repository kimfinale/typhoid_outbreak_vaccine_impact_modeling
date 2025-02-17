# draw_parameter_samples
# draw samples of parameter sets with size of nruns
# First create sobol low discrepancy sequence of random numbers
# Second, adjust the sequence based on the known distribution for the parameter
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
  # dur_vcamp_lb <- dat[dat$Parameter == "Duration_Campaign",]$Lower_95
  # dur_vcamp_ub <- dat[dat$Parameter == "Duration_Campaign",]$Upper_95
  dur_vcamp_min <- dat[dat$Parameter == "Duration_Campaign",]$Min
  dur_vcamp_max <- dat[dat$Parameter == "Duration_Campaign",]$Max
  dur_vcamp_sd <- dat[dat$Parameter == "Duration_Campaign",]$SD
  # dur_vcamp_sd_approx <- (dur_vcamp_ub - dur_vcamp_est) / qnorm(0.975)

  dur_vcamp_max <- Inf

  dur_vcamp_sample <- truncnorm::qtruncnorm(sobol_seq$camp_dur, a=dur_vcamp_min,
                                            b=dur_vcamp_max, mean=dur_vcamp_est,
                                            sd=dur_vcamp_sd)


  # calculate the standard deviation given mean and p, which is the value you want to have at the lower alpha percentile.
  # assume a truncated normal distribution
  # ve_delay_est <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Value
  # ve_delay_lb <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Lower_95
  ve_delay_min <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Min
  ve_delay_max <- dat[dat$Parameter == "Delay_Vacc_Eff",]$Max
  # ve_delay_sd_approx <- (ve_delay_est - ve_delay_lb) / qnorm(0.975)
  # this is an approximation max and min of the sample captured with +- 3 sd
  # ve_delay_sd_approx <- ((ve_delay_max - ve_delay_min)*.997)/6
  # # compare the following two samples
  # xx <- rnorm(1e5, mean=3, sd=4)
  # xx2 <- truncnorm::rtruncnorm(1e5, a=min(xx),
  # b=max(xx), mean=mean(xx), sd=((max(xx) - min(xx))*0.997)/6)
  # xx2 <- truncnorm::rtruncnorm(1e5, a=quantile(xx, probs=0.015), b=quantile(xx, probs=0.985), mean=mean(xx), sd=((max(xx) - min(xx))*0.997)/6)

  # ve_delay_sample <- truncnorm::qtruncnorm(sobol_seq$ve_delay, a=ve_delay_min,
  #                                          b=ve_delay_max, mean=ve_delay_est,
  #                                          sd=ceiling(ve_delay_sd_approx))

  ve_delay_sample <- qunif(sobol_seq$ve_delay,
                           min=ve_delay_min, max=ve_delay_max)

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
  ive_data <- readRDS("outputs/ive_yrep_20241010.rds")
  n_ive_sample <- nrow(ive_data$yrep)

  p_trans <- list(run_id = 1:nruns,
                  dur_vcamp = dur_vcamp_sample,
                  ve_delay = ve_delay_sample,
                  dve_U5 = dve_sample[,1],
                  dve_5plus = dve_sample[,2],
                  ive_id = round(sobol_seq$ive_id * n_ive_sample))

  params_transformed <- as.data.frame(do.call('cbind', p_trans))
  names(params_transformed) <- param_names

  return(params_transformed)
}


# compute the vaccine impact for an outbreak with weekly time series
# data refer to a single outbreak
vacc_impact_outbreak_weekly <- function(data = NULL,
                                        # age_dist = NULL,
                                        vacc_week = NULL,
                                        case_trigger = NULL,
                                        vacc_cov = NULL,
                                        dve = NULL,
                                        ive_data = NULL,
                                        ive_rownum = NULL,
                                        week_delay = NULL) {

  # yr_ch <- as.character(data.table::year(data$start_date[1]))
  # # prop under 5 varies by country and year
  # prop_U5 <-
  #   age_dist[(age_dist$ISO3 == data$country[1] &
  #               as.character(age_dist$year) == yr_ch), "0-4"]
  increment <- 0.02 # increment used in posterior predictive values
  ive <- as.numeric(ive_data[ive_rownum, (round(vacc_cov/increment)+1)])

  ive <- 0.7
  # output dataframe
  df <- data[, c("study_id", "time_since_outbreak", "id", "suspected",
                 "confirmed", "no_of_patients")]

  df$week_vaccination <- vacc_week
  df$week_delay_to_vacc_eff <- week_delay
  df$case_trigger <- case_trigger
  df$vacc_cov <- vacc_cov
  df$ive <- ive

  # vaccination happens after outbreak ended
  # if (week_delay > nrow(data)) {
  if (week_delay > max(data$time_since_outbreak)) {
    df$suspected_averted <- 0
  }
  else {
    # fraction to be averted via vaccination
    fa <- 1 - (vacc_cov*(1-dve)*(1-ive)+(1-vacc_cov)*(1-ive))
    # apply delay to vaccine effectiveness
    frac_averted <- rep(fa, nrow(data))

    if (week_delay >= 1) { # if week_delay <= 0, pre-emptive vacc scenario
      frac_averted[1:week_delay] <- 0
    }
    # update output dataframe
    df$suspected_averted <- data$suspected * frac_averted
  }

  return(df)
}


# run_vacc_impact_outbreak_weekly
# compute the impact of vaccine - static model
# This function computes vaccine impact from Week 1 as though the vaccine
# recipients have already acquired vaccine-derived immunity
run_vacc_impact_outbreak_weekly <- function(outbreak_data = NULL, # time series
                                            ive_data = NULL,
                                            parameters = NULL,
                                            age_dist = NULL,
                                            vacc_week = NULL,
                                            case_trigger = NULL,
                                            vacc_cov = NULL,
                                            runid = NULL) {

  outbreak_ids <- unique(outbreak_data$id)
  # lists to store simulation results
  lst <- vector("list", length(outbreak_ids))
  lst2 <- vector("list", length(vacc_cov))
  lst3 <- vector("list", length(vacc_week))

  p <- parameters[runid, c("dur_campaign", "delay_vacc_effect",
                           "vacc_efficacy_direct",
                           "vacc_efficacy_indirect_id")]
  # indirect effectiveness was pre-calculated and
  # here we only randomly select a column indicated a sample
  # ive_col <- parameters$vacc_eff_indirect_id[runid] + 1
  # col_name <- paste0("X", ive_col)
  dve <- p$vacc_effect_direct

  if (!is.null(vacc_week) & is.null(case_trigger)) {
    for (k in 1:length(vacc_week)) {
      time_to_vaccination <- (vacc_week[k] - 1)*7 + 3.5 # vaccinations start at mid-day
      week_delay <- round((time_to_vaccination + p$dur_campaign/2 + p$delay_vacc_effect)/7)
      for (j in 1:length(vacc_cov)) {
        for (i in 1:length(outbreak_ids)) {
          d <- vacc_impact_outbreak_weekly(
            data = outbreak_data[outbreak_data$id == outbreak_ids[i],],
            # age_dist = age_dist,
            vacc_week = vacc_week[k],
            case_trigger = case_trigger,
            vacc_cov = vacc_cov[j],
            dve = dve,
            ive_data = ive_data,
            ive_rownum = p$vacc_efficacy_indirect_id,
            week_delay = week_delay)

          lst[[i]] <- d
        }
        lst2[[j]] <- rbindlist(lst)
      }
      lst3[[k]] <- rbindlist(lst2)
    }
  }
  else if (!is.null(case_trigger) & is.null(vacc_week)) {
    for (k in 1:length(case_trigger)) {
      for (j in 1:length(vacc_cov)) {
        for (i in 1:length(outbreak_ids)) {
          data <- outbreak_data[outbreak_data$id == outbreak_ids[i],]
          cum_case <- cumsum(data$sCh)
          vacc_week <- min(which(cum_case >= case_trigger[k]))
          time_to_vaccination <- (vacc_week - 1)*7 + 3.5 # vaccinations start at mid-day
          week_delay <- round((time_to_vaccination + p$dur_campaign/2 + p$delay_vacc_effect)/7)
          d <- vacc_impact_outbreak_weekly(
            data = data,
            # age_dist = age_dist,
            vacc_week = vacc_week,
            case_trigger = case_trigger[k],
            vacc_cov = vacc_cov[j],
            dve = dve,
            ive_data = ive_data,
            ive_rownum = p$vacc_efficacy_indirect_id,
            week_delay = week_delay)

          lst[[i]] <- d
        }
        lst2[[j]] <- rbindlist(lst)
      }
      lst3[[k]] <- rbindlist(lst2)
    }
  }
  res <- rbindlist(lst3)
  res$runid <- runid

  return (as.data.frame(res))
}

# compute_yld
# compute years of life lost due to disabilities
compute_yld <- function(cases, parms=NULL, option=1){
  if (is.null(parms)) {
    stop("Parameters to compute YLD must be provided")
  }
  # pm <- parms[parms$Parameter == "Prop_Moderate",]$Value
  # ps <- parms[parms$Parameter == "Prop_Severe",]$Value
  # ps_gi <- parms[parms$Parameter == "Prop_Severe_GI_Bleeding",]$Value
  # ps_perf <- parms[parms$Parameter == "Prop_Severe_Other_Complications",]$Value
  #
  # durm <- parms[parms$Parameter == "Dur_Illness_Moderate",]$Value
  # durs <- parms[parms$Parameter == "Dur_Illness_Severe",]$Value
  # durs_gi <- parms[parms$Parameter == "Dur_Illness_Severe_GI_Bleeding",]$Value
  # durs_perf <- parms[parms$Parameter == "Dur_Illness_Severe_Other_Complications",]$Value
  #
  # wm <- parms[parms$Parameter == "Disability_Weight_Moderate",]$Value
  # ws <- parms[parms$Parameter == "Disability_Weight_Severe",]$Value
  # ws_gi <- parms[parms$Parameter == "Disability_Weight_Severe_GI_Bleeding",]$Value
  # ws_perf <- parms[parms$Parameter == "Disability_Weight_Severe_Ileal_Perforation",]$Value
  #
  # prop_tot <- pm + ps + ps_gi + ps_perf
  # dis_wt_tot <- pm/(prop_tot)*wm + ps/(prop_tot)*ws + ps_gi/(prop_tot)*ws_gi + ps_perf/(prop_tot)*ws_perf
  #
  # yld_weighted_mean <- durm*pm/(prop_tot)*wm + durs*ps/(prop_tot)*ws +
  #   durs_gi*ps_gi/(prop_tot)*ws_gi + durs_perf*ps_perf/(prop_tot)*ws_perf

  # You might execute the codes above to get the weighted mean YLD
  yld_weighted_mean <- parms[parms$Parameter == "YLD_Weighted_Mean",]$Value
  yld_tot <- cases * yld_weighted_mean

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



# extract admin names from id
get_adm <- function(x) {
  # d <- grep("[a-zA-Z+]", strsplit(x, "::|-")[[1]], value=TRUE)
  d <- lapply(x, function(z) grep("[a-zA-Z+]", strsplit(z, "::|-")[[1]], value=TRUE))
  return (d)
}


# extract outbreak ids that may not occur because of the lasting immunity from vaccination occurred during the previous outbreak
# lb and ub is the lower and upper bounds of the duration of the immunity
# we only examine outbreaks that occur in the same admin 2 locations
find_ghost_outbreaks <- function(x, lb=0, ub=1000) {
  mean_dur_campaign <- 6 # days duration of the vaccination campaign
  mean_delay_vacc_eff <- 11 # days, mean delay before the immunity kicks in
  # time to vaccination since the start of the outbreak
  # it assumes that the vaccination starts at mid-day; hence, 3.5 days
  # time_to_vaccination <- (week_vaccination - 1) * 7 + 3.5 # days
  # date immunity: time when the population in this locality acquire immunity
  # full duration of a campaign is added because it assumes that preventing
  # next outbreak is possible when the entire target population is covered
  x |>
    mutate(date_immunity = start_date + (week_vaccination - 1) * 7 + 3.5 +
             mean_dur_campaign + mean_delay_vacc_eff) -> x1

  x1 |>
    group_by(admin0, admin1, admin2) |>
    reframe(min_date_immunity = min(date_immunity)) -> x2

  x1 <- left_join(x1, x2, by=c("admin0", "admin1", "admin2"))
  x1 <- mutate(x1, date_diff = date_immunity - min_date_immunity)

  # identify the first outbreak in the region (date_diff = 0) and the outbreaks
  # that could disappear (ghost outbreaks, lb <= date_diff <= ub)
  ghost_id <- c()
  ghost_id <- c(ghost_id, x1[x1$date_diff > lb & x1$date_diff <= ub,]$data_id)

  # outbreaks intervened (i.e., date_diff=0) and the outbreaks that could have
  # been protected (i.e, lb <= date_diff <= ub) are removed

  x1 <- x1[x1$date_diff > ub, ] # extract outbreaks that would still occur
  # repeat the process of selecting the first outbreak in the region and
  # identify the ghost outbreaks
  while (!is.null(x1) & nrow(x1) > 1) {
    x1 |>
      group_by(admin0, admin1, admin2) |>
      reframe(min_date_immunity = min(date_immunity)) -> x2

    # x1 <- left_join(x1[, names(x1)[!names(x1) %in% c("min_date_immunity")]], x2,
    #                 by=c("admin0", "admin1", "admin2"))
    x1 <- left_join(subset(x1, select=-c(min_date_immunity)), x2,
                    by=c("admin0", "admin1", "admin2"))
    x1 <- mutate(x1, date_diff = date_immunity - min_date_immunity)
    ghost_id <- c(ghost_id, x1[x1$date_diff > lb & x1$date_diff <= ub,]$data_id)

    x1 <- x1[x1$date_diff > ub, ]
  }

  return (ghost_id)
}


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

# final epidemic size
# the simplest SIR model
# R = 1 - exp(-R0*R) where R is the final epidemic size (or R(\infty) for the SIR model)
final_epidemic_size <- function(R0 = 2) {
  y = function(x) x - 1 + exp(-R0*x)
  final_size <- uniroot(y, interval=c(1e-6,1-1e-6))$root

  return(final_size)

}

# print parameter values
print_params <- function(params){
  n <- names(params)
  for(i in seq_along(n)){
    cat(paste0(n[i], "=", params[n[i]]), ", ")
  }
}

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) {
  # if(x>=1 | x<=0){
  #   stop("The function can take the numbers between 0 and 1 with both exclusive")
  # }
  log(x/(1-x))
}

get_output_days <- function(x) {
  unitdays = gsub("^([0-9]+).*", "\\1", x$date_range)
  l = length(unique(unitdays))
  if (l != 1){
    stop("Output days have more than one kind")
  }
  return(as.double(unique(unitdays)))
}
#
my_discrete_colors <-
  c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
    "#FF7F00","black", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
    "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")


# tstamp
# create the time stamp yearmonthday by default and hour, minute, and second can be added
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

print_params <- function(params){
  n <- names(params)
  str <- paste0(n, "=", params[n])
  print(str)
}

# calculates alpha and beta parameters for the Beta distribution
calc_beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


# calculates alpha and beta parameters for the Beta distribution
calc_lognorm_params <- function(mean, sd){
  v <- sd*sd
  m <- mean
  phi <- sqrt(v + m*m);
  mu <- log(m*m/phi);                # mean of log(Y)
  sigma <- sqrt(log(phi*phi/(m*m))); # std dev of log(Y)

  return(list(mu = mu, sigma = sigma))
}

calc_sd_lognorm <- function(mu, p, alpha=0.025, max=10) {
  eq <- function(x) exp(mu + x * qnorm(alpha)) - p
  uniroot(eq, interval=c(0,max))$root
}

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
quan_func <- function(x) exp(mu + sqrt(2*sigma*sigma) / erf(2*p-1))
# MLE grid search -  1D ---------------------------------------------------------

## copied from https://stackoverflow.com/questions/29067916/error-function-erfz
## if you want the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

MLE_check <- function(p_name = "local_rep_prop", theta_tab,nn=1e3){

  # theta_tab <- seq(0.001,0.01,0.001)
  store_lik <- NULL

  for(ii in 1:length(theta_tab)){

    theta[[p_name]] <- theta_tab[ii]

    # Run SMC and output likelihoods
    output_smc <- smc_model(theta,
                            nn=1e3 # number of particles
    )
    store_lik <- rbind(store_lik,c(theta_tab[ii],output_smc$lik))

  }

  colnames(store_lik) <- c("param","lik")
  store_lik <- as_tibble(store_lik)

}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_2D <- function(p1_name = "local_rep_prop", p2_name = "confirmed_prop",
                         theta_tab1, theta_tab2,nn=1e3,
                         filename = 1){

  # p1_name = "local_rep_prop"; p2_name = "confirmed_prop"; theta_tab1 = seq(0.01,0.05,0.01); theta_tab2 = seq(0.3,1,0.1)

  store_lik <- NULL

  for(ii in 1:length(theta_tab1)){

    for(jj in 1:length(theta_tab2)){

      theta[[p1_name]] <- theta_tab1[ii]
      theta[[p2_name]] <- theta_tab2[jj]

      # Run SMC and output likelihooda
      output_smc <- smc_model(theta,
                              nn=1e3 # number of particles
      )
      store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],output_smc$lik))

    }
  }

  colnames(store_lik) <- c("param1","param2","lik")
  store_lik <- as_tibble(store_lik)

  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))

}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_3D <- function(p1_name = "local_rep_prop",
                         p2_name = "confirmed_prop",
                         p3_name = "betavol",
                         theta_tab1,
                         theta_tab2,
                         theta_tab3,
                         nn=1e3,
                         filename = 1){

  # p1_name = "local_rep_prop"; p2_name = "confirmed_prop"; p3_name = "betavol"; theta_tab1 = seq(0.01,0.05,0.02); theta_tab2 = seq(0.6,1,0.2); theta_tab3 = seq(0.1,0.3,0.1)

  store_lik <- NULL

  out_fit <- foreach(ii = 1:length(theta_tab1)) %dopar% {
    #for(ii in 1:length(theta_tab1)){
    for(jj in 1:length(theta_tab2)){
      for(kk in 1:length(theta_tab3)){
        theta[[p1_name]] <- theta_tab1[ii]
        theta[[p2_name]] <- theta_tab2[jj]
        theta[[p3_name]] <- theta_tab3[kk]
        # Run SMC and output likelihooda
        output_smc <- smc_model(theta,
                                nn=1e3 # number of particles
        )
        store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],theta_tab3[kk],output_smc$lik))
      }
    }
    store_lik
  }
  # Collate results
  store_lik <- NULL
  for(ii in 1:length(theta_tab1)){

    store_lik <- rbind(store_lik,out_fit[[ii]])

  }

  colnames(store_lik) <- c("param1","param2","param3","lik")
  store_lik <- as_tibble(store_lik)

  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))

}

# Compute acceptance probability ------------------------------------------



figure_size <- data.frame(journal=c("Nature","Elsevier","Lancet"),
                          single=c(89,90,75),
                          double=c(183,190,154),
                          unit2=c("mm","mm","mm"))




seir_jl <- function(u, p, t){

  S <- u[1];
  E <- u[2];
  # the number of compartments to model the duration of infectiousness
  I <- u[3];
  R <- u[4];
  C <- u[5];

  # population size
  pop <- S+E+I+R

  epsilon <- p[1] # 1/latent period
  gamma <- p[2] # 1/duration of infectiousness
  beta <- p[3] # transmission rate
  omega <- p[4] # 1/omega = duration of natural immunity
  # force of infection
  foi <- beta*I/pop

  muEI <- epsilon
  muIR <- gamma
  muRS <- omega

  # differential equations
  dS <- - foi*S + muRS*R
  dE <- foi*S - muEI*E
  dI <- muEI*E - muIR*I
  dR <- muIR*I - muRS*R

  dC <- muEI*E

  return(c(dS,dE,dI,dR,dC))
}

run_seir_jl <- function(model=NULL,
                        epsilon=1/4,
                        gamma=1/7,
                        omega=1/(4*365),
                        R0=2,
                        tend=100,
                        saveat=1,
                        nms=c("t","S","E","I","R","C")){
  # library(diffeqr)
  de <- diffeqr::diffeq_setup()
  # simulation
  # R0 <- 2 # basic reproduction number
  # epsilon <- 1/4 # 1/epsilon = incubation period
  # gamma <- 1/7 # 1/gamma = duration of infectiousness
  beta <- R0*gamma # instantaneous transmission rate
  # omega <- 1/(4*365) # natural immunity waning rate
  # parameters
  params <- c(epsilon=epsilon, gamma=gamma, beta=beta, omega=omega)
  u0 <- c(0.99, 0, 0.01, 0, 0)
  tend <- 100 #
  tspan <- c(0.0, tend)

  prob <- de$ODEProblem(model, u0, tspan, params)
  sol <- de$solve(prob, de$Tsit5(), saveat=saveat)
  mat <- sapply(sol$u, identity)
  udf <- as.data.frame(t(mat))
  out <- cbind(data.frame(t=sol$t), udf)
  names(out) <- nms

  return(out)
}



# ggplot2 themes


function (base_size = 11, base_family = "", base_line_size = base_size/22,
          base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family,
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid = element_line(colour = "grey92"),
          panel.grid.minor = element_line(linewidth = rel(0.5)),
          strip.background = element_rect(fill = "grey85",
                                          colour = "grey20"),
          legend.key = element_rect(fill = "white",
                                    colour = NA), complete = TRUE)
}

theme_pub <- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title=element_text(face="bold",
                                    size=rel(1.2),
                                    hjust=0.5,
                                    margin=margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face="plain", size=rel(1)),
            axis.title.y = element_text(angle=90, vjust=2),
            axis.title.x = element_text(vjust=-0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(3,"mm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="plain"),
            plot.margin=unit(c(3,3,3,3),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="plain")
    ))
}

library(ggplot2)
library(grid)

# define consistent ggplot theme to apply to all figures
theme_ms <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1),
            panel.grid=element_blank()
      ))
}

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)

}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)

}


### Dark theme for ggplot plots

theme_dark_grey <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = 'grey20'),
            plot.background = element_rect(colour = NA, fill = '#262626'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'grey70'),
            axis.line.x = element_line(colour="grey70"),
            axis.line.y = element_line(colour="grey70"),
            axis.ticks = element_line(colour="grey70"),
            panel.grid.major = element_line(colour="#262626"),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill ='#262626'),
            legend.text = element_text(color = 'white'),
            legend.key = element_rect(colour = NA, fill = '#262626'),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'white'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
            strip.text = element_text(face="bold", colour = 'white')
    ))
}

scale_fill_Publication_dark <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)

}

scale_colour_Publication_dark <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#f2f2f2")), ...)

}


# theme_transparent <- function(base_size=14, base_family="sans") {
#    library(grid)
#    library(ggthemes)
#    (theme_foundation(base_size=base_size, base_family=base_family)
#       + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
#                                         size = rel(1.2), hjust = 0.5),
#               text = element_text(),
#               panel.background = element_rect(colour = NA, fill = 'transparent'),
#               plot.background = element_rect(colour = NA, fill = 'transparent'),
#               panel.border = element_rect(colour = NA),
#               axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
#               axis.title.y = element_text(angle=90,vjust =2),
#               axis.title.x = element_text(vjust = -0.2),
#               axis.text = element_text(colour = 'grey70'),
#               axis.line.x = element_line(colour="grey70"),
#               axis.line.y = element_line(colour="grey70"),
#               axis.ticks = element_line(colour="grey70"),
#               panel.grid.major = element_line(colour="#262626"),
#               panel.grid.minor = element_blank(),
#               legend.background = element_rect(fill = 'transparent'),
#               legend.text = element_text(color = 'white'),
#               legend.key = element_rect(colour = NA, fill = 'grey20'),
#               legend.position = "bottom",
#               legend.direction = "horizontal",
#               legend.box = "vetical",
#               legend.key.size= unit(0.5, "cm"),
#               #legend.margin = unit(0, "cm"),
#               legend.title = element_text(face="italic", colour = 'white'),
#               plot.margin=unit(c(10,5,5,5),"mm"),
#               strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
#               strip.text = element_text(face="bold", colour = 'white')
#       ))
# }

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

theme_dark_blue <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA, fill = '#282C33'),
            plot.background = element_rect(colour = NA, fill = '#282C33'),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(colour = 'grey70'),
            axis.line.x = element_line(colour="grey70"),
            axis.line.y = element_line(colour="grey70"),
            axis.ticks = element_line(colour="grey70"),
            panel.grid.major = element_line(colour="#343840"),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill ='#282C33'),
            legend.text = element_text(color = 'white'),
            legend.key = element_rect(colour = NA, fill = '#282C33'),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic", colour = 'white'),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#2D3A4C",fill="#2D3A4C"),
            strip.text = element_text(face="bold", colour = 'white')
    ))
}



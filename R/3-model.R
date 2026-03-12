############ sobol sequence simualtion
library(TruncExpFam)
library(randtoolbox)
library(truncnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

n_samples <- 200

sobol_matrix <- sobol(n = n_samples, dim = 3, scrambling = 0, seed = 123)

# scale <- (2^2) / 8 # 2
# rate <- 8 / 2^2
# shape_cd <- (mean_cd / sd_cd)^2  # 16

param_sets <- data.frame(
  dve         = qtruncnorm(sobol_matrix[, 1], a = 0, b = 1, mean = 0.83,  sd = 0.087),
  campaign_duration = qtruncgamma(sobol_matrix[, 2], a = 3, b = 21, shape = 16, scale = 0.5, rate = 2),
  immuno_delay = qunif(sobol_matrix[, 3], min = 14, max = 28)
)
# summary(param_sets$campaign_duration)
# hist(param_sets$campaign_duration, breaks = 50, main = "Campaign Duration", col = "skyblue")
# 
# param_sets_long <- param_sets %>%
#   pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")
# 
# ggplot(param_sets_long, aes(x = value)) +
#   geom_histogram(bins = 50, fill = "steelblue", color = "white") +
#   facet_wrap(~parameter, scales = "free") +
#   theme_minimal()
# 
# 
# summary(param_sets)
# 
# ive <- 0.17
# ?sobol_matrix

########### weekly outbreaks model

ts_weekly <- ts1 %>%
  filter(time_unit == "weekly") %>%
  group_by(study_id) %>%
  mutate(week = row_number(),
         date = date_mid)

outbreak_ids <- unique(ts_weekly$study_id)
# 
# vacc_covs <- seq(0.6, 0.95, by = 0.05)
# vacc_weeks <- seq(1, 52, by = 1) 

vacc_covs <- 0.80
vacc_weeks <- 8

vacc_impact_outbreak_weekly <- function(data = NULL,
                                        age_dist = NULL,
                                        vacc_week = NULL,
                                        vacc_cov = NULL,
                                        week_delay = NULL,
                                        dve = NULL,
                                        ive = NULL, 
                                        vacc_supply = NULL,
                                        runid = NULL) {
  
  yr_ch <- as.character(data.table::year(data$date[1]))
  
  filtered_data <- age_dist[age_dist$ISO3 == data$country_iso[1] & 
                              as.character(age_dist$Year) == yr_ch, ]
  filtered_data <- na.omit(filtered_data)
  # total_population <- filtered_data$total_population[1] 
  # prop_u2 <- sum(filtered_data[, c("0", "1")], na.rm = TRUE) / total_population
  # prop_2_to_4 <- sum(filtered_data[, c("2", "3", "4")], na.rm = TRUE) / total_population
  # prop_5_to_16 <- sum(filtered_data[, c("5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")],
  #                     na.rm = TRUE) / total_population
  # prop_over16 <- sum(filtered_data[, as.character(17:100)], na.rm = TRUE) / total_population  
  # 
  df <- data[, c("country_iso", "week", "study_id", "date",
                 "suspected_cases", "confirmed_cases",
                 "outbreak_start_date", "outbreak_end_date", "deaths",
                 "Vaccine_Delivery", "coi_per_patient", "tot_tested_amr", 
                 "tot_amr_cases", "resistant", "mdr", "fqns", "xdr", 
                 "prop_amr_cases", "prop_resistant_cases", "prop_mdr_cases", 
                 "prop_fqns_cases", "prop_xdr_cases")]
  
  df$vacc_week <- vacc_week
  
  week_delay_to_vacc_effect <- week_delay
  
  # df$case_trigger <- case_trigger
  df$vacc_cov <- vacc_cov
  # df$prop_u2 <- prop_u2
  # df$prop_2_to_4 <- prop_2_to_4
  # df$prop_5_to_16 <- prop_5_to_16
  # df$prop_over16 <- prop_over16
  # 
  df$ori_occurred <- FALSE
  
  nrow <- nrow(data) - 3 # outbreak missed if it ends 3 weeks prior 
  
  # vaccination happens after outbreak ended
  if (week_delay_to_vacc_effect > nrow) {
    
    df$s_ch_averted <- 0
    df$c_ch_averted <- 0
    df$death_averted <- 0
    
    df$runid <- runid
  }
  else {
    df$ori_occurred <- TRUE

    fa <- 1 - (vacc_cov*(1-dve)*(1-ive)+(1-vacc_cov)*(1-ive))
    
    # apply delay to vaccine effectiveness
    frac_averted <- rep(fa, nrow(data))
    
    if (week_delay_to_vacc_effect >= 1) { # if week_delay <= 0, pre-emptive vacc scenario
      frac_averted[1:week_delay_to_vacc_effect] <- 0
    }
    
    df$s_ch_averted <- data$suspected_cases * frac_averted
    
    df$c_ch_averted <- data$confirmed_cases * frac_averted
    
    df$death_averted <- data$deaths * frac_averted
    
    df$s_ch_averted_amr <- df$s_ch_averted * data$prop_amr_cases
    df$s_ch_averted_resistant <- df$s_ch_averted * data$prop_resistant_cases
    df$s_ch_averted_mdr <- df$s_ch_averted * data$prop_mdr_cases
    df$s_ch_averted_fqns <- df$s_ch_averted * data$prop_fqns_cases
    df$s_ch_averted_xdr <- df$s_ch_averted * data$prop_xdr_cases
    
    df$runid <- runid
  }
  return(df)
}

results_list <- list()
# ive <- seq(0, 50, by = 5)/100 
ive <- 0.17
for (study in outbreak_ids) {
  outbreak_data <- ts_weekly[ts_weekly$study_id == study, ]
  
  # Loop over vaccination weeks and coverage levels
  for (w in vacc_weeks) {
    for (vc in vacc_covs) {
      
      sample_results <- purrr::map_dfr(1:nrow(param_sets), function(i) {
        p <- param_sets[i, ]
        
        vacc_impact_outbreak_weekly(
          data = outbreak_data,
          age_dist = age_dist,
          vacc_week = w,
          vacc_cov = vc,
          week_delay = vacc_weeks[w] + p$immuno_delay/7 + (round(p$campaign_duration/7 / 2)), # vaccine takes effect halfway through campaign duration
          dve = p$dve,
          ive = ive, # fixed value based on Khanam study for now
          runid = i
        )
      })
      
      sample_results$vacc_week <- w
      sample_results$vacc_cov <- vc
      sample_results$study_id <- study
      
      list_name <- paste0("study", study, "_w", w, "_vc", vc * 100)
      results_list[[list_name]] <- sample_results
    }
  }
}

svim <- data.table::rbindlist(results_list, fill = TRUE)

population <- dat %>%
  select(study_id, population)
population <- as.data.frame(population)

svim <- left_join(svim, population, by = c("study_id"))

country <- dat %>% select(study_id, country_iso) %>% rename(country1 = country_iso)
svim <- left_join(svim, country, by = c("study_id"))

svim <- svim %>% select(-country_iso)
svim <- svim %>%
  rename(country_iso = country1)

year <- ts1 %>%
  group_by(study_id) %>%
  summarize(year = first(year(as.Date(start_date))))

svim <- left_join(svim, year, by = c("study_id"))

dates <- dat %>% select(study_id, outbreak_start_date, outbreak_end_date)

svim <- svim %>% select(-outbreak_end_date, -outbreak_start_date)
svim <- left_join(svim, dates, by = c("study_id"))

svim <- svim %>% select(-Vaccine_Delivery)
svim <- left_join(svim, costs2, by = c("country_iso" = "ISO3"))

svim_sum <- svim_sum_across_week(svim)
svim_sum1 <- left_join(svim_sum, gdp_long, by = c("country", "year"))
svim_sum1 <- left_join(svim_sum1, workforce_long, by = c("country", "year"))

svim_sum1 <- left_join(svim_sum1,
                  life_exp_data[, c("country","year", as.character(avg_age_inf))], 
                  by = c("country","year"))
svim_sum1 <- svim_sum1 %>% rename(life_exp = `13`)
svim_sum1$life_exp <- as.numeric(svim_sum1$life_exp)

# svim_sum <- add_cea_variables(svim_sum)
svim_sum1 <- add_cea_results(svim_sum1)

svim_fin <- fin_median(svim_sum1)

svim_fin <- fread("weekly_output_n200_fin.csv")


############# monthly outbreaks model
ts_monthly <- ts1 %>%
  filter(time_unit == "monthly") %>%
  group_by(study_id) %>%
  mutate(month = row_number(),
         date = date_mid)

vacc_covs <- seq(0.6, 0.95, by = 0.05)
vacc_month <- seq(1, 12, 1)

vacc_impact_outbreak_monthly <- function(data = NULL,
                                         age_dist = NULL,
                                         vacc_month = NULL,
                                         vacc_cov = NULL,
                                         week_delay = NULL,
                                         dve = NULL,
                                         ive = NULL, 
                                         runid = NULL)  {
  yr_ch <- as.character(data.table::year(data$date[1]))
  
  filtered_data <- age_dist[age_dist$ISO3 == data$country_iso[1] & 
                              as.character(age_dist$Year) == yr_ch, ]
  filtered_data <- na.omit(filtered_data)
  # total_population <- filtered_data$total_population[1] 
  # prop_u2 <- sum(filtered_data[, c("0", "1")], na.rm = TRUE) / total_population
  # prop_2_to_4 <- sum(filtered_data[, c("2", "3", "4")], na.rm = TRUE) / total_population
  # prop_5_to_16 <- sum(filtered_data[, c("5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16")],
  #                     na.rm = TRUE) / total_population
  # prop_over16 <- sum(filtered_data[, as.character(17:100)], na.rm = TRUE) / total_population  
  
  df <- data[, c("country_iso", "month", "study_id", "date",
                 "suspected_cases", "confirmed_cases",
                 "outbreak_start_date", "outbreak_end_date", "deaths",
                 "Vaccine_Delivery", "coi_per_patient", "tot_tested_amr", 
                 "tot_amr_cases", "resistant", "mdr", "fqns", "xdr", 
                 "prop_amr_cases", "prop_resistant_cases", "prop_mdr_cases", 
                 "prop_fqns_cases", "prop_xdr_cases")]
  
  df$vacc_month <- vacc_month
  # df$vacc_week <- vacc_month * 4
  
  # df$case_trigger <- case_trigger
  df$vacc_cov <- vacc_cov
 
  df$ori_occurred <- FALSE
  
  week_delay_to_vacc_effect <- week_delay / 4
  # df$case_trigger <- case_trigger
  df$vacc_cov <- vacc_cov

  df$ori_occurred <- FALSE
  
  nrow <- nrow(data) - 3
  
  # vaccination happens after outbreak ended
  if (week_delay_to_vacc_effect > nrow) {

    df$s_ch_averted <- 0
    df$c_ch_averted <- 0
    df$death_averted <- 0
    
    df$runid <- runid
  }
  else {
    df$ori_occurred <- TRUE
    fa <- 1 - (vacc_cov*(1-dve)*(1-ive)+(1-vacc_cov)*(1-ive))

     # apply delay to vaccine effectiveness
    frac_averted <- rep(fa, nrow(data))

    if (week_delay_to_vacc_effect >= 1) { # if week_delay <= 0, pre-emptive vacc scenario
      frac_averted[1:week_delay_to_vacc_effect] <- 0
    }
    

    df$s_ch_averted <- data$suspected_cases * frac_averted 
    
    df$c_ch_averted <- data$confirmed_cases * frac_averted 
    
    df$death_averted <- data$deaths * frac_averted 
    
    df$s_ch_averted_amr <- df$s_ch_averted * data$prop_amr_cases
    df$s_ch_averted_resistant <- df$s_ch_averted * data$prop_resistant_cases
    df$s_ch_averted_mdr <- df$s_ch_averted * data$prop_mdr_cases
    df$s_ch_averted_fqns <- df$s_ch_averted * data$prop_fqns_cases
    df$s_ch_averted_xdr <- df$s_ch_averted * data$prop_xdr_cases
    
    df$runid <- runid
  }
  return(df)
}

outbreak_ids <- unique(ts_monthly$study_id)

vacc_covs <- seq(0.6, 0.95, by = 0.05)
vacc_month <- seq(1, 12, 1)
results_list <- list()

vacc_covs <- c(0.80)
vacc_month <- 2

for (study in outbreak_ids) {
  outbreak_data <- ts_monthly[ts_monthly$study_id == study, ]
  
  for (w_idx in seq_along(vacc_month)) {
    w <- vacc_month[w_idx]
    
    for (vc in vacc_covs) {
      
      sample_results <- purrr::map_dfr(1:nrow(param_sets), function(i) {
        p <- param_sets[i, ]
        
        vacc_impact_outbreak_monthly(
          data = outbreak_data,
          age_dist = age_dist,
          vacc_month = w,
          vacc_cov = vc,
          week_delay = round((w * 4) + p$immuno_delay/7 + (p$campaign_duration/ 7 / 2)),
          dve = p$dve,
          ive = ive,
          runid = i
        )
      })
      
      sample_results$vacc_month <- w
      sample_results$vacc_cov <- vc
      sample_results$study_id <- study
      
      list_name <- paste0("study", study,"_w", w, "_vc", vc * 100)
      results_list[[list_name]] <- sample_results
    }
  }}

svim_month <- data.table::rbindlist(results_list, fill = TRUE)

svim_month <- left_join(svim_month, population, by = c("study_id"))
svim_month$vacc_week <- svim_month$vacc_month * 4

svim_month <- svim_month %>%
  mutate(country_iso = ifelse(study_id == "Lutterloh 2012", "MWI", country_iso))

svim_month <- svim_month %>%
  select(-Vaccine_Delivery, -coi_per_patient)
svim_month <- left_join(svim_month, costs2, by = c("country_iso" = "ISO3"))

year <- ts1 %>%
  group_by(study_id) %>%
  summarize(year = first(year(as.Date(start_date))))

svim_month <- left_join(svim_month, year, by = c("study_id"))

dates <- dat %>% select(study_id, outbreak_start_date, outbreak_end_date)

svim_month <- svim_month %>% select(-outbreak_end_date, -outbreak_start_date)
svim_month <- left_join(svim_month, dates, by = c("study_id"))

svim_month_sum <- svim_sum_across_monthly(svim_month)

# svim_month_sum1 <- add_cea_variables(svim_month_sum)
svim_month_sum1 <- add_cea_results(svim_month_sum1)
svim_month_fin <- fin_median(svim_month_sum1)

fwrite(svim_month_fin, "monthly_output_n200_fin_ive.csv")


#### model for sensitivity analyses for IVE

vacc_covs <- c(0.80)        
vacc_week <- 8      
ive_vals <- seq(0, 50, by = 5)/100

results_list <- list()

for (study in outbreak_ids) {
  outbreak_data <- ts_weekly[ts_weekly$study_id == study, ]
  
  for (vc in vacc_covs) {
    for (iv in ive_vals) {
      
      sample_results <- purrr::map_dfr(1:nrow(param_sets), function(i) {
        p <- param_sets[i, ]
        
        vacc_impact_outbreak_weekly(
          data = outbreak_data,
          age_dist = age_dist,
          vacc_week = vacc_week,  
          vacc_cov = vc,
          week_delay = vacc_week + p$immuno_delay/7 + (round(p$campaign_duration/7 / 2)), 
          dve = p$dve,
          ive = iv,
          runid = i
        )
      })
      
      sample_results$vacc_week <- vacc_week
      sample_results$vacc_cov <- vc
      sample_results$study_id <- study
      sample_results$ive <- iv
      
      list_name <- paste0("study", study, "_w", vacc_week, "_vc", vc * 100, "_ive", iv*100)
      results_list[[list_name]] <- sample_results
    }
  }
}

svim <- data.table::rbindlist(results_list, fill = TRUE)

fwrite(svim, "weekly_output_n200_ive.csv")

population <- dat %>%
  select(study_id, population)
population <- as.data.frame(population)
svim_sum <- svim_sum %>% select(-population.x, -population.y)
svim <- left_join(svim, population, by = c("study_id"))

country <- dat %>% select(study_id, country_iso) %>% rename(country1 = country_iso)
svim <- left_join(svim, country, by = c("study_id"))

svim <- svim %>% select(-country_iso)
svim <- svim %>%
  rename(country_iso = country1)

year <- ts1 %>%
  group_by(study_id) %>%
  summarize(year = first(year(as.Date(start_date))))

svim <- left_join(svim, year, by = c("study_id"))

dates <- dat %>% select(study_id, outbreak_start_date, outbreak_end_date)

svim <- svim %>% select(-outbreak_end_date, -outbreak_start_date)
svim <- left_join(svim, dates, by = c("study_id"))

svim <- svim %>% select(-Vaccine_Delivery, -coi_per_patient)
svim <- left_join(svim, costs2, by = c("country_iso" = "ISO3"))

svim_sum <- svim_sum_across_week(svim)
svim_sum1 <- left_join(svim_sum, gdp_long, by = c("country", "year"))
svim_sum1 <- left_join(svim_sum1, workforce_long, by = c("country", "year"))

svim_sum1 <- left_join(svim_sum1,
                       life_exp_data[, c("country","year", as.character(avg_age_inf))], 
                       by = c("country","year"))
svim_sum1 <- svim_sum1 %>% rename(life_exp = `13`)
svim_sum1$life_exp <- as.numeric(svim_sum1$life_exp)

# svim_sum <- add_cea_variables(svim_sum)
svim_sum1 <- add_cea_results(svim_sum1)

svim_fin <- fin_median(svim_sum1)

fwrite(svim_fin, "weekly_output_n200_ive_fin.csv")


ts_monthly <- ts1 %>%
  filter(time_unit == "monthly") %>%
  group_by(study_id) %>%
  mutate(month = row_number(),
         date = date_mid)
outbreak_ids <- unique(ts_monthly$study_id)

vacc_covs <- c(0.80)        
vacc_month <- 2                # fixed at month 2
ive_vals <- seq(0, 50, by = 5)/100

results_list <- list()

for (study in outbreak_ids) {
  outbreak_data <- ts_monthly[ts_monthly$study_id == study, ]
  
  for (vc in vacc_covs) {
    for (iv in ive_vals) {
      
      sample_results <- purrr::map_dfr(1:nrow(param_sets), function(i) {
        p <- param_sets[i, ]
        
        vacc_impact_outbreak_monthly(
          data = outbreak_data,
          age_dist = age_dist,
          vacc_month = vacc_month,  
          vacc_cov = vc,
          week_delay = round((vacc_month * 4) + p$immuno_delay/7 + (p$campaign_duration/7 / 2)), 
          dve = p$dve,
          ive = iv,   # loop over ive values
          runid = i
        )
      })
      
      sample_results$vacc_month <- vacc_month
      sample_results$vacc_cov <- vc
      sample_results$study_id <- study
      sample_results$ive <- iv
      
      list_name <- paste0("study", study,"_m", vacc_month, "_vc", vc * 100, "_ive", iv*100)
      results_list[[list_name]] <- sample_results
    }
  }
}

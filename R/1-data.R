library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)

age_dist <- as.data.frame(read_xls("WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES1.xls")) # [United Nations World Population Prospects 2024 - Population]
age_dist$`100+` <- as.numeric(age_dist$`100+`)
age_dist <- age_dist %>%
  mutate(across(as.character(0:99), as.numeric)) %>%
  rename(`100` = `100+`) %>%
  mutate(total_population = rowSums(across(as.character(0:100)), na.rm = TRUE))

life_exp_data <- as.data.frame(fread("wpp_life_expectancy_20241022.csv")) # [United Nations World Population Prospects 2024 - Mortality]
workforce_data <- read_xls("Workforce_Worldbank.xls") # [World Bank Data, Labor force participation rate, total (% of total population ages 15-64) (modeled ILO estimate)]
cpi_percent <- read_xls("imf-dm-export-20250216.xls")
gdp_data <- read_xls("API_NY.GDP.PCAP.xls")

parms <- as.data.frame(fread("parameters.csv"))
parms <- parms %>%
  filter(Disease == "Typhoid") %>%
  select(2:9)

costs <- as.data.frame(read_xlsx("country_costs.xlsx"))

pr_moderate <- parms[parms$Parameter == "Prop_Moderate", "Value"]
pr_severe <- parms[parms$Parameter == "Prop_Severe", "Value"]
pr_gi_bleed <- parms[parms$Parameter == "Prop_GI_Bleed", "Value"]
pr_intestinal_perf <- parms[parms$Parameter == "Prop_Intestinal_Perf", "Value"]
dur_moderate <- parms[parms$Parameter == "Duration_Illness_Moderate", "Value"]
dur_severe <- parms[parms$Parameter == "Duration_Illness_Severe", "Value"]

wt_moderate <- parms[parms$Parameter == "Disability_Weight_Moderate", "Value"]
wt_severe <- parms[parms$Parameter == "Disability_Weight_Severe", "Value"]
wt_gi_bleed <- parms[parms$Parameter == "Disability_Weight_GI_Bleed", "Value"]
wt_intestinal_perf <- parms[parms$Parameter == "Disability_Weight_Intestinal_Perf", "Value"]

vacc_cost_per_dose <- parms[parms$Parameter == "Vaccine_Cost", "Value"]
vacc_equipment_cost <- parms[parms$Parameter == "Vaccine_Equipment_Cost", "Value"]

# costs include indirect costs like productivity lost already

avg_age_inf <- parms[parms$Parameter == "Mean_Age_Infection", "Value"]

dose_regimen <- 1

pr_tot <- pr_moderate + pr_severe + pr_gi_bleed + pr_intestinal_perf

cpi_percent <- cpi_percent %>%
  filter(`Inflation rate, average consumer prices (Annual percent change)` == "World") %>%
  pivot_longer(
    cols = `2017`:`2025`,         
    names_to = "year",             
    values_to = "value"            
  ) %>%
  select(year, value)

cpi_percent$year <- as.numeric(cpi_percent$year)
cpi_percent$value <- as.numeric(cpi_percent$value)

# year starts at the t + 1 year of reported cost, so if 2015 cost then should use
# % change CPI starting 2016
xi <- (cpi_percent$value/100 + 1)
cumulative_inflation_factor <- xi[1] * xi[2] * xi[3] * xi[4] * xi[5] * xi[6] * xi[7] * xi[8] * xi[9]
cumulative_inflation_factor2 <- xi[7] * xi[8] * xi[9]
costs1 <- costs %>%
  mutate(Value_New = ifelse(!(ISO3 %in% c("DEU", "CHN", "USA", "MEX", "FJI")),
                            Value * cumulative_inflation_factor, 
                            Value * cumulative_inflation_factor2))

vacc_equipment_cost_adj <- vacc_equipment_cost * cumulative_inflation_factor

costs2 <- costs1 %>%
  distinct(ISO3, Estimate, Value_New) %>%
  pivot_wider(
    names_from = Estimate,
    values_from = Value_New
  )

costs2$coi_per_patient <- 
  (pr_moderate / pr_tot * costs2$Cost_Outpatient) + 
  (pr_severe / pr_tot * costs2$Cost_Inpatient) +
  (pr_gi_bleed / pr_tot * costs2$Cost_Inpatient) +
  (pr_intestinal_perf / pr_tot * costs2$Cost_Inpatient)

names(life_exp_data)[names(life_exp_data) == "ISO3 Alpha-code"] <- "country"
names(life_exp_data)[names(life_exp_data) == "Year"] <- "year"
names(workforce_data)[names(workforce_data) == "Country Code"] <- "country"

workforce_long <- 
  workforce_data |> 
  select(c("country", `2001`:`2023`)) |>
  pivot_longer(cols=`2001`:`2023`, names_to="year", values_to="pct_workforce")
workforce_long$year <- as.integer(workforce_long$year)

names(gdp_data)[names(gdp_data) == "Country Code"] <- "country"
gdp_long <- 
  gdp_data |> 
  select(c("country", `2001`:`2023`)) |>
  pivot_longer(cols=`2001`:`2023`, names_to="year", values_to="gdp")
gdp_long$year <- as.integer(gdp_long$year)


#################################

dat <- fread("/Users/monicaduong/Desktop/IVI/Typhoid/R/typhoid_outbreak_summary_v1.csv")

ts <- fread("/Users/monicaduong/Desktop/IVI/Typhoid/R/typhoid_outbreak_ts_v2.csv")
dat$population <- as.numeric(dat$population)
dat$tot_deaths <- as.numeric(dat$tot_deaths)

### exclude outbreaks based on eligibility criteria
dat <- dat %>% 
  filter(!study_id %in% c("Cherian 2015", "Wang 2022 (2)", "Imanishi 2014 (total)", 
                          "Polonsky 2014 (total)", "Limpitikul 2014", "Limpitikul 2014 (1)", "Limpitikul 2014 (2)",
                          "Poncin 2022", "Scobie 2014", "Muehlen 2007", "Bayram 2011", "Keddy 2011", "Michel 2005",
                          "Burnsed 2018", "Hu 2022", "Bano-Zaidi 2018", "Wang 2022 (1)", "Makungo 2020", "Al-Sanouri 2008", 
                          "Hechaichi 2023", "Holt 2011", "Roy 2016", "Nimonkar 2022", "Srinivasan 2022"))

dat <- dat %>% 
  group_by(study_id) |> 
  mutate(duration_weeks = as.integer((as.Date(end_date) - 
                                        as.Date(start_date))/7))
dat <- dat %>% 
  mutate(intervention_date = ifelse(intervention_date == "-", NA, intervention_date))
dat$intervention_date <- as.Date(dat$intervention_date)
dat$start_date <- as.Date(dat$start_date)
dat$intervention_timing <- as.numeric(dat$intervention_date - dat$start_date)/7
dat$attack_rate <- as.numeric(dat$attack_rate)
dat <- dat %>%
  rename(outbreak_start_date = start_date,
         outbreak_end_date = end_date) %>%
  select(study_id, country_iso, location, WHO_region, outbreak_start_date, outbreak_end_date, 
         time_unit, tot_deaths, population, attack_rate, intervention_timing)

amr <- fread("amr-data.csv") # 29 outbreaks with tot_tested_amr not NA
amr$tot_amr_cases <- as.integer(amr$tot_amr_cases)
amr$prop_amr_cases <- amr$tot_amr_cases / amr$tot_tested_amr
amr$prop_resistant_cases <- amr$resistant / amr$tot_tested_amr
amr$prop_mdr_cases <- amr$mdr / amr$tot_tested_amr
amr$prop_fqns_cases <- amr$fqns / amr$tot_tested_amr
amr$prop_xdr_cases <- amr$xdr / amr$tot_tested_amr

amr1 <- amr %>%
  select(study_id, amr_status, tot_tested_amr, tot_amr_cases, resistant, mdr, fqns, xdr,
         prop_amr_cases, prop_resistant_cases, prop_mdr_cases, prop_fqns_cases, prop_xdr_cases)

dat <- left_join(dat, amr1, by = c("study_id"))
ts <- ts %>% filter(study_id %in% dat$study_id)

ts <- left_join(ts, dat, by = c("study_id"))

ts <- ts %>%
  group_by(study_id) %>%
  mutate(cum_cases = sum(total_cases),
         prop_cases = total_cases / cum_cases,
         deaths = tot_deaths * prop_cases) # use distribution of cases time series for deaths

ts <- ts %>%
  select(-V10, -V11, -V12, -V13, -V14) # clean

ts <- left_join(ts, costs2, by = c("country_iso" = "ISO3"))


# ts <- left_join(ts, duration, by = c("study_id"))

# rm(a, b, c, d, amr, costs, costs1, costs2, cpi_percent, gdp_data, gdp_long, life_exp_data)
# rm(workforce_data, workforce_long)
# gc()

ts <- ts%>%
  rename(suspected_cases2 = suspected_cases,
         suspected_cases = total_cases
  ) %>%
  ungroup()

################## rescale daily outbreaks to weekly and fortnightly to monthly

ts <- ts %>%
  mutate(end_date = as.Date(end_date))

# ts <- left_join(ts, dat, by = c("study_id"))

# 1. Process daily outbreaks: rescale to weekly and assign midpoint

ts_daily <- ts %>% 
  filter(time_unit == "daily") %>%
  mutate(
    week_start = floor_date(end_date, unit = "week", week_start = 1)
  ) %>%
  group_by(study_id, week_start) %>%
  summarise(
    suspected_cases = sum(suspected_cases, na.rm = TRUE),
    time_unit = "weekly",
    date_mid = first(week_start + days(3)),
    WHO_region = first(WHO_region),
    country_iso = first(country_iso),
    start_date = first(start_date),
    end_date = first(end_date),
    .groups = "drop"
  )

ts_biweekly <- ts %>%
  filter(time_unit == "fortnightly") %>%
  mutate(
    month_start = floor_date(end_date, unit = "month")
  ) %>%
  group_by(study_id, month_start) %>%
  summarise(
    suspected_cases = sum(suspected_cases, na.rm = TRUE),
    time_unit = "monthly",
    date_mid = first(month_start + days(15)),
    WHO_region = first(WHO_region),
    country_iso = first(country_iso),
    start_date = first(start_date),
    end_date = first(end_date),
    .groups = "drop"
  )

# Handle weekly and monthly as-is
ts_other <- ts %>%
  filter(time_unit %in% c("weekly", "monthly")) %>%
  mutate(date_mid = case_when(
    time_unit == "weekly"  ~ end_date - days(3),
    time_unit == "monthly" ~ end_date - days(15),
    TRUE                   ~ end_date
  ))

ts1 <- bind_rows(ts_other, ts_biweekly, ts_daily)

dat |>
  group_by(study_id) |>
  summarize(duration_weeks = as.integer((as.Date(outbreak_end_date) -
                                           as.Date(outbreak_start_date))/7)) -> duration

ts1 <- ts1 %>%
  group_by(study_id) %>%
  mutate(cum_cases = sum(suspected_cases))

ts1 <- left_join(ts1, duration, by = c("study_id"))



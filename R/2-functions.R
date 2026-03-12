compute_yld <- function(cases, parms = NULL, option = 1, amr_multiplier = 1) {
  if (is.null(parms)) stop("Parameters to compute YLD must be provided")
  
  pm <- parms[parms$Parameter == "Prop_Moderate", "Value"]
  ps <- parms[parms$Parameter == "Prop_Severe", "Value"]
  pgi <- parms[parms$Parameter == "Prop_GI_Bleed", "Value"]
  pip <- parms[parms$Parameter == "Prop_Intestinal_Perf", "Value"]
  
  dur_moderate <- parms[parms$Parameter == "Duration_Illness_Moderate", "Value"]
  dur_severe <- parms[parms$Parameter == "Duration_Illness_Severe", "Value"]
  
  wm <- parms[parms$Parameter == "Disability_Weight_Moderate", "Value"]
  ws <- parms[parms$Parameter == "Disability_Weight_Severe", "Value"]
  wgi <- parms[parms$Parameter == "Disability_Weight_GI_Bleed", "Value"]
  wip <- parms[parms$Parameter == "Disability_Weight_Intestinal_Perf", "Value"]
  
  dis_wt_mod <- pm/(pm + ps + pgi + pip) * wm
  dis_wt_severe <- ps/(pm + ps + pgi + pip) * ws +
    pgi/(pm + ps + pgi + pip) * wgi +
    pip/(pm + ps + pgi + pip) * wip
  
  yld_per_capita <- (dur_moderate/365 * dis_wt_mod) +
    (dur_severe/365 * dis_wt_severe)
  
  yld_tot <- cases * yld_per_capita * amr_multiplier # multiplier is 2 for AMR cases
  return(yld_tot)
}

get_life_expectancy <- function(life_exp_data, age, data) {
  
  country_iso <- unique(data$country_iso)
  year <- unique(data.table::year(data$date))
  
  life_exp <- life_exp_data |> 
    dplyr::filter(`ISO3 Alpha-code` %in% country_iso, 
                  Year %in% as.character(year)) |> 
    dplyr::select(as.character(age)) |> 
    unlist() |> 
    as.numeric()
  
  return(life_exp)
}


compute_yll <- function(deaths,
                        life_exps=NULL,
                        parms=NULL) {
  
  dr <- parms[parms$Parameter == "Rate_Discount", "Value"]
  yll <- deaths * (1/dr) * (1 - exp(-dr * life_exps))
  
  return(yll)
}

svim_sum_across_week <- function(svim, case_trigger = FALSE) {

  by_cols <- if (case_trigger) {
    c("runid", "id_outbreak", "case_trigger", "vacc_cov")
  } else {
    c("study_id", "vacc_week", "vacc_cov", "runid", "ive")
  }
  
  svim %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      country = first(country_iso),  # Take the first value
      start_date = as.Date(first(outbreak_start_date)),
      end_date = as.Date(first(outbreak_end_date)),
      
      year = first(year),
      population = first(population),
      # week_delay_to_vacc_effect = first(week_delay_to_vacc_effect),
      c_ch_tot = sum(confirmed_cases, na.rm = TRUE),
      s_ch_averted_tot = sum(s_ch_averted, na.rm = TRUE),
      s_ch_tot = sum(suspected_cases, na.rm = TRUE),
      
      s_ch_averted_tot = sum(s_ch_averted, na.rm = TRUE),
      s_ch_averted_amr_tot = sum(s_ch_averted_amr, na.rm = TRUE),
      s_ch_averted_resistant_tot = sum(s_ch_averted_resistant, na.rm = TRUE),
      s_ch_averted_mdr_tot = sum(s_ch_averted_mdr, na.rm = TRUE),
      s_ch_averted_fqns_tot = sum(s_ch_averted_fqns, na.rm = TRUE),
      s_ch_averted_xdr_tot = sum(s_ch_averted_xdr, na.rm = TRUE),
      death_tot = sum(deaths, na.rm = TRUE),
      death_averted_tot = sum(death_averted, na.rm = TRUE),
      pct_reduc_case = 100 * sum(s_ch_averted, na.rm = TRUE) / sum(suspected_cases, na.rm = TRUE),
      pct_reduc_case_amr = 100 * sum(s_ch_averted_amr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_amr_cases)),
      pct_reduc_case_mdr = 100 * sum(s_ch_averted_mdr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_mdr_cases)),
      pct_reduc_case_fqns = 100 * sum(s_ch_averted_fqns, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_fqns_cases)),
      pct_reduc_case_xdr = 100 * sum(s_ch_averted_xdr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_xdr_cases)),
      ori_occurred = first(ori_occurred),
      
      vaccine_delivery_cost = first(Vaccine_Delivery),
      coi_per_patient = first(coi_per_patient),
      .groups = "drop"
    )
}


svim_sum_across_monthly <- function(svim, case_trigger = FALSE) {

  by_cols <- if (case_trigger) {
    c("runid", "id_outbreak", "case_trigger", "vacc_cov")
  } else {
    c("study_id", "vacc_week", "vacc_cov", "runid", "ive")
  }
  
  svim %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      country = first(country_iso),
      year = first(year),
      start_date = as.Date(first(outbreak_start_date)),
      end_date = as.Date(first(outbreak_end_date)),
      population = first(population),
      # week_delay_to_vacc_effect = first(week_delay_to_vacc_effect),
      c_ch_tot = sum(confirmed_cases, na.rm = TRUE),
      s_ch_averted_tot = sum(s_ch_averted, na.rm = TRUE),
      s_ch_tot = sum(suspected_cases, na.rm = TRUE),
      s_ch_tot_amr = sum(tot_amr_cases, na.rm = TRUE),
      s_ch_tot_resistant = sum(resistant, na.rm = TRUE),
      s_ch_tot_mdr = sum(mdr, na.rm=TRUE),
      s_ch_tot_fqns = sum(fqns, na.rm = TRUE),
      s_ch_tot_xdr = sum(xdr, na.rm = TRUE),
      s_ch_averted_tot = sum(s_ch_averted, na.rm = TRUE),
      s_ch_averted_amr_tot = sum(s_ch_averted_amr, na.rm = TRUE),
      s_ch_averted_resistant_tot = sum(s_ch_averted_resistant, na.rm = TRUE),
      s_ch_averted_mdr_tot = sum(s_ch_averted_mdr, na.rm = TRUE),
      s_ch_averted_fqns_tot = sum(s_ch_averted_fqns, na.rm = TRUE),
      s_ch_averted_xdr_tot = sum(s_ch_averted_xdr, na.rm = TRUE),
      death_tot = sum(deaths, na.rm = TRUE),
      death_averted_tot = sum(death_averted, na.rm = TRUE),
      pct_reduc_case = 100 * sum(s_ch_averted, na.rm = TRUE) / sum(suspected_cases, na.rm = TRUE),
      ori_occurred = first(ori_occurred),
      pct_reduc_case = 100 * sum(s_ch_averted, na.rm = TRUE) / sum(suspected_cases, na.rm = TRUE),
      pct_reduc_case_amr = 100 * sum(s_ch_averted_amr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_amr_cases)),
      pct_reduc_case_mdr = 100 * sum(s_ch_averted_mdr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_mdr_cases)),
      pct_reduc_case_fqns = 100 * sum(s_ch_averted_fqns, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_fqns_cases)),
      pct_reduc_case_xdr = 100 * sum(s_ch_averted_xdr, na.rm = TRUE) / (sum(suspected_cases, na.rm = TRUE) * first(prop_xdr_cases)),
      ori_occurred = first(ori_occurred),
      vaccine_delivery_cost = first(Vaccine_Delivery),
      coi_per_patient = first(coi_per_patient),
      .groups = "drop"
    )
}

add_cea_variables <- function(svim) {
  # svim$year <- data.table::year(as.Date(svim$start_date))
  svim <- left_join(svim,
                    life_exp_data[, c("country","year", as.character(avg_age_inf))], 
                    by = c("country","year"))
  names(svim)[names(svim) == as.character(avg_age_inf)] <- "life_exp"
  svim$life_exp <- as.numeric(svim$life_exp)
  svim <- left_join(svim,
                    gdp_long[, c("country", "year", "gdp")],
                    by = c("country", "year"))
  
  svim <- left_join(svim, 
                    workforce_long[, c("country", "year", "pct_workforce")],
                    by = c("country", "year"))
  return(svim)
}

add_cea_results <- function(svim) {
  svim |> mutate(
    
    # Sum of all predicted AMR cases averted
    s_ch_averted_amr_total = s_ch_averted_amr_tot,
    
    # Compute burden for AMR and non-AMR averted cases
    yld_averted_amr = compute_yld(s_ch_averted_amr_total, parms = parms, amr_multiplier = 2),
    yld_averted_non_amr = compute_yld(s_ch_averted_tot - s_ch_averted_amr_total, parms = parms),
    yld_averted = yld_averted_amr + yld_averted_non_amr,
    
    yll_averted = compute_yll(death_averted_tot, life_exp, parms),
    daly_averted = yld_averted + yll_averted,
    
    coi_averted_amr = s_ch_averted_amr_total * coi_per_patient * 2,
    coi_averted_non_amr = (s_ch_averted_tot - s_ch_averted_amr_total) * coi_per_patient,
    coi_averted = coi_averted_amr + coi_averted_non_amr,
    
    cod_averted = death_averted_tot * gdp * life_exp,
    
    vacc_dose = population * vacc_cov * dose_regimen,
    vacc_cost = (vacc_cost_per_dose + vaccine_delivery_cost + vacc_equipment_cost_adj) * vacc_dose,
    net_cost_averted = vacc_cost - coi_averted - cod_averted,

    cost_per_daly_averted = ifelse(daly_averted > 0, net_cost_averted / daly_averted, 0),
    cost_per_case_averted = ifelse(s_ch_averted_tot > 0, net_cost_averted / s_ch_averted_tot, 0),
    cost_per_death_averted = ifelse(death_averted_tot > 0, net_cost_averted / death_averted_tot, 0),

    case_averted_per_1000_OCV = ifelse(vacc_dose > 0, 1000 * s_ch_averted_tot / vacc_dose, 0),
    amr_case_averted_per_1000_OCV = ifelse(vacc_dose > 0, 1000 * s_ch_averted_amr_total / vacc_dose, 0),
    death_averted_per_1000_OCV = ifelse(vacc_dose > 0, 1000 * death_averted_tot / vacc_dose, 0),
    mdr_case_averted_per_1000_OCV =
      ifelse(vacc_dose > 0, 1000 * s_ch_averted_mdr_tot / vacc_dose, 0),
    fqns_case_averted_per_1000_OCV =
      ifelse(vacc_dose > 0, 1000 * s_ch_averted_fqns_tot / vacc_dose, 0),
    xdr_case_averted_per_1000_OCV =
      ifelse(vacc_dose > 0, 1000 * s_ch_averted_xdr_tot / vacc_dose, 0),

    cost_eff_threshold = 3 * gdp,
    ratio_cost_per_daly_averted_to_gdp = ifelse(gdp > 0, cost_per_daly_averted / gdp, 0),
    cost_effective = ifelse(cost_per_daly_averted > cost_eff_threshold, FALSE, TRUE)
  ) -> svim_cea
  
  return(svim_cea)
}


fin_median <- function(svim) {
  svim |>
    group_by(study_id,vacc_week, vacc_cov, ive) |>
    summarise(country = first(country), 
              start_date = first(start_date), 
              population = first(population),
              c_ch_tot = first(c_ch_tot), 
              s_ch_tot = first(s_ch_tot),
              death_tot = first(death_tot),
              ori_occurred = first(ori_occurred),
              
              s_ch_averted_median = median(s_ch_averted_tot, na.rm = TRUE),
              s_ch_averted_low = quantile(s_ch_averted_tot, 0.025, na.rm = TRUE),
              s_ch_averted_high = quantile(s_ch_averted_tot, 0.975, na.rm = TRUE),
              s_ch_averted_amr_median = median(s_ch_averted_amr_tot, na.rm = TRUE),
              s_ch_averted_amr_low = quantile(s_ch_averted_amr_tot, 0.025, na.rm = TRUE),
              s_ch_averted_amr_high = quantile(s_ch_averted_amr_tot, 0.975, na.rm = TRUE),
              
              s_ch_averted_resistant_median = median(s_ch_averted_resistant_tot, na.rm = TRUE),
              
              # s_ch_averted_resistant_iqr_low = quantile(s_ch_averted_resistant_tot, 0.25, na.rm = TRUE),
              # s_ch_averted_resistant_iqr_high = quantile(s_ch_averted_resistant_tot, 0.75, na.rm = TRUE),
              # 
              s_ch_averted_mdr_median = median(s_ch_averted_mdr_tot, na.rm = TRUE),
              # 
              # s_ch_averted_mdr_iqr_low = quantile(s_ch_averted_mdr_tot, 0.25, na.rm = TRUE),
              # s_ch_averted_mdr_iqr_high = quantile(s_ch_averted_mdr_tot, 0.75, na.rm = TRUE),
              # 
              s_ch_averted_fqns_median = median(s_ch_averted_fqns_tot, na.rm = TRUE),
              # 
              # s_ch_averted_fqns_iqr_low = quantile(s_ch_averted_fqns_tot, 0.25, na.rm = TRUE),
              # s_ch_averted_fqns_iqr_high = quantile(s_ch_averted_fqns_tot, 0.75, na.rm = TRUE),
              # 
              s_ch_averted_xdr_median = median(s_ch_averted_xdr_tot, na.rm = TRUE),
              # 
              # s_ch_averted_xdr_iqr_low = quantile(s_ch_averted_xdr_tot, 0.25, na.rm = TRUE),
              # s_ch_averted_xdr_iqr_high = quantile(s_ch_averted_xdr_tot, 0.75, na.rm = TRUE),
              # 
              death_averted_tot_median = median(death_averted_tot, na.rm = TRUE),
              death_averted_tot_low = quantile(death_averted_tot, 0.025, na.rm = TRUE),
              death_averted_tot_high = quantile(death_averted_tot, 0.975, na.rm = TRUE),
              
              pct_reduc_case_median = median(pct_reduc_case, na.rm = TRUE), 
              pct_reduc_case_low = quantile(pct_reduc_case, 0.025, na.rm = TRUE),
              pct_reduc_case_high = quantile(pct_reduc_case, 0.975, na.rm = TRUE),
              
              pct_reduc_case_median_amr = median(pct_reduc_case_amr, na.rm = TRUE), 
              pct_reduc_case_median_mdr = median(pct_reduc_case_mdr, na.rm = TRUE), 
              pct_reduc_case_median_fqns = median(pct_reduc_case_fqns, na.rm = TRUE), 
              pct_reduc_case_median_xdr = median(pct_reduc_case_xdr, na.rm = TRUE), 
              cost_per_daly_averted_median = median(cost_per_daly_averted, na.rm = TRUE), 
              cost_per_daly_averted_low = quantile(cost_per_daly_averted, 0.025, na.rm = TRUE),
              cost_per_daly_averted_high = quantile(cost_per_daly_averted, 0.975, na.rm = TRUE),
              
              cost_per_case_averted_median = median(cost_per_case_averted, na.rm = TRUE), 
              cost_per_case_averted_low = quantile(cost_per_case_averted, 0.025, na.rm = TRUE),
              cost_per_case_averted_high = quantile(cost_per_case_averted, 0.975, na.rm = TRUE),
              
              cost_per_death_averted_median = median(cost_per_death_averted, na.rm = TRUE), 
              cost_per_death_averted_low = quantile(cost_per_death_averted, 0.025, na.rm = TRUE),
              cost_per_death_averted_high = quantile(cost_per_death_averted, 0.975, na.rm = TRUE),
              
              case_averted_per_1000_median = median(case_averted_per_1000_OCV, na.rm = TRUE), 
              case_averted_per_1000_low = quantile(case_averted_per_1000_OCV, 0.025, na.rm = TRUE),
              case_averted_per_1000_high = quantile(case_averted_per_1000_OCV, 0.975, na.rm = TRUE),
              
              amr_case_averted_per_1000_median = median(amr_case_averted_per_1000_OCV, na.rm = TRUE), 
              fqns_case_averted_per_1000_median = median(fqns_case_averted_per_1000_OCV, na.rm = TRUE), 
              mdr_case_averted_per_1000_median = median(mdr_case_averted_per_1000_OCV, na.rm = TRUE), 
              xdr_case_averted_per_1000_median = median(xdr_case_averted_per_1000_OCV, na.rm = TRUE), 
              
              death_averted_per_1000_median = median(death_averted_per_1000_OCV, na.rm = TRUE), 
              # death_averted_per_1000_iqr_low = quantile(death_averted_per_1000_OCV, 0.25, na.rm = TRUE),
              # death_averted_per_1000_iqr_high = quantile(death_averted_per_1000_OCV, 0.75, na.rm = TRUE),
              # 
              ratio_cost_per_daly_averted_to_gdp_median = median(ratio_cost_per_daly_averted_to_gdp, na.rm = TRUE)
              # ratio_cost_per_daly_averted_to_gdp_iqr_low = quantile(ratio_cost_per_daly_averted_to_gdp, 0.25, na.rm = TRUE),
              # ratio_cost_per_daly_averted_to_gdp_iqr_high = quantile(ratio_cost_per_daly_averted_to_gdp, 0.75, na.rm = TRUE)) -> svim_cea
              
    )-> svim_cea
  return(svim_cea)
}

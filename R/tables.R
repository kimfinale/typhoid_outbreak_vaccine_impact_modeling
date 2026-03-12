library(dplyr)
library(tidyr)

# Supplementary Table 2

svim_fin <- fread("weekly_output_n200_fin.csv")

df <- svim_fin %>%
  filter(vacc_week %in% c(1, 4, 8, 16), vacc_cov %in% c(0.65, 0.8, 0.9)) %>%
  select(study_id, vacc_week, vacc_cov,
         s_ch_averted_median, s_ch_averted_low, s_ch_averted_high,
         pct_reduc_case_median, pct_reduc_case_low, pct_reduc_case_high,
         cost_per_daly_averted_median, cost_per_daly_averted_low, cost_per_daly_averted_high,
         cost_per_case_averted_median, cost_per_case_averted_low, cost_per_case_averted_high,
         cost_per_death_averted_median, cost_per_death_averted_low, cost_per_death_averted_high,
         case_averted_per_1000_median, case_averted_per_1000_low, case_averted_per_1000_high)

df_summary <- df %>%
  # filter(ori_occurred == TRUE) %>%
  group_by(vacc_week, vacc_cov) %>%
  summarise(cases_averted_median = median(s_ch_averted_median),
            cases_averted_lb = quantile(s_ch_averted_median, probs = 0.025),
            cases_averted_ub = quantile(s_ch_averted_median, probs = 0.975), 
            pct_cases_averted_median = median(pct_reduc_case_median),
            pct_cases_averted_lb = quantile(pct_reduc_case_median, probs = 0.025),
            pct_cases_averted_ub = quantile(pct_reduc_case_median, probs = 0.975),
            cost_per_daly_averted_median2 = median(cost_per_daly_averted_median, na.rm = TRUE),
            cost_per_daly_averted_lb = quantile(cost_per_daly_averted_median, probs = 0.025, na.rm = TRUE),
            cost_per_daly_averted_ub = quantile(cost_per_daly_averted_median, probs = 0.975, na.rm = TRUE),
            case_averted_per_1000_median2 = median(case_averted_per_1000_median, na.rm = TRUE),
            case_averted_per_1000_median_lb = quantile(case_averted_per_1000_median, probs = 0.025, na.rm = TRUE),
            case_averted_per_1000_median_ub = quantile(case_averted_per_1000_median, probs = 0.975, na.rm = TRUE),
  )

df_combined <- df_summary %>%
  mutate(
    s_ch_averted = paste0(round(s_ch_averted_median, 2), 
                          " (", round(s_ch_averted_low, 2), ", ", round(s_ch_averted_high, 2), ")"),
    pct_reduc_case = paste0(round(pct_reduc_case_median, 2), 
                            " (", round(pct_reduc_case_low, 2), ", ", round(pct_reduc_case_high, 2), ")"),
    cost_per_daly_averted = paste0(round(cost_per_daly_averted_median, 2), 
                                   " (", round(cost_per_daly_averted_low, 2), ", ", round(cost_per_daly_averted_high, 2), ")"),
    cost_per_case_averted = paste0(round(cost_per_case_averted_median, 2), 
                                   " (", round(cost_per_case_averted_low, 2), ", ", round(cost_per_case_averted_high, 2), ")"),
    cost_per_death_averted = paste0(round(cost_per_death_averted_median, 2), 
                                    " (", round(cost_per_death_averted_low, 2), ", ", round(cost_per_death_averted_high, 2), ")"),
    case_averted_per_1000 = paste0(round(case_averted_per_1000_median, 2), 
                                   " (", round(case_averted_per_1000_low, 2), ", ", round(case_averted_per_1000_high, 2), ")")
  ) %>%
  select(vacc_week, vacc_cov,
         s_ch_averted, pct_reduc_case, 
         cost_per_daly_averted, cost_per_case_averted, 
         cost_per_death_averted, case_averted_per_1000)


df_combined <- df %>%
  mutate(
    s_ch_averted = paste0(round(s_ch_averted_median, 2), 
                          " (", round(s_ch_averted_low, 2), ", ", round(s_ch_averted_high, 2), ")"),
    pct_reduc_case = paste0(round(pct_reduc_case_median, 2), 
                            " (", round(pct_reduc_case_low, 2), ", ", round(pct_reduc_case_high, 2), ")"),
    cost_per_daly_averted = paste0(round(cost_per_daly_averted_median, 2), 
                                   " (", round(cost_per_daly_averted_low, 2), ", ", round(cost_per_daly_averted_high, 2), ")"),
    cost_per_case_averted = paste0(round(cost_per_case_averted_median, 2), 
                                   " (", round(cost_per_case_averted_low, 2), ", ", round(cost_per_case_averted_high, 2), ")"),
    cost_per_death_averted = paste0(round(cost_per_death_averted_median, 2), 
                                    " (", round(cost_per_death_averted_low, 2), ", ", round(cost_per_death_averted_high, 2), ")"),
    case_averted_per_1000 = paste0(round(case_averted_per_1000_median, 2), 
                                   " (", round(case_averted_per_1000_low, 2), ", ", round(case_averted_per_1000_high, 2), ")")
  ) %>%
  select(study_id, vacc_week, vacc_cov,
         s_ch_averted, pct_reduc_case, 
         cost_per_daly_averted, cost_per_case_averted, 
         cost_per_death_averted, case_averted_per_1000)

# Supplementary Table 3

df <- fread("monthly_output_n200_fin.csv")


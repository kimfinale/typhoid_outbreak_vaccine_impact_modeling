# compute the impact of vaccine in a static way
# account for age (< 5 yo vs. 5+ yo)
#
run_vacc_impact <- function(outbreak_data = NULL, # time series
                            outbreak_data2 = NULL, # summary data
                            parameters = NULL,
                            age_dist = NULL,
                            vacc_cov = seq(0.1, 0.9, by=0.1),
                            runid = 1) {

  outbreak_ids <- unique(outbreak_data$id)
  # lists to store simulation results
  lst <- vector("list", length(outbreak_ids))
  lst2 <- vector("list", length(vacc_coverage))

  p <- parameters[runid, c("dur_campaign", "delay_vacc_eff",
                       "vacc_eff_direct_U5", "vacc_eff_direct_5plus")]
  # indirect effectiveness was pre-calculated and
  # here we only randomly select a column indicated a sample
  ive_col <- parameters$vacc_eff_indirect_id[runid] + 1
  col_name <- paste0("X", ive_col)
  dve <- c(p$vacc_eff_direct_U5, p$vacc_eff_direct_5plus)
  # individuals vaccinated at the mid-point, i.e., dur_campaign/2
  week_delay <- round((p$dur_campaign/2 + p$delay_vacc_eff)/7)
  # cat(sprintf("id = %d, dve[1] = %.2f, dve[2] = %.2f, week_delay = %d\n",
  #             runid, dve[1], dve[2], week_delay))

  for (j in 1:length(vacc_coverage)) {
    for (i in 1:length(outbreak_ids)) {
      d <- outbreak_data[outbreak_data$id == outbreak_ids[i],]
      # to account for that children under 5 have lower efficacy
      yr_ch <- strsplit(as.character(d$date[1]),"-")[[1]][1]
      # prop under 5 according to year and country
      # this is a data.table way and works fine.
      prop_U5 <- age_dist[ISO3 == d$country[1] & as.character(year) == yr_ch]$`0-4`
      # categorize by confirmed cases
      d$confirmed <- outbreak_data2[outbreak_data2$ID_outbreak == outbreak_ids[i], ]$total_confirmed_cases
      d$prop_U5 <- prop_U5
      # reported cases per each week
      d$case_wk_tot <- d$sCh
      d$case_wk_U5 <- d$sCh * prop_U5
      d$case_wk_5up <- d$sCh * (1-prop_U5)
      # total cases per outbreak  - stays constant across week
      d$case_total_U5 <- d$case_total * prop_U5
      d$case_total_5up <- d$case_total * (1-prop_U5)
      # cases to occur during the remainder of the outbreak
      # i.e., decreases over the week
      d$case_rem_U5 <- d$case_rem * prop_U5
      d$case_rem_5up <- d$case_rem * (1-prop_U5)

      # cat(d$id[1], ", propU5 =", prop_U5, "\n")
      vc <- vacc_coverage[j]
      d$vacc_cov <- vc
      # this only accounts for vaccine-induced immunity but not existing in the population
      # this could be another reason for underestimation for the vaccine impact
      # variables that may different across simulations
      vacc_cov_eff <- vc * (prop_U5 * dve[1] + (1-prop_U5) * dve[2])
      d$vacc_cov_eff <- vacc_cov_eff
      ive <- as.numeric(ive_pred[X1 == round(vacc_cov_eff, digits=2), ..col_name])
      d$ive <- ive
      if (i == 1) {
        cat(sprintf("vacc_cov_eff = %.2f, ive = %.2f\n", vacc_cov_eff, ive))
      }
      # # weekly case
      frac_averted_U5 <- 1 - (vc*(1-dve[1])*(1-ive)+(1-vc)*(1-ive))
      frac_averted_5up <- 1 - (vc*(1-dve[2])*(1-ive)+(1-vc)*(1-ive))

      d$case_wk_averted_U5 <- d$case_wk_U5 * frac_averted_U5
      d$case_wk_averted_5up <- d$case_wk_5up * frac_averted_5up
      d$case_wk_averted_tot <- d$case_wk_averted_U5 + d$case_wk_averted_5up

      d$case_rem_averted_U5 <- d$case_rem_U5 * frac_averted_U5
      d$case_rem_averted_5up <- d$case_rem_5up * frac_averted_5up
      d$case_rem_averted_tot <- d$case_rem_averted_U5 + d$case_rem_averted_5up

      d$death_averted_U5 <- d$death_total * prop_U5 * frac_averted_U5
      d$death_averted_5up <- d$death_total * (1-prop_U5) * frac_averted_5up
      d$death_averted <- d$death_averted_U5 + d$death_averted_5up
      # set vaccine impact variables at zero before the vaccine takes effect
      d[i:week_delay, c("case_wk_averted_U5", "case_wk_averted_5up",
                        "case_wk_averted_tot",  "case_rem_averted_U5", "case_rem_averted_5up",
                        "case_rem_averted_tot")] <- 0
      # deaths are given as NA for some cases
      # leave as NA if it is NA, making to zero will create errors

      for (xx in 1:week_delay) {
        if(!is.na(d[xx, "death_averted_U5"])){
          d[xx, "death_averted_U5"] <- 0
        }
        if(!is.na(d[xx, "death_averted_5up"])){
          d[xx, "death_averted_5up"] <- 0
        }
        if(!is.na(d[xx, "death_averted"])){
          d[xx, "death_averted"] <- 0
        }
      }
      d$pct_case_averted_U5 = 100 * d$case_rem_averted_U5 / d$case_total_U5
      d$pct_case_averted_5up = 100 * d$case_rem_averted_5up / d$case_total_5up
      d$pct_case_averted_tot = 100 * d$case_rem_averted_tot / d$case_total

      d$pct_gain_case_averted_U5 = c(-diff(d$pct_case_averted_U5),NA)
      d$pct_gain_case_averted_5up = c(-diff(d$pct_case_averted_5up),NA)
      d$pct_gain_case_averted_tot = c(-diff(d$pct_case_averted_tot),NA)

      lst[[i]] <- d
    }
    lst2[[j]] <- rbindlist(lst)
  }
  res <- rbindlist(lst2)
  res$runid <- runid

  return (res)
}

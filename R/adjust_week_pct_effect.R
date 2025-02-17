
# The following columns are NA's at the las row and there don't have to be set to zero
# c("pct_gain_case_averted_U5",
#   "pct_gain_case_averted_5up",
#   "pct_gain_case_averted_tot")


# the function adds (sort of) emptry rows for the outbreaks that ended before
# the week_final parameter such that summary statistics correctly account for
# weeks after the outbreak

adjust_week_pct_effect <- function(dat, week_final) {
  # these columns will be zero when new rows are added
  columns_zero <- c("sCh", "case_wk_tot", "case_wk_U5", "case_wk_5up", "case_rem",
                    "case_rem_U5", "case_rem_5up", "case_wk_averted_U5" ,
                    "case_wk_averted_5up", "case_wk_averted_tot", "case_rem_averted_U5",
                    "case_rem_averted_5up",  "case_rem_averted_tot",
                    # "death_averted_U5", "death_averted_5up", "death_averted",
                    "pct_case_averted_U5", "pct_case_averted_5up",
                    "pct_case_averted_tot")

  idx <- unique(dat$data_id)
  vc <- unique(dat$vacc_cov)
  lst2 <- vector("list", length(vc))
  for (j in 1:length(vc)) {
    cat("Filling the rows for the vaccine coverage =", vc[j], "\n")
    lst <- vector("list", length(idx))
    for (i in 1:length(idx)) {
      d <- dat[data_id == idx[i] & vacc_cov == vc[j]]
      maxwk <- max(d$week)
      if (maxwk < week_final) {
        nr <- week_final - maxwk
        for (k in 1:nr) {
          last_row <- d[nrow(d),]
          if (k == 1){
            # last_row[,c(6,12,16:37)] <- 0
            last_row[, columns_zero] <- 0
            if(!is.na(last_row[, "death_averted_U5"])){
                last_row[, "death_averted_U5"] <- 0
            }
            if(!is.na(last_row[, "death_averted_5up"])){
              last_row[, "death_averted_5up"] <- 0
            }
            if(!is.na(last_row[, "death_averted"])){
              last_row[, "death_averted"] <- 0
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

  return(newdata)
}

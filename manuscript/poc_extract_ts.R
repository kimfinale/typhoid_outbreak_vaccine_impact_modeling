setwd(Sys.getenv("RENEWAL_ROOT", "."))
ts <- read.csv("data/Typhoid_Outbreak_Time_Series_2000_2022_Timeseries.csv", stringsAsFactors = FALSE)
names(ts)[1:9] <- c("study", "start", "end", "n_pat", "suspected", "confirmed", "probable", "onset_tf", "conf_tf")
ts$study <- trimws(ts$study)
for (s in c("Kabwama 2017", "Neil 2012", "Aye 2004")) {
  d <- ts[ts$study == s, c("start", "end", "n_pat", "suspected", "confirmed")]
  cat("\n=== ", s, " (", nrow(d), " bins) ===\n", sep = "")
  print(d, row.names = FALSE)
  cat("sum suspected:", sum(as.numeric(d$suspected), na.rm = TRUE),
      " sum confirmed:", sum(as.numeric(d$confirmed), na.rm = TRUE), "\n")
}

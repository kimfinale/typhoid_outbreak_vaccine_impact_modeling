# Consolidate the background-incidence / outbreak-recurrence extractions.
#   Rscript manuscript/consolidate_bg.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
dir <- "manuscript/data/bg_extract"
files <- list.files(dir, pattern = "^b_\\d+\\.csv$", full.names = TRUE)
rows <- lapply(files, function(f) tryCatch(read.csv(f, colClasses = "character", check.names = FALSE),
                                           error = function(e) NULL))
d <- do.call(rbind, Filter(Negate(is.null), rows))
d <- d[order(as.integer(d$idx)), ]
d$has_bg  <- !(is.na(d$background_incidence) | d$background_incidence %in% c("NA", ""))
d$has_freq <- !(is.na(d$outbreak_frequency) | d$outbreak_frequency %in% c("NA", ""))
write.csv(d, "manuscript/data/tf_background_recurrence.csv", row.names = FALSE)
cat("=== extracted:", nrow(d), "papers | with background_incidence:", sum(d$has_bg),
    "| with recurrence info:", sum(d$has_freq), "===\n\n")
for (i in seq_len(nrow(d))) {
  cat(sprintf("• %s [%s]\n", d$study[i], d$region[i]))
  if (d$has_bg[i])   cat("    bg  :", d$background_incidence[i], "\n")
  if (d$has_freq[i]) cat("    freq:", d$outbreak_frequency[i], "\n")
}

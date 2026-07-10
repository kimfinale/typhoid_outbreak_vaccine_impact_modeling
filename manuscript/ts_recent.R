# Build the ORI time-series candidate pool restricted to studies published >= 2023
# (the tf outbreak dataset already covers through 2022). Extracts publication year
# from the MEDLINE DP field, joins to the parsed records + TS prefilter signal.
#   Rscript manuscript/ts_recent.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
lines <- readLines("manuscript/data/pubmed-typhoidtia-set.txt", warn = FALSE, encoding = "UTF-8")

pmid <- NA_character_; yr <- NA_integer_; res <- list()
for (ln in lines) {
  if (grepl("^PMID- ", ln)) {
    if (!is.na(pmid)) res[[length(res) + 1]] <- c(pmid, yr)
    pmid <- trimws(sub("^PMID- ", "", ln)); yr <- NA_integer_
  } else if (grepl("^DP  - ", ln) && is.na(yr)) {
    y <- sub("^DP  - ", "", ln); yr <- suppressWarnings(as.integer(sub("^.*?(\\d{4}).*$", "\\1", y)))
  }
}
if (!is.na(pmid)) res[[length(res) + 1]] <- c(pmid, yr)
yrs <- data.frame(pmid = vapply(res, `[`, "", 1), year = as.integer(vapply(res, `[`, "", 2)), stringsAsFactors = FALSE)

parsed <- read.csv("manuscript/data/pubmed_parsed.csv", colClasses = "character")
parsed$pmid <- trimws(parsed$pmid)
cand   <- read.csv("manuscript/data/ts_candidates.csv", colClasses = "character")  # has ts_signal
m <- merge(parsed, yrs, by = "pmid", all.x = TRUE)
m$ts_signal <- cand$ts_signal[match(m$pmid, trimws(cand$pmid))]
m$ts_signal[is.na(m$ts_signal)] <- "none"

cat("year coverage: parsed", nrow(parsed), " with-year", sum(!is.na(m$year)), "\n")
cat("year>=2023 (all):", sum(m$year >= 2023, na.rm = TRUE), "\n")
cat("year>=2023 AND ts_signal strong:", sum(m$year >= 2023 & m$ts_signal == "strong", na.rm = TRUE), "\n")
cat("year>=2023 AND ts_signal strong|weak:", sum(m$year >= 2023 & m$ts_signal %in% c("strong","weak"), na.rm = TRUE), "\n")
cat("year distribution (2020+):\n"); print(table(m$year[m$year >= 2020]))

# pool = 2023+ records carrying ANY outbreak/temporal signal (strong or weak), max recall on the recent set
pool <- m[!is.na(m$year) & m$year >= 2023 & m$ts_signal %in% c("strong", "weak"),
          c("pmid", "title", "abstract", "language", "pubtype", "year", "ts_signal")]
pool <- pool[order(-as.integer(pool$ts_signal == "strong"), pool$pmid), ]
write.csv(pool, "manuscript/data/ts_recent_pool.csv", row.names = FALSE)
cat("\nwrote manuscript/data/ts_recent_pool.csv (", nrow(pool), "records )\n")

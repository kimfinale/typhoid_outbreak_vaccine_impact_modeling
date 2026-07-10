# Consolidate the TS-confirmation outputs into ori_timeseries_candidates.csv.
# All records are >=2023 by construction, so they are new relative to the tf
# outbreak dataset (2000-2022) -- no further de-dup vs tf needed.
#   Rscript manuscript/ts_consolidate.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
rdir <- "manuscript/data/ts_batches"
cols <- c("pmid", "timeseries", "ts_type", "period", "reason")
files <- list.files(rdir, pattern = "_screened\\.csv$", full.names = TRUE)

L <- lapply(files, function(f) {
  d <- read.csv(f, colClasses = "character", check.names = FALSE)
  for (k in cols) if (!k %in% names(d)) d[[k]] <- NA_character_
  d[, cols]
})
S <- do.call(rbind, L); S$pmid <- trimws(S$pmid); S <- S[!duplicated(S$pmid), ]
S$timeseries <- tolower(trimws(S$timeseries))

pool <- read.csv("manuscript/data/ts_recent_pool.csv", colClasses = "character"); pool$pmid <- trimws(pool$pmid)
missing <- setdiff(pool$pmid, S$pmid); extra <- setdiff(S$pmid, pool$pmid)
cat(sprintf("pool=%d screened=%d missing=%d extra=%d\n", nrow(pool), nrow(S), length(missing), length(extra)))
cat("\n== timeseries decision (all 458) ==\n"); print(table(S$timeseries, useNA = "ifany"))

out <- merge(pool[, c("pmid", "title", "year", "language", "ts_signal")], S, by = "pmid", all.x = TRUE)
cand <- out[out$timeseries %in% c("yes", "maybe"), c("pmid", "timeseries", "ts_type", "period", "year", "title", "reason")]
cand <- cand[order(match(cand$timeseries, c("yes", "maybe")), cand$ts_type, cand$pmid), ]
write.csv(cand, "manuscript/data/ori_timeseries_candidates.csv", row.names = FALSE)

cat("\n== ts_type x decision (yes/maybe) ==\n"); print(table(ts_type = cand$ts_type, decision = cand$timeseries))
cat("\nwrote manuscript/data/ori_timeseries_candidates.csv (", nrow(cand), "rows )\n")
if (length(missing)) writeLines(missing, file.path(rdir, "_ts_missing.txt"))

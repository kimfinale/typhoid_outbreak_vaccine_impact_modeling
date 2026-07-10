setwd(Sys.getenv("RENEWAL_ROOT", "."))
m <- read.csv("manuscript/data/ts_yes_manifest.csv", colClasses = "character")
l <- read.csv("manuscript/data/pdfs/timeseries/_fetch_log.csv", colClasses = "character")
got <- l$pmid[l$source != "FAILED"]
m$got <- m$pmid %in% got
cat("=== OBTAINED (", sum(m$got), ") by ts_type ===\n")
g <- m[m$got, ]
g <- g[order(match(g$ts_type, c("outbreak_curve", "dated_counts", "surveillance_trend"))), ]
for (i in seq_len(nrow(g)))
  cat(sprintf("  %s [%s | %s] %s\n", g$pmid[i], g$ts_type[i], g$period[i], substr(g$title[i], 1, 62)))
cat("\n=== FAILED (", sum(!m$got), ") ===\n")
f <- m[!m$got, ]
for (i in seq_len(nrow(f)))
  cat(sprintf("  %s [%s] %s\n", f$pmid[i], f$ts_type[i], substr(f$title[i], 1, 62)))

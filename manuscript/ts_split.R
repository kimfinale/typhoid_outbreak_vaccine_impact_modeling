# Split ts_recent_pool.csv into <=120-record chunks for the TS-confirmation pass.
#   Rscript manuscript/ts_split.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
d <- read.csv("manuscript/data/ts_recent_pool.csv", colClasses = "character")
dir.create("manuscript/data/ts_batches", showWarnings = FALSE)
sz <- 115; np <- ceiling(nrow(d) / sz)
for (p in seq_len(np)) {
  idx <- ((p - 1) * sz + 1):min(p * sz, nrow(d))
  write.csv(d[idx, ], sprintf("manuscript/data/ts_batches/ts_%02d.csv", p), row.names = FALSE)
  cat(sprintf("ts_%02d.csv : %d records\n", p, length(idx)))
}
cat("wrote", np, "TS chunks\n")

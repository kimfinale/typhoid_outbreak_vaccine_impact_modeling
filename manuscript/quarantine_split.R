# Split a filter-blocked re-screen sub-file into 50-record quarantine chunks.
#   Rscript manuscript/quarantine_split.R rs_13_2
setwd(Sys.getenv("RENEWAL_ROOT", "."))
stub <- commandArgs(trailingOnly = TRUE)[1]
rdir <- "manuscript/data/batches/rescreen"
d <- read.csv(file.path(rdir, paste0(stub, ".csv")), colClasses = "character")
sz <- 50; np <- ceiling(nrow(d) / sz)
for (p in seq_len(np)) {
  idx <- ((p - 1) * sz + 1):min(p * sz, nrow(d))
  write.csv(d[idx, ], file.path(rdir, sprintf("q_%s_%d.csv", stub, p)), row.names = FALSE)
}
cat("wrote", np, "quarantine chunks of <=", sz, "for", stub, "\n")

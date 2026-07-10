# Split the FLAGGED batches into small (<=250) sub-files for a clean re-screen.
# Reads the canonical, uncontaminated repo inputs batch_XX.csv (NOT scratchpad).
#   Rscript manuscript/split_rescreen.R 03 09 13 15
setwd(Sys.getenv("RENEWAL_ROOT", "."))
bids <- commandArgs(trailingOnly = TRUE)
if (!length(bids)) stop("pass batch ids, e.g. 03 09 13 15")
bdir <- "manuscript/data/batches"
outdir <- file.path(bdir, "rescreen")
dir.create(outdir, showWarnings = FALSE)
sz <- 250; man <- character(0)
for (bid in bids) {
  inf <- file.path(bdir, sprintf("batch_%s.csv", bid))
  d <- read.csv(inf, colClasses = "character")
  n <- nrow(d); np <- ceiling(n / sz)
  for (p in seq_len(np)) {
    idx <- ((p - 1) * sz + 1):min(p * sz, n)
    f <- file.path(outdir, sprintf("rs_%s_%d.csv", bid, p))
    write.csv(d[idx, ], f, row.names = FALSE)
    man <- c(man, sprintf("rs_%s_%d.csv : %d records", bid, p, length(idx)))
  }
}
cat(paste(man, collapse = "\n"), "\n")
cat("wrote", length(man), "sub-files to", outdir, "\n")

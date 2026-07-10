# Assemble the re-screen sub-file outputs (rescreen/rs_<bid>_<p>_screened.csv) into
# per-batch batch_<bid>_screened.csv, so the STRICT consolidation accepts them.
# Verifies each assembled batch's PMID set == canonical batch_<bid>.csv before writing.
#   Rscript manuscript/assemble_rescreen.R 03 09 13 15
setwd(Sys.getenv("RENEWAL_ROOT", "."))
bids <- commandArgs(trailingOnly = TRUE)
if (!length(bids)) stop("pass batch ids, e.g. 03 09 13 15")
bdir <- "manuscript/data/batches"; rdir <- file.path(bdir, "rescreen")
keep <- c("pmid", "decision", "layer", "reason", "citation_chase", "language_exclude")

for (bid in bids) {
  parts <- sort(Sys.glob(file.path(rdir, sprintf("rs_%s_*_screened.csv", bid))))
  if (!length(parts)) { cat(sprintf("batch %s: NO screened sub-files yet\n", bid)); next }
  d <- do.call(rbind, lapply(parts, function(f) {
    x <- read.csv(f, colClasses = "character", check.names = FALSE)
    for (k in keep) if (!k %in% names(x)) x[[k]] <- NA_character_
    x[, keep]
  }))
  d$pmid <- trimws(d$pmid); d <- d[!duplicated(d$pmid), ]
  inpm <- trimws(read.csv(file.path(bdir, sprintf("batch_%s.csv", bid)), colClasses = "character")$pmid)
  miss <- setdiff(inpm, d$pmid); extra <- setdiff(d$pmid, inpm)
  cat(sprintf("batch %s: %d parts, %d rows; missing=%d extra=%d\n",
              bid, length(parts), nrow(d), length(miss), length(extra)))
  if (length(miss) == 0 && length(extra) == 0) {
    d <- d[match(inpm, d$pmid), ]
    write.csv(d, file.path(bdir, sprintf("batch_%s_screened.csv", bid)), row.names = FALSE)
    cat(sprintf("  -> wrote batch_%s_screened.csv (VERIFIED)\n", bid))
  } else {
    if (length(miss)) writeLines(miss, file.path(rdir, sprintf("_still_missing_%s.txt", bid)))
    cat(sprintf("  -> NOT assembled (coverage gap); see _still_missing_%s.txt\n", bid))
  }
}

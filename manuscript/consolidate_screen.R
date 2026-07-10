# STRICT consolidation of the title/abstract screen.
# Accept a batch's screened output ONLY if manuscript/data/batches/batch_XX_screened.csv
# exists and its PMID set EXACTLY equals the canonical batch_XX.csv input set. Any batch
# that is missing / unreadable / mismatched (e.g. nested batches whose parents never wrote
# a verified merged file, or scratchpad contamination) is FLAGGED and its records are dumped
# — sourced from the master parse, never scratchpad — to _rescreen_input.csv for a clean pass.
# Only when all 15 batches are OK does it write the final screening_ta_pubmed.csv.
#   Rscript manuscript/consolidate_screen.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
bdir <- "manuscript/data/batches"
keep <- c("pmid", "decision", "layer", "reason", "citation_chase", "language_exclude")

read_screened <- function(f) {
  l1 <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) character(0))
  if (!length(l1)) return(NULL)
  has_header <- grepl("(^|,)\\s*pmid\\s*(,|$)", l1, ignore.case = TRUE) && grepl("decision", l1, ignore.case = TRUE)
  d <- tryCatch({
    if (has_header) read.csv(f, colClasses = "character", check.names = FALSE)
    else { x <- read.csv(f, header = FALSE, colClasses = "character", check.names = FALSE)
           if (ncol(x) >= 6) names(x)[1:6] <- keep; x }
  }, error = function(e) NULL)
  if (is.null(d) || !all(c("pmid", "decision") %in% names(d))) return(NULL)
  for (k in keep) if (!k %in% names(d)) d[[k]] <- NA_character_
  d <- d[, keep]; d$pmid <- trimws(d$pmid)
  d[grepl("^[0-9]+$", d$pmid), ]
}

master <- read.csv("manuscript/data/pubmed_parsed.csv", colClasses = "character")
master$pmid <- trimws(master$pmid)

inputs <- sort(Sys.glob(file.path(bdir, "batch_[0-9][0-9].csv")))
clean <- list(); flagged <- character(0); report <- character(0)
for (inf in inputs) {
  bid  <- sub(".*batch_([0-9]+)\\.csv$", "\\1", inf)
  inpm <- trimws(read.csv(inf, colClasses = "character")$pmid)
  scf  <- file.path(bdir, sprintf("batch_%s_screened.csv", bid))
  status <- "MISSING-FILE"
  if (file.exists(scf)) {
    d <- read_screened(scf)
    if (is.null(d)) status <- "UNREADABLE"
    else {
      d <- d[!duplicated(d$pmid), ]
      if (setequal(d$pmid, inpm)) { clean[[bid]] <- d; status <- sprintf("OK (%d)", nrow(d)) }
      else status <- sprintf("MISMATCH have=%d/%d missing=%d extra=%d",
                             length(intersect(d$pmid, inpm)), length(inpm),
                             length(setdiff(inpm, d$pmid)), length(setdiff(d$pmid, inpm)))
    }
  }
  if (!grepl("^OK", status)) flagged <- c(flagged, bid)
  report <- c(report, sprintf("  batch %s: %s", bid, status))
}
cat("=== per-batch verification (screened PMID-set == canonical input) ===\n")
cat(paste(report, collapse = "\n"), "\n")
cat(sprintf("\nCLEAN=%d  FLAGGED=%s\n", length(clean), if (length(flagged)) paste(flagged, collapse = ",") else "none"))

if (length(flagged)) {
  fp <- unique(unlist(lapply(flagged, function(bid)
    trimws(read.csv(file.path(bdir, sprintf("batch_%s.csv", bid)), colClasses = "character")$pmid))))
  mm <- master[master$pmid %in% fp, c("pmid", "title", "abstract", "language", "pubtype")]
  write.csv(mm, file.path(bdir, "_rescreen_input.csv"), row.names = FALSE)
  cat(sprintf("\nNOT FINAL. Wrote %d records (batches %s) -> _rescreen_input.csv. Re-screen cleanly, then re-run.\n",
              nrow(mm), paste(flagged, collapse = ",")))
  quit(save = "no")
}

# all 15 clean -> finalise
S <- do.call(rbind, clean); S$decision <- tolower(trimws(S$decision)); S <- S[!duplicated(S$pmid), ]
missing <- setdiff(master$pmid, S$pmid); extra <- setdiff(S$pmid, master$pmid)
cat(sprintf("\nmaster=%d  screened=%d  missing=%d  extra=%d\n", nrow(master), nrow(S), length(missing), length(extra)))

out <- merge(master[, c("pmid", "title", "language")], S, by = "pmid", all.x = TRUE)
out$decision[is.na(out$decision)] <- "UNSCREENED"
out$layer[is.na(out$layer) | out$layer == ""] <- "none"
out$citation_chase   <- tolower(trimws(out$citation_chase))
out$language_exclude <- tolower(trimws(out$language_exclude))
out <- out[order(match(out$decision, c("include", "maybe", "exclude", "UNSCREENED")), out$pmid), ]
write.csv(out, "manuscript/data/screening_ta_pubmed.csv", row.names = FALSE)

cat("\n== decision counts ==\n"); print(table(out$decision, useNA = "ifany"))
adv <- out[out$decision %in% c("include", "maybe"), ]
cat("\n== include/maybe by layer ==\n"); print(table(layer = adv$layer, decision = adv$decision))
cat("\ncitation_chase=yes:", sum(out$citation_chase == "yes", na.rm = TRUE),
    " | language_exclude=yes:", sum(out$language_exclude == "yes", na.rm = TRUE), "\n")
cat("\nFINAL: wrote manuscript/data/screening_ta_pubmed.csv (", nrow(out), "rows )\n")

# Build a download manifest for the 32 yes-rated TS articles: pmid, doi, year, title.
# Extracts DOI from the MEDLINE AID/LID [doi] fields.
#   Rscript manuscript/ts_manifest.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
lines <- readLines("manuscript/data/pubmed-typhoidtia-set.txt", warn = FALSE, encoding = "UTF-8")
doitag <- "[doi]"
pmid <- NA_character_; doi <- NA_character_; res <- list()
for (ln in lines) {
  if (startsWith(ln, "PMID- ")) {
    if (!is.na(pmid)) res[[length(res) + 1]] <- c(pmid, doi)
    pmid <- trimws(sub("^PMID- ", "", ln)); doi <- NA_character_
  } else if (is.na(doi) && grepl(doitag, ln, fixed = TRUE) &&
             (startsWith(ln, "LID - ") || startsWith(ln, "AID - "))) {
    v <- sub("^(LID|AID)\\s*-\\s*", "", ln)
    doi <- trimws(sub("\\s*\\[doi\\].*$", "", v))
  }
}
if (!is.na(pmid)) res[[length(res) + 1]] <- c(pmid, doi)
d <- data.frame(pmid = vapply(res, `[`, "", 1), doi = vapply(res, `[`, "", 2), stringsAsFactors = FALSE)

yes <- read.csv("manuscript/data/ori_timeseries_candidates.csv", colClasses = "character")
yes <- yes[yes$timeseries == "yes", ]
m <- merge(yes[, c("pmid", "year", "title", "ts_type", "period")], d, by = "pmid", all.x = TRUE)
m <- m[order(m$ts_type, m$pmid), c("pmid", "doi", "year", "ts_type", "period", "title")]
write.csv(m, "manuscript/data/ts_yes_manifest.csv", row.names = FALSE)
cat("TS yes articles:", nrow(m), "| with DOI:", sum(!is.na(m$doi) & m$doi != ""), "\n")

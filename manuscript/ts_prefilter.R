# Prefilter the parsed PubMed set for TYPHOID TIME-SERIES candidates (for ORI modeling).
# Broad, recall-first regex over title+abstract; LLM confirms the `timeseries` tag later.
#   Rscript manuscript/ts_prefilter.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
df <- read.csv("manuscript/data/pubmed_parsed.csv", colClasses = "character")
ta <- paste(df$title, df$abstract)

typh <- grepl("typhoid|enteric fever|paratyph|typhi\\b|salmonella enterica", ta, ignore.case = TRUE)

# outbreak / temporal signal (a digitizable epidemic curve or dated case counts)
strong <- grepl(paste0(
  "outbreak|epidemic|epi[- ]?curve|epidemic curve|point[- ]?source|food[- ]?borne|",
  "water[- ]?borne|water[- ]?related|attack rate|time[- ]?series|",
  "\\bweekly\\b|\\bmonthly\\b|\\bdaily\\b|by (week|month|day)|per (week|month|day)|",
  "date of onset|onset date|index case|cluster of|notified cases|reported cases"),
  ta, ignore.case = TRUE)
# weaker temporal / surveillance signal (endemic trend -> likely 'maybe')
weak <- grepl("surveillance|incidence|seasonal|temporal trend|over time|between 20|between 19|from 20|from 19",
              ta, ignore.case = TRUE)

df$ts_signal <- ifelse(strong, "strong", ifelse(weak, "weak", "none"))
df$ts_candidate <- typh & (strong | weak)

cat("total:", nrow(df), " typhoid-context:", sum(typh), "\n")
cat("ts_candidate:", sum(df$ts_candidate),
    " (strong:", sum(typh & strong), " weak-only:", sum(typh & weak & !strong), ")\n")

cand <- df[df$ts_candidate, c("pmid", "title", "abstract", "language", "pubtype", "ts_signal")]
write.csv(cand, "manuscript/data/ts_candidates.csv", row.names = FALSE)
cat("wrote manuscript/data/ts_candidates.csv (", nrow(cand), "rows )\n")

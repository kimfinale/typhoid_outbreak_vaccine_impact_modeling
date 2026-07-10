# Parse the PubMed MEDLINE export into a tidy table + screening batches.
#   Rscript manuscript/parse_pubmed.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
path <- "manuscript/data/pubmed-typhoidtia-set.txt"
lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

recs <- list(); cur <- list(); curtag <- NA_character_
add <- function(cur, tag, v) { cur[[tag]] <- if (is.null(cur[[tag]])) v else paste(cur[[tag]], v); cur }
for (ln in lines) {
  if (grepl("^[A-Z]{2,4}\\s*- ", ln)) {
    curtag <- sub("^([A-Z]{2,4})\\s*-.*$", "\\1", ln)
    v <- sub("^[A-Z]{2,4}\\s*- ", "", ln)
    if (curtag == "PMID") { if (length(cur)) recs[[length(recs) + 1]] <- cur; cur <- list() }
    cur <- add(cur, curtag, v)
  } else if (grepl("^\\s{4,}\\S", ln) && !is.na(curtag) && length(cur)) {
    cur <- add(cur, curtag, trimws(ln))
  }
}
if (length(cur)) recs[[length(recs) + 1]] <- cur

g <- function(r, k) { x <- r[[k]]; if (is.null(x)) NA_character_ else x }
df <- data.frame(
  pmid     = vapply(recs, g, "", "PMID"),
  title    = vapply(recs, g, "", "TI"),
  abstract = vapply(recs, g, "", "AB"),
  language = vapply(recs, g, "", "LA"),
  pubtype  = vapply(recs, g, "", "PT"),
  stringsAsFactors = FALSE)
df <- df[!is.na(df$pmid) & df$pmid != "", ]
cat("parsed records:", nrow(df), "| with abstract:", sum(!is.na(df$abstract) & df$abstract != ""), "\n")
cat("languages:", paste(names(sort(table(df$language), decreasing = TRUE))[1:5], collapse = ", "), "\n")
write.csv(df, "manuscript/data/pubmed_parsed.csv", row.names = FALSE)

# batches of 500 for parallel screening
dir.create("manuscript/data/batches", showWarnings = FALSE)
bs <- 500; nb <- ceiling(nrow(df) / bs)
for (i in seq_len(nb)) {
  idx <- ((i - 1) * bs + 1):min(i * bs, nrow(df))
  write.csv(df[idx, ], sprintf("manuscript/data/batches/batch_%02d.csv", i), row.names = FALSE)
}
cat("wrote", nb, "batches of", bs, "to manuscript/data/batches/\n")

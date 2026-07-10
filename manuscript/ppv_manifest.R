# Build a PPV extraction manifest from candidate_pool.csv: resolve each study to its
# Zotero-stored PDF (keys are in the provenance column) where available.
#   Rscript manuscript/ppv_manifest.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
zstore <- "C:/Users/jonghoon.kim/Zotero/storage"
cp <- read.csv("manuscript/data/candidate_pool.csv", colClasses = "character", check.names = FALSE)

find_pdf <- function(keys) {
  for (k in keys) {
    d <- file.path(zstore, k)
    if (dir.exists(d)) {
      pdfs <- list.files(d, pattern = "\\.pdf$", full.names = TRUE, ignore.case = TRUE)
      if (length(pdfs)) return(pdfs[which.max(file.info(pdfs)$size)])
    }
  }
  NA_character_
}

rows <- lapply(seq_len(nrow(cp)), function(i) {
  prov <- cp$provenance[i]
  keys <- if (grepl("zotero", prov, ignore.case = TRUE))
    unique(regmatches(prov, gregexpr("[A-Z0-9]{8}", prov))[[1]]) else character(0)
  pdf <- if (length(keys)) find_pdf(keys) else NA_character_
  data.frame(study = cp$study[i], layer = cp$layer[i], setting_label = cp$setting_label[i],
             zotero_keys = paste(keys, collapse = "|"),
             pdf_path = pdf, has_pdf = !is.na(pdf),
             data_recoverable = cp$data_recoverable[i], stringsAsFactors = FALSE)
})
m <- do.call(rbind, rows)
write.csv(m, "manuscript/data/ppv_manifest.csv", row.names = FALSE)
cat("PPV studies:", nrow(m), "| local Zotero PDF found:", sum(m$has_pdf),
    "| no PDF (web/none):", sum(!m$has_pdf), "\n\n")
cat("=== studies WITHOUT a local PDF (need web fetch) ===\n")
cat(paste0("  ", m$study[!m$has_pdf]), sep = "\n")

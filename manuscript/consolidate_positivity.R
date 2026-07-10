# Consolidate the per-PDF positivity extractions into one clean table + usability class.
#   Rscript manuscript/consolidate_positivity.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
dir <- "manuscript/data/ppv_extract"
files <- list.files(dir, pattern = "^e_\\d+\\.csv$", full.names = TRUE)
rows <- lapply(files, function(f) tryCatch(read.csv(f, colClasses = "character", check.names = FALSE),
                                           error = function(e) NULL))
d <- do.call(rbind, Filter(Negate(is.null), rows))
num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("NA", "", "na"), NA, x)))
d$suspected <- num(d$suspected); d$tested_bc <- num(d$tested_bc); d$confirmed_bc <- num(d$confirmed_bc)
d$positivity <- ifelse(!is.na(d$tested_bc) & d$tested_bc > 0, round(d$confirmed_bc / d$tested_bc, 3), NA)
d$flags[is.na(d$flags)] <- ""

classify <- function(r) {
  if (grepl("MISFILED", r$notes)) return(c("EXCLUDE", "misfiled/unreadable PDF"))
  if (is.na(r$confirmed_bc)) return(c("EXCLUDE", "no confirmed count"))
  if (is.na(r$tested_bc) || r$tested_bc == 0) return(c("EXCLUDE", "no blood-culture-tested denominator"))
  if (r$tested_bc < 20) return(c("EXCLUDE", "too few tested (<20)"))
  if (!is.na(r$positivity) && r$positivity >= 0.8) return(c("EXCLUDE", "positivity>=0.8 (confirmed-only/selective)"))
  if (grepl("high_income", r$flags)) return(c("EXCLUDE", "high-income point-source"))
  if (grepl("surveillance_program", r$flags) || grepl("surveillance", r$setting))
    return(c("USABLE_surveillance", "all-tested surveillance (endemic stratum)"))
  c("USABLE_outbreak", "community-outbreak; selective testing")
}
cl <- t(sapply(seq_len(nrow(d)), function(i) classify(d[i, ])))
d$usability <- cl[, 1]; d$reason <- cl[, 2]
d <- d[order(d$usability, d$study), ]
keep <- c("idx","study","country","outbreak_year","setting","suspected","tested_bc","confirmed_bc",
          "positivity","usability","reason","flags","case_definition","notes")
write.csv(d[, keep], "manuscript/data/ppv_positivity_table.csv", row.names = FALSE)

cat("=== extracted:", nrow(d), "studies ===\n")
print(table(d$usability))
cat("\n=== USABLE_outbreak (match existing community-syndromic estimand) ===\n")
o <- d[d$usability == "USABLE_outbreak", c("study","suspected","tested_bc","confirmed_bc","positivity")]
print(o, row.names = FALSE)
cat("\n=== USABLE_surveillance (endemic all-tested stratum) ===\n")
s <- d[d$usability == "USABLE_surveillance", c("study","suspected","tested_bc","confirmed_bc","positivity")]
print(s, row.names = FALSE)
cat("\n(existing model community set: Aye 22/3, Neil 63/25, Kabwama 364/51)\n")

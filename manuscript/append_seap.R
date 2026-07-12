# Add the SEAP outpatient cohort (Aiemjoy 2020, user-provided PDF) to the positivity
# table as a large endemic-surveillance point. Enter COUNTS, not the reported "PPV"
# (which is blood-culture-referenced = positivity = pi*Se_BC, not our latent-true pi).
#   Rscript manuscript/append_seap.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
d <- read.csv("manuscript/data/ppv_positivity_table.csv", stringsAsFactors = FALSE)
if (any(grepl("SEAP", d$study))) d <- d[!grepl("Aiemjoy|SEAP 2020", d$study), ]  # idempotent
new <- d[1, ]; new[, ] <- NA
new$idx <- "SEAP"
new$study <- "Aiemjoy/SEAP 2020 (Bangladesh/Nepal/Pakistan)"
new$country <- "multi-country Asia"; new$outbreak_year <- 2016
new$setting <- "population-surveillance"
new$suspected <- 20899; new$tested_bc <- 20899; new$confirmed_bc <- 2116   # S. Typhi (Paratyphi +297)
new$positivity <- round(2116 / 20899, 3)
new$usability <- "USABLE_surveillance"
new$reason <- "all-tested outpatient surveillance; S. Typhi only (Paratyphi +297)"
new$flags <- "surveillance_program"
new$case_definition <- ">=3d fever; all blood-cultured"
new$notes <- "reported PPV 10.5pct is BC-referenced = positivity; user-provided Aiemjoy 2020 PDF"
new$stratum <- "endemic surveillance (all-tested)"
new$implied_pi <- round(min((2116 / 20899) / 0.619, 1), 3)
d <- rbind(d, new)
write.csv(d, "manuscript/data/ppv_positivity_table.csv", row.names = FALSE)
cat("added SEAP: 20899 tested / 2116 S. Typhi -> positivity", new$positivity,
    " implied pi ~", new$implied_pi, "\n")

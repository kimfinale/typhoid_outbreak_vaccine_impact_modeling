# =============================================================================
# run_case_definitions.R — synthesize the case-definition typology: type and
# poolability tallies, a poolability figure, and render the report. The
# extraction itself is the curated CSVs (full-text reading of the Zotero PDFs).
#   Rscript case_definitions/run_case_definitions.R
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(123); cat("seed: 123 (no stochastic computation)\n")
dir.create("case_definitions/figures", showWarnings = FALSE)

D <- read.csv("case_definitions/tab_case_definitions.csv", stringsAsFactors = FALSE, check.names = FALSE)
F_ <- read.csv("case_definitions/tab_fitting_design.csv", stringsAsFactors = FALSE, check.names = FALSE)

ren <- D[D$renewal_set == "Y", ]
cat("\n=== Observation-process type among the 13 renewal outbreaks ===\n")
print(table(ren$type))
cat("\n=== Invariance status (renewal 13) ===\n")
print(table(F_$invariance_status[F_$renewal_set == "Y"]))

# Pool-membership tally for percent-reduction (renewal 13).
fr <- F_[F_$renewal_set == "Y", ]
pool <- fr %>% mutate(pool_class = case_when(
  grepl("EXCLUDE|exclude", recommendation) ~ "exclude (Type D perforation)",
  grepl("caution", recommendation) ~ "include w/ caution (invariance)",
  grepl("true-case", recommendation) ~ "Type B (harmonize absolute)",
  grepl("suspected-series pool", recommendation) ~ "suspected-series pool (A+C; rho via alpha)",
  TRUE ~ "other")) %>% count(pool_class)
cat("\n=== Recommendation classes (renewal 13) ===\n"); print(as.data.frame(pool), row.names = FALSE)

# Figure: type distribution colored by pool recommendation.
dd <- ren %>% count(type) %>% mutate(type = factor(type, levels = c("A","B","C","D","F")))
g <- ggplot(dd, aes(type, n, fill = type)) +
  geom_col() + geom_text(aes(label = n), vjust = -0.3) +
  scale_fill_manual(values = c(A = "#0072B2", B = "#2CA02C", C = "#9467BD",
    D = "#D55E00", F = "grey60"), guide = "none") +
  scale_x_discrete(labels = c(A = "A\nsuspected", B = "B\nconfirmed-only",
    C = "C\nmixed", D = "D\nperforation", F = "F\nambiguous")) +
  labs(x = NULL, y = "outbreaks (of 13)",
       title = "Case-definition observation types (13 renewal outbreaks)",
       subtitle = "D (Muyembe, perforation series) should leave the renewal/theta pool; B (Qamar, Yousafzai) need inverse harmonization") +
  theme_minimal(base_size = 11)
ggsave("case_definitions/figures/fig_type_distribution.png", g, width = 8, height = 5, dpi = 300)

if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  system2("quarto", c("render", shQuote("case_definitions/report_case_definitions.qmd")),
          stdout = FALSE, stderr = FALSE)
  cat("Rendered report_case_definitions.html\n")
}
cat("=== Case-definition synthesis complete ===\n")

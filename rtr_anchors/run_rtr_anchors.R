# =============================================================================
# run_rtr_anchors.R — synthesize the literature R_tr anchors: coverage stats,
# the transparent prop_secondary -> R_tr_household mapping, a coverage figure,
# and render the report. The extraction itself is the curated CSVs (from full-
# text reading of the Zotero PDFs + the Vellore household study).
#   Rscript rtr_anchors/run_rtr_anchors.R
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(dplyr); library(ggplot2) })
set.seed(123); cat("seed: 123 (no stochastic computation)\n")

dir.create("rtr_anchors/figures", showWarnings = FALSE)
A <- read.csv("rtr_anchors/tab_Rtr_anchors.csv", stringsAsFactors = FALSE, check.names = FALSE)
P <- read.csv("rtr_anchors/tab_Rtr_priors_for_model.csv", stringsAsFactors = FALSE, check.names = FALSE)

ren <- A[A$renewal_set == "Y", ]                       # the 13 renewal outbreaks
cov <- ren %>% count(anchor_type) %>% mutate(pct = round(100 * n / sum(n)))
cat("\n=== Coverage among the 13 renewal outbreaks ===\n")
print(as.data.frame(cov), row.names = FALSE)
n_indep <- sum(ren$anchor_type %in% c("independent_level", "independent_ordinal"))
n_level <- sum(ren$anchor_type == "independent_level")
n_pooled <- sum(P$constraint_type[P$study_id %in% ren$study_id] == "pooled_fallback")
cat(sprintf("\nIndependent p2p quantity: %d/13 | level-anchoring (denominator): %d/13 | pooled fallback: %d/13\n",
            n_indep, n_level, n_pooled))

# Transparent mapping: a within-population secondary FRACTION p_sec maps to a
# one-generation household transmission number R_tr_household ~ p_sec/(1-p_sec).
p_sec_polonsky <- 0.053; p_sec_vellore <- 0.125
rtr_hh <- function(p) round(p / (1 - p), 3)
cat(sprintf("\nMapping prop_secondary -> R_tr_household (one-generation): Polonsky 5.3%% -> %.3f; Vellore pooled 12.5%% -> %.3f\n",
            rtr_hh(p_sec_polonsky), rtr_hh(p_sec_vellore)))

# Coverage figure.
cov$anchor_type <- factor(cov$anchor_type,
  levels = c("independent_level", "independent_ordinal", "none"))
g <- ggplot(cov, aes(anchor_type, n, fill = anchor_type)) +
  geom_col() + geom_text(aes(label = n), vjust = -0.3) +
  scale_fill_manual(values = c(independent_level = "#0072B2",
    independent_ordinal = "#9467BD", none = "grey70"), guide = "none") +
  scale_x_discrete(labels = c(independent_level = "independent\n(level: denominator)",
    independent_ordinal = "independent\n(ordinal: OR only)", none = "none")) +
  labs(x = NULL, y = "outbreaks (of 13)",
       title = "Independent R_tr anchors are thin",
       subtitle = "Only 2/13 give a level-anchoring quantity (Polonsky prop-secondary); 3/13 ordinal contact ORs; 8/13 fall back to a pooled prior") +
  theme_minimal(base_size = 11)
ggsave("rtr_anchors/figures/fig_anchor_coverage.png", g, width = 8, height = 5, dpi = 300)

if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  system2("quarto", c("render", shQuote("rtr_anchors/report_Rtr_anchors.qmd")),
          stdout = FALSE, stderr = FALSE)
  cat("Rendered report_Rtr_anchors.html\n")
}
cat("=== R_tr anchors synthesis complete ===\n")

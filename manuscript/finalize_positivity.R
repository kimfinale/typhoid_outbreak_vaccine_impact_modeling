# Finalize the positivity table into 2 clean strata + implied-pi gradient figure.
#   Rscript manuscript/finalize_positivity.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages(library(ggplot2))
d <- read.csv("manuscript/data/ppv_positivity_table.csv", stringsAsFactors = FALSE)
d$positivity <- ifelse(!is.na(d$tested_bc) & d$tested_bc > 0, d$confirmed_bc / d$tested_bc, NA)

# --- clean re-classification by DENOMINATOR type (labels were noisy) ---------
d$stratum <- "excluded"
usable <- !is.na(d$tested_bc) & d$tested_bc >= 20 & !is.na(d$confirmed_bc) &
          !grepl("high_income", d$flags) & (is.na(d$positivity) | d$positivity < 0.8)
ratio <- ifelse(!is.na(d$suspected) & d$suspected > 0, d$tested_bc / d$suspected, NA)
d$stratum[usable & !is.na(ratio) & ratio < 0.7]  <- "outbreak (selective testing)"
d$stratum[usable & (is.na(ratio) | ratio >= 0.7)] <- "endemic surveillance (all-tested)"
# N'Cho usable subset is a hospital review during the outbreak -> annotate
d$reason[grepl("N'Cho", d$study)] <- "hospital case-mgmt review subset (286 tested); whole outbreak 3187/191 lacks tested denom"

# --- implied community PPV: pi ~ positivity / Se_BC (ref 5 mL median 0.62) ----
Se_ref <- 0.619
d$implied_pi <- pmin(d$positivity / Se_ref, 1)
write.csv(d, "manuscript/data/ppv_positivity_table.csv", row.names = FALSE)

cat("=== strata ===\n"); print(table(d$stratum))
u <- d[d$stratum != "excluded", c("study","suspected","tested_bc","confirmed_bc","positivity","implied_pi","stratum")]
u <- u[order(u$stratum, u$implied_pi), ]; print(u, row.names = FALSE)

# --- figure: implied pi by study, by stratum --------------------------------
u$lab <- sub(" \\(.*", "", u$study)
u$lab <- factor(u$lab, levels = u$lab[order(u$implied_pi)])
inmodel <- c("Aye 2004","Neil 2012","Kabwama 2017")
u$inmodel <- ifelse(sub(" .*","",u$study) %in% sub(" .*","",inmodel), "in current model", "newly extracted")
p <- ggplot(u, aes(implied_pi, lab, colour = stratum, shape = inmodel)) +
  geom_point(size = 3.2) +
  geom_vline(xintercept = c(0.23, 0.65), linetype = 3, colour = "#7570B3") +
  scale_colour_manual(values = c("outbreak (selective testing)" = "#B85C1A",
                                 "endemic surveillance (all-tested)" = "#0E6E78")) +
  scale_shape_manual(values = c("in current model" = 17, "newly extracted" = 16)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(title = "Implied community PPV (π) by setting — 40-PDF positivity extraction",
       subtitle = "π ≈ positivity / Se_BC(5 mL); dotted lines = current model outbreak range (Aye–Kabwama–Neil)",
       x = "Implied PPV  π = P(true typhoid | suspected)", y = NULL, colour = "Stratum", shape = NULL) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12.5), legend.position = "right")
ggsave("manuscript/figures/fig1d_ppv_extraction.png", p, width = 9.5, height = 5, dpi = 130)
cat("\nwrote manuscript/figures/fig1d_ppv_extraction.png\n")

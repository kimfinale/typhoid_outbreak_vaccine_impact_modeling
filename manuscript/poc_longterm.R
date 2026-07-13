# POC: long-term (beyond-outbreak) TCV impact per site, vs the one-time outbreak-averted.
# Uses the background-incidence / recurrence parameters extracted from the outbreak
# papers (BACKGROUND_RECURRENCE.md) + Se_BC from the latent-class fit.
#   Rscript manuscript/poc_longterm.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages(library(ggplot2)); suppressMessages(library(tidyr)); suppressMessages(library(dplyr))
source("renewal/R/longterm.R")
se_bc <- read.csv("latent_class_ppv/tables/final_parameters.csv")
Se <- se_bc$med[se_bc$param == "Se_BC_5mL"]                      # 0.62

# --- per-site parameters (grounded in the extraction; GBD-approx where noted) ----
# endemic_inc = CONFIRMED endemic incidence /100k/yr; recur = outbreaks/yr;
# mean_ob_true = true cases per recurrent outbreak in the campaign cohort.
sites <- data.frame(
  site = c("Kampala (Kabwama)", "Harare (Poncin campaign)", "Vellore (SEFI)", "Fiji Northern (Scobie)", "Songkhla (Limpitikul)"),
  Npop = c(500000, 373027, 150000, 130000, 120000),
  endemic_inc = c(150, 50, 1977, 42, 1.8),                      # /100k/yr confirmed
  recur = c(0.2, 1.0, 0.0, 0.2, 0.1),                            # outbreaks/yr
  mean_ob_true = c(1000, 1300, 0, 200, 300),                     # true cases/recurrent outbreak
  outbreak_averted = c(2058, 1300, 80, 300, 400),               # one-time (renewal/proxy), true
  inc_source = c("GBD-approx Uganda", "recurrence-dominated (Poncin: annual since 2010)",
                 "SEFI 1977/100k (measured)", "Scobie 33-52/100k (measured)", "Limpitikul 1.8/100k (measured)"),
  stringsAsFactors = FALSE)

# base scenario: VE0 0.80, coverage 0.85, tau 15y (TCV slow waning), horizon 7y
lt <- lapply(seq_len(nrow(sites)), function(i) with(sites[i, ],
  longterm_averted(Npop, endemic_inc, recur, mean_ob_true, VE0 = .80, coverage = .85,
                   tau_years = 15, horizon_years = 7, se_bc = Se)))
sites$lt_total   <- sapply(lt, `[[`, "averted")
sites$lt_endemic <- sapply(lt, `[[`, "averted_endemic")
sites$lt_recur   <- sapply(lt, `[[`, "averted_recur")
sites$ratio_lt_vs_outbreak <- round(sites$lt_total / sites$outbreak_averted, 1)
write.csv(sites, "manuscript/data/longterm_sites.csv", row.names = FALSE)

cat("=== Long-term vs one-time outbreak-averted (base: VE0 0.80, cov 0.85, tau 15y, H 7y) ===\n")
for (i in seq_len(nrow(sites))) cat(sprintf(
  "  %-26s outbreak=%5.0f  long-term(7y)=%7.0f  (endemic %5.0f + recur %5.0f)  = %4.1fx outbreak\n",
  sites$site[i], sites$outbreak_averted[i], sites$lt_total[i],
  sites$lt_endemic[i], sites$lt_recur[i], sites$ratio_lt_vs_outbreak[i]))

# --- sensitivity: protection horizon H = 3,5,7,10 years -------------------------
Hs <- c(3, 5, 7, 10)
sens <- do.call(rbind, lapply(Hs, function(H) {
  a <- sapply(seq_len(nrow(sites)), function(i) with(sites[i, ],
    longterm_averted(Npop, endemic_inc, recur, mean_ob_true, tau_years = 15, horizon_years = H, se_bc = Se)$averted))
  data.frame(site = sites$site, H = H, lt = a)
}))

# --- figures -------------------------------------------------------------------
p1 <- sites %>% select(site, Outbreak = outbreak_averted, `Endemic (7y)` = lt_endemic, `Recurrent outbreaks (7y)` = lt_recur) %>%
  pivot_longer(-site, names_to = "component", values_to = "cases") %>%
  mutate(component = factor(component, levels = c("Outbreak", "Endemic (7y)", "Recurrent outbreaks (7y)"))) %>%
  ggplot(aes(reorder(site, cases), cases, fill = component)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("Outbreak" = "#9E9E9E", "Endemic (7y)" = "#0E6E78", "Recurrent outbreaks (7y)" = "#B85C1A")) +
  labs(title = "TCV impact: one-time outbreak vs long-term (7-yr protection)",
       subtitle = "True typhoid cases averted per site; long-term = endemic + averted future outbreaks",
       x = NULL, y = "True typhoid cases averted", fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))
ggsave("manuscript/figures/fig3a_longterm_vs_outbreak.png", p1, width = 9, height = 5, dpi = 130)

p2 <- ggplot(sens, aes(H, lt, colour = site)) + geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Long-term averted vs protection window", subtitle = "Grows with the duration of TCV protection",
       x = "Protection horizon H (years)", y = "True typhoid cases averted (long-term)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))
ggsave("manuscript/figures/fig3b_longterm_horizon.png", p2, width = 9, height = 5, dpi = 130)
cat("\nwrote manuscript/data/longterm_sites.csv + fig3a/fig3b\n")

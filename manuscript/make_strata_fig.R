# Figure: community PPV by endemicity stratum (two-stratum fit).
#   Rscript manuscript/make_strata_fig.R -> manuscript/figures/fig1e_ppv_strata.png
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages(library(ggplot2))
res <- read.csv("latent_class_ppv/tables/community_ppv_strata.csv", stringsAsFactors = FALSE)
par <- read.csv("latent_class_ppv/tables/community_ppv_strata_params.csv", stringsAsFactors = FALSE)
g <- function(p) par[par$param == p, c("med", "lo", "hi")]
surv <- g("pi_surv"); out <- g("pi_out"); or <- g("OR_out_vs_surv")

res$lab <- sub(" \\(.*", "", res$study)
res <- res[order(res$stratum, res$pi_med), ]
res$lab <- factor(res$lab, levels = res$lab)
cols <- c("outbreak (selective testing)" = "#B85C1A", "endemic surveillance (all-tested)" = "#0E6E78")

p <- ggplot(res, aes(pi_med, lab, colour = stratum)) +
  annotate("rect", xmin = surv$lo, xmax = surv$hi, ymin = -Inf, ymax = Inf, fill = "#0E6E78", alpha = .10) +
  annotate("rect", xmin = out$lo,  xmax = out$hi,  ymin = -Inf, ymax = Inf, fill = "#B85C1A", alpha = .10) +
  geom_vline(xintercept = surv$med, colour = "#0E6E78", linetype = 2) +
  geom_vline(xintercept = out$med,  colour = "#B85C1A", linetype = 2) +
  geom_errorbarh(aes(xmin = pi_lo, xmax = pi_hi), height = .25, linewidth = .6) +
  geom_point(size = 3) +
  annotate("text", x = surv$med, y = 0.6, label = sprintf("surveillance\nπ=%.2f", surv$med), colour = "#0E6E78", size = 3, vjust = 0) +
  annotate("text", x = out$med,  y = 0.6, label = sprintf("outbreak\nπ=%.2f", out$med), colour = "#B85C1A", size = 3, vjust = 0) +
  scale_colour_manual(values = cols) + scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, .8, .2)) +
  labs(title = sprintf("Community PPV by endemicity: outbreak vs surveillance (OR %.1f [%.1f–%.1f])",
                       or$med, or$lo, or$hi),
       subtitle = "Two-stratum partial-pooling fit on the 40-PDF positivity extraction; bands = stratum-mean 95% CrI",
       x = "PPV  π = P(true typhoid | suspected)", y = NULL, colour = "Stratum") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 12))
ggsave("manuscript/figures/fig1e_ppv_strata.png", p, width = 9, height = 5, dpi = 130)
cat("wrote manuscript/figures/fig1e_ppv_strata.png\n")

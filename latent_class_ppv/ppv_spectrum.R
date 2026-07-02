# =============================================================================
# PPV across the spectrum of clinical case definitions x prevalence.
# PPV = prev*Se_clin / [prev*Se_clin + (1-prev)*(1-Sp_clin)]   (Bayes)
# Se_clin/Sp_clin per definition from the literature (data/clinical_definition_accuracy.csv);
# overlays real anchor points (Pemba WHO defs; our outbreak PPV estimates).
#   Rscript latent_class_ppv/ppv_spectrum.R
# =============================================================================
ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(ggplot2); library(dplyr) })
dir.create("latent_class_ppv/figures", showWarnings = FALSE)
dir.create("latent_class_ppv/tables",  showWarnings = FALSE)

defs <- read.csv("latent_class_ppv/data/clinical_definition_accuracy.csv", check.names = FALSE, stringsAsFactors = FALSE)
defs <- defs[order(defs$order), ]
ppv  <- function(prev, Se, Sp) prev*Se / (prev*Se + (1-prev)*(1-Sp))

# --- PPV curves over prevalence (typhoid among febrile presenters) -----------
prev <- 10^seq(log10(0.01), log10(0.7), length.out = 200)
curve <- do.call(rbind, lapply(seq_len(nrow(defs)), function(i)
  data.frame(variant = defs$variant[i], tier = defs$tier[i], prev = prev,
             PPV = ppv(prev, defs$Se_clin[i], defs$Sp_clin[i]))))
curve$variant <- factor(curve$variant, levels = defs$variant)

# --- summary table: PPV at reference prevalences -----------------------------
refp <- c(`endemic 0.03`=0.03, `moderate 0.10`=0.10, `outbreak 0.30`=0.30, `high 0.50`=0.50)
tab <- do.call(rbind, lapply(seq_len(nrow(defs)), function(i)
  data.frame(variant = defs$variant[i], Se = defs$Se_clin[i], Sp = defs$Sp_clin[i],
             t(round(ppv(refp, defs$Se_clin[i], defs$Sp_clin[i]), 2)))))
names(tab)[4:7] <- names(refp)
write.csv(tab, "latent_class_ppv/tables/ppv_spectrum_table.csv", row.names = FALSE)
cat("=== PPV by definition x prevalence (typhoid among febrile) ===\n"); print(tab, row.names = FALSE)

# --- real anchor points ------------------------------------------------------
# Pemba (Thriemer 2012): WHO suspected PPV 2.7% and WHO probable PPV 73.1% at prev ~0.021
pemba <- data.frame(label = c("Pemba WHO-suspected","Pemba WHO-probable"),
                    prev = c(0.021, 0.021), PPV = c(0.027, 0.731))
# our outbreak PPV estimates (syndromic definitions): back out implied prevalence
# from the WHO-suspected curve (Se 0.83, Sp 0.36) to place them consistently
Se0 <- 0.83; Sp0 <- 0.36
implied_prev <- function(pi) pi*(1-Sp0) / (Se0*(1-pi) + pi*(1-Sp0))
ours <- data.frame(label = c("Aye 2004","Kabwama 2017","Neil 2012"),
                   pi = c(0.26, 0.23, 0.65))
ours$prev <- implied_prev(ours$pi); ours$PPV <- ours$pi
cat("\nOur outbreak PPVs placed on the WHO-suspected curve -> implied typhoid-among-febrile prevalence:\n")
print(data.frame(study=ours$label, PPV=ours$pi, implied_prev=round(ours$prev,2)), row.names=FALSE)

# --- figure ------------------------------------------------------------------
gg <- ggplot(curve, aes(prev, PPV, colour = variant)) +
  annotate("rect", xmin = 0.02, xmax = 0.08, ymin = 0, ymax = 1, alpha = 0.06, fill = "grey40") +
  annotate("text", x = 0.045, y = 0.97, label = "endemic\nsurveillance", size = 3, colour = "grey40") +
  annotate("rect", xmin = 0.2, xmax = 0.5, ymin = 0, ymax = 1, alpha = 0.06, fill = "#D55E00") +
  annotate("text", x = 0.31, y = 0.05, label = "outbreak", size = 3, colour = "#D55E00") +
  geom_line(linewidth = 1) +
  geom_point(data = pemba, aes(prev, PPV), inherit.aes = FALSE, shape = 17, size = 3) +
  geom_text(data = pemba, aes(prev, PPV, label = label), inherit.aes = FALSE,
            hjust = -0.08, vjust = 0.4, size = 3) +
  geom_point(data = ours, aes(prev, PPV), inherit.aes = FALSE, shape = 21, fill = "white", size = 3) +
  geom_text(data = ours, aes(prev, PPV, label = label), inherit.aes = FALSE,
            vjust = -0.8, size = 3) +
  scale_x_log10(labels = scales::percent, breaks = c(0.01,0.02,0.05,0.1,0.2,0.5)) +
  scale_colour_viridis_d(option = "D", end = 0.92, name = "case definition (broad -> strict)") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "typhoid prevalence among febrile presenters (log scale)",
       y = "PPV of the case definition",
       title = "PPV across the spectrum of typhoid case definitions",
       subtitle = "syndromic-only definitions (Sp ~0.2-0.36) cap PPV; adding a serologic tier (Sp ~0.85-1.0) lifts it steeply\ntriangles = Pemba WHO defs (Thriemer 2012); circles = our outbreak PPV estimates on the syndromic curve") +
  theme_minimal(base_size = 11) + theme(legend.position = "right")
ggsave("latent_class_ppv/figures/fig_ppv_spectrum.png", gg, width = 9.5, height = 5.5, dpi = 300)
cat("\n=== PPV-spectrum figure written ===\n")

# Allocation-scenario figure for the example outbreak (Kabwama 2017).
#   Rscript revision/make_allocation_figure.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("revision/allocation_scenarios.R")
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
OB <- Sys.getenv("OUTBREAK", "Kabwama 2017")

r <- run_allocation_scenarios(OB)
cat(sprintf("%s  PPV=%.2f  protect wk %d/%d\n", OB, r$pi, r$t_eff, r$Tn))
print(r$summ, row.names=FALSE)

cols <- c(Multiplicative="#0072B2", Additive="#D55E00", `Time-varying PPV`="#009E73",
          `Renewal-consistent`="#CC79A7", `Reporting-shock`="#E69F00")
th <- theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey93"), legend.position="top", legend.title=element_blank(),
        plot.title=element_text(face="bold", size=12.5), plot.subtitle=element_text(colour="grey40", size=9.5))

# Panel A: the five true-typhoid allocations (same cumulative area), over suspected
susp <- unique(r$df[,c("week","S")])
gA <- ggplot() +
  geom_area(data=susp, aes(week, S), fill="grey88") +
  geom_line(data=r$df, aes(week, I, colour=alloc), linewidth=1) +
  geom_vline(xintercept=r$t_eff, linetype="dotted", colour="grey45") +
  scale_colour_manual(values=cols) +
  labs(x=NULL, y="Cases / week",
       title=sprintf("A  Five weekly allocations of true typhoid - %s (PPV %.2f)", OB, r$pi),
       subtitle="Grey = reported suspected. Each coloured curve has the same cumulative total (pi x suspected) but a different temporal shape; dotted = protection onset") + th

# Panel B: reduction under each allocation (true = point, suspected = point, gap = PPV dilution)
s <- r$summ; s$alloc <- factor(s$alloc, levels=rev(names(cols)))
gB <- ggplot(s) +
  geom_segment(aes(y=alloc, yend=alloc, x=susp_red, xend=true_red), colour="grey75", linewidth=1) +
  geom_point(aes(y=alloc, x=true_red, colour=alloc), size=4) +
  geom_point(aes(y=alloc, x=susp_red), colour="#8A8A8A", size=3.5) +
  geom_text(aes(y=alloc, x=true_red, label=sprintf("true %.0f%%", true_red)), hjust=-0.2, size=3.1, colour="grey20") +
  geom_text(aes(y=alloc, x=susp_red, label=sprintf("%.0f%%", susp_red)), hjust=1.25, size=3.1, colour="grey45") +
  scale_colour_manual(values=cols, guide="none") +
  scale_x_continuous(limits=c(0, max(s$true_red)*1.18), labels=function(x)paste0(x,"%")) +
  labs(x="Reduction after one vaccination campaign", y=NULL,
       title="B  Impact depends on allocation; the surveillance dilution does not",
       subtitle=sprintf("Coloured = true typhoid averted; grey = visible in suspected surveillance. Dilution is exactly PPV (%.2f) in every scenario.", r$pi)) +
  th + theme(panel.grid.major.y=element_blank())

ok <- requireNamespace("patchwork", quietly=TRUE)
fig <- if (ok) { library(patchwork); (gA/gB) + patchwork::plot_layout(heights=c(1.15,1)) } else gA
ggsave(file.path(outdir, "fig_allocation_scenarios.png"), fig, width=8.5, height=8.5, dpi=300)
cat("\nWrote", file.path(outdir,"fig_allocation_scenarios.png"),
    sprintf("\nTrue-reduction range across allocations: %.0f%%-%.0f%%  (suspected %.0f%%-%.0f%%)\n",
            min(r$summ$true_red), max(r$summ$true_red), min(r$summ$susp_red), max(r$summ$susp_red)))

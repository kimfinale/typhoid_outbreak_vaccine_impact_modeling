# Schematic of the latent-class PPV model (latent_class_ppv/stan/model_final.stan).
#   Rscript manuscript/make_model_schematic.R  -> manuscript/figures/fig0_latent_class_model.png
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages(library(ggplot2))
figdir <- "manuscript/figures"; dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

box <- function(x, y, w, h, fill, lab, size = 3.0, face = "plain", col = "#1A1D1E")
  list(rect = data.frame(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2, fill = fill),
       txt  = data.frame(x = x, y = y, lab = lab, size = size, face = face, col = col))
N <- list(
  # shared transportable test-accuracy parameters (top band)
  box(3.10, 8.4, 4.3, 1.15, "#DCEBEC", "Se_BC(V) = inv.logit( α0 + α1·log V + τ·u )\nvolume-anchored  ·  α1 ~ N(0.36, 0.15)  [Antillón 2018]", 3.0, "bold"),
  box(6.55, 8.4, 2.0, 1.15, "#DCEBEC", "Se_BM ~ Beta(12,2)\n≈ 0.90", 3.0, "bold"),
  box(8.95, 8.4, 2.3, 1.15, "#DCEBEC", "Sp ~ Beta(200,1)\n≈ 0.995  (pinned)", 3.0, "bold"),
  # LEFT layer : paired culture (identifies Se, phi)
  box(2.7, 5.7, 4.4, 1.0, "#E6F0EA", "LATENT: true typhoid  ·  hospital PPV  φ_s   [LOCAL]", 3.0, "bold", "#1B7A4B"),
  box(2.7, 3.7, 4.4, 1.3, "#FFFFFF", "blood-culture ±  &  bone-marrow ±\n2×2 cell counts ~ Multinomial\n(Hui–Walter conditional-independence LCM)", 2.8),
  box(2.7, 1.35, 4.4, 0.95, "#F3EFE6", "DATA: 10 paired blood + bone-marrow\nculture studies  (Mogasale 2016; hospital)", 2.8, "italic"),
  # RIGHT layer : outbreak positivity (identifies pi)
  box(8.1, 5.7, 4.6, 1.0, "#EEEAF4", "LATENT: true typhoid  ·  community PPV  π_o   [LOCAL]", 3.0, "bold", "#5A4B9E"),
  box(8.1, 3.7, 4.6, 1.3, "#FFFFFF", "positives k_o among n_o tested\n~ Binomial( n_o ,  π_o·Se_BC(V) + (1−π_o)(1−Sp) )", 2.8),
  box(8.1, 1.35, 4.6, 0.95, "#F3EFE6", "DATA: outbreak suspected → confirmed\n(Aye 2004, Neil 2012, Kabwama 2017)", 2.8, "italic")
)
rects <- do.call(rbind, lapply(N, `[[`, "rect"))
txts  <- do.call(rbind, lapply(N, `[[`, "txt"))
edge <- function(x1, y1, x2, y2, col = "#5B6566") data.frame(x = x1, y = y1, xend = x2, yend = y2, col = col)
E <- rbind(
  edge(3.1, 7.82, 2.9, 4.36), edge(6.55, 7.82, 3.4, 4.36), edge(8.95, 7.82, 3.9, 4.36),  # params -> left lik
  edge(3.1, 7.82, 7.6, 4.36), edge(8.95, 7.82, 8.5, 4.36),                                # Se_BC,Sp -> right lik
  edge(2.7, 5.2, 2.7, 4.36, "#1B7A4B"), edge(8.1, 5.2, 8.1, 4.36, "#5A4B9E"),             # latent -> lik
  edge(2.7, 1.83, 2.7, 3.05, "#B08A3E"), edge(8.1, 1.83, 8.1, 3.05, "#B08A3E")            # data -> lik
)
p <- ggplot() +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
            colour = "#C9C4B8", linewidth = .4) +
  geom_segment(data = E, aes(x = x, y = y, xend = xend, yend = yend, colour = I(col)),
               arrow = arrow(length = unit(.14, "cm"), type = "closed"), linewidth = .45) +
  geom_text(data = txts, aes(x, y, label = lab, fontface = face, size = I(size), colour = I(col)), lineheight = .92) +
  annotate("text", x = 3.1, y = 9.35, label = "TRANSPORTABLE  —  test accuracy (reused across settings)",
           fontface = "bold", size = 3.2, colour = "#0E6E78") +
  annotate("segment", x = 0.5, xend = 10.3, y = 9.05, yend = 9.05, colour = "#0E6E78", linewidth = .3) +
  annotate("text", x = 2.7, y = 6.45, label = "—— paired-culture layer ——", size = 2.9, colour = "#1B7A4B", fontface = "bold") +
  annotate("text", x = 8.1, y = 6.45, label = "—— outbreak-positivity layer ——", size = 2.9, colour = "#5A4B9E", fontface = "bold") +
  annotate("text", x = 5.5, y = 0.35,
           label = "LOCAL population parameters: hospital φ_s  vs  community π_o  (the PPV of the clinical case definition, by setting)",
           size = 2.9, fontface = "italic", colour = "#5B6566") +
  scale_x_continuous(limits = c(0.3, 10.5)) + scale_y_continuous(limits = c(0, 9.7)) +
  labs(title = "Latent-class model for typhoid diagnostic accuracy and case-definition PPV",
       subtitle = "model_final.stan · two data layers share transportable test parameters; each identifies a local PPV") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 10.5, colour = "#5B6566", hjust = 0),
        plot.margin = margin(14, 14, 10, 14))
ggsave(file.path(figdir, "fig0_latent_class_model.png"), p, width = 10.5, height = 6.6, dpi = 135)
cat("wrote", file.path(figdir, "fig0_latent_class_model.png"), "\n")

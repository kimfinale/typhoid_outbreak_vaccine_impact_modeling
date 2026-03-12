library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readxl)
library(patchwork)
library(gridExtra)
library(ggpattern)
library(ggtext)

svim_fin <- fread("weekly_output_n200_ive_fin.csv")
svim_month_fin <- fread("monthly_output_n200_fin_ive.csv")

my_discrete_colors <-
  c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
    "#FF7F00", "gold1", "skyblue2", "palegreen2", "#FDBF6F",
    "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
format_num <- function(x, digits=0) {
  format(round(x, digits=digits), big.mark=",", trim=TRUE)
}

# Figure 2
ts_single_outbreak <- ts1 %>% 
  filter(study_id == "Polonsky 2014, Dzivaresekwa") %>%
  select(study_id, suspected_cases) %>%
  mutate(week = row_number())

pct_reduc <- svim_fin %>% 
  filter(study_id == "Polonsky 2014, Dzivaresekwa") %>%
  filter(vacc_cov == 0.85, vacc_week == 8) %>%
  pull(pct_reduc_case_median)

single_outbreak <- ts_single_outbreak %>%
  mutate(cases = ifelse(week < 8, suspected_cases, suspected_cases * pct_reduc / 100))

ggplot(single_outbreak, aes(x = week, y = suspected_cases, group = 1)) +
  geom_rect_pattern(aes(xmin=week, xmax=lead(week), ymin=0, ymax=suspected_cases),
                    pattern = 'stripe',
                    fill            = NA,
                    colour          = 'black', 
                    pattern_density = 0.35, 
                    pattern_fill    = 'steelblue',
                    pattern_colour  = NA) + 
  geom_rect(
    aes(xmin=week, xmax=lead(week), ymin=0, ymax=counterfactual_cases_week8), 
    fill="grey", alpha=0.8, inherit.aes = FALSE) +
  geom_segment(data = single_outbreak,
               aes(x = week, xend = lead(week), 
                   y = counterfactual_cases_week8, yend = counterfactual_cases_week8),
               linetype = "dashed", colour = "black") +
  geom_step()+
  # 
  # geom_ribbon(aes(ymin = suspected_cases, ymax = counterfactual_cases_week8), 
  #             fill = "red", alpha = 0.3) +
  # geom_ribbon(aes(ymin = 0, ymax = counterfactual_cases_week8), 
  #             fill = "steelblue", alpha = 0.3) +
  ylim(0, 170) +
  scale_x_continuous(breaks = seq(1, 24, by = 1), labels = as.character(seq(1, 24, by = 1))) +
  
  # geom_segment(aes(x=2,y=25,xend=3,yend=11),
  #              arrow=arrow(angle=10, length=unit(3,"mm")),
  #              inherit.aes = FALSE)+
  annotate("text", x=16.75, y=109,
           label= "% Case Reduction = 1 - ")+
  geom_segment(aes(x=19,
                   xend=19.8,
                   y=108.5, yend=108.5), inherit.aes = FALSE)+
  geom_rect(aes(xmin=19.15, xmax=19.65, 
                ymin=110, ymax=117), color='black', fill="grey", alpha=0.5)+
  geom_rect_pattern(
    aes(xmin=19.15,xmax=19.65, 
        ymin=100, ymax=107),
    pattern = 'stripe',
    fill            = NA,
    colour          = 'black', 
    pattern_density = 0.35, 
    pattern_fill    = 'steelblue',
    pattern_colour  = NA)+
  annotate("text", x=21, y=109,
           label= "x 100%")+
  # geom_segment(aes(x=17.35,
  #                  xend=18.15,
  #                  y=108.5, yend=108.5), inherit.aes = FALSE)+
 annotate("text", x=2.7, y=139, label= "Outbreak detected", hjust=0.8, angle = 90, size = 3.35) +
  # geom_segment(aes(x=4,y=62,xend=5,yend=48),
  #              arrow=arrow(angle=10, length=unit(3,"mm")),
  #              inherit.aes = FALSE)+
  # geom_segment(aes(x=8,y=120,xend=8,yend=79),
  #              arrow=arrow(angle=10, length=unit(3,"mm")),
  #              inherit.aes = FALSE)+
  geom_vline(xintercept = 3, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 8, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 12, linetype = "dashed", alpha = 0.5) +
  annotate("text", x=7.7, y=138.5, label= "Vaccination", hjust=0.8, angle = 90, size = 3.35) +
  
  annotate("text", x=11.7, y=138.5, label= "Vaccine takes effect", hjust=0.8, angle = 90, size = 3.35) +
  # geom_segment(aes(x = 3.05, y = 150, xend = 7.95, yend = 150),
  #              arrow = arrow(ends = "both", length = unit(1, "mm")), size = 0.2,
  #              inherit.aes = FALSE) +
  annotate("text", x = 5.5, y = 159, size = 3.5, label = "Outbreak\nresponse delay", hjust = 0.5) +
  
  annotate("text", x = 10, y = 162, size = 3.5, label = "Vaccine\nimmunogenicity\nrandom delay", hjust = 0.5) +
  # geom_segment(aes(x = 8.05, y = 150, xend = 23.95, yend = 150),
  #              arrow = arrow(ends = "both", length = unit(1, "mm")), size = 0.2,
  #              inherit.aes = FALSE) +
  annotate("text", x = 18, y = 157, size = 3.8, label = "Vaccine impact", hjust = 0.5) +
  # geom_segment(aes(x = 1, y = 150, xend = 2.95, yend = 150),
  #              arrow = arrow(ends = "both", length = unit(1, "mm")), size = 0.2,
  #              inherit.aes = FALSE) +
  annotate("segment", x = 17.35, xend = 18.15, y = 108.5, yend = 108.5) +
  annotate("segment", x = 3.05, xend = 7.95, y = 150, yend = 150,
           arrow = arrow(ends = "both", length = unit(1, "mm")), linewidth = 0.2) +
  annotate("segment", x = 8.05, xend = 23.95, y = 150, yend = 150,
           arrow = arrow(ends = "both", length = unit(1, "mm")), linewidth = 0.2) +
  annotate("segment", x = 1, xend = 2.95, y = 150, yend = 150,
           arrow = arrow(ends = "both", length = unit(1, "mm")), linewidth = 0.2) +

  annotate("text", x = 1.9, y = 159, size = 3.5, label = "Outbreak\nstart", hjust = 0.5) +
  # 
  # geom_line(aes(x = vacc_week, y = counterfactual_cases_week8), linetype = "dashed", color = "red") +
  # geom_segment(aes(x = 12, y = 32.5, xend = 12, yend = 48),
  #              arrow = arrow(ends = "both", length = unit(1.5, "mm")), linetype = "dashed", size = 0.3,
  #              inherit.aes = FALSE, color = "red") +
  # geom_segment(aes(x = 12.5, y = 39, xend = 14, yend = 55),
  #              arrow = arrow(length = unit(1, "mm")), size = 0.2,
  #              inherit.aes = FALSE) +
  # annotate("text", x = 14, y = 61, size = 3.35, label = "Averted cases due to\nvaccination", 
  #          hjust = 0.5) +
  
labs(x = "", y = "Number of Cases", title = "Hypothetical counterfactual outbreak scenario") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.945, size = 11.5, vjust = -8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) -> p1


p1

single_outbreak1 <- ts1 %>% filter(study_id == "Nimonkar 2022")
single_outbreak2 <- ts1 %>% filter(study_id == "Aye 2004")
single_outbreak3 <- ts1 %>% filter(study_id == "Davis 2018")
single_outbreak4 <- ts1 %>% filter(study_id == "Muti 2014")

ggplot() +
  geom_rect_pattern(data = single_outbreak3,
                    aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = suspected_cases),
                    pattern = 'stripe',
                    fill = 'white',
                    colour = NA,
                    pattern_density = 0.35,
                    pattern_fill = 'thistle',
                    pattern_colour = NA,
                    alpha = 0.8) +
  
  geom_rect(data = single_outbreak3,
            aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = cases), 
            fill = "darkorchid", alpha = 0.65, inherit.aes = FALSE) +
  
  geom_segment(data = filter(single_outbreak3, vacc_week %in% c(11:24)),
               aes(x = vacc_week, xend = vacc_week,
                   y = lag(suspected_cases), yend = suspected_cases),
               colour = "thistle") +
  geom_segment(data = filter(single_outbreak3, vacc_week %in% c(11:24)),
               aes(x = vacc_week, xend = lead(vacc_week),
                   y = suspected_cases, yend = suspected_cases),
               colour = "thistle") +
  
  geom_segment(data = filter(single_outbreak3, vacc_week %in% c(11:24)),
               aes(x = lead(vacc_week), xend = lead(vacc_week),
                   y = suspected_cases, yend = lead(suspected_cases)),
               colour = "thistle", size = 0.5) +
  
  geom_segment(data = filter(single_outbreak3, vacc_week == 24),
               aes(x = lead(vacc_week), xend = lead(vacc_week),
                   y = 0, yend = suspected_cases),
               colour = "thistle", size = 0.5) +
  
  geom_rect_pattern(data = single_outbreak4,
                    aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = suspected_cases),
                    pattern = 'stripe',
                    fill = 'white',
                    colour = NA,
                    pattern_density = 0.35,
                    pattern_fill = 'darkseagreen2',
                    pattern_colour = NA,
                    alpha = 1, size = 0.5) + 
  
  geom_rect(data = single_outbreak4,
            aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = cases), 
            fill = "seagreen4", alpha = 0.65, inherit.aes = FALSE) +
  
  geom_segment(data = filter(single_outbreak4, vacc_week %in% c(11:24)),
               aes(x = vacc_week, xend = vacc_week,
                   y = lag(suspected_cases), yend = suspected_cases),
               colour = "darkseagreen2") +
  geom_segment(data = filter(single_outbreak4, vacc_week %in% c(11:24)),
               aes(x = vacc_week, xend = lead(vacc_week),
                   y = suspected_cases, yend = suspected_cases),
               colour = "darkseagreen2") +
  
  geom_segment(data = filter(single_outbreak4, vacc_week %in% c(11:24)),
               aes(x = lead(vacc_week), xend = lead(vacc_week),
                   y = suspected_cases, yend = lead(suspected_cases)),
               colour = "darkseagreen2", size = 0.5) +
  
  geom_segment(data = filter(single_outbreak4, vacc_week == 24),
               aes(x = lead(vacc_week), xend = lead(vacc_week),
                   y = 0, yend = suspected_cases),
               colour = "darkseagreen2", size = 0.5) +
  geom_rect(data = single_outbreak1,
            aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = cases), 
            fill = "sienna3", alpha = 0.85, inherit.aes = FALSE) +
  
  geom_rect(data = single_outbreak2,
            aes(xmin = vacc_week, xmax = lead(vacc_week), ymin = 0, ymax = cases), 
            fill = "skyblue3", alpha = 0.85, inherit.aes = FALSE) +
  
  geom_vline(xintercept = c(8, 11), linetype = "dashed", alpha = 0.5) +
  
  scale_x_continuous(breaks = seq(1, 24, by = 1), labels = as.character(seq(1, 24, by = 1))) +
  
  geom_segment(aes(x = 7.65, y = 2.2, xend = 4.8, yend = 50),
               arrow = arrow(length = unit(1, "mm")), size = 0.25, inherit.aes = FALSE) +
  annotate("text", x = 4, y = 62, size = 2.75, label = "Outbreak ends before\neffective immunity\n(no averted cases)", hjust = 0.5) +
  
  geom_segment(aes(x = 2.5, y = 2.5, xend = 2, yend = 25),
               arrow = arrow(length = unit(1, "mm")), size = 0.25, inherit.aes = FALSE) +
  annotate("text", x = 1.8, y = 37, size = 2.75, label = "Outbreak ends\nbefore vaccination\n(no averted cases)", hjust = 0.5) +
  
  geom_segment(aes(x = 13.5, y = 38.25, xend = 13.5, yend = 85),
               arrow = arrow(ends = "both", length = unit(1, "mm")), size = 0.5, inherit.aes = FALSE, color = "red") +
  geom_segment(aes(x = 19.5, y = 27.9, xend = 19.5, yend = 62),
               arrow = arrow(ends = "both", length = unit(1, "mm")), size = 0.5, inherit.aes = FALSE, color = "red") +
  geom_segment(aes(x = 19.4, y = 45, xend = 18, yend = 75),
               arrow = arrow(length = unit(1, "mm")), size = 0.25, inherit.aes = FALSE) +
  geom_segment(aes(x = 13.7, y = 55, xend = 16, yend = 75),
               arrow = arrow(length = unit(1, "mm")), size = 0.25, inherit.aes = FALSE) +
  annotate("text", x = 17, y = 80, size = 2.75, label = "Averted cases\ndue to ORI", hjust = 0.5) +
  
  labs(x = "Week", y = "Number of Cases", title = "Range of individual outbreak scenarios") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.945, vjust = -9, size = 11.5),
    legend.title = element_blank()
  ) -> p2

p2

p1 <- p1 + theme(plot.margin = margin(0, -14, 0, -1))  
p2 <- p2 + theme(plot.margin = margin(-15, 10, 0, 0))  

grid <- grid.arrange(p1, p2, nrow = 2, heights = c(1.5,1))

ggsave("fig2.jpg", grid, width = 10, height = 10, dpi = 300)

# Figure 3

# Week 7 markers (start + 6*7 days) for intervention timing
week7_lines <- ts1 %>%
  group_by(study_id) %>%
  summarise(
    week7_date = min(date_mid) + days(6 * 7),
    .groups = "drop"
  )

ts1$intervention_timing2 <- as.integer(round(ts1$intervention_timing))

week_lines <- ts1 %>%
  group_by(id_outbreak) %>%
  summarise(
    intervention_date = min(date_mid) + days(7 * intervention_timing2[1]),
    .groups = "drop"
  )
week_lines <- week_lines %>%
  mutate(id_outbreak = factor(id_outbreak, levels = outbreak_start_order))

chars_dat$location[chars_dat$location == "MalawiÐMozambique"] <- "Malawi-Mozambique"
# chars_dat$location[chars_dat$location == "Amritsar, India."] <- "Amritsar, India"

# Label with location, add * if AMR positive
labels_df <- chars_dat %>%
  mutate(
    label = ifelse(is.na(amr_positive), location, 
                   ifelse(amr_positive == 1, paste0(location, "*"), location)),
    .groups = "drop"
  )
labels_df <- labels_df %>%
  mutate(
    label = ifelse(is.na(amr_positive), location,
                   ifelse(amr_positive == 1,
                          paste0(location, "<span style='font-size:14pt'>*</span>"),
                          location))
  )

# Add AMR status as a factor with labels
ts1 <- ts1 %>%
  mutate(
    amr_status = factor(ifelse(is.na(amr_positive),"AMR-Negative",
                               ifelse(amr_positive == 1, "AMR-Positive", "AMR-Negative"))
    ))

xlims <- as.Date(c("1997-01-01", "2020-01-31"))

plt_a <- ggplot(ts1, aes(x = date_mid, y = suspected_cases)) +
  geom_col(aes(fill = who_region), width = 25) +
  facet_wrap(~id_outbreak, ncol = 1, strip.position = "right", scales = "fixed") +
  geom_vline(
    data = week_lines,
    aes(xintercept = as.numeric(intervention_date), linetype = "Intervention Timing"),
    color = "black",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_richtext(
    data = labels_df,
    aes(x = as.Date("1997-01-01"), y = 10000, label = label),
    hjust = 0, vjust = 1,
    size = 3.5,       
    fill = NA, label.color = NA,  # removes background and box
    inherit.aes = FALSE
  ) +
  scale_y_log10() +
  scale_x_date(
    limits = xlims,
    breaks = seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by = "2 years"),
    date_labels = "%Y",
    expand = c(0.01, 0)
  ) +
  geom_point(
    aes(shape = amr_status),
    size = 0,     
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = c("red3", "blue2", "green4")) +
  scale_linetype_manual(
    name = NULL,
    values = c("Intervention Timing" = "dashed")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("AMR-Positive" = 8),  # 8 = star symbol
    labels = c("AMR-Negative" = "", "AMR-Positive" = "AMR-Positive")
  ) +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 3, override.aes = list(size = 2, alpha = 1))
  ) +
  labs(
    x = "Time",
    y = "Monthly Suspected Cases (log)",
    fill = "WHO Region"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text.y.right = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(0.5, "mm"),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.caption = element_text(hjust = 0, size = 10),
    legend.spacing.y = unit(-2, "mm"), 
    axis.text.x = element_text(size = 10)
  )
plt_a


plt_b <- ggplot(chars_dat, aes(x = id_outbreak, y = duration_weeks)) +
  geom_col(fill = "grey20") +
  coord_flip() +
  scale_x_discrete(breaks = NULL, labels = NULL) +
  labs(y = "Duration (weeks)", x = "") +
  theme_minimal(base_size = 9) +
  theme(plot.margin = margin(0, 0, 0,0),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7))
plt_b

# boxplot within each barplot - decided to remove in final figure

# duration_box_grob <- ggplotGrob(
#   ggplot(chars_dat, aes(y = duration_weeks)) +
#     geom_boxplot(width = 0.4, fill = "grey80", outlier.shape = NA) +
#     labs(y = "Duration (weeks)") +
#     theme_minimal(base_size = 9) +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.title = element_text(size = 9, hjust = 0.5),
#       plot.margin = margin(5, 5, 5, 5)
#     )
# )
# 
# plt_b2 <- plt_b +
#   annotation_custom(
#     grob = duration_box_grob,
#     xmin = 8, xmax = 13,
#     ymin = 50,
#     ymax = 150
#   )
# plt_b2

plt_c <- ggplot(chars_dat, aes(x = id_outbreak, y = cum_cases)) +
  geom_col(fill = "grey20") +
  coord_flip() +
  scale_x_discrete(breaks = NULL, labels = NULL) +
  labs(y = "Cumulative cases", x = "") +
  theme_minimal(base_size = 9) +
  scale_y_continuous(labels = function(x) {
    ifelse(x == 0, "0", paste0(x / 1000, "k"))
  }) +
  theme(plot.margin = margin(0, 0, 0, 0), 
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 7))

plt_c

# cum_cases_box_grob <- ggplotGrob(
#   ggplot(chars_dat, aes(y = cum_cases)) +
#     geom_boxplot(width = 0.4, fill = "grey80", outlier.shape = NA) +
#     labs(y = "Cumulative cases") +
#     theme_minimal(base_size = 9) +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.title = element_text(size = 9, hjust = 0.5),
#       plot.margin = margin(5, 5, 5, 5)
#     )
# )
# 
# plt_c2 <- plt_c +
#   annotation_custom(
#     grob = cum_cases_box_grob,
#     xmin = 8, xmax = 13,
#     ymin = 3000,
#     ymax = 10000
#   )

plt_d <- ggplot(chars_dat, aes(x = id_outbreak, y = attack_rate)) +
  geom_col(fill = "grey20") +
  coord_flip() +
  scale_x_discrete(breaks = NULL, labels = NULL) +
  labs(y = "Attack rate", x = "") +
  theme_minimal(base_size = 9) +
  theme(plot.margin = margin(0, 0, 0, 0),
        
        axis.text.x = element_text(size =7),
        axis.title.x = element_text(size = 7))
plt_d
# 
# ar_box_grob <- ggplotGrob(
#   ggplot(chars_dat, aes(y = attack_rate)) +
#     geom_boxplot(width = 0.4, fill = "grey80", outlier.shape = NA) +
#     labs(y = "Attack rate") +
#     theme_minimal(base_size = 9) +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       plot.title = element_text(size = 9, hjust = 0.5),
#       plot.margin = margin(5, 5, 5, 5)
#     )
# )

# plt_d2 <- plt_d +
#   annotation_custom(
#     grob = ar_box_grob,
#     xmin = 8, xmax = 13,
#     ymin = 2,
#     ymax = 7.2
#   )
# 
# plt_b2 <- plt_b2 + theme(legend.position = "none")
# plt_c2 <- plt_c2 + theme(legend.position = "none")
# plt_d2 <- plt_d2 + theme(legend.position = "none")

plt_combined <- (plt_a | plt_b2 | plt_c2 | plt_d2) +
  plot_layout(ncol = 4, widths = c(1.3, 0.3, 0.3, 0.3))

plt_combined
ggsave("figure3.jpg", plt_combined, width = 11.5, height = 6, dpi = 500)


# Figure 4
df <- svim_fin %>% filter(vacc_cov %in% c(0.65, 0.80, 0.90), vacc_week %in% c(1, 2, 4, 6 ,8 ,10, 12))

# df <- svim_month_fin %>% filter(vacc_cov %in% c(0.65, 0.80, 0.90), vacc_week %in% c(4, 8, 12, 16, 20,24))

ggplot(df, aes(x = as.factor(vacc_week), y = pct_reduc_case_median, color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 85)) +
  scale_color_manual(
    name = "Vaccine coverage",
    values = my_discrete_colors[1:3],
    labels = c("65%", "80%", "90%")
  ) +
  labs(
    x = "",
    y = "Case Reduction Percent",
    title = ""
  ) +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    legend.direction = "horizontal"
  ) -> p1

p1

ggplot(df, aes(x = as.factor(vacc_week), y = s_ch_averted_median +1,
               color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  # geom_point(aes(color = type),
  #          position = position_jitter(width = 0.2, height = 0),
  #          size = 1.8, alpha = 0.7) +
  
  # scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 80)) +
  # 
  # ylim(0, 100) +
  scale_y_log10() +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  labs(
    x = "Week of vaccination",
    y = "Case Reduction Number (log)"
  ) +
  
  theme_light() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> p2
p2


dt_filtered2 <- df %>% filter(ori_occurred == TRUE)
ggplot(dt_filtered2, aes(x = as.factor(vacc_week), y = case_averted_per_1000_median,
                         color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  scale_y_log10() +
  labs(
    x = "Delay in ORI (weeks)",
    y = "Cases Averted per 1,000 Doses",
    title = ""
  ) +
  
  theme_light() +
  theme(
    legend.position = "none"
  ) -> p3
p3


ggplot(df, aes(x = as.factor(vacc_week), y = cost_per_daly_averted_median,
               color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  labs(
    x = "Delay in ORI (weeks)",
    y = "Cost per DALY Averted (Log)"
  ) +
  scale_y_log10() +
  theme_light() +
  theme(
    legend.position = "none"
  ) -> p4
p4

legend <- cowplot::get_legend(p1)

grid <- grid.arrange(
  legend,            # Legend on top
  arrangeGrob(p1 + theme(legend.position = "none", plot.margin = margin(-25, 5, 18, 5)), 
              p2 + theme(legend.position = "none", plot.margin = margin(-7, 5, 18, 7)), 
              p3 + theme(legend.position = "none", plot.margin = margin(-33, 5, 9, 5.5)), 
              p4 + theme(legend.position = "none", plot.margin = margin(-16, 5, 8, 1.5)), 
              nrow = 2),
  ncol = 1,
  heights = c(0.1, 1))


ggsave("figure4.jpg", grid, width = 10, height = 8, dpi = 300)
# ggsave("fig4-monthly.jpg", grid, width = 10, height = 8, dpi = 300)


# Figure 5
chars2 <- dat %>% 
  filter(study_id %in% ts_weekly$study_id) %>%
  select(study_id, amr_positive, cum_cases, duration_weeks)

size_thresh <- quantile(chars2$cum_cases, 0.75, na.rm = TRUE)
duration_thresh <- quantile(chars2$duration_weeks, 0.75, na.rm = TRUE)
# 
# size <- chars2 %>% filter(cum_cases >= size_thresh)
# dur <- chars2 %>% filter(duration_weeks >= duration_thresh)
# 
# summary(size)
# summary(dur)

df <- left_join(svim_fin, chars2, by = c("study_id"))
# cor(chars2$duration_weeks, chars2$cum_cases)

df_cat <- df %>%
  mutate(
    category_all = "All",
    category_size = ifelse(cum_cases >= size_thresh, ">75th Size", NA),
    category_duration = ifelse(duration_weeks >= duration_thresh, ">75th Duration", NA),
    category_amr = amr_positive
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("category_"),
    names_to = "category_type",
    values_to = "category"
  ) %>%
  filter(!is.na(category))

df_cat <- df_cat %>% filter(vacc_cov == 0.8, vacc_week == 8)

# remove AMR-negative outbreaks due to no impact
df_cat2 <- df_cat %>% filter(category != "Negative")

df_cat2 <- df_cat2 %>%
  mutate(category = factor(category,
                           levels = c("All", ">75th Duration", ">75th Size", "Positive", "Not tested")
  ))

category_colors <- c(
  "All" = "lightblue",
  ">75th Duration" = "#A4A4D5",
  ">75th Size" = "#DDB5D5",
  "Positive" = "#FFE0B2",
  "Not tested" = "#E65100"
)

common_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

df_cat3 <- df_cat2 %>%
  mutate(
    category2 = ifelse(category == "Positive", "AMR-Positive",
                       ifelse(category == "Not tested", "AMR-Not tested", category)),
    category2 = factor(category2,
                       levels = c("All", ">75th Duration", ">75th Size",
                                  "AMR-Positive", "AMR-Not tested"))
  )

p1 <- ggplot(df_cat3, aes(x = category2, y = s_ch_averted_median)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(color = category), width = 0.2, size = 2, alpha = 0.7) +
  
  scale_color_manual(values = category_colors) +
  labs(y = "Cases Averted", x = "") +
  common_theme
p1

p2 <- ggplot(df_cat3, aes(x = category2, y = pct_reduc_case_median)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(color = category), width = 0.2, size = 2, alpha = 0.7) +
  
  scale_color_manual(values = category_colors) +
  #ylim(0, 100) +
  labs(y = "% Cases Averted", x = "") +
  common_theme

p2

p3 <- ggplot(df_cat3, aes(x = category2, y = cost_per_daly_averted_median)) +
  
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(color = category), width = 0.2, size = 2, alpha = 0.7) +
  
  scale_color_manual(values = category_colors) +
  #ylim(0, 700000) +
  labs(y = "Cost per DALY averted", x = "") +
  common_theme
p3

df_cat3_filter <- df_cat3 %>% filter(ori_occurred == TRUE)
p4 <- ggplot(df_cat3_filter, aes(x = category2, y = case_averted_per_1000_median)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(color = category), width = 0.2, size = 2, alpha = 0.7) +
  
  scale_color_manual(values = category_colors) +
  ylim(0, 10.5) +
  labs(y = "Cases Averted per 1,000 Doses", x = "") +
  common_theme
p4

(p1 | p2) / (p3 | p4) ->p
p
ggsave("figure4.jpg", p, height = 6, width = 8)


# Supplementary Figure 1
df <- svim_month_fin %>% filter(vacc_cov %in% c(0.65, 0.80, 0.90), vacc_week %in% c(4, 8, 12, 16, 20,24))

ggplot(df, aes(x = as.factor(vacc_week), y = pct_reduc_case_median, color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 85)) +
  scale_color_manual(
    name = "Vaccine coverage",
    values = my_discrete_colors[1:3],
    labels = c("65%", "80%", "90%")
  ) +
  labs(
    x = "",
    y = "Case Reduction Percent",
    title = ""
  ) +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "right",
    legend.direction = "horizontal"
  ) -> p1

p1

ggplot(df, aes(x = as.factor(vacc_week), y = s_ch_averted_median +1,
               color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  # geom_point(aes(color = type),
  #          position = position_jitter(width = 0.2, height = 0),
  #          size = 1.8, alpha = 0.7) +
  
  # scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 80)) +
  # 
  # ylim(0, 100) +
  scale_y_log10() +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  labs(
    x = "Week of vaccination",
    y = "Case Reduction Number (log)"
  ) +
  
  theme_light() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> p2
p2


dt_filtered2 <- df %>% filter(ori_occurred == TRUE)
ggplot(dt_filtered2, aes(x = as.factor(vacc_week), y = case_averted_per_1000_median,
                         color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  scale_y_log10() +
  labs(
    x = "Delay in ORI (weeks)",
    y = "Cases Averted per 1,000 Doses",
    title = ""
  ) +
  
  theme_light() +
  theme(
    legend.position = "none"
  ) -> p3
p3


ggplot(df, aes(x = as.factor(vacc_week), y = cost_per_daly_averted_median,
               color = as.factor(vacc_cov))) +
  geom_boxplot(width = 0.6, outlier.size = 0.6) +
  scale_color_manual(
    name = "Vaccine coverage %",
    values = my_discrete_colors[1:3]
  ) +
  labs(
    x = "Delay in ORI (weeks)",
    y = "Cost per DALY Averted (Log)"
  ) +
  scale_y_log10() +
  theme_light() +
  theme(
    legend.position = "none"
  ) -> p4
p4

legend <- cowplot::get_legend(p1)

grid <- grid.arrange(
  legend,       
  arrangeGrob(p1 + theme(legend.position = "none", plot.margin = margin(-25, 5, 18, 5)), 
              p2 + theme(legend.position = "none", plot.margin = margin(-7, 5, 18, 7)), 
              p3 + theme(legend.position = "none", plot.margin = margin(-33, 5, 9, 5.5)), 
              p4 + theme(legend.position = "none", plot.margin = margin(-16, 5, 8, 1.5)), 
              nrow = 2),
  ncol = 1,
  heights = c(0.1, 1))


ggsave("supp-fig1.jpg", grid, width = 10, height = 8, dpi = 300)

# Supplementary Figure 2
d_ive <- svim_fin %>%
  #filter(ori_occurred == TRUE) %>%
  group_by(ive) %>%
  summarise(cases_averted_median = median(s_ch_averted_median),
            cases_averted_lb = quantile(s_ch_averted_median, probs = 0.025),
            cases_averted_ub = quantile(s_ch_averted_median, probs = 0.975), 
            pct_cases_averted_median = median(pct_reduc_case_median),
            pct_cases_averted_lb = quantile(pct_reduc_case_median, probs = 0.025),
            pct_cases_averted_ub = quantile(pct_reduc_case_median, probs = 0.975),
            cost_per_daly_averted_median2 = median(cost_per_daly_averted_median, na.rm = TRUE),
            cost_per_daly_averted_lb = quantile(cost_per_daly_averted_median, probs = 0.025, na.rm = TRUE),
            cost_per_daly_averted_ub = quantile(cost_per_daly_averted_median, probs = 0.975, na.rm = TRUE),
            case_averted_per_1000_median2 = median(case_averted_per_1000_median, na.rm = TRUE),
            case_averted_per_1000_median_lb = quantile(case_averted_per_1000_median, probs = 0.025, na.rm = TRUE),
            case_averted_per_1000_median_ub = quantile(case_averted_per_1000_median, probs = 0.975, na.rm = TRUE),
  )

ggplot(d_ive, aes(x = ive, y = pct_cases_averted_median)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = pct_cases_averted_lb, ymax = pct_cases_averted_ub), 
                width = 0.02, color = "steelblue") +
  
  # highlight point with legend
  geom_point(aes(x = 0.17, y = 30.0, shape = "Main model", color = "Main model"), size = 4) +
  geom_errorbar(aes(x = 0.17, ymin = 0, ymax = 68.5, color = "Main model"), width = 0.02) +
  
  scale_shape_manual(name = "", values = c("Main model" = 17)) +
  scale_color_manual(name = "", values = c("Main model" = "red")) +
  ylim(0, 85) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) + 
  labs(
    x = "Indirect Vaccine Effectiveness (IVE)",
    y = "% Cases Averted", subtitle = "Outbreaks reporting weekly incidence"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = c(0.05, 0.80),
    legend.justification = c(0.05, -0.5) 
  ) -> p1
p1

d_ive_month <- svim_month_fin %>%
  #filter(ori_occurred == TRUE) %>%
  group_by(ive) %>%
  summarise(cases_averted_median = median(s_ch_averted_median),
            cases_averted_lb = quantile(s_ch_averted_median, probs = 0.025),
            cases_averted_ub = quantile(s_ch_averted_median, probs = 0.975), 
            pct_cases_averted_median = median(pct_reduc_case_median),
            pct_cases_averted_lb = quantile(pct_reduc_case_median, probs = 0.025),
            pct_cases_averted_ub = quantile(pct_reduc_case_median, probs = 0.975),
            cost_per_daly_averted_median2 = median(cost_per_daly_averted_median, na.rm = TRUE),
            cost_per_daly_averted_lb = quantile(cost_per_daly_averted_median, probs = 0.025, na.rm = TRUE),
            cost_per_daly_averted_ub = quantile(cost_per_daly_averted_median, probs = 0.975, na.rm = TRUE),
            case_averted_per_1000_median2 = median(case_averted_per_1000_median, na.rm = TRUE),
            case_averted_per_1000_median_lb = quantile(case_averted_per_1000_median, probs = 0.025, na.rm = TRUE),
            case_averted_per_1000_median_ub = quantile(case_averted_per_1000_median, probs = 0.975, na.rm = TRUE),
  )

ggplot(d_ive_month, aes(x = ive, y = pct_cases_averted_median)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = pct_cases_averted_lb, ymax = pct_cases_averted_ub), 
                width = 0.02, color = "steelblue") +
  
  geom_point(aes(x = 0.17, y = 69, shape = "Main model", color = "Main model"), size = 4) +
  geom_errorbar(aes(x = 0.17, ymin = 35, ymax = 71.2, color = "Main model"), width = 0.02) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) + 
  scale_shape_manual(name = "", values = c("Main model" = 17)) +
  scale_color_manual(name = "", values = c("Main model" = "red")) +
  ylim(0, 85) +
  labs(
    x = "Indirect Vaccine Effectiveness (IVE)",
    y = NULL, subtitle = "Outbreaks reporting monthly incidence"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = c(0.05, 0.8), 
    legend.justification = c(0.05, -0.5) 
  ) -> p2
p2

p <- p1 + p2
p
ggsave("supp3.jpg", p, height = 5, width = 8)


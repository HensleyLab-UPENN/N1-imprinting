# Load packages 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)

# Load NAI results
h1n1_df <- read_csv("Ort-et-al_N1-imprinting/Data/H1N1_Titers.csv")
h3n2_df <- read_csv("Ort-et-al_N1-imprinting/Data/H3N2_Titers.csv")


# For H1N1, determine acute and convalescent titers and add extra columns for plotting purposes
h1n1_df_acute <- h1n1_df[h1n1_df$Visit==1 & h1n1_df$DPSO < 7,]
h1n1_df_acute$Timepoint = "Acute"
h1n1_df_acute <- h1n1_df_acute[, c("Participant", "Timepoint", "CA09", "WI22", "DC24", "BC24", "TL22")]

h1n1_df_conv <- h1n1_df[h1n1_df$Participant %in% h1n1_df_acute$Participant & h1n1_df$Visit > 1,]
h1n1_df_conv <- h1n1_df_conv %>%
  group_by(Participant) %>%
  mutate(CA09 = max(CA09)) %>%
  mutate(WI22 = max(WI22)) %>%
  mutate(DC24 = max(DC24)) %>%
  mutate(BC24 = max(BC24)) %>%
  mutate(TL22 = max(TL22)) %>%
  mutate(Timepoint = "Convalescent") %>%
  ungroup()
h1n1_df_conv <- h1n1_df_conv[!duplicated(h1n1_df_conv$Participant), c("Participant", "Timepoint", "CA09", "WI22", "DC24", "BC24", "TL22")]

h1n1_df_combined <- rbind(h1n1_df_acute, h1n1_df_conv)
h1n1_df_combined

h1n1_df_combined_long <- pivot_longer(h1n1_df_combined, cols = c(3:7),
                                      names_to = "Virus",
                                      values_to = "Titer")
h1n1_df_combined_long <- h1n1_df_combined_long %>%
  mutate(Strain_num = as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))),
         x_pos = case_when(
           Timepoint == "Acute" ~ Strain_num - 0.15,
           Timepoint == "Convalescent" ~ Strain_num + 0.15
         ))


# For H3N2, determine acute and convalescent titers and add extra columns for plotting purposes
h3n2_df_acute <- h3n2_df[h3n2_df$Visit==1,]
h3n2_df_acute$Timepoint = "Acute"
h3n2_df_acute <- h3n2_df_acute[, c("Participant", "Timepoint", "CA09", "WI22", "DC24", "BC24", "TL22")]

h3n2_df_conv <- h3n2_df[h3n2_df$Participant %in% h3n2_df_acute$Participant & h3n2_df$Visit > 1,]
h3n2_df_conv <- h3n2_df_conv %>%
  group_by(Participant) %>%
  mutate(CA09 = max(CA09)) %>%
  mutate(WI22 = max(WI22)) %>%
  mutate(DC24 = max(DC24)) %>%
  mutate(BC24 = max(BC24)) %>%
  mutate(TL22 = max(TL22)) %>%
  mutate(Timepoint = "Convalescent") %>%
  ungroup()
h3n2_df_conv <- h3n2_df_conv[!duplicated(h3n2_df_conv$Participant), c("Participant", "Timepoint", "CA09", "WI22", "DC24", "BC24", "TL22")]

h3n2_df_combined <- rbind(h3n2_df_acute, h3n2_df_conv)
h3n2_df_combined

h3n2_df_combined_long <- pivot_longer(h3n2_df_combined, cols = c(3:7),
                                      names_to = "Virus",
                                      values_to = "Titer")
h3n2_df_combined_long <- h3n2_df_combined_long %>%
  mutate(Strain_num = as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))),
         x_pos = case_when(
           Timepoint == "Acute" ~ Strain_num - 0.15,
           Timepoint == "Convalescent" ~ Strain_num + 0.15
         ))


# For H1N1, run paired t-tests and prepare significance data
h1n1_pvals <- h1n1_df_combined_long %>%
  group_by(Virus) %>%
  summarize(
    p = t.test(
      log2(Titer[Timepoint == "Acute"]),
      log2(Titer[Timepoint == "Convalescent"]),
      paired = TRUE
    )$p.value
  ) %>%
  mutate(
    y.position = case_when(
      Virus == "CA09" ~ log2(max(h1n1_df_combined_long[h1n1_df_combined_long$Virus == "CA09",]$Titer)) + 1,
      Virus == "WI22" ~ log2(max(h1n1_df_combined_long[h1n1_df_combined_long$Virus == "WI22",]$Titer)) + 1,
      Virus == "DC24" ~ log2(max(h1n1_df_combined_long[h1n1_df_combined_long$Virus == "DC24",]$Titer)) + 1,
      Virus == "BC24" ~ log2(max(h1n1_df_combined_long[h1n1_df_combined_long$Virus == "BC24",]$Titer)) + 1,
      Virus == "TL22" ~ log2(max(h1n1_df_combined_long[h1n1_df_combined_long$Virus == "TL22",]$Titer)) + 1
    ),
    #c(15, 15, 15, 15), # Adjust to avoid overlapping data
    group1 = case_when(
      Virus == "CA09" ~ 0.85,
      Virus == "WI22" ~ 1.85,
      Virus == "DC24" ~ 2.85,
      Virus == "BC24" ~ 3.85,
      Virus == "TL22" ~ 4.85
    ),
    group2 = group1 + 0.3,
    p.label = ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**",
                            ifelse(p < 0.05, "*", "ns"))),
    p.adj = p.adjust(p, method = "holm"),
    p.adj.round = ifelse(p.adj < 0.001, scientific(p.adj, digits=2),
                         ifelse(p.adj > 0.1, format(round(p.adj, 2)), format(round(p.adj, 4))))
  )

h1n1_summary_stats <- h1n1_df_combined_long %>%
  group_by(Virus, Timepoint) %>%
  summarize(
    n = n(),
    geo_mean = mean(log2(Titer)),
    se = sd(log2(Titer)) / sqrt(n),
    ci_lower = mean(log2(Titer) - qt(0.975, df = n - 1) * se),
    ci_upper = mean(log2(Titer) + qt(0.975, df = n - 1) * se)
  ) %>%
  ungroup()

h1n1_summary_stats <- h1n1_summary_stats %>%
  mutate(x_pos = case_when(
    Timepoint == "Acute" ~ as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))) - 0.15,
    Timepoint == "Convalescent" ~ as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))) + 0.15
  ))


# For H3N2, run paired t-tests and prepare significance data
h3n2_pvals <- h3n2_df_combined_long %>%
  group_by(Virus) %>%
  summarize(
    p = t.test(
      log2(Titer[Timepoint == "Acute"]),
      log2(Titer[Timepoint == "Convalescent"]),
      paired = TRUE
    )$p.value
  ) %>%
  mutate(
    y.position = case_when(
      Virus == "CA09" ~ log2(max(h3n2_df_combined_long[h3n2_df_combined_long$Virus == "CA09",]$Titer)) + 1,
      Virus == "WI22" ~ log2(max(h3n2_df_combined_long[h3n2_df_combined_long$Virus == "WI22",]$Titer)) + 1,
      Virus == "DC24" ~ log2(max(h3n2_df_combined_long[h3n2_df_combined_long$Virus == "DC24",]$Titer)) + 1,
      Virus == "BC24" ~ log2(max(h3n2_df_combined_long[h3n2_df_combined_long$Virus == "BC24",]$Titer)) + 1,
      Virus == "TL22" ~ log2(max(h3n2_df_combined_long[h3n2_df_combined_long$Virus == "TL22",]$Titer)) + 1
    ),
    #c(15, 15, 15, 15), # Adjust to avoid overlapping data
    group1 = case_when(
      Virus == "CA09" ~ 0.85,
      Virus == "WI22" ~ 1.85,
      Virus == "DC24" ~ 2.85,
      Virus == "BC24" ~ 3.85,
      Virus == "TL22" ~ 4.85
    ),
    group2 = group1 + 0.3,
    p.label = ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**",
                            ifelse(p < 0.05, "*", "ns"))),
    p.adj = p.adjust(p, method = "holm"),
    p.adj.round = ifelse(p.adj < 0.001, scientific(p.adj, digits=2),
                         ifelse(p.adj > 0.1, format(round(p.adj, 2)), format(round(p.adj, 4))))
  )

h3n2_summary_stats <- h3n2_df_combined_long %>%
  group_by(Virus, Timepoint) %>%
  summarize(
    n = n(),
    geo_mean = mean(log2(Titer)),
    se = sd(log2(Titer)) / sqrt(n),
    ci_lower = mean(log2(Titer) - qt(0.975, df = n - 1) * se),
    ci_upper = mean(log2(Titer) + qt(0.975, df = n - 1) * se)
  ) %>%
  ungroup()

h3n2_summary_stats <- h3n2_summary_stats %>%
  mutate(x_pos = case_when(
    Timepoint == "Acute" ~ as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))) - 0.15,
    Timepoint == "Convalescent" ~ as.numeric(factor(Virus, levels = c("CA09", "WI22", "DC24", "BC24", "TL22"))) + 0.15
  ))


# plot H1N1 data
plot = ggplot(h1n1_df_combined_long, aes(x = x_pos, y = Titer)) +
  geom_hline(yintercept=20, color = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_line(aes(group = interaction(Participant, Virus)), color = "grey80", linewidth = 0.3) +
  geom_jitter(aes(color = Timepoint), size = 0.5, width = 0, height = 0.1) +
  scale_color_manual(values = c("Acute" = "black", "Convalescent" = "#704961")) +
  scale_x_continuous(
    breaks = 1:5,
    labels = case_match(c("CA09", "WI22", "DC24", "BC24", "TL22"),
                        "CA09" ~ "A/California/07/2009",
                        "WI22" ~ "A/Wisconsin/50/2022",
                        "DC24" ~ "A/Dairy Cow/Texas/24-008749-002-v/2024",
                        "BC24" ~ "A/British Columbia/PHL-2032/2024",
                        "TL22" ~ "A/Thailand/8/2022")
  ) +
  scale_y_continuous(transform = "log2",
                     # limits = c(8, 131072),
                     breaks = c(20, 80, 320, 1280, 5120, 20480, 81920),
                     labels = c(20, 80, 320, "1,280", "5,120", "20,480", "81,920"),
                     minor_breaks = c(10, 40, 160, 640, 2560, 10240, 40960)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    rect = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 7),
    axis.text = element_text(size = 5), #change font size of axis text
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5), #change font size of axis titles
    panel.border = element_rect(linewidth = 0.2, color = "black", fill = NA),
    axis.line = element_line(linewidth = 0.1),
    strip.text = element_text(size = 12), # Adjust facet label front size
    plot.margin = margin(10, 30, 10, 10),
    axis.minor.ticks.length = rel(0.5),
    axis.ticks = element_line(linewidth = 0.25)
  ) +
  labs(
    y = expression("Geometric mean "~NAI[50]~" titer"),
    title = "H1N1 infections (n=13)\n"
  ) +
  stat_pvalue_manual(
    h1n1_pvals,
    # label = "p.label",
    label = "p.adj.round",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 1.75,
    na.rm = TRUE,
    vjust = -0.5
  ) +
  geom_errorbar(
    data = h1n1_summary_stats,
    aes(x = x_pos, y = 2^geo_mean, ymin = 2^ci_lower, ymax = 2^ci_upper, color = factor(Timepoint)),
    size = 0.3,
    width = 0.125,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = h1n1_summary_stats,
    aes(
      x = x_pos - 0.125, xend = x_pos + 0.125,
      y = 2^geo_mean, yend = 2^geo_mean,
      color = factor(Timepoint)
    ),
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  guides(
    y = guide_axis(minor.ticks = TRUE)
  ) +
  geom_vline(xintercept = 2.5, linewidth = 0.2, color = "black") +
  geom_vline(xintercept = 4.5, linewidth = 0.2, color = "black") +
  annotate("text", label = "H1N1", x=1.5, y=400000, size = 2) +
  annotate("text", label = "H5N1", x=3.5, y=400000, size = 2) +
  annotate("text", label = "H3N2", x=5, y=400000, size = 2) +
  coord_cartesian(ylim = c(8, 131072), clip="off")
  
  ggsave("H1N1_infected.pdf", path = "Ort-et-al_N1-imprinting/Figures/Figure2", width = 1.5, height = 1.5, dpi = 300, scale = 2, bg = "transparent")


# plot H3N2 data
plot = ggplot(h3n2_df_combined_long, aes(x = x_pos, y = Titer)) +
  geom_hline(yintercept=20, color = "grey80", linewidth = 0.3, linetype = "dashed") +
  geom_line(aes(group = interaction(Participant, Virus)), color = "grey80", linewidth = 0.3) +
  geom_jitter(aes(color = Timepoint), size = 0.5, width = 0, height = 0.1) +
  scale_color_manual(values = c("Acute" = "black", "Convalescent" = "#863f20")) +
  scale_x_continuous(
    breaks = 1:5,
    labels = case_match(c("CA09", "WI22", "DC24", "BC24", "TL22"),
                        "CA09" ~ "A/California/07/2009",
                        "WI22" ~ "A/Wisconsin/50/2022",
                        "DC24" ~ "A/Dairy Cow/Texas/24-008749-002-v/2024",
                        "BC24" ~ "A/British Columbia/PHL-2032/2024",
                        "TL22" ~ "A/Thailand/8/2022")
  ) +
  scale_y_continuous(transform = "log2",
                     # limits = c(8, 131072),
                     breaks = c(20, 80, 320, 1280, 5120, 20480, 81920),
                     labels = c(20, 80, 320, "1,280", "5,120", "20,480", "81,920"),
                     minor_breaks = c(10, 40, 160, 640, 2560, 10240, 40960)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    rect = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 7),
    axis.text = element_text(size = 5), #change font size of axis text
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5), #change font size of axis titles
    panel.border = element_rect(linewidth = 0.2, color = "black", fill = NA),
    axis.line = element_line(linewidth = 0.1),
    strip.text = element_text(size = 12), # Adjust facet label front size
    plot.margin = margin(10, 30, 10, 10),
    axis.minor.ticks.length = rel(0.5),
    axis.ticks = element_line(linewidth = 0.25)
  ) +
  labs(
    y = expression("Geometric mean "~NAI[50]~" titer"),
    title = "H3N2 infections (n=12)\n"
  ) +
  stat_pvalue_manual(
    h3n2_pvals,
    # label = "p.label",
    label = "p.adj.round",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.02,
    size = 1.75,
    na.rm = TRUE,
    vjust = -0.5
  ) +
  geom_errorbar(
    data = h3n2_summary_stats,
    aes(x = x_pos, y = 2^geo_mean, ymin = 2^ci_lower, ymax = 2^ci_upper, color = factor(Timepoint)),
    size = 0.3,
    width = 0.125,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = h3n2_summary_stats,
    aes(
      x = x_pos - 0.125, xend = x_pos + 0.125,
      y = 2^geo_mean, yend = 2^geo_mean,
      color = factor(Timepoint)
    ),
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  guides(
    y = guide_axis(minor.ticks = TRUE)
  ) +
  geom_vline(xintercept = 2.5, linewidth = 0.2, color = "black") +
  geom_vline(xintercept = 4.5, linewidth = 0.2, color = "black") +
  annotate("text", label = "H1N1", x=1.5, y=400000, size = 2) +
  annotate("text", label = "H5N1", x=3.5, y=400000, size = 2) +
  annotate("text", label = "H3N2", x=5, y=400000, size = 2) + 
  coord_cartesian(ylim = c(8, 131072), clip="off")
  
  ggsave("H3N2_infected.pdf", path = "Ort-et-al_N1-imprinting/Figures/Figure2", width = 1.5, height = 1.5, dpi = 300, scale = 2, bg = "transparent")
  
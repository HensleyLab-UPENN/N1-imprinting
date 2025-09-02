# Load packages 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(ggpmisc)


# Load NAI results
wide_df_2005 <- read_csv("Ort-et-al_N1-imprinting/Data/2005_Titers.csv")
wide_df_2017 <- read_csv("Ort-et-al_N1-imprinting/Data/2017_Titers.csv")


# CA09 vs each H5N1
# VN04
for (virus in list("VN04", "DC24", "BC24")) {
  plot = ggplot(data=NULL, aes(x=CA09, y=!!sym(virus))) +
    geom_jitter(
      data = wide_df_2005,
      shape = 22,
      color = "grey40",
      size = 1.5,
      fill = "#506E82",
      width = 0.15,
      height = 0.15,
      alpha = 0.6,
      stroke = 0.1
    ) +
    geom_jitter(
      data = wide_df_2017,
      shape = 21,
      color = "grey40",
      size = 1.5,
      fill = "#768655",
      width = 0.15,
      height = 0.15,
      alpha = 0.6,
      stroke = 0.1
    ) +
    geom_smooth(
      data = wide_df_2005,
      method = "lm",
      col = "#506E82",
      fill = "#506E82",
      alpha = 0.4,
      span = .6
    ) +
    stat_correlation(
      data = wide_df_2005,
      geom = "text",
      use_label("R", "p"),
      coef.keep.zeros = T,
      size = 1.5,
      label.y = 12.91,
      label.x = 4,
      color = "#506E82"
    ) +
    geom_smooth(
      data = wide_df_2017,
      method = "lm",
      col = "#768655",
      fill = "#768655",
      alpha = 0.4,
      span = .6
    ) +
    stat_correlation(
      data = wide_df_2017,
      geom = "text",
      use_label("R", "p"),
      coef.keep.zeros = T,
      size = 1.5,
      label.y = 12.41,
      label.x = 4,
      color = "#768655"
    ) +
    coord_cartesian(xlim = c(8, 8192), ylim = c(8, 8192)) +
    scale_x_continuous(transform = "log2", breaks = c(10,20,40,80,160,320,640,1280,2560,5120)) +
    scale_y_continuous(transform = "log2", breaks = c(10,20,40,80,160,320,640,1280,2560,5120)) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      rect = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 7),
      axis.text = element_text(size = 5),
      axis.title = element_text(size = 5),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(linewidth = 0.1),
      strip.text = element_text(size = 12),
      plot.margin = margin(10, 30, 10, 10),
      axis.ticks = element_line(linewidth = 0.25)
    ) +
    geom_abline(linetype ="dashed",
                color = "black",
                linewidth = 0.2) +
    labs(
      y = "H5N1 A/Vietnam/1203/2004\ngeometric mean NAI50 titer",
      x = "H1N1 A/California/07/2009\ngeometric mean NAI50 titer"
    )
  
  ggsave(paste0(virus, ".pdf"), path = "Ort-et-al_N1-imprinting/Figures/Extended_Figure1", width = 1.68, height = 1.5, dpi = 300, scale = 2, bg = "transparent")
}
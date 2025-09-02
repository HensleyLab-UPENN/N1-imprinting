# Load packages 
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)


# Load NAI results
wide_df_2017 <- read_csv("Ort-et-al_N1-imprinting/Data/2017_Titers.csv")
long_df_2017 <- pivot_longer(wide_df_2017, cols = c(2:5),
                             names_to = "virus",
                             values_to = "titer")

# plot 2017 titers vs YOB
# fit with loess
for (virus in list("CA09", "VN04", "DC24", "BC24")) {
    plot = ggplot(long_df_2017[long_df_2017[, "virus"] == virus,], aes(x = yob, y = titer)) +
    geom_vline(xintercept = c(1957, 1968, 1977, 2009), linetype ="dashed", 
               color = "black", linewidth = 0.2) +
    geom_jitter(shape = 21, color = "grey40", size = 1.5,
                fill = "#768655",
                width = 0.1, height = 0.1,
                alpha = 0.6, stroke = 0.1) +
    geom_smooth(col = "#768655", fill = "#768655",
                alpha = 0.4, span = 0.5) +
    annotate("text", x = 1960.5, y = 7680, label = "1957",
             size = 1.5, color = "black") +
    annotate("text", x = 1971.5, y = 7680, label = "1968",
             size = 1.5, color = "black") +
    annotate("text", x = 1980.5, y = 7680, label = "1977",
             size = 1.5, color = "black") +
    annotate("text", x = 2012.5, y = 7680, label = "2009",
             size = 1.5, color = "black") +
    scale_x_continuous(breaks = seq(1920, 2020, 10)) +
    coord_cartesian(ylim = c(8, 8192)) + 
    scale_y_continuous(trans = "log2",
                       labels = c(10,20,40,80,160,320,640,"1,280","2,560","5,120"),
                       breaks = c(10,20,40,80,160,320,640,1280,2560,5120)) + 
    theme_classic() +
    theme(panel.grid.major = element_blank(),
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
    labs(y = "Geometric mean NAI50 titer",
         x = "Year of birth",
         title = case_match(
           virus,
           "CA09" ~ "H1N1 A/California/07/2009",
           "VN04" ~ "H5N1 A/Vietnam/1203/2004",
           "DC24" ~ "H5N1 A/Dairy Cow/Texas/24-008749-002-v/2024",
           "BC24" ~ "H5N1 A/British Columbia/PHL-2032/2024"
          )
    )
    
    ggsave(paste0(virus, ".pdf"), path = "Ort-et-al_N1-imprinting/Figures/Figure1", width = 1.5, height = 1, dpi = 300, scale = 2, bg = "transparent")
}


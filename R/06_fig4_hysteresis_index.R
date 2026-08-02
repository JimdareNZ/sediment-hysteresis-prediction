# ══════════════════════════════════════════════════════════════════════════════
# 06_fig4_hysteresis_index.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 4: Violin plots of continuous hysteresis index h across sites.
#
# Input:  data/Final_Hysteresis_Dataset_Spatial_Ante.csv
# Output: figures/Fig4_Hysteresis_Index_Distribution.png
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(tidyverse)

df <- read_csv("data/Final_Hysteresis_Dataset_Spatial_Ante.csv")

df <- df %>%
  mutate(
    Site_Code = case_when(
      Site == "Tuangi" ~ "PKO-D",  Site == "Inanga" ~ "PKO-UW",
      Site == "Kokopu" ~ "PKO-UE", Site == "Kokako" ~ "MGT",
      Site == "Pipi"   ~ "PUN",    Site == "Kotuku" ~ "WHR-D",
      Site == "Matata" ~ "PNG-M",  Site == "Tuna"   ~ "PNG-TR1"
    ),
    Sub_Catchment = case_when(
      Site %in% c("Tuangi", "Inanga", "Kokopu") ~ "Pokopoko",
      Site == "Kokako" ~ "Mangatoetoe", Site == "Pipi" ~ "Puanene",
      Site == "Kotuku" ~ "Wharere", Site %in% c("Matata", "Tuna") ~ "Pongakawa"
    )
  ) %>%
  filter(!Zuec_class %in% c(5, 6, 7, 8))

site_order <- c("PKO-D", "PKO-UW", "PKO-UE", "MGT", "PUN", "WHR-D", "PNG-M", "PNG-TR1")
df <- df %>% mutate(Site_Code = factor(Site_Code, levels = site_order))

sc_colours <- c("Pokopoko" = "#f4a582", "Mangatoetoe" = "#ca0020",
                "Puanene" = "#9970ab", "Wharere" = "#404040", "Pongakawa" = "#92c5de")

fig4 <- ggplot(df %>% filter(!is.na(Zuec_h), !is.na(Site_Code)),
               aes(x = Site_Code, y = Zuec_h, fill = Sub_Catchment)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_violin(aes(colour = Sub_Catchment), alpha = 0.6, linewidth = 0.4, trim = TRUE, scale = "width") +
  geom_boxplot(aes(colour = Sub_Catchment), width = 0.15, fill = "white",
               outlier.shape = NA, linewidth = 0.5) +
  scale_fill_manual(values = sc_colours, name = "Sub-catchment") +
  scale_colour_manual(values = sc_colours, name = "Sub-catchment") +
  annotate("text", x = 0.5, y = 0.55, label = "Clockwise",
           hjust = 0, size = 3, colour = "grey40", fontface = "italic") +
  annotate("text", x = 0.5, y = -0.55, label = "Anticlockwise",
           hjust = 0, size = 3, colour = "grey40", fontface = "italic") +
  scale_y_continuous(name = "Hysteresis index (h)", limits = c(-0.75, 0.75),
                     breaks = seq(-0.75, 0.75, 0.25)) +
  labs(x = "Monitoring site") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9), axis.title = element_text(size = 10),
        legend.position = "right", legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(fill = guide_legend(override.aes = list(colour = unname(sc_colours), alpha = 0.8)),
         colour = "none")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig4_Hysteresis_Index_Distribution.png", fig4,
       dpi = 300, width = 10, height = 6, units = "in", bg = "white")

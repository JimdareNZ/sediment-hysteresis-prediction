# ══════════════════════════════════════════════════════════════════════════════
# 06_fig3_class_distribution.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 3: Percentage distribution of hysteresis classes across all sites.
#
# Input:  data/Final_Hysteresis_Dataset_Spatial_Ante.csv
# Output: figures/Fig3_Hysteresis_Class_Distribution.png
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(tidyverse)

df <- read_csv("data/Final_Hysteresis_Dataset_Spatial_Ante.csv")

df <- df %>%
  mutate(Site_Code = case_when(
    Site == "Tuangi" ~ "PKO-D",  Site == "Inanga" ~ "PKO-UW",
    Site == "Kokopu" ~ "PKO-UE", Site == "Kokako" ~ "MGT",
    Site == "Pipi"   ~ "PUN",    Site == "Kotuku" ~ "WHR-D",
    Site == "Matata" ~ "PNG-M",  Site == "Tuna"   ~ "PNG-TR1"
  ))

site_order <- c("PKO-D", "PKO-UW", "PKO-UE", "MGT", "PUN", "WHR-D", "PNG-M", "PNG-TR1")
df <- df %>% mutate(Site_Code = factor(Site_Code, levels = site_order))

class_labels <- c(
  "1" = "Class I (Clockwise)",          "2" = "Class II (Mixed, CW dominant)",
  "3" = "Class III (Mixed, ACW dominant)", "4" = "Class IV (Anticlockwise)",
  "5" = "Class V (Dilution, CW)",       "6" = "Class VI (Dilution, Mixed CW)",
  "7" = "Class VII (Dilution, Mixed ACW)", "8" = "Class VIII (Dilution, ACW)"
)

class_colours <- c(
  "1" = "#d73027", "2" = "#f4a582", "3" = "#92c5de", "4" = "#2166ac",
  "5" = "#b8b8b8", "6" = "#969696", "7" = "#737373", "8" = "#525252"
)

df_summary <- df %>%
  filter(!is.na(Zuec_class), !is.na(Site_Code)) %>%
  mutate(Zuec_class = as.character(Zuec_class)) %>%
  group_by(Site_Code, Zuec_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Site_Code) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup() %>%
  mutate(Zuec_class = factor(Zuec_class, levels = as.character(1:8)))

site_n <- df_summary %>% group_by(Site_Code) %>% summarise(total = first(total))

fig3 <- ggplot(df_summary, aes(x = Site_Code, y = prop, fill = Zuec_class)) +
  geom_bar(stat = "identity", width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(data = site_n, aes(x = Site_Code, y = 1.04, label = paste0("n=", total)),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  scale_fill_manual(values = class_colours, labels = class_labels, name = "Hysteresis class") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Monitoring site", y = "Percentage of events (%)") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9), axis.title = element_text(size = 10),
        legend.position = "right", legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8), legend.key.size = unit(0.5, "cm"),
        plot.background = element_rect(fill = "white", colour = NA))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig3_Hysteresis_Class_Distribution.png", fig3,
       dpi = 300, width = 10, height = 6, units = "in", bg = "white")

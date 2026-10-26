# ==========================================
# 1. SETUP & LIBRARIES
# ==========================================
library(tidyverse)
library(readxl)
library(patchwork)
library(ggpubr)
library(svglite)
library(ggtext)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

MET_MIC  <- file.path(INPUT_DIR, "inputMetMic.xlsx")
METADATA <- read_excel(MET_MIC, sheet = "METADATA")
Met      <- read_excel(MET_MIC, sheet = "Met")

# ==========================================
# 2. DATA PROCESSING & MAPPING
# ==========================================
METADATA <- METADATA %>%
  mutate(Treatments_Clean = factor(case_when(
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched C. asiatica",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "C. asiatica",
    TRUE ~ "Placebo"
  ), levels = c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")))

color_palette <- c(
  "Placebo" = "#999999",
  "C. asiatica" = "#e6c60d",
  "Fiber-enriched C. asiatica" = "#5ca08e"
)

# ==========================================
# 3. PLOTTING
# ==========================================
FONT <- "Arial"
class_cols <- list("Total SCFAs" = c(4, 6, 7, 9), "Total CEs" = c(18, 19, 20),
                   "Total FAs" = 10:16, "Total UTs" = 40)

df_plot <- bind_rows(lapply(names(class_cols), function(cl) {
  Met %>%
    mutate(Total_Val = rowSums(select(., all_of(colnames(Met)[class_cols[[cl]]])), na.rm = TRUE)) %>%
    select(ID, Total_Val) %>%
    mutate(Class = cl)
})) %>%
  left_join(METADATA %>% select(ID, Weeks, Treatments_Clean), by = "ID") %>%
  mutate(Weeks = factor(paste0("W", Weeks), levels = c("W0", "W2", "W4")),
         Class = factor(Class, levels = names(class_cols))) %>%
  filter(!is.na(Treatments_Clean))

facet_labels <- c(
  "Placebo"                    = "Placebo",
  "C. asiatica"                = "*C. asiatica*",
  "Fiber-enriched C. asiatica" = "Fiber-<br>enriched<br>*C. asiatica*"
)

plot_class <- function(cl) {
  ggplot(filter(df_plot, Class == cl), aes(x = Weeks, y = Total_Val, fill = Treatments_Clean)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, color = "black", linewidth = 0.25) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.25, color = "black") +
    facet_grid(. ~ Treatments_Clean, labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = color_palette) +
    labs(title = cl, y = NULL, x = NULL) +
    theme_bw(base_size = 7, base_family = FONT) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_markdown(size = 7, halign = 0.5, lineheight = 1.1, margin = margin(1, 0.5, 1, 0.5)),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(0, 0, 1, 0)),
      plot.title.position = "plot",
      axis.text = element_text(size = 7, colour = "black"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      panel.spacing.x = unit(0.6, "mm"),
      plot.margin = margin(0.5, 1, 0.5, 1, "mm")
    )
}

final_plot <- wrap_plots(lapply(names(class_cols), plot_class), nrow = 1)

# ==========================================
# 4. SAVE
# ==========================================
# Panel size as placed in Figures_EDITED pptx (13-inch-high slides) at 180 mm print width
W <- 180; H <- 36
ggsave(file.path(RESULT_DIR, "Figure4A.png"), plot = final_plot, width = W, height = H, units = "mm", dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure4A.svg"), plot = final_plot, width = W, height = H, units = "mm",
       device = svglite, system_fonts = list(Arial = "Arial"))

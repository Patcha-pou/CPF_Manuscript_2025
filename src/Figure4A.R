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
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched Gotu kola juice (F-GKJ)",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "Gotu kola juice-only (GKJ)",
    TRUE ~ "Placebo (CON)"
  ), levels = c("Placebo (CON)", "Gotu kola juice-only (GKJ)", "Fiber-enriched Gotu kola juice (F-GKJ)")))

color_palette <- c(
  "Placebo (CON)" = "#999999",
  "Gotu kola juice-only (GKJ)" = "#e6c60d",
  "Fiber-enriched Gotu kola juice (F-GKJ)" = "#5ca08e"
)

# ==========================================
# 3. PLOTTING
# ==========================================
FONT <- "Arial"
class_cols <- list(
  "Total SCFAs (feces)" = c("Acetic acid", "Propionic acid", "Butanoic acid", "Pentanoic acid",
                            "2-Methylpropionic acid", "3-Methylbutanoic acid", "2-Methylbutanoic acid"),
  "Total Cholesterol and CEs (plasma)" = c("total CE", "Cholesterol"),
  "Total FAs (plasma)"  = "Total Lipid",
  "Total UTs (plasma)"  = c("p-Cresolsulfate", "Phenysulfate", "Indoxyl sulfate", "2-Methoxyphenylsulfate",
                            "Hippuric acid", "Cinnamoylglycine", "m-Hydroxyhippuric acid",
                            "Phenylacetyl-L-glutamine", "Trimethylamine N-oxide", "Dimethylglycine"))
class_units <- c("Total SCFAs (feces)" = "µmol/g\n(wet weight)", "Total Cholesterol and CEs (plasma)" = "µM",
                 "Total FAs (plasma)" = "µM", "Total UTs (plasma)" = "µM")

df_plot <- bind_rows(lapply(names(class_cols), function(cl) {
  vals <- Met %>% select(all_of(class_cols[[cl]])) %>% mutate(across(everything(), as.numeric))
  Met %>%
    mutate(Total_Val = ifelse(rowSums(!is.na(vals)) == 0, NA, rowSums(vals, na.rm = TRUE))) %>%
    select(ID, Total_Val) %>%
    mutate(Class = cl)
})) %>%
  left_join(METADATA %>% select(ID, Weeks, Treatments_Clean), by = "ID") %>%
  mutate(Weeks = factor(paste0("W", Weeks), levels = c("W0", "W2", "W4")),
         Class = factor(Class, levels = names(class_cols))) %>%
  filter(!is.na(Treatments_Clean))

facet_labels <- c(
  "Placebo (CON)"                    = "CON",
  "Gotu kola juice-only (GKJ)"                = "GKJ",
  "Fiber-enriched Gotu kola juice (F-GKJ)" = "F-GKJ"
)

plot_class <- function(cl) {
  ggplot(filter(df_plot, Class == cl, !is.na(Total_Val)), aes(x = Weeks, y = Total_Val, fill = Treatments_Clean)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, color = "black", linewidth = 0.25) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.25, color = "black") +
    facet_grid(. ~ Treatments_Clean, labeller = as_labeller(facet_labels)) +
    scale_fill_manual(values = color_palette) +
    (if (cl == "Total UTs (plasma)") scale_y_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000"))
     else scale_y_continuous()) +
    labs(title = cl, y = class_units[[cl]], x = NULL) +
    theme_bw(base_size = 7, base_family = FONT) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_markdown(size = 7, halign = 0.5, lineheight = 1.1, margin = margin(1, 0.5, 1, 0.5)),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold", margin = margin(0, 0, 1, 0)),
      plot.title.position = "plot",
      axis.title.y = element_text(size = 8, face = "bold", margin = margin(0, 0.5, 0, 0, "mm")),
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

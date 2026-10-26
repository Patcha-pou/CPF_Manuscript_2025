# ==========================================
# 1. SETUP & LIBRARIES
# ==========================================
library(tidyverse)
library(readxl)
library(missForest)
library(patchwork)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

MET_MIC  <- file.path(INPUT_DIR, "inputMetMic.xlsx")
METADATA <- read_excel(MET_MIC, sheet = "METADATA")
Mic      <- read_excel(MET_MIC, sheet = "Mic")

METADATA <- METADATA %>%
  mutate(Treatments = factor(case_when(
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched C. asiatica",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "C. asiatica",
    TRUE ~ "Placebo"
  ), levels = c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")))

# ==========================================
# 2. IMPUTATION (missForest)
# ==========================================
Mic_ids <- Mic$ID
input_data_raw <- Mic %>% dplyr::select(-ID) %>% as.data.frame()
input_data_raw[] <- lapply(input_data_raw, function(x) as.numeric(as.character(x)))

set.seed(123)
imputed_result <- missForest(input_data_raw)
input_data_raw <- imputed_result$ximp
rownames(input_data_raw) <- Mic_ids

# ==========================================
# 3. FEATURE SELECTION & RENAMING
# ==========================================
target_columns <- c(
  "genus103:Ruminococcus",
  "genus63:Fusicatenibacter",
  "genus100:Ruminococcaceae_Incertae_Sedis",
  "genus73:Lachnospiraceae_UCG-004"
)

rename_map <- c(
  "genus103:Ruminococcus" = "Ruminococcus",
  "genus63:Fusicatenibacter" = "Fusicatenibacter",
  "genus100:Ruminococcaceae_Incertae_Sedis" = "Ruminococcaceae\nIncertae Sedis",
  "genus73:Lachnospiraceae_UCG-004" = "Lachnospiraceae UCG-004"
)

# ==========================================
# 4. PLOTTING FUNCTION
# ==========================================
color_palette <- c(
  "Placebo" = "#999999", 
  "C. asiatica" = "#e6c60d", 
  "Fiber-enriched C. asiatica" = "#5ca08e"
)

data_long <- input_data_raw %>%
  dplyr::select(any_of(target_columns)) %>%
  rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "Taxa", values_to = "Value") %>%
  left_join(METADATA %>% dplyr::select(ID, Subject, Weeks, Treatments), by = "ID") %>%
  mutate(
    Weeks = as.numeric(as.character(Weeks)),
    Taxa = factor(rename_map[Taxa], levels = rename_map)
  )

FONT <- "Arial"
p_average <- ggplot(data_long, aes(x = Weeks, y = Value, color = Treatments, group = Treatments)) +
  stat_summary(fun = median, fun.min = function(x) quantile(x, 0.25),
               fun.max = function(x) quantile(x, 0.75), geom = "errorbar",
               width = 0.25, linewidth = 0.3, position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "line", linewidth = 0.6, position = position_dodge(0.4)) +
  stat_summary(fun = median, geom = "point", size = 0.9, position = position_dodge(0.4)) +
  facet_wrap(~ Taxa, scales = "free_y", ncol = 1) +
  
  scale_color_manual(
    values = color_palette,
    labels = c(
      "Placebo" = "Placebo",
      "C. asiatica" = expression(italic("C. asiatica")),
      "Fiber-enriched C. asiatica" = expression(paste("Fiber-enriched ", italic("C. asiatica")))
    )
  ) +
  
  scale_x_continuous(breaks = c(0, 2, 4)) +
  theme_bw(base_size = 7, base_family = FONT) +
  labs(y = "Relative Abundance", 
       x = "Weeks") +
  guides(color = guide_legend(ncol = 1)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7), 
    legend.key.size = unit(3, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    strip.text = element_text(face = "bold.italic", size = 7, margin = margin(1, 1, 1, 1)),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7, colour = "black"),
    plot.title = element_text(size = 9, face = "bold", hjust = 0),
    plot.margin = margin(1, 1, 1, 1, "mm"),
    panel.grid.minor = element_blank()
  )

# ==========================================
# 5. SAVE FILES
# ==========================================
# Panel size as placed in Figures_EDITED pptx at 180 mm print width
W <- 50.6; H <- 122
ggsave(file.path(RESULT_DIR, "Figure3B.png"), plot = p_average, width = W, height = H, units = "mm", dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure3B.svg"), plot = p_average, width = W, height = H, units = "mm",
       device = svglite::svglite, system_fonts = list(Arial = "Arial"))

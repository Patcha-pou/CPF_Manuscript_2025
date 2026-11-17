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
Met      <- read_excel(MET_MIC, sheet = "Met")

METADATA <- METADATA %>%
  mutate(Treatments = factor(case_when(
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched Gotu kola juice (F-GKJ)",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "Gotu kola juice-only (GKJ)",
    TRUE ~ "Placebo (CON)"
  ), levels = c("Placebo (CON)", "Gotu kola juice-only (GKJ)", "Fiber-enriched Gotu kola juice (F-GKJ)"))) %>%
  as.data.frame()
rownames(METADATA) <- METADATA$ID

# ==========================================
# 2. IMPUTATION (missForest)
# ==========================================
met_ids <- Met$ID
input_data_raw <- Met %>% dplyr::select(-ID) %>% as.data.frame()
input_data_raw[] <- lapply(input_data_raw, function(x) as.numeric(as.character(x)))

set.seed(123)
imputed_result <- missForest(input_data_raw)
input_data_raw <- imputed_result$ximp
rownames(input_data_raw) <- met_ids

# ==========================================
# 3. FEATURE RENAMING (Mapping)
# ==========================================
rename_map <- c(
  "total.Lipid" = "Total Lipid",
  "ω6.ω3.ratio" = "omega6/omega3 ratio",
  "CE.20.5." = "CE(20:5)",
  "MUFA.SFA.ratio" = "MUFA/SFA ratio",
  "PUFA.SFA.ratio" = "PUFA/SFA ratio",
  "CE.16.0." = "CE(16:0)",
  "CE.20.4." = "CE(20:4)",
  "CE.17.0." = "CE(17:0)",
  "CE.18.0." = "CE(18:0)",
  "CE.20.0." = "CE(20:0)",
  "CE.20.2." = "CE(20:2)",
  "CE.18.1." = "CE(18:1)",
  "Acetic.acid" = "Acetic acid",
  "Butanoic.acid" = "Butanoic acid",
  "total.CE" = "Total CE",
  "CE.18.2." = "CE(18:2)",
  "Pentanoic.acid" = "Pentanoic acid",
  "CE.16.1." = "CE(16:1)",
  "X2.Methylhexanoic.acid" = "2-Methylhexanoic acid",
  "CE.15.0." = "CE(15:0)",
  "X2.Methypropionic.acid" = "2-Methypropionic acid",
  "X3.Methylbutanoic.acid" = "3-Methylbutanoic acid",
  "X2.Methylbutanoic.acid" = "2-Methylbutanoic acid",
  "p.Cresolsulfate" = "p-Cresolsulfate",
  "Indoxyl.sulfate" = "Indoxyl sulfate",
  "X2.Methoxyphenylsulfate" = "2-Methoxyphenylsulfate",
  "Hippuric.acid" = "Hippuric acid",
  "Phenylacetyl.L.glutamine" = "Phenylacetyl-L-glutamine",
  "Trimethylamine.N.oxide" = "Trimethylamine N-oxide"
)

# ==========================================
# 4. PLOTTING FUNCTION
# ==========================================
color_palette <- c(
  "Placebo (CON)" = "#999999", 
  "Gotu kola juice-only (GKJ)" = "#e6c60d", 
  "Fiber-enriched Gotu kola juice (F-GKJ)" = "#5ca08e"
)

FONT <- "Arial"
strip_breaks <- c("2-Methypropionic acid" = "2-Methyl-\npropionic\nacid", "2-Methylpropionic acid" = "2-Methyl-\npropionic\nacid", "3-Methylbutanoic acid" = "3-Methyl-\nbutanoic\nacid",
                  "2-Methylbutanoic acid" = "2-Methyl-\nbutanoic\nacid", "2-Methylhexanoic acid" = "2-Methyl-\nhexanoic\nacid",
                  "MUFA/SFA ratio" = "MUFA/SFA\nratio", "PUFA/SFA ratio" = "PUFA/SFA\nratio",
                  "Dimethylglycine" = "Dimethyl-\nglycine")
strip_label <- function(x) ifelse(x %in% names(strip_breaks), strip_breaks[x], sub(" ratio$", "\nratio", x))
tick_label <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 3))
PANEL_W <- unit(16, "mm"); PANEL_H <- unit(10, "mm")
plot_metabolite <- function(column_index, ylabel = NULL, x_title = NULL, title = NULL) {
  raw_name <- colnames(Met)[column_index]
  met_name <- ifelse(raw_name %in% names(rename_map), rename_map[raw_name], raw_name)

  data_long <- input_data_raw %>%
    dplyr::select(all_of(met_name)) %>%
    rownames_to_column("ID") %>%
    pivot_longer(-ID, names_to = "Metabolite", values_to = "Value") %>%
    left_join(METADATA %>% dplyr::select(ID, Subject, Weeks, Treatments), by = "ID") %>%
    mutate(Weeks = as.numeric(as.character(Weeks)))

  ggplot(data_long, aes(x = Weeks, y = Value, color = Treatments, group = Treatments)) +
    stat_summary(fun = median, fun.min = function(x) quantile(x, 0.25),
                 fun.max = function(x) quantile(x, 0.75), geom = "errorbar",
                 width = 0.25, linewidth = 0.25, position = position_dodge(0.4)) +
    stat_summary(fun = median, geom = "line", linewidth = 0.45, position = position_dodge(0.4)) +
    stat_summary(fun = median, geom = "point", size = 0.6, position = position_dodge(0.4)) +
    facet_wrap(~ Metabolite, labeller = as_labeller(strip_label)) +
    ggh4x::force_panelsizes(rows = PANEL_H, cols = PANEL_W) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 2), labels = tick_label) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(breaks = c(0, 2, 4)) +
    theme_bw(base_size = 7, base_family = FONT) +
    labs(title = title, y = ylabel, x = x_title) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 7, lineheight = 0.9, margin = margin(0.5, 0.3, 0.5, 0.3)),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7, colour = "black"),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.title.position = "plot",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm")
    )
}

# ==========================================
# 5. ASSEMBLY (SCFA, CE, FA, ratio, Other rows, in the heatmap order; one 4-column grid)
# ==========================================
plots <- c(
  Map(function(i, k) plot_metabolite(i, if (k == 1) "\u00b5mol/g\n(wet weight)", title = if (k == 1) "Short-chain fatty acids (SCFAs)"),
      c(4, 6, 7, 9), 1:4),
  list(plot_metabolite(11, "\u00b5M", title = "Cholesterol and CEs")),
  Map(function(i, k) plot_metabolite(i, if (k == 1) "\u00b5M", title = if (k == 1) "Fatty acids (FAs)"), 25:27, 1:3),
  Map(function(i, k) plot_metabolite(i, if (k == 1) "Ratio"), 28:30, 1:3),
  list(plot_metabolite(40, "\u00b5M", x_title = "Weeks", title = "Uremic toxins (UTs)"))
)
design <- "
ABCD
E###
FGH#
IJK#
L###
"
final_plot <- wrap_plots(plots, design = design, guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0)) &
  guides(color = guide_legend(ncol = 1, title.position = "left"))

# ==========================================
# 6. SAVE FILES
# ==========================================
# Panel size as placed in Figures_EDITED pptx (13-inch-high slides) at 180 mm print width
W <- 99; H <- 134
ggsave(file.path(RESULT_DIR, "Figure4C.png"), plot = final_plot, width = W, height = H, units = "mm", dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure4C.svg"), plot = final_plot, width = W, height = H, units = "mm",
       device = svglite::svglite, system_fonts = list(Arial = "Arial"))

# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)
library(svglite)

# --- 2. Clear Graphics ---
while (!is.null(dev.list())) dev.off()

# --- 3. Read & Process Data ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

TAXA <- read_excel(file.path(INPUT_DIR, "Downstream_analysis.xlsx"), sheet = "nonRedundantTAXA_filtered")
TAXA <- TAXA[,1:2]
colnames(TAXA) <- c("feature","Shortform")
TAXA$Shortform <- sub("^genus.*?:", "", TAXA$Shortform)

process_pval <- function(path, treat_name) {
  if(!file.exists(path)) return(NULL)
  df <- read.delim(path, check.names = FALSE)
  df$Time2_Time1 <- paste0(df$value, "/W0")
  df$Treatments <- treat_name
  return(df %>% select(feature, pval, qval, coef, Time2_Time1, Treatments))
}

MAASLIN_DIR <- file.path(OUTPUT_DIR, "Mic")

Placebo_pval <- process_pval(file.path(MAASLIN_DIR, "output_Placebo/all_results.tsv"), "Placebo")
Asiatica_pval <- process_pval(file.path(MAASLIN_DIR, "output_Asiatica/all_results.tsv"), "Asiatica")
Fiber_added_pval <- process_pval(file.path(MAASLIN_DIR, "output_Fiber_added/all_results.tsv"), "Fiber_added")
stopifnot(!is.null(Placebo_pval), !is.null(Asiatica_pval), !is.null(Fiber_added_pval))
pval_total <- rbind(Placebo_pval, Asiatica_pval, Fiber_added_pval)

# --- 4. Matrix Construction & Clustering ---
target_cols <- c("Placebo_W2/W0", "Placebo_W4/W0", "Asiatica_W2/W0", "Asiatica_W4/W0", "Fiber_added_W2/W0", "Fiber_added_W4/W0")

pval_matrix_data <- pval_total %>%
  mutate(Treatments_Weeks = paste0(Treatments, "_", Time2_Time1),
         coef = as.numeric(coef),
         sign_asso = (-log10(as.numeric(pval))) * sign(coef)) %>%
  filter(Treatments_Weeks %in% target_cols) %>%
  group_by(feature, Treatments_Weeks) %>%
  summarise(coef = mean(coef, na.rm = TRUE), sign_asso = mean(sign_asso, na.rm = TRUE),
            pval = min(pval, na.rm = TRUE),
            qval = min(qval, na.rm = TRUE), .groups = "drop")

sig_features <- pval_matrix_data %>% filter(-log10(pval) > 1) %>% pull(feature) %>% unique()

dup_feat <- sum(duplicated(TAXA$feature))
if (dup_feat > 0) {
  message("TAXA: dropped ", dup_feat, " duplicate 'feature' row(s)")
  TAXA <- TAXA[!duplicated(TAXA$feature), ]
}

wide_of <- function(value_col, fill) {
  pval_matrix_data %>%
    filter(feature %in% sig_features) %>%
    select(feature, Treatments_Weeks, all_of(value_col)) %>%
    pivot_wider(names_from = Treatments_Weeks, values_from = all_of(value_col),
                values_fill = fill) %>%
    as.data.frame()
}

wide_asso <- wide_of("coef", 0)
wide_sign <- wide_of("sign_asso", 0)
wide_pval <- wide_of("pval", 1)
wide_qval <- wide_of("qval", 1)
stopifnot(identical(wide_asso$feature, wide_pval$feature),
          identical(wide_asso$feature, wide_qval$feature))

row_labels <- TAXA$Shortform[match(wide_asso$feature, TAXA$feature)]
n_unmatched <- sum(is.na(row_labels))
if (n_unmatched > 0)
  message("TAXA: ", n_unmatched, " of ", nrow(wide_asso),
          " features had no Shortform — falling back to the raw feature name ",
          "(check that MaAsLin2 feature names match the sheet, e.g. ':' vs '.')")
row_labels[is.na(row_labels)] <- wide_asso$feature[is.na(row_labels)]
row_labels <- make.unique(row_labels, sep = " ")

matrix_asso <- as.matrix(wide_asso[, target_cols])
rownames(matrix_asso) <- row_labels

matrix_stars <- as.matrix(wide_qval[, target_cols])
matrix_stars <- ifelse(matrix_stars < 0.01, "**", ifelse(matrix_stars < 0.05, "*", ""))
rownames(matrix_stars) <- row_labels

has_star <- rowSums(as.matrix(wide_pval[, target_cols]) < 0.05) > 0
final_asso <- matrix_asso[has_star, ]
final_stars <- matrix_stars[has_star, ]

# --- Hierarchical clustering ---
final_sign <- as.matrix(wide_sign[, target_cols])[has_star, ]
row.clus <- hclust(dist(final_sign))

# --- 5. Highlighting ---
target_highlights <- c("Ruminococcus", "Fusicatenibacter", "Ruminococcaceae_Incertae_Sedis", "Lachnospiraceae_UCG-004")
highlight_features <- rownames(final_asso)[rownames(final_asso) %in% target_highlights]

# --- 6. Styling & Legend ---
FONT <- "Arial"
lim_fc  <- 1.5
col_fun = colorRamp2(c(-lim_fc, 0, lim_fc), c("blue", "white", "red"))

lgd_asso = Legend(col_fun = col_fun, title = "log2 fold change (vs W0)", at = c(-lim_fc, 0, lim_fc), labels = c("\u2264 -1.5", "0", "\u2265 1.5"), 
                  direction = "horizontal", title_position = "leftcenter", title_gap = unit(4, "mm"),
                  legend_width = unit(25, "mm"), grid_height = unit(2.5, "mm"),
                  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT),
                  labels_gp = gpar(fontsize = 7, fontfamily = FONT))

# --- 7. Plot Object Definition ---
ht <- Heatmap(final_asso, 
              name = "Association", 
              col = col_fun, 
              show_heatmap_legend = FALSE, 
              
              cluster_rows = row.clus,
              show_row_dend = TRUE,
              row_dend_side = "left",
              row_dend_width = unit(8, "mm"),

              cluster_columns = FALSE,
              border = TRUE,
              
              column_split = factor(rep(c("Placebo", "C. asiatica", "Fiber-added"), each = 2),
                                    levels = c("Placebo", "C. asiatica", "Fiber-added")),
              column_labels = rep(c("W2", "W4"), 3),
              column_names_rot = 0, 
              column_names_gp = gpar(fontsize = 7, fontfamily = FONT),
              column_names_centered = TRUE,
              column_names_side = "bottom", 
              column_gap = unit(2, "mm"),
              column_title = NULL,
              
              rect_gp = gpar(col = "white", lwd = 0.5), 
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(final_stars[i, j], x, y, gp = gpar(fontsize = 8, fontfamily = FONT))
              },
              
              row_names_side = "right", 
              row_names_gp = gpar(fontsize = 7, fontface = "italic", fontfamily = FONT,
                                  fill = ifelse(rownames(final_asso) %in% highlight_features, "#90EE90", NA)),
              
              top_annotation = HeatmapAnnotation(
                title = anno_block(
                  align_to = list("Placebo" = 1:2, "C. asiatica" = 3:4, "Fiber-added" = 5:6),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = "white", col = "black", lwd = 0.5))
                    if(nm == "Placebo") {
                      grid.text("Placebo", gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT))
                    } else if(nm == "C. asiatica") {
                      grid.text(expression(bold(italic("C. asiatica"))), gp = gpar(fontsize = 7, fontfamily = FONT))
                    } else {
                      grid.text("Fiber-enriched", y = unit(0.72, "npc"), gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT))
                      grid.text(expression(bold(italic("C. asiatica"))), y = unit(0.28, "npc"), gp = gpar(fontsize = 7, fontfamily = FONT))
                    }
                    
                  }, height = unit(6, "mm")
                ),
                show_annotation_name = FALSE
              ),
              
              bottom_annotation = HeatmapAnnotation(
                spacer = anno_empty(border = FALSE, height = unit(0.5, "mm"))
              )
)

# --- 8. Save Files ---
# Panel size as placed in Figures_EDITED pptx at 180 mm print width
W <- 126.4; H <- 120
draw_fig <- function() draw(ht, padding = unit(c(0.5, 0.5, 0.5, 1), "mm"), annotation_legend_list = list(lgd_asso),
                            heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
png(file.path(RESULT_DIR, "Figure3A.png"), width = W, height = H, units = "mm", res = 600)
draw_fig()
dev.off()

svglite(file.path(RESULT_DIR, "Figure3A.svg"), width = W / 25.4, height = H / 25.4,
        system_fonts = list(Arial = "Arial"))
draw_fig()
dev.off()

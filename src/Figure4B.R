# --- 1. Load Libraries ---
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)

# --- 2. Clear Graphics ---
while (!is.null(dev.list())) dev.off()
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 3. Read & Process Data ---
process_pval <- function(path, treat_name) {
  if(!file.exists(path)) return(NULL)
  df <- read.delim(path, check.names = FALSE)
  df$Time2_Time1 <- paste0(df$value, "/W0")
  df$Treatments <- treat_name
  return(df %>% select(feature, pval, qval, coef, Time2_Time1, Treatments))
}

MAASLIN_DIR <- file.path(OUTPUT_DIR, "Met")

Placebo_pval <- process_pval(file.path(MAASLIN_DIR, "output_Placebo/all_results.tsv"), "Placebo")
Asiatica_pval <- process_pval(file.path(MAASLIN_DIR, "output_Asiatica/all_results.tsv"), "Asiatica")
Fiber_added_pval <- process_pval(file.path(MAASLIN_DIR, "output_Fiber_added/all_results.tsv"), "Fiber_added")
stopifnot(!is.null(Placebo_pval), !is.null(Asiatica_pval), !is.null(Fiber_added_pval))
pval_total <- rbind(Placebo_pval, Asiatica_pval, Fiber_added_pval)

Change_feature_name <- data.frame(
  feature = c("total.Lipid", "ω6.ω3.ratio", "CE.20.5.", "MUFA.SFA.ratio", "PUFA.SFA.ratio", "CE.16.0.", 
              "CE.20.4.", "CE.17.0.", "CE.18.0.", "CE.20.0.", "CE.20.2.", "CE.18.1.", "Acetic.acid", 
              "Butanoic.acid", "total.CE", "CE.18.2.", "Pentanoic.acid", "CE.16.1.", "X2.Methylhexanoic.acid", 
              "CE.15.0.", "X2.Methypropionic.acid", "X3.Methylbutanoic.acid", "X2.Methylbutanoic.acid",
              "p.Cresolsulfate","Phenysulfate","Indoxyl.sulfate","X2.Methoxyphenylsulfate",
              "Hippuric.acid","Cinnamoylglycine"," m.Hydroxyhippuric.acid",
              "Phenylacetyl.L.glutamine","Trimethylamine.N.oxide","Dimethylglycine"),
  new_feature = c("Total Lipid", "omega6/omega3 ratio", "CE(20:5)", "MUFA/SFA ratio", "PUFA/SFA ratio", "CE(16:0)", 
                  "CE(20:4)", "CE(17:0)", "CE(18:0)", "CE(20:0)", "CE(20:2)", "CE(18:1)", "Acetic acid", 
                  "Butanoic acid", "Total CE", "CE(18:2)", "Propionic acid", "CE(16:1)", "2-Methylhexanoic acid", 
                  "CE(15:0)", "2-Methylpropionic acid", "3-Methylbutanoic acid", "2-Methylbutanoic acid",
                  "p-Cresolsulfate","Phenysulfate","Indoxyl sulfate","2-Methoxyphenylsulfate",
                  "Hippuric acid","Cinnamoylglycine", " m-Hydroxyhippuric acid",
                  "Phenylacetyl-L-glutamine","Trimethylamine N-oxide","Dimethylglycine")
)

pval_total <- pval_total %>%
  left_join(Change_feature_name, by = "feature") %>%
  mutate(feature = ifelse(is.na(new_feature), as.character(feature), as.character(new_feature)))

# --- 4. Custom Order ---
custom_order <- c(
  "Acetic acid", "Propionic acid", "2-Methylpropionic acid", "Butanoic acid", "3-Methylbutanoic acid", "2-Methylbutanoic acid", "2-Methylhexanoic acid", 
  "Total CE", "Cholesterol", "CE(15:0)", "CE(16:0)", "CE(16:1)", "CE(17:0)", "CE(18:0)", "CE(18:1)", "CE(18:2)", "CE(20:0)", "CE(20:2)", "CE(20:4)", "CE(20:5)", 
  "Total Lipid", "SFA", "MUFA", "PUFA", "MUFA/SFA ratio", "PUFA/SFA ratio", "omega6/omega3 ratio", 
  "p-Cresolsulfate","Phenysulfate","Indoxyl sulfate","2-Methoxyphenylsulfate","Hippuric acid","Cinnamoylglycine", "Phenylacetyl-L-glutamine","Trimethylamine N-oxide","Dimethylglycine"
)

# --- 5. Matrix Construction ---
target_cols <- c("Placebo_W2/W0", "Placebo_W4/W0", "Asiatica_W2/W0", "Asiatica_W4/W0", "Fiber_added_W2/W0", "Fiber_added_W4/W0")
pval_matrix_data <- pval_total %>%
  mutate(Treatments_Weeks = paste0(Treatments, "_", Time2_Time1),
         coef = as.numeric(coef)) %>%
  filter(Treatments_Weeks %in% target_cols) %>%
  group_by(feature, Treatments_Weeks) %>%
  summarise(coef = mean(coef), pval = min(pval), qval = min(qval), .groups = "drop")

matrix_asso <- pval_matrix_data %>%
  select(feature, Treatments_Weeks, coef) %>%
  pivot_wider(names_from = Treatments_Weeks, values_from = coef, values_fill = 0) %>%
  filter(feature %in% custom_order) %>%
  arrange(match(feature, custom_order)) %>%
  as.data.frame() %>%
  { .rows <- .$feature; .mat <- as.matrix(.[, target_cols]); rownames(.mat) <- .rows; .mat }

matrix_stars <- pval_matrix_data %>%
  select(feature, Treatments_Weeks, qval) %>%
  pivot_wider(names_from = Treatments_Weeks, values_from = qval) %>%
  filter(feature %in% custom_order) %>%
  arrange(match(feature, custom_order)) %>%
  as.data.frame() %>%
  { .mat <- as.matrix(.[, target_cols]); ifelse(.mat < 0.01, "**", ifelse(.mat < 0.05, "*", "")) }

missing_feat <- setdiff(custom_order, rownames(matrix_asso))
if (length(missing_feat) > 0)
  stop("custom_order names not found after renaming: ",
       paste(missing_feat, collapse = ", "),
       "\nCheck them against Change_feature_name$new_feature.")

# --- 6. Grouping Logic ---
n_scfa <- 7
n_ce   <- 13
n_fa   <- 7
n_ut   <- 9
total_rows <- nrow(matrix_asso)

minor_split <- factor(c(rep("SCFAs", n_scfa), rep("CEs", n_ce), rep("FAs", n_fa), rep("UTs", n_ut)), 
                      levels = c("SCFAs", "CEs", "FAs", "UTs"))

highlight_features <- c(
  "2-Methylhexanoic acid", "2-Methylpropionic acid", "3-Methylbutanoic acid", 
  "2-Methylbutanoic acid", "Cholesterol", "Total Lipid", "SFA", "MUFA", 
  "PUFA", "MUFA/SFA ratio", "PUFA/SFA ratio", "omega6/omega3 ratio", "Dimethylglycine"
)

FONT <- "Arial"
lim_fc  <- 1.5
col_fun = colorRamp2(c(-lim_fc, 0, lim_fc), c("blue", "white", "red"))
lgd_asso = Legend(col_fun = col_fun, 
                  title = "log2 fold change (vs W0)", 
                  at = c(-lim_fc, 0, lim_fc), labels = c("\u2264 -1.5", "0", "\u2265 1.5"), 
                  direction = "horizontal", title_position = "leftcenter", title_gap = unit(3, "mm"),
                  legend_width = unit(25, "mm"), grid_height = unit(2.5, "mm"),
                  title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT),
                  labels_gp = gpar(fontsize = 7, fontfamily = FONT))

# --- 8. Final Plotting ---

ht <- Heatmap(matrix_asso,
              name = "Association",
              col = col_fun,
              show_heatmap_legend = FALSE, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = minor_split,
              row_title = NULL,
              row_gap = unit(1.2, "mm"),
              border = TRUE,
              
              column_split = factor(rep(c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica"), each = 2),
                                    levels = c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")),
              column_labels = rep(c("W2", "W4"), 3),
              column_names_rot = 0,
              column_names_gp = gpar(fontsize = 7, fontfamily = FONT),
              column_names_centered = TRUE,
              column_names_side = "bottom",
              column_gap = unit(1.2, "mm"),
              show_column_names = TRUE,
              column_title = NULL,
              
              bottom_annotation = HeatmapAnnotation(
                spacer = anno_empty(border = FALSE, height = unit(0.5, "mm"))
              ),
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(matrix_stars[i, j], x, y, gp = gpar(fontsize = 8, fontfamily = FONT))
              },
              
              rect_gp = gpar(col = "white", lwd = 0.5),
              
              row_names_side = "right", 
              row_names_gp = gpar(
                fontsize = 7, fontfamily = FONT,
                fill = ifelse(rownames(matrix_asso) %in% highlight_features, "#90EE90", NA)
              ),
              
              left_annotation = rowAnnotation(
                Major = anno_block(
                  align_to = list("Feces" = 1:n_scfa, "Plasma" = (n_scfa + 1):total_rows),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.75))
                    grid.text(nm, rot = 90, gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT))
                  }, width = unit(3.5, "mm")
                ),
                Minor = anno_block(
                  align_to = list("SCFAs" = 1:n_scfa, "CEs" = (n_scfa + 1):(n_scfa + n_ce),
                                  "FAs" = (n_scfa + n_ce + 1):(n_scfa + n_ce + n_fa),
                                  "UTs" = (n_scfa + n_ce + n_fa + 1):total_rows),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.75))
                    grid.text(nm, rot = 90, gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT))
                  }, width = unit(3.5, "mm")
                )
              ),
              
              top_annotation = HeatmapAnnotation(
                title = anno_block(
                  align_to = list("Placebo" = 1:2, "C. asiatica" = 3:4, "Fiber-enriched C. asiatica" = 5:6),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = "white", col = "black", lwd = 0.5))
                    gp <- gpar(fontsize = 7, fontface = "bold", fontfamily = FONT)
                    if(nm == "Placebo") {
                      grid.text("Placebo", gp = gp)
                    } else if(nm == "C. asiatica") {
                      grid.text(expression(bold(italic("C. asiatica"))), gp = gp)
                    } else {
                      grid.text("Fiber-", y = unit(0.8, "npc"), gp = gp)
                      grid.text("enriched", y = unit(0.5, "npc"), gp = gp)
                      grid.text(expression(bold(italic("C. asiatica"))), y = unit(0.2, "npc"), gp = gp)
                    }
                  }, height = unit(9, "mm")
                )
              )
)


# --- 9. Save --------------------------------------------------------------
# Panel size as placed in Figures_EDITED pptx (13-inch-high slides) at 180 mm print width
W <- 80; H <- 134
draw_fig <- function() draw(ht, padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"), annotation_legend_list = list(lgd_asso),
                            heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
png(file.path(RESULT_DIR, "Figure4B.png"), width = W, height = H, units = "mm", res = 600)
draw_fig()
dev.off()

svglite(file.path(RESULT_DIR, "Figure4B.svg"), width = W / 25.4, height = H / 25.4,
        system_fonts = list(Arial = "Arial"))
draw_fig()
dev.off()

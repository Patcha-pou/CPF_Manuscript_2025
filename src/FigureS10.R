# ==============================================================================
# SCRIPT: Figure S10 - All-sample mixed-model microbe x metabolite associations
# ==============================================================================
pacman::p_load(openxlsx, dplyr, ComplexHeatmap, circlize, grid, ggplot2, patchwork, svglite)

# --- 1. Paths -----------------------------------------------------------------
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
TABLE_S5   <- file.path(OUTPUT_DIR, "Table_S5.xlsx")

P_CUT <- 0.05
Q_CUT <- 0.05
FONT  <- "Arial"
gp_txt <- function(size, face = "plain") gpar(fontsize = size, fontface = face, fontfamily = FONT)

# --- 2. Pairs with p < 0.05 across all samples (Table S5) ---------------------
rm_all <- read.xlsx(TABLE_S5, sheet = "Table_S5_All_samples")
sel <- rm_all %>% filter(p < P_CUT) %>% arrange(Metabolite, desc(beta))
stopifnot(nrow(sel) > 0)

mat_r <- sel %>%
  transmute(All = beta, CON = beta_CON, GKJ = beta_GKJ, `F-GKJ` = beta_FGKJ,
            `W0-W2` = beta_W2, `W0-W4` = beta_W4) %>%
  as.matrix()
rownames(mat_r) <- paste(sel$Microbe, sel$Metabolite, sep = " | ")
mat_star <- cbind(ifelse(sel$p.adj < Q_CUT, "*", ""), matrix("", nrow(sel), 5))

lim <- ceiling(max(abs(mat_r), na.rm = TRUE) * 10) / 10
col_fun <- colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "Standardized β", at = c(-lim, 0, lim),
              direction = "horizontal", legend_width = unit(25, "mm"), grid_height = unit(2.5, "mm"),
              title_position = "topcenter", title_gp = gp_txt(8, "bold"), labels_gp = gp_txt(7))

# --- 3. Heatmap ---------------------------------------------------------------
met_lv <- unique(sel$Metabolite)
wrap_met <- function(x) sub("Methyl", "Methyl-\n", x)
ht <- Heatmap(mat_r, name = "beta", col = col_fun, show_heatmap_legend = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,
              row_split = factor(wrap_met(sel$Metabolite), levels = wrap_met(met_lv)), row_gap = unit(1, "mm"),
              row_title_rot = 0, row_title_side = "left", row_title_gp = gp_txt(7, "bold"),
              column_split = factor(c("All", rep("Group", 3), rep("Time", 2)), levels = c("All", "Group", "Time")),
              column_title = NULL, column_gap = unit(1.5, "mm"),
              column_names_side = "top", column_names_rot = 90,
              column_names_gp = gp_txt(7, "bold"),
              rect_gp = gpar(col = "white", lwd = 0.5),
              cell_fun = function(j, i, x, y, width, height, fill)
                grid.text(mat_star[i, j], x, y - unit(0.4, "mm"), gp = gp_txt(9, "bold")),
              row_labels = sel$Microbe, row_names_side = "right", row_names_gp = gp_txt(7, "italic"),
              width = unit(6 * 6, "mm"), height = unit(nrow(mat_r) * 2.9, "mm"))
p <- wrap_elements(full = grid.grabExpr(
  draw(ht, annotation_legend_list = list(lgd), annotation_legend_side = "bottom",
       padding = unit(c(2, 2, 2, 2), "mm"))))

# --- 4. Save ------------------------------------------------------------------
W <- 180 / 25.4
H <- (nrow(mat_r) * 2.9 + length(met_lv) * 1 + 30) / 25.4
ggsave(file.path(RESULT_DIR, "Fig_S10.png"), p, width = W, height = H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S10.svg"), p, width = W, height = H, device = "svg")
message("Saved Fig_S10.{png,svg}")

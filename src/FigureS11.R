# ==============================================================================
# SCRIPT: Figure S11 - Between-group differences in microbe x metabolite associations
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

# --- 2. Pairs with interaction p < 0.05 (Table S5) ----------------------------
int_all <- read.xlsx(TABLE_S5, sheet = "Table_S5_Interaction")
sel <- int_all %>% filter(p_interaction < P_CUT) %>% arrange(Metabolite, p_interaction)
stopifnot(nrow(sel) > 0)

mat_d <- sel %>%
  transmute(`GKJ vs CON` = diff_GKJ_vs_CON, `F-GKJ vs CON` = diff_FGKJ_vs_CON,
            `F-GKJ vs GKJ` = diff_FGKJ_vs_GKJ) %>%
  as.matrix()
mat_q <- sel %>% transmute(q_GKJ_vs_CON, q_FGKJ_vs_CON, q_FGKJ_vs_GKJ) %>% as.matrix()
mat_star <- ifelse(mat_q < Q_CUT, "*", "")

lim <- ceiling(max(abs(mat_d), na.rm = TRUE) * 10) / 10
col_fun <- colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "Difference in standardized β", at = c(-lim, 0, lim),
              direction = "horizontal", legend_width = unit(25, "mm"), grid_height = unit(2.5, "mm"),
              title_position = "topcenter", title_gp = gp_txt(8, "bold"), labels_gp = gp_txt(7))

# --- 3. Heatmap ---------------------------------------------------------------
met_lv <- unique(sel$Metabolite)
wrap_met <- function(x) sub("Methyl", "Methyl-\n", x)
ht <- Heatmap(mat_d, name = "diff", col = col_fun, show_heatmap_legend = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,
              row_split = factor(wrap_met(sel$Metabolite), levels = wrap_met(met_lv)), row_gap = unit(1, "mm"),
              row_title_rot = 0, row_title_side = "left", row_title_gp = gp_txt(7, "bold"),
              column_names_side = "top", column_names_rot = 90,
              column_names_gp = gp_txt(7, "bold"),
              rect_gp = gpar(col = "white", lwd = 0.5),
              cell_fun = function(j, i, x, y, width, height, fill)
                grid.text(mat_star[i, j], x, y - unit(0.4, "mm"), gp = gp_txt(9, "bold")),
              row_labels = sel$Microbe, row_names_side = "right", row_names_gp = gp_txt(7, "italic"),
              width = unit(3 * 7, "mm"), height = unit(nrow(mat_d) * 2.9, "mm"))
p <- wrap_elements(full = grid.grabExpr(
  draw(ht, annotation_legend_list = list(lgd), annotation_legend_side = "bottom",
       padding = unit(c(2, 2, 2, 2), "mm"))))

# --- 4. Save ------------------------------------------------------------------
W <- 180 / 25.4
H <- (nrow(mat_d) * 2.9 + length(met_lv) * 1 + 40) / 25.4
ggsave(file.path(RESULT_DIR, "Fig_S11.png"), p, width = W, height = H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S11.svg"), p, width = W, height = H, device = "svg")
message("Saved Fig_S11.{png,svg}")

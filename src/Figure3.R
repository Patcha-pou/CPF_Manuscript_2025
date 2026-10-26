# ==============================================================================
# SCRIPT: Figure 3 (A + B assembled) - genus-level response to the intervention
# (A) Within-arm MaAsLin2 heatmap (as Figure3A.R); (B) medians of highlighted genera (as Figure3B.R)
# ==============================================================================
pacman::p_load(dplyr, tidyr, tibble, readxl, missForest, ComplexHeatmap, circlize,
               grid, ggplot2, patchwork, svglite, ragg)

# --- 1. Paths -----------------------------------------------------------------
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

FONT   <- "Arial"
gp_txt <- function(size, face = "plain", ...) gpar(fontsize = size, fontface = face, fontfamily = FONT, ...)
PAL <- c("Placebo" = "#999999", "C. asiatica" = "#e6c60d",
         "Fiber-enriched C. asiatica" = "#5ca08e")

# --- 2. Panel A: within-arm MaAsLin2 heatmap ---------------------------------
TAXA <- read_excel(file.path(INPUT_DIR, "Downstream_analysis.xlsx"),
                   sheet = "nonRedundantTAXA_filtered")[, 1:2]
colnames(TAXA) <- c("feature", "Shortform")
TAXA$Shortform <- sub("^genus.*?:", "", TAXA$Shortform)
TAXA <- TAXA[!duplicated(TAXA$feature), ]

read_arm <- function(arm) {
  read.delim(file.path(OUTPUT_DIR, "Mic", paste0("output_", arm), "all_results.tsv"),
             check.names = FALSE) %>%
    mutate(Treatments_Weeks = paste0(arm, "_", value, "/W0")) %>%
    select(feature, Treatments_Weeks, pval, qval, coef)
}
target_cols <- c("Placebo_W2/W0", "Placebo_W4/W0", "Asiatica_W2/W0", "Asiatica_W4/W0",
                 "Fiber_added_W2/W0", "Fiber_added_W4/W0")
res <- bind_rows(lapply(c("Placebo", "Asiatica", "Fiber_added"), read_arm)) %>%
  filter(Treatments_Weeks %in% target_cols) %>%
  mutate(sign_asso = -log10(pval) * sign(coef)) %>%
  group_by(feature, Treatments_Weeks) %>%
  summarise(coef = mean(coef), sign_asso = mean(sign_asso),
            pval = min(pval), qval = min(qval), .groups = "drop")

sig_features <- res %>% filter(pval < 0.05) %>% pull(feature) %>% unique()
wide_of <- function(v, fill) {
  res %>% filter(feature %in% sig_features) %>%
    select(feature, Treatments_Weeks, all_of(v)) %>%
    pivot_wider(names_from = Treatments_Weeks, values_from = all_of(v), values_fill = fill) %>%
    column_to_rownames("feature") %>% as.matrix() %>% .[, target_cols]
}
mat_fc   <- wide_of("coef", 0)
mat_sign <- wide_of("sign_asso", 0)
mat_q    <- wide_of("qval", 1)
labs_row <- TAXA$Shortform[match(rownames(mat_fc), TAXA$feature)]
labs_row[is.na(labs_row)] <- rownames(mat_fc)[is.na(labs_row)]
rownames(mat_fc) <- rownames(mat_sign) <- rownames(mat_q) <- make.unique(labs_row, sep = " ")
mat_star <- ifelse(mat_q < 0.01, "**", ifelse(mat_q < 0.05, "*", ""))
mat_star[is.na(mat_star)] <- ""

row_clus <- hclust(dist(mat_sign))
G_SIZE   <- c(G1 = 12, G2 = 9, G3 = 13)
stopifnot(sum(G_SIZE) == nrow(mat_sign))
G_END    <- cumsum(G_SIZE)
HIGHLIGHT <- c("Ruminococcus", "Fusicatenibacter", "Ruminococcaceae_Incertae_Sedis",
               "Lachnospiraceae_UCG-004")

lim_fc  <- 1.5
col_fun <- colorRamp2(c(-lim_fc, 0, lim_fc), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "log2 fold change\n(vs W0)",
              at = c(-lim_fc, 0, lim_fc), labels = c("≤ -1.5", "0", "≥ 1.5"),
              direction = "horizontal", legend_width = unit(30, "mm"), grid_height = unit(3, "mm"),
              title_position = "leftcenter",
              title_gp = gp_txt(8, "bold"), labels_gp = gp_txt(7))

ARM_LAB <- c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")
ht <- Heatmap(mat_fc, name = "log2FC", col = col_fun, show_heatmap_legend = FALSE,
              cluster_rows = row_clus, show_row_dend = TRUE, row_dend_side = "left",
              row_dend_width = unit(8, "mm"),
              cluster_columns = FALSE, border = FALSE,
              column_split = factor(rep(ARM_LAB, each = 2), levels = ARM_LAB),
              column_title = NULL, column_gap = unit(1.5, "mm"),
              column_labels = rep(c("W2 vs W0", "W4 vs W0"), 3),
              column_names_rot = 45, column_names_centered = TRUE,
              column_names_side = "bottom", column_names_gp = gp_txt(7),
              rect_gp = gpar(col = "white", lwd = 0.5),
              cell_fun = function(j, i, x, y, width, height, fill)
                grid.text(mat_star[i, j], x, y, gp = gp_txt(8)),
              row_names_side = "right",
              row_labels = paste0(rownames(mat_fc), "   "),
              row_names_gp = gp_txt(7, "italic",
                                    fill = ifelse(rownames(mat_fc) %in% HIGHLIGHT, "#90EE90", NA)),
              width = unit(6 * 10, "mm"), height = unit(nrow(mat_fc) * 3.4, "mm"),
              left_annotation = rowAnnotation(
                grp = anno_empty(border = FALSE, width = unit(6, "mm")),
                show_annotation_name = FALSE),
              top_annotation = HeatmapAnnotation(
                arm = anno_block(
                  align_to = list("Placebo" = 1:2, "C. asiatica" = 3:4,
                                  "Fiber-enriched C. asiatica" = 5:6),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = "white", col = "black"))
                    if (nm == "Placebo") {
                      grid.text("Placebo", gp = gp_txt(7, "bold"))
                    } else if (nm == "C. asiatica") {
                      grid.text("C. asiatica", gp = gp_txt(7, "bold.italic"))
                    } else {
                      grid.text("Fiber-enriched", y = unit(0.68, "npc"), gp = gp_txt(7, "bold"))
                      grid.text("C. asiatica", y = unit(0.3, "npc"), gp = gp_txt(7, "bold.italic"))
                    }
                  }, height = unit(7, "mm")),
                show_annotation_name = FALSE))
n_row  <- nrow(mat_fc)
g_y    <- function(k) {
  top <- if (k == 1) 0 else G_END[k - 1]
  c(bottom = 1 - G_END[k] / n_row, top = 1 - top / n_row)
}
draw_groups <- function() {
  for (cs in 1:3) decorate_heatmap_body("log2FC", column_slice = cs, {
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1))
    for (k in 1:2) grid.lines(c(0, 1), rep(1 - G_END[k] / n_row, 2),
                              gp = gpar(col = "black", lwd = 1))
  })
  decorate_annotation("grp", {
    for (k in seq_along(G_SIZE)) {
      yy <- g_y(k)
      grid.rect(y = unit(mean(yy), "npc"), height = unit(diff(rev(yy)), "npc"),
                gp = gpar(fill = "black", col = "white", lwd = 1))
      grid.text(names(G_SIZE)[k], y = unit(mean(yy), "npc"), gp = gp_txt(8, "bold", col = "white"))
    }
  })
}
# --- 3. Panel B: median relative abundance of highlighted genera -------------
MET_MIC  <- file.path(INPUT_DIR, "inputMetMic.xlsx")
METADATA <- read_excel(MET_MIC, sheet = "METADATA") %>%
  mutate(Treatments = factor(case_when(
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched C. asiatica",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "C. asiatica",
    TRUE ~ "Placebo"), levels = names(PAL)))
Mic <- read_excel(MET_MIC, sheet = "Mic")
mic <- Mic %>% select(-ID) %>% as.data.frame()
mic[] <- lapply(mic, function(x) as.numeric(as.character(x)))
set.seed(123)
mic <- missForest(mic)$ximp
rownames(mic) <- Mic$ID

TARGET <- c("genus103:Ruminococcus" = "Ruminococcus",
            "genus63:Fusicatenibacter" = "Fusicatenibacter",
            "genus100:Ruminococcaceae_Incertae_Sedis" = "Ruminococcaceae Incertae Sedis",
            "genus73:Lachnospiraceae_UCG-004" = "Lachnospiraceae UCG-004")
dB <- mic %>% select(all_of(names(TARGET))) %>%
  rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "Taxa", values_to = "Value") %>%
  left_join(METADATA %>% select(ID, Weeks, Treatments), by = "ID") %>%
  mutate(Weeks = as.numeric(as.character(Weeks)),
         Taxa = factor(TARGET[Taxa], levels = TARGET))

plot_one <- function(tx, show_x) {
  ggplot(dB %>% filter(Taxa == tx), aes(Weeks, Value, colour = Treatments, group = Treatments)) +
    stat_summary(fun = median, geom = "line", linewidth = 0.7) +
    stat_summary(fun = median, geom = "point", size = 1.3) +
    facet_wrap(~ Taxa, labeller = label_wrap_gen(width = 18)) +
    scale_colour_manual(values = PAL) +
    scale_x_continuous(breaks = c(0, 2, 4)) +
    labs(x = if (show_x) "Weeks" else NULL, y = NULL) +
    theme_bw(base_size = 7, base_family = FONT) +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold.italic", size = 7, margin = margin(1, 1, 1, 1)),
          axis.title = element_text(size = 8, face = "bold"),
          axis.text = element_text(size = 7, colour = "black"),
          axis.text.x = if (show_x) element_text(size = 7, colour = "black") else element_blank(),
          panel.grid.minor = element_blank(), plot.margin = margin(1, 2, 1, 1))
}
pB <- setNames(lapply(seq_along(TARGET), function(k) plot_one(TARGET[k], k == length(TARGET))),
               TARGET)
get_legend_B <- function() {
  g <- ggplotGrob(ggplot(dB, aes(Weeks, Value, colour = Treatments)) + geom_line() +
                    scale_colour_manual(values = PAL, labels = c(
                      "Placebo" = "Placebo",
                      "C. asiatica" = expression(italic("C. asiatica")),
                      "Fiber-enriched C. asiatica" = expression("Fiber-enriched"~italic("C. asiatica")))) +
                    guides(colour = guide_legend(ncol = 1)) +
                    theme_bw(base_size = 7, base_family = FONT) +
                    theme(legend.title = element_blank(), legend.text = element_text(size = 7),
                          legend.key.height = unit(3, "mm")))
  idx <- grep("guide-box", g$layout$name)
  Filter(function(x) !inherits(x, "zeroGrob"), g$grobs[idx])[[1]]
}

# --- 4. Page layout (mm) with dashed connectors -------------------------------
W <- 180; H <- 160
X_A <- c(2, 118); X_B <- c(134, 179)
B_TOP <- 14; B_BOT <- 132
vp_mm <- function(x, y_top, w, h) viewport(x = unit(x, "mm"), y = unit(H - y_top, "mm"),
                                           width = unit(w, "mm"), height = unit(h, "mm"),
                                           just = c("left", "top"))
loc_in <- function(x, y) { l <- deviceLoc(x, y); c(x = as.numeric(l$x), y = as.numeric(l$y)) }
lab_ord <- rownames(mat_fc)[row_clus$order]

draw_page <- function() {
  grid.newpage()
  pushViewport(vp_mm(X_A[1], 4, diff(X_A), H - 6))
  draw(ht, newpage = FALSE, annotation_legend_list = list(lgd),
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       padding = unit(c(1, 1, 1, 1), "mm"))
  draw_groups()
  src <- list()
  decorate_row_names("log2FC", {
    for (tx in names(TARGET)) {
      nm <- TARGET[[tx]]; key <- gsub(" ", "_", nm)
      i  <- match(key, lab_ord)
      src[[nm]] <<- rbind(top    = loc_in(unit(1, "npc"), unit(1 - (i - 1) / n_row, "npc")),
                          bottom = loc_in(unit(1, "npc"), unit(1 - i / n_row, "npc")))
    }
  })
  upViewport()
  # B: title, four panels, legend
  grid.text("Dynamic change of microbe\nabundance (median)",
            x = unit(mean(X_B), "mm"), y = unit(H - 5, "mm"), just = "top",
            gp = gp_txt(8, "bold"))
  grid.text("Relative abundance", x = unit(X_B[1] - 5, "mm"),
            y = unit(H - (B_TOP + B_BOT) / 2, "mm"), rot = 90, gp = gp_txt(8, "bold"))
  hk <- (B_BOT - B_TOP) / length(pB)
  tgt <- list()
  for (k in seq_along(pB)) {
    h <- hk + if (k == length(pB)) 6 else 0
    pushViewport(vp_mm(X_B[1], B_TOP + (k - 1) * hk, diff(X_B), h))
    print(pB[[k]], newpage = FALSE)
    tgt[[names(pB)[k]]] <- rbind(
      top    = loc_in(unit(0, "npc"), unit(1, "npc")),
      bottom = loc_in(unit(0, "npc"), unit(1, "npc") - unit(hk, "mm")))
    upViewport()
  }
  pushViewport(vp_mm(X_B[1], B_BOT + 8, diff(X_B), 16)); grid.draw(get_legend_B()); upViewport()
  for (nm in names(tgt)) for (e in c("top", "bottom"))
    grid.lines(unit(c(src[[nm]][e, "x"], tgt[[nm]][e, "x"]), "inches"),
               unit(c(src[[nm]][e, "y"], tgt[[nm]][e, "y"]), "inches"),
               gp = gpar(lty = "22", col = "grey45", lwd = 0.5))
  grid.text("A", x = unit(1, "mm"), y = unit(H - 2, "mm"), just = c("left", "top"), gp = gp_txt(9, "bold"))
  grid.text("B", x = unit(X_B[1] - 8, "mm"), y = unit(H - 2, "mm"), just = c("left", "top"), gp = gp_txt(9, "bold"))
}

agg_png(file.path(RESULT_DIR, "Figure3.png"), width = W, height = H, units = "mm", res = 600)
draw_page(); invisible(dev.off())
svglite(file.path(RESULT_DIR, "Figure3.svg"), width = W / 25.4, height = H / 25.4)
draw_page(); invisible(dev.off())
message("Saved Figure3.{png,svg} (", W, " x ", H, " mm)")

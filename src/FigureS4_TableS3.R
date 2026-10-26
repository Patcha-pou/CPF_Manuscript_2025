# ==============================================================================
# SCRIPT: Figure S4 + Table S3 - Genus-level Group x Time MaAsLin2 models (DiD)
# ==============================================================================
pacman::p_load(readxl, dplyr, tidyr, tibble, Maaslin2, ComplexHeatmap, circlize,
               grid, ggplot2, patchwork, svglite, openxlsx)

# --- 1. Paths -----------------------------------------------------------------
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

P_SHOW <- 0.05
Q_CUT  <- 0.05

GROUP_LV  <- c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")
GROUP_MAP <- c(Placebo = "Placebo", Asiatica = "C. asiatica",
               Fiber_added = "Fiber-enriched C. asiatica")
ARM_ID <- c("Placebo" = "Pla", "C. asiatica" = "CA", "Fiber-enriched C. asiatica" = "FCA")
PAL <- c("Placebo" = "#999999", "C. asiatica" = "#e6c60d",
         "Fiber-enriched C. asiatica" = "#5ca08e")
COL_LV <- c("CA-Pla W2", "CA-Pla W4", "FCA-Pla W2", "FCA-Pla W4", "FCA-CA W2", "FCA-CA W4")

# --- 2. Metadata (covariates as in Maaslin2.R: Age, BMI, Gender, Energy) -----
meta <- read.csv(file.path(INPUT_DIR, "input_metadata.csv"), row.names = 1)
mf0 <- meta %>%
  rownames_to_column("SampleID") %>%
  transmute(SampleID, Subject, Group = GROUP_MAP[Treatments], Time = Weeks,
            Age, BMI, Gender, Energy)
stopifnot(complete.cases(mf0))

# --- 3. Genus table ----------------------------------------------------------
rel <- read.csv(file.path(INPUT_DIR, "input_data_rel.csv"), row.names = 1, check.names = FALSE)
rel <- rel[mf0$SampleID, ]

TAXA <- read_excel(file.path(INPUT_DIR, "Downstream_analysis.xlsx"),
                   sheet = "nonRedundantTAXA_filtered") %>%
  transmute(feature = Feature,
            label   = sub("^genus.*?:", "", ShortForm),
            Phylum  = sub("^p__", "", Phylum))

# --- 4. MaAsLin2 with Group x Time dummies -----------------------------------
run_did <- function(data, ref, layer) {
  arms <- setdiff(GROUP_LV, ref)
  md <- data.frame(mf0[, setdiff(names(mf0), "SampleID")], row.names = mf0$SampleID)
  int_cols <- character(0)
  for (a in arms) for (wk in c("W2", "W4")) {
    nm <- paste0(ARM_ID[a], "_", wk)
    md[[nm]] <- as.numeric(md$Group == a & md$Time == wk)
    int_cols <- c(int_cols, nm)
  }
  out <- file.path(sub("/+$", "", OUTPUT_DIR), layer, paste0("ref_", ARM_ID[ref]))
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  fit <- Maaslin2(
    data.frame(data, check.names = FALSE), md, out,
    transform = "LOG", normalization = "NONE", standardize = FALSE,
    fixed_effects = c("Group", "Time", int_cols, "Age", "BMI", "Gender", "Energy"),
    random_effects = "Subject",
    reference = c(paste0("Group,", ref), "Time,W0"),
    plot_heatmap = FALSE, plot_scatter = FALSE)
  fit$results %>%
    filter(metadata %in% int_cols) %>%
    transmute(feature, contrast = paste0(sub("_.*", "", metadata), "-", ARM_ID[ref]),
              week = sub(".*_", "", metadata), coef, stderr, pval, qval, N)
}

did_all <- bind_rows(
  run_did(rel, "Placebo", "Mic_DiD"),
  run_did(rel, "C. asiatica", "Mic_DiD") %>% filter(contrast == "FCA-CA")) %>%
  mutate(lower.CL = coef - qnorm(0.975) * stderr, upper.CL = coef + qnorm(0.975) * stderr) %>%
  left_join(TAXA, by = "feature") %>%
  relocate(label, Phylum, .after = feature)

# --- Table S3 ---
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", halign = "center",
                  fgFill = "#EDEDED", border = "bottom")
add_sheet <- function(wb, name, df) {
  df <- as.data.frame(df)
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs, withFilter = TRUE)
  freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet(wb, "Genus_DiD", did_all %>%
  arrange(contrast, week, pval) %>%
  transmute(Genus = label, Phylum, Feature = feature, Contrast = contrast, Week = week,
            DiD_log2FC = coef, SE = stderr, CI95_lower = lower.CL, CI95_upper = upper.CL,
            p_value = pval, q_value = qval, N))
saveWorkbook(wb, file.path(OUTPUT_DIR, "Table_S3.xlsx"), overwrite = TRUE)
cat("Wrote", file.path(OUTPUT_DIR, "Table_S3.xlsx"), "\n")
message("Genera tested: ", n_distinct(did_all$feature),
        " | p<", P_SHOW, " in any contrast: ", n_distinct(did_all$feature[did_all$pval < P_SHOW]),
        " | q<", Q_CUT, ": ", n_distinct(did_all$feature[did_all$qval < Q_CUT]))

# --- 5. Figure (A: Figure 3A heatmap style, B: Figure 3B line style) --------
FONT <- "Arial"
gp_txt <- function(size, face = "plain") gpar(fontsize = size, fontface = face, fontfamily = FONT)
CONTRAST_LAB <- c("CA vs\nPlacebo", "FCA vs\nPlacebo", "FCA vs\nCA")

sel <- did_all %>% filter(pval < P_SHOW) %>% distinct(feature) %>% pull(feature)
stopifnot(length(sel) > 0)
wide_of <- function(v) {
  did_all %>% filter(feature %in% sel) %>%
    mutate(col_id = paste(contrast, week)) %>%
    select(label, col_id, all_of(v)) %>%
    pivot_wider(names_from = col_id, values_from = all_of(v)) %>%
    column_to_rownames("label") %>% as.matrix() %>% .[, COL_LV, drop = FALSE]
}
mat_coef <- wide_of("coef")
mat_p    <- wide_of("pval")
mat_q    <- wide_of("qval")
mat_star <- ifelse(mat_q < 0.01, "**", ifelse(mat_q < 0.05, "*", ""))
lim <- ceiling(max(abs(mat_coef)))
col_fun <- colorRamp2(c(-lim, 0, lim), c("blue", "white", "red"))

lgd <- Legend(col_fun = col_fun, title = "DiD (log2 FC)", at = c(-lim, 0, lim),
              legend_height = unit(25, "mm"), grid_width = unit(3, "mm"),
              title_gp = gp_txt(8, "bold"), labels_gp = gp_txt(7))
split <- factor(rep(CONTRAST_LAB, each = 2), levels = CONTRAST_LAB)
row_clus <- hclust(dist(mat_coef))
ht <- Heatmap(mat_coef, name = "DiD", col = col_fun, show_heatmap_legend = FALSE,
              cluster_rows = row_clus, show_row_dend = TRUE,
              row_dend_side = "left", row_dend_width = unit(10, "mm"),
              cluster_columns = FALSE, border = TRUE,
              column_split = split, column_title = NULL, column_gap = unit(1.5, "mm"),
              column_labels = rep(c("W2 vs W0", "W4 vs W0"), 3),
              column_names_rot = 45, column_names_centered = TRUE,
              column_names_side = "bottom", column_names_gp = gp_txt(7),
              rect_gp = gpar(col = "white", lwd = 0.5),
              cell_fun = function(j, i, x, y, width, height, fill)
                grid.text(mat_star[i, j], x, y, gp = gp_txt(9)),
              row_names_side = "right", row_names_gp = gp_txt(7, "italic"),
              width = unit(6 * 11, "mm"), height = unit(nrow(mat_coef) * 3.2, "mm"),
              top_annotation = HeatmapAnnotation(
                title = anno_block(
                  align_to = setNames(list(1:2, 3:4, 5:6), CONTRAST_LAB),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = "white", col = "black"))
                    grid.text(nm, gp = gp_txt(8, "bold"))
                  }, height = unit(8, "mm")),
                show_annotation_name = FALSE))
pA <- wrap_elements(full = suppressWarnings(grid.grabExpr(
  draw(ht, annotation_legend_list = list(lgd), heatmap_legend_side = "left",
       annotation_legend_side = "left", padding = unit(c(2, 2, 2, 2), "mm")))))

wrap_strip <- function(x) {
  x <- sub("_", "\n", x)
  vapply(strsplit(x, "\n"), function(w) {
    w <- ifelse(nchar(w) > 13, paste0(substr(w, 1, 10), "-\n", substring(w, 11)), w)
    paste(w, collapse = "\n")
  }, character(1))
}

lab_ord <- rownames(mat_coef)[row_clus$order]
log_maaslin <- function(x) { y <- x; y[y == 0] <- min(x[x > 0]) / 2; log2(y) }
traj <- rel[, sel, drop = FALSE] %>%
  mutate(across(everything(), log_maaslin)) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "feature", values_to = "value") %>%
  left_join(mf0 %>% select(SampleID, Group, Time), by = "SampleID") %>%
  group_by(feature, Group, Time) %>%
  summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop") %>%
  left_join(TAXA, by = "feature") %>%
  mutate(label = factor(label, levels = lab_ord),
         Group = factor(Group, levels = GROUP_LV),
         Weeks = as.numeric(sub("W", "", Time)))

pB <- ggplot(traj, aes(Weeks, mean, colour = Group, group = Group)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.4,
                position = position_dodge(0.6), linewidth = 0.3) +
  geom_line(position = position_dodge(0.6), linewidth = 0.6) +
  geom_point(position = position_dodge(0.6), size = 1.2) +
  facet_wrap(~ label, scales = "free_y", ncol = 6,
             labeller = labeller(label = wrap_strip)) +
  scale_colour_manual(values = PAL, labels = c(
    "Placebo" = "Placebo",
    "C. asiatica" = expression(italic("C. asiatica")),
    "Fiber-enriched C. asiatica" = expression("Fiber-enriched"~italic("C. asiatica")))) +
  scale_x_continuous(breaks = c(0, 2, 4)) +
  labs(x = "Weeks", y = expression(bold(log[2]~"relative abundance (mean ± SE)"))) +
  theme_bw(base_size = 7, base_family = FONT) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7),
        strip.text = element_text(face = "bold.italic", size = 7, margin = margin(1, 1, 1, 1)),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7, colour = "black"),
        panel.grid.minor = element_blank(), panel.spacing = unit(1.5, "mm"))

n_row_B <- ceiling(length(sel) / 6)
hA <- nrow(mat_coef) * 3.2 + 32
hB <- n_row_B * 22 + 15
p <- pA / pB + plot_layout(heights = c(hA, hB)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))

W <- 180 / 25.4
H <- (hA + hB + 5) / 25.4
ggsave(file.path(RESULT_DIR, "Fig_S4.png"), p, width = W, height = H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S4.svg"), p, width = W, height = H, device = "svg")
message("Saved Fig_S4.{png,svg}")

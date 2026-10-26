# ==============================================================================
# SCRIPT: Figure 4 (A + B + C assembled) - fecal and plasma metabolite responses
# (A) class totals (as Figure4A.R); (B) within-arm MaAsLin2 heatmap (as Figure4B.R);
# (C) median trajectories (as Figure4C.R)
# ==============================================================================
pacman::p_load(dplyr, tidyr, tibble, readxl, missForest, ComplexHeatmap, circlize,
               grid, ggplot2, ggtext, patchwork, svglite, ragg)

# --- 1. Paths -----------------------------------------------------------------
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

FONT <- "Arial"
FS_S <- 5.5
FS_M <- 6
FS_L <- 6.5
FS_T <- 8
gp_txt <- function(size, face = "plain", ...) gpar(fontsize = size, fontface = face, fontfamily = FONT, ...)
GROUP_LV <- c("Placebo", "C. asiatica", "Fiber-enriched C. asiatica")
PAL <- setNames(c("#999999", "#e6c60d", "#5ca08e"), GROUP_LV)
LEG_LAB <- c("Placebo" = "Placebo",
             "C. asiatica" = expression(italic("C. asiatica")),
             "Fiber-enriched C. asiatica" = expression("Fiber-enriched"~italic("C. asiatica")))
theme_fig <- theme_bw(base_size = FS_S, base_family = FONT) +
  theme(axis.title = element_text(size = FS_M, face = "bold"),
        axis.text = element_text(size = FS_S, colour = "black"),
        axis.ticks = element_line(linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1.5, 1, 1.5))

MET_MIC  <- file.path(INPUT_DIR, "inputMetMic.xlsx")
METADATA <- read_excel(MET_MIC, sheet = "METADATA") %>%
  mutate(Treatments = factor(case_when(
    grepl("Fiber", Treatments, ignore.case = TRUE) ~ "Fiber-enriched C. asiatica",
    grepl("Asiatica", Treatments, ignore.case = TRUE) ~ "C. asiatica",
    TRUE ~ "Placebo"), levels = GROUP_LV))
Met <- read_excel(MET_MIC, sheet = "Met")

# --- 2. Panel A: class totals over weeks (column indexes as in Figure4A.R) ---
FACET_LAB <- c("Placebo" = "Placebo", "C. asiatica" = "*C. asiatica*",
               "Fiber-enriched C. asiatica" = "Fiber-<br>enriched<br>*C. asiatica*")
plot_total <- function(cols_idx, title_text) {
  d <- Met %>%
    mutate(Total_Val = rowSums(select(., all_of(colnames(Met)[cols_idx])), na.rm = TRUE)) %>%
    select(ID, Total_Val) %>%
    left_join(METADATA %>% select(ID, Weeks, Treatments), by = "ID") %>%
    mutate(Weeks = factor(paste0("W", Weeks), levels = c("W0", "W2", "W4"))) %>%
    filter(!is.na(Treatments))
  ggplot(d, aes(Weeks, Total_Val, fill = Treatments)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, colour = "black", linewidth = 0.2) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.2, colour = "black") +
    facet_grid(. ~ Treatments, labeller = as_labeller(FACET_LAB)) +
    scale_fill_manual(values = PAL) +
    labs(subtitle = title_text, x = NULL, y = NULL) +
    theme_fig +
    theme(panel.grid.major = element_blank(), legend.position = "none",
          panel.spacing = unit(0.8, "mm"),
          strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.25),
          strip.text = element_markdown(face = "bold", size = FS_S, halign = 0.5,
                                        lineheight = 0.9, margin = margin(0.5, 0.5, 0.5, 0.5)),
          plot.subtitle = element_text(hjust = 0.5, size = FS_M, face = "bold",
                                       margin = margin(0, 0, 0.5, 0)),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}
pA <- plot_total(c(4, 6, 7, 9), "Total SCFAs") / plot_total(c(18, 19, 20), "Total CEs") /
  plot_total(10:16, "Total FAs") / plot_total(40, "Total UTs") +
  plot_annotation(title = "Metabolome data over weeks",
                  theme = theme(plot.title = element_text(size = FS_L, face = "bold", hjust = 0.5,
                                                          family = FONT)))

# --- 3. Panel B: within-arm MaAsLin2 heatmap (as Figure4B.R) -----------------
read_arm <- function(arm) {
  read.delim(file.path(OUTPUT_DIR, "Met", paste0("output_", arm), "all_results.tsv"),
             check.names = FALSE) %>%
    mutate(Treatments_Weeks = paste0(arm, "_", value, "/W0")) %>%
    select(feature, Treatments_Weeks, pval, qval, coef)
}
RENAME <- c("total.Lipid" = "Total Lipid", "ω6.ω3.ratio" = "omega6/omega3 ratio",
            "CE.20.5." = "CE(20:5)", "MUFA.SFA.ratio" = "MUFA/SFA ratio",
            "PUFA.SFA.ratio" = "PUFA/SFA ratio", "CE.16.0." = "CE(16:0)", "CE.20.4." = "CE(20:4)",
            "CE.17.0." = "CE(17:0)", "CE.18.0." = "CE(18:0)", "CE.20.0." = "CE(20:0)",
            "CE.20.2." = "CE(20:2)", "CE.18.1." = "CE(18:1)", "Acetic.acid" = "Acetic acid",
            "Butanoic.acid" = "Butanoic acid", "total.CE" = "Total CE", "CE.18.2." = "CE(18:2)",
            "Pentanoic.acid" = "Propionic acid", "CE.16.1." = "CE(16:1)",
            "X2.Methylhexanoic.acid" = "2-Methylhexanoic acid", "CE.15.0." = "CE(15:0)",
            "X2.Methypropionic.acid" = "2-Methylpropionic acid",
            "X3.Methylbutanoic.acid" = "3-Methylbutanoic acid",
            "X2.Methylbutanoic.acid" = "2-Methylbutanoic acid",
            "p.Cresolsulfate" = "p-Cresolsulfate", "Indoxyl.sulfate" = "Indoxyl sulfate",
            "X2.Methoxyphenylsulfate" = "2-Methoxyphenylsulfate",
            "Hippuric.acid" = "Hippuric acid",
            "Phenylacetyl.L.glutamine" = "Phenylacetyl-L-glutamine",
            "Trimethylamine.N.oxide" = "Trimethylamine N-oxide")
ROW_ORDER <- list(
  SCFAs = c("Acetic acid", "Propionic acid", "2-Methylpropionic acid", "Butanoic acid",
            "3-Methylbutanoic acid", "2-Methylbutanoic acid", "2-Methylhexanoic acid"),
  CEs   = c("Total CE", "Cholesterol", "CE(15:0)", "CE(16:0)", "CE(16:1)", "CE(17:0)",
            "CE(18:0)", "CE(18:1)", "CE(18:2)", "CE(20:0)", "CE(20:2)", "CE(20:4)", "CE(20:5)"),
  FAs   = c("Total Lipid", "SFA", "MUFA", "PUFA", "MUFA/SFA ratio", "PUFA/SFA ratio",
            "omega6/omega3 ratio"),
  UTs   = c("p-Cresolsulfate", "Phenysulfate", "Indoxyl sulfate", "2-Methoxyphenylsulfate",
            "Hippuric acid", "Cinnamoylglycine", "Phenylacetyl-L-glutamine",
            "Trimethylamine N-oxide", "Dimethylglycine"))
HIGHLIGHT <- c("2-Methylhexanoic acid", "2-Methylpropionic acid", "3-Methylbutanoic acid",
               "2-Methylbutanoic acid", "Cholesterol", "Total Lipid", "SFA", "MUFA", "PUFA",
               "MUFA/SFA ratio", "PUFA/SFA ratio", "omega6/omega3 ratio", "Dimethylglycine")
target_cols <- c("Placebo_W2/W0", "Placebo_W4/W0", "Asiatica_W2/W0", "Asiatica_W4/W0",
                 "Fiber_added_W2/W0", "Fiber_added_W4/W0")

res <- bind_rows(lapply(c("Placebo", "Asiatica", "Fiber_added"), read_arm)) %>%
  mutate(feature = ifelse(feature %in% names(RENAME), RENAME[feature], feature)) %>%
  filter(Treatments_Weeks %in% target_cols) %>%
  group_by(feature, Treatments_Weeks) %>%
  summarise(coef = mean(coef), qval = min(qval), .groups = "drop")
row_order <- unlist(ROW_ORDER, use.names = FALSE)
wide_of <- function(v, fill) {
  m <- res %>% select(feature, Treatments_Weeks, all_of(v)) %>%
    pivot_wider(names_from = Treatments_Weeks, values_from = all_of(v), values_fill = fill) %>%
    column_to_rownames("feature") %>% as.matrix()
  miss <- setdiff(row_order, rownames(m))
  if (length(miss)) stop("metabolites not found: ", paste(miss, collapse = ", "))
  m[row_order, target_cols]
}
mat_fc   <- wide_of("coef", 0)
mat_q    <- wide_of("qval", 1)
mat_star <- ifelse(mat_q < 0.01, "**", ifelse(mat_q < 0.05, "*", ""))
cls <- factor(rep(names(ROW_ORDER), lengths(ROW_ORDER)), levels = names(ROW_ORDER))

lim_fc  <- 1.5
col_fun <- colorRamp2(c(-lim_fc, 0, lim_fc), c("blue", "white", "red"))
lgd <- Legend(col_fun = col_fun, title = "log2 fold change (vs W0)",
              at = c(-lim_fc, 0, lim_fc), labels = c("≤ -1.5", "0", "≥ 1.5"),
              direction = "horizontal", legend_width = unit(22, "mm"), grid_height = unit(2, "mm"),
              title_position = "topcenter",
              title_gp = gp_txt(FS_M, "bold"), labels_gp = gp_txt(FS_S))
n_scfa <- length(ROW_ORDER$SCFAs)
block <- function(idx, lwd) anno_block(
  align_to = idx,
  panel_fun = function(index, nm) {
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = lwd))
    grid.text(nm, rot = 90, gp = gp_txt(FS_M, "bold"))
  }, width = unit(3.5, "mm"))
ROW_H <- 2.5
ht <- Heatmap(mat_fc, name = "log2FC", col = col_fun, show_heatmap_legend = FALSE,
              cluster_rows = FALSE, cluster_columns = FALSE,
              border = TRUE, border_gp = gpar(lwd = 0.5),
              row_split = cls, row_title = NULL, row_gap = unit(1, "mm"),
              column_split = factor(rep(GROUP_LV, each = 2), levels = GROUP_LV),
              column_title = NULL, column_gap = unit(1, "mm"),
              column_labels = rep(c("W2 vs W0", "W4 vs W0"), 3),
              column_names_rot = 45, column_names_centered = TRUE,
              column_names_side = "bottom", column_names_gp = gp_txt(FS_S),
              rect_gp = gpar(col = "white", lwd = 0.3),
              cell_fun = function(j, i, x, y, width, height, fill)
                grid.text(mat_star[i, j], x, y, gp = gp_txt(FS_M)),
              row_names_side = "right",
              row_names_gp = gp_txt(FS_S, fill = ifelse(rownames(mat_fc) %in% HIGHLIGHT, "#90EE90", NA)),
              width = unit(6 * 6.5, "mm"), height = unit(nrow(mat_fc) * ROW_H, "mm"),
              left_annotation = rowAnnotation(
                Major = block(list("Feces" = 1:n_scfa, "Plasma" = (n_scfa + 1):nrow(mat_fc)), 0.5),
                Minor = block(split(seq_along(cls), cls), 0.7),
                show_annotation_name = FALSE),
              top_annotation = HeatmapAnnotation(
                arm = anno_block(
                  align_to = setNames(list(1:2, 3:4, 5:6), GROUP_LV),
                  panel_fun = function(index, nm) {
                    grid.rect(gp = gpar(fill = "white", col = "black", lwd = 0.5))
                    if (nm == "Placebo") {
                      grid.text("Placebo", gp = gp_txt(FS_S, "bold"))
                    } else if (nm == "C. asiatica") {
                      grid.text("C. asiatica", gp = gp_txt(FS_S, "bold.italic"))
                    } else {
                      grid.text("Fiber-enriched", y = unit(0.7, "npc"), gp = gp_txt(FS_S, "bold"))
                      grid.text("C. asiatica", y = unit(0.28, "npc"), gp = gp_txt(FS_S, "bold.italic"))
                    }
                  }, height = unit(5, "mm")),
                show_annotation_name = FALSE))

# --- 4. Panel C: median trajectories (columns and imputation as Figure4C.R) --
met <- Met %>% select(-ID) %>% as.data.frame()
met[] <- lapply(met, function(x) as.numeric(as.character(x)))
set.seed(123)
met <- missForest(met)$ximp
rownames(met) <- Met$ID

wrap_strip <- function(x, width = 12) {
  x <- sub("(methyl)(?=[a-z])", "\\1- ", x, perl = TRUE, ignore.case = TRUE)
  vapply(strwrap(x, width = width, simplify = FALSE), function(ln) {
    ln <- vapply(ln, function(w) {
      if (nchar(w) <= width) return(w)
      h <- gregexpr("[-/]", w)[[1]]; h <- h[h > 2 & h < nchar(w) - 2]
      if (length(h)) paste0(substr(w, 1, max(h)), "\n", substring(w, max(h) + 1))
      else paste0(substr(w, 1, width - 1), "-\n", substring(w, width))
    }, character(1))
    paste(ln, collapse = "\n")
  }, character(1))
}
plot_traj <- function(cols_idx, title, ylab, ncol, show_x = TRUE, strip_w = 12) {
  nm <- colnames(Met)[cols_idx]
  d <- met %>% select(all_of(nm)) %>%
    rownames_to_column("ID") %>%
    pivot_longer(-ID, names_to = "Metabolite", values_to = "Value") %>%
    left_join(METADATA %>% select(ID, Weeks, Treatments), by = "ID") %>%
    mutate(Weeks = as.numeric(as.character(Weeks)),
           Metabolite = factor(Metabolite, levels = nm))
  ggplot(d, aes(Weeks, Value, colour = Treatments, group = Treatments)) +
    stat_summary(fun = median, geom = "line", linewidth = 0.4) +
    stat_summary(fun = median, geom = "point", size = 0.5) +
    facet_wrap(~ Metabolite, scales = "free_y", ncol = ncol,
               labeller = labeller(Metabolite = function(x) wrap_strip(x, strip_w))) +
    scale_colour_manual(values = PAL, labels = LEG_LAB, name = NULL) +
    scale_x_continuous(breaks = c(0, 2, 4)) +
    labs(title = title, x = if (show_x) "Weeks" else NULL, y = ylab) +
    theme_fig +
    theme(legend.position = "none", panel.spacing = unit(1, "mm"),
          strip.text = element_text(face = "bold", size = FS_S, lineheight = 0.9,
                                    margin = margin(0.5, 0.5, 0.5, 0.5)),
          plot.title = element_text(size = FS_M, face = "bold", hjust = 0.5,
                                    margin = margin(0, 0, 0.5, 0)))
}
pC <- list(
  SCFAs = plot_traj(c(4, 6, 7, 9), "Short-chain fatty acids (SCFAs)", "µmol/g", 4),
  CEs   = plot_traj(11, "Cholesteryl esters (CEs)", "µM", 1),
  FAs   = plot_traj(24:30, "Fatty acids (FAs)", "µM", 4, strip_w = 6),
  UTs   = plot_traj(40, "Uremic toxins (UTs)", "µM", 1))
get_legend_C <- function() {
  g <- ggplotGrob(pC$SCFAs + theme(legend.position = "bottom",
                                   legend.text = element_text(size = FS_S),
                                   legend.key.width = unit(3, "mm"),
                                   legend.key.height = unit(2, "mm")))
  idx <- grep("guide-box", g$layout$name)
  Filter(function(x) !inherits(x, "zeroGrob"), g$grobs[idx])[[1]]
}

# --- 5. Page layout (mm) -------------------------------------------------------
W <- 180; H <- 130
X_A <- c(0, 45); X_B <- c(47, 112); X_C <- c(119, 180)
C_ROWS <- list(SCFAs = c(4, 31), CEs = c(34, 56), FAs = c(59, 101), UTs = c(104, 125))
C_WIDTH_SINGLE <- 0.36

vp_mm <- function(x, y_top, w, h) viewport(x = unit(x, "mm"), y = unit(H - y_top, "mm"),
                                           width = unit(w, "mm"), height = unit(h, "mm"),
                                           just = c("left", "top"))
loc_in <- function(x, y) { l <- deviceLoc(x, y); c(x = as.numeric(l$x), y = as.numeric(l$y)) }

draw_page <- function() {
  grid.newpage()
  ctr <- list()

  pushViewport(vp_mm(X_A[1] + 3, 2, diff(X_A) - 3, H - 4)); print(pA, newpage = FALSE); upViewport()

  pushViewport(vp_mm(X_B[1] + 2, 2, diff(X_B) - 2, H - 4))
  draw(ht, newpage = FALSE, annotation_legend_list = list(lgd),
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       padding = unit(c(1, 0, 0, 0), "mm"))
  row_anchor <- function(first, last) {
    k  <- as.character(cls[match(first, rownames(mat_fc))])
    in_k <- rownames(mat_fc)[cls == k]; nk <- length(in_k)
    i1 <- match(first, in_k); i2 <- match(last, in_k)
    out <- NULL
    decorate_row_names("log2FC", slice = which(levels(cls) == k), {
      out <<- rbind(top    = loc_in(unit(1, "npc"), unit(1 - (i1 - 1) / nk, "npc")),
                    bottom = loc_in(unit(1, "npc"), unit(1 - i2 / nk, "npc")))
    })
    out
  }
  ctr$SCFAs <- row_anchor("2-Methylpropionic acid", "2-Methylhexanoic acid")
  ctr$CEs   <- row_anchor("Cholesterol", "Cholesterol")
  ctr$FAs   <- row_anchor("Total Lipid", "omega6/omega3 ratio")
  ctr$UTs   <- row_anchor("Dimethylglycine", "Dimethylglycine")
  upViewport()

  tgt <- list()
  for (k in names(pC)) {
    w <- diff(X_C) * if (k %in% c("CEs", "UTs")) C_WIDTH_SINGLE else 1
    h <- diff(C_ROWS[[k]])
    pushViewport(vp_mm(X_C[1], C_ROWS[[k]][1], w, h))
    print(pC[[k]], newpage = FALSE)
    tgt[[k]] <- rbind(top    = loc_in(unit(0, "npc"), unit(1, "npc") - unit(3, "mm")),
                      bottom = loc_in(unit(0, "npc"), unit(0, "npc") + unit(3, "mm")))
    upViewport()
  }
  pushViewport(vp_mm(X_C[1], 125.5, diff(X_C), 4.5)); grid.draw(get_legend_C()); upViewport()

  for (k in names(ctr)) for (e in c("top", "bottom"))
    grid.lines(unit(c(ctr[[k]][e, "x"], tgt[[k]][e, "x"]), "inches"),
               unit(c(ctr[[k]][e, "y"], tgt[[k]][e, "y"]), "inches"),
               gp = gpar(lty = "22", col = "grey45", lwd = 0.5))

  for (tg in list(c("A", X_A[1]), c("B", X_B[1]), c("C", X_C[1] - 4)))
    grid.text(tg[1], x = unit(as.numeric(tg[2]) + 1, "mm"), y = unit(H - 3, "mm"),
              just = c("left", "top"), gp = gp_txt(FS_T, "bold"))
}

agg_png(file.path(RESULT_DIR, "Figure4.png"), width = W, height = H, units = "mm", res = 600)
draw_page(); invisible(dev.off())
svglite(file.path(RESULT_DIR, "Figure4.svg"), width = W / 25.4, height = H / 25.4)
draw_page(); invisible(dev.off())
message("Saved Figure4.{png,svg} (", W, " x ", H, " mm)")

# ==============================================================================
# SCRIPT: Figure S7 + Figure S8 + Table S5 - Microbe x Metabolite delta-correlation heatmaps
# ==============================================================================
pacman::p_load(readxl, dplyr, tidyr, tibble, purrr, pheatmap, RColorBrewer,
               grid, ggplot2, cowplot, svglite, openxlsx)

# --- 1. Path Configuration ---
BASE_DIR <- "//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025"
if (dir.exists(BASE_DIR)) setwd(BASE_DIR)

INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
FILE_PATH  <- file.path(INPUT_DIR, "CPF_Downstream_analysis.xlsx")
TABLE_S5   <- file.path(OUTPUT_DIR, "Table_S5.xlsx")
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

microbe_df <- read_excel(FILE_PATH, sheet = "Microbes (Rarefied)")
clinical   <- read_excel(FILE_PATH, sheet = "Clinical")
fecal_met  <- read_excel(FILE_PATH, sheet = "Fecal_Metabolites (uM)")

GROUPS <- c(A = "Placebo", B = "Asiatica", C = "Fiber_added")
P_CUT  <- 0.05
TOP_N  <- 25
MIN_N  <- 6

# --- 2. Delta data for one group / one timepoint ------------------------------
get_delta_data <- function(group_name, time_target) {
  df_group      <- clinical %>% filter(Treatments == group_name)
  group_samples <- df_group$SampleID
  
  mic_group_data <- microbe_df %>%
    group_by(Genus) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup()
  
  mic_temp  <- mic_group_data %>% select(Genus, all_of(group_samples))
  presence  <- rowSums(mic_temp %>% select(-Genus) > 0)
  threshold <- 0.5 * ncol(mic_temp %>% select(-Genus))
  valid_genera <- mic_temp$Genus[presence >= threshold]
  
  common_sub <- sort(intersect(df_group$Subject[df_group$Weeks == 0],
                               df_group$Subject[df_group$Weeks == time_target]))
  if(length(common_sub) < MIN_N) return(NULL)
  
  map_ids <- df_group %>%
    filter(Subject %in% common_sub, Weeks %in% c(0, time_target)) %>%
    arrange(Subject, Weeks)
  
  mic_final <- mic_group_data %>%
    filter(Genus %in% valid_genera) %>%
    column_to_rownames("Genus") %>%
    select(all_of(map_ids$SampleID))
  
  mic_clr <- apply(sweep(mic_final, 2, colSums(mic_final), "/") + 0.1, 2,
                   function(x) log(x) - mean(log(x)))
  
  met_data <- fecal_met %>% filter(ID %in% map_ids$SampleID) %>%
    column_to_rownames("ID") %>% t()
  met_data   <- met_data[, map_ids$SampleID]
  met_scaled <- t(apply(log1p(met_data), 1, scale))
  colnames(met_scaled) <- colnames(met_data)
  
  calc_delta <- function(mat) {
    as.matrix(sapply(common_sub, function(s) {
      mat[, map_ids$SampleID[map_ids$Subject == s & map_ids$Weeks == time_target]] -
        mat[, map_ids$SampleID[map_ids$Subject == s & map_ids$Weeks == 0]]
    }))
  }
  
  list(m = calc_delta(mic_clr), mt = calc_delta(met_scaled))
}

# --- 3. Top-25 correlation matrix --------------------------------------------
build_matrix <- function(dat) {
  cor_res <- map_dfr(rownames(dat$mt), function(met) {
    map_dfr(rownames(dat$m), function(mic) {
      x  <- dat$m[mic, ]
      y  <- dat$mt[met, ]
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < MIN_N) return(NULL)
      test <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      data.frame(from = mic, to = met, rho = test$estimate, p_value = test$p.value)
    })
  })
  
  full_stats <- cor_res %>%
    mutate(p.adj = p.adjust(p_value, method = "BH"),
           abs_rho = abs(rho)) %>%
    transmute(Microbe = from, Metabolite = to, rho, p = p_value, p.adj, abs_rho)
  rownames(full_stats) <- NULL
  
  top <- cor_res %>%
    slice_max(order_by = abs(rho), n = TOP_N, with_ties = FALSE) %>%
    mutate(rho_shown = ifelse(p_value < P_CUT, rho, 0))
  
  mat <- top %>%
    select(from, to, rho_shown) %>%
    pivot_wider(names_from = to, values_from = rho_shown, values_fill = 0) %>%
    column_to_rownames("from") %>%
    as.matrix()
  
  keep_r <- rowSums(mat != 0) > 0
  keep_c <- colSums(mat != 0) > 0
  message("  dropped ", sum(!keep_r), " all-zero taxa, ", sum(!keep_c), " all-zero metabolites")
  mat <- mat[keep_r, keep_c, drop = FALSE]
  
  list(mat = mat, stats = top, full = full_stats)
}

# --- 4. Draw & save -----------------------------------------------------------
PAL    <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
BREAKS <- seq(-1, 1, length.out = 101)

CELL_W  <- 52
CELL_H  <- 20
FS      <- 15
TITLE_H <- 0.55

make_heatmap <- function(mat, title) {
  ph <- pheatmap(mat,
                 cluster_rows     = nrow(mat) >= 2,
                 cluster_cols     = ncol(mat) >= 2,
                 color            = PAL,
                 breaks           = BREAKS,
                 display_numbers  = TRUE,
                 number_format    = "%.2f",
                 number_color     = "black",
                 cellwidth        = CELL_W,
                 cellheight       = CELL_H,
                 fontsize         = FS,
                 fontsize_row     = FS,
                 fontsize_col     = FS,
                 fontsize_number  = 12,
                 border_color     = "grey70",
                 angle_col        = 90,
                 legend           = FALSE,
                 main             = "",
                 silent           = TRUE)
  
  chr_in     <- FS * 0.55 / 72
  row_lab_in <- max(nchar(rownames(mat))) * chr_in
  col_lab_in <- max(nchar(colnames(mat))) * chr_in
  
  w      <- ncol(mat) * CELL_W / 72 + row_lab_in + 1.5
  body_h <- nrow(mat) * CELL_H / 72 + col_lab_in + 1.2
  
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold", size = FS * 1.3)
  
  list(grob = cowplot::plot_grid(title_grob, ph$gtable, ncol = 1,
                                 rel_heights = c(TITLE_H, body_h)),
       w    = w,
       h    = TITLE_H + body_h)
}

save_figure <- function(panels, stem) {
  grobs <- lapply(panels, `[[`, "grob")
  hs    <- vapply(panels, `[[`, numeric(1), "h")
  w     <- max(vapply(panels, `[[`, numeric(1), "w"))
  
  combined <- cowplot::plot_grid(plotlist = grobs, ncol = 1,
                                 rel_heights = hs,
                                 labels = names(panels),
                                 label_size = 26, label_fontface = "bold",
                                 hjust = 0, vjust = 1)
  
  ggsave(file.path(RESULT_DIR, paste0(stem, ".png")), combined,
         width = w, height = sum(hs), dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(file.path(RESULT_DIR, paste0(stem, ".svg")), combined,
         width = w, height = sum(hs), device = "svg", bg = "white", limitsize = FALSE)
  
  message("wrote ", stem, ".{png,svg}  — panels ", paste(names(panels), collapse = "/"))
}

# --- 5. Execution -------------------------------------------------------------
make_title <- function(group_name, week)
  paste0("Top 25 Delta correlations (rho) – ", group_name,
         " (Week ", week, " vs Week 0)")

figs <- list(S7 = 2, S8 = 4)

full_tabs <- list()

for (fig in names(figs)) {
  wk <- figs[[fig]]
  panels <- list()
  
  for (panel in names(GROUPS)) {
    grp <- GROUPS[[panel]]
    message("Figure ", fig, " panel ", panel, " - ", grp)
    
    dat <- get_delta_data(grp, wk)
    if (is.null(dat)) {
      message("  (<6 paired subjects — skipped)")
      next
    }
    
    res <- build_matrix(dat)
    panels[[panel]] <- make_heatmap(res$mat, make_title(grp, wk))
    
    full_tabs[[paste0("Table_S5_Full_W", wk, "_", grp)]] <- res$full
  }
  
  if (length(panels) > 0) save_figure(panels, paste0("Fig_", fig))
}

# --- 5b. Repeated-measures correlation across all samples --------------------
rm_cor <- function(x, y, subj) {
  fit <- lm(y ~ factor(subj) + x)
  a   <- anova(fit)
  ss_x <- a["x", "Sum Sq"]; ss_r <- a["Residuals", "Sum Sq"]
  data.frame(r = sign(coef(fit)[["x"]]) * sqrt(ss_x / (ss_x + ss_r)),
             df = a["Residuals", "Df"], p = a["x", "Pr(>F)"])
}

all_ids   <- clinical$SampleID[clinical$SampleID %in% names(microbe_df) &
                                 clinical$SampleID %in% fecal_met$ID]
mic_all   <- microbe_df %>% group_by(Genus) %>%
  summarise(across(where(is.numeric), sum)) %>% ungroup() %>%
  column_to_rownames("Genus") %>% select(all_of(all_ids))
mic_all   <- mic_all[rowSums(mic_all > 0) >= 0.5 * ncol(mic_all), ]
clr_all   <- apply(sweep(mic_all, 2, colSums(mic_all), "/") + 0.1, 2,
                   function(x) log(x) - mean(log(x)))
met_all   <- fecal_met %>% filter(ID %in% all_ids) %>% column_to_rownames("ID") %>% t()
met_all   <- t(apply(log1p(met_all[, all_ids]), 1, scale))
colnames(met_all) <- all_ids
meta_all  <- clinical[match(all_ids, clinical$SampleID), ]

rm_all <- map_dfr(rownames(met_all), function(met) {
  map_dfr(rownames(clr_all), function(mic) {
    arm <- lapply(GROUPS, function(g) {
      k <- meta_all$Treatments == g
      rm_cor(clr_all[mic, k], met_all[met, k], meta_all$Subject[k])
    })
    cbind(data.frame(Microbe = mic, Metabolite = met),
          rm_cor(clr_all[mic, ], met_all[met, ], meta_all$Subject),
          r_Placebo = arm$A$r, p_Placebo = arm$A$p,
          r_Asiatica = arm$B$r, p_Asiatica = arm$B$p,
          r_Fiber_added = arm$C$r, p_Fiber_added = arm$C$p)
  })
}) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  relocate(p.adj, .after = p) %>%
  arrange(p)
full_tabs[["Table_S5_All_samples"]] <- rm_all
message("All-sample rmcorr: ", nrow(rm_all), " pairs, ", sum(rm_all$p < 0.05),
        " with p < 0.05, ", sum(rm_all$p.adj < 0.05), " with q < 0.05")

# --- 6. Append full correlation stats to Table_S5.xlsx ------------------------
if (length(full_tabs) > 0) {
  wb <- createWorkbook()
  
  hs <- createStyle(textDecoration = "bold", halign = "center",
                    fgFill = "#EDEDED", border = "bottom")
  
  for (nm in names(full_tabs)) {
    df <- as.data.frame(full_tabs[[nm]])
    if (nm %in% names(wb)) removeWorksheet(wb, nm)
    addWorksheet(wb, nm)
    writeData(wb, nm, df, headerStyle = hs, withFilter = TRUE)
    freezePane(wb, nm, firstActiveRow = 2)
    setColWidths(wb, nm, cols = seq_len(ncol(df)), widths = "auto")
    message("  sheet ", nm, " (", nrow(df), " rows)")
  }
  
  saveWorkbook(wb, TABLE_S5, overwrite = TRUE)
  cat("Wrote", TABLE_S5, "\n")
}
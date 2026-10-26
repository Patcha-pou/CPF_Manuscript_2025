# ==============================================================================
# SCRIPT: Figure 5C + 5D (Delta Networks, right column of the Figure 5 layout)
# ==============================================================================
pacman::p_load(readxl, dplyr, tidyr, igraph, ggraph, tidygraph, tibble, purrr,
               ggplot2, svglite, openxlsx)

# --- 1. Path Configuration ---
BASE_DIR <- "//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025"
if (dir.exists(BASE_DIR)) setwd(BASE_DIR)

INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
FILE_PATH  <- file.path(INPUT_DIR, "CPF_Downstream_analysis.xlsx")
TABLE_S5   <- file.path(OUTPUT_DIR, "Table_S5.xlsx")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

GROUP <- "Fiber_added"
TOP_N <- 25
MIN_N <- 6

microbe_df <- read_excel(FILE_PATH, sheet = "Microbes (Rarefied)")
clinical   <- read_excel(FILE_PATH, sheet = "Clinical")
fecal_met  <- read_excel(FILE_PATH, sheet = "Fecal_Metabolites (uM)")

# --- 2. Data Processing Function ---
get_delta_data <- function(time_target) {
  df_group <- clinical %>% filter(Treatments == GROUP)
  group_samples <- df_group$SampleID

  mic_group_data <- microbe_df %>%
    group_by(Genus) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup()

  mic_temp <- mic_group_data %>% select(Genus, all_of(group_samples))
  presence <- rowSums(mic_temp %>% select(-Genus) > 0)
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

  met_data <- fecal_met %>% filter(ID %in% map_ids$SampleID) %>% column_to_rownames("ID") %>% t()
  met_data <- met_data[, map_ids$SampleID]
  met_scaled <- t(apply(log1p(met_data), 1, scale))
  colnames(met_scaled) <- colnames(met_data)

  calc_delta <- function(mat) {
    res <- sapply(common_sub, function(s) {
      mat[, map_ids$SampleID[map_ids$Subject == s & map_ids$Weeks == time_target]] -
        mat[, map_ids$SampleID[map_ids$Subject == s & map_ids$Weeks == 0]]
    })
    return(as.matrix(res))
  }

  return(list(m = calc_delta(mic_clr), mt = calc_delta(met_scaled)))
}

# --- 2b. Edge statistics: read back from Table S5 -----------------------------
read_cor_stats <- function(week) {
  sheet <- paste0("Table_S5_Full_W", week, "_", GROUP)
  if (!file.exists(TABLE_S5))
    stop("missing ", TABLE_S5, " — run FigureS7S8_TableS5.R first", call. = FALSE)
  if (!sheet %in% excel_sheets(TABLE_S5))
    stop("missing sheet '", sheet, "' in ", TABLE_S5,
         " — run FigureS7S8_TableS5.R first", call. = FALSE)

  read_excel(TABLE_S5, sheet = sheet) %>%
    transmute(from = Microbe, to = Metabolite, rho, p_value = p, q_value = p.adj)
}

# --- 3. Plotting Function ---
process_and_save <- function(time_val, week_label, name) {
  dat <- get_delta_data(time_val)
  if(is.null(dat)) return(NULL)

  cor_res <- read_cor_stats(time_val)

  final_edges <- cor_res %>%
    slice_max(order_by = abs(rho), n = TOP_N, with_ties = FALSE)

  node_names <- unique(c(final_edges$from, final_edges$to))
  combined_delta_mat <- rbind(dat$m, dat$mt)

  node_info <- data.frame(name = node_names) %>%
    mutate(
      Type = ifelse(name %in% rownames(dat$m), "Microbe", "Metabolite"),
      Type = factor(Type, levels = c("Microbe", "Metabolite")),
      abs_delta = abs(rowMeans(combined_delta_mat[name, , drop=FALSE]))
    )


  g <- tbl_graph(nodes = node_info,
                 edges = final_edges %>% mutate(abs_rho = abs(rho),
                                                FDR = ifelse(q_value < 0.05, "q < 0.05", "q \u2265 0.05")),
                 directed = FALSE)

  set.seed(123)
  p <- ggraph(g, layout = 'fr', weights = abs_rho) +
    geom_edge_link(aes(color = rho, linetype = FDR), edge_width = 1.5, alpha = 0.8) +
    geom_node_point(aes(fill = Type, shape = Type, size = abs_delta), color = "black", stroke = 0.5) +
    geom_node_text(aes(label = name,
                       fontface = ifelse(Type == "Microbe", "italic", "plain")),
                   repel = TRUE, size = 5) +

    scale_fill_manual(values = c("Microbe" = "#089099", "Metabolite" = "#F4A261"),
                      breaks = c("Microbe", "Metabolite")) +
    scale_shape_manual(values = c("Microbe" = 21, "Metabolite" = 22),
                       breaks = c("Microbe", "Metabolite")) +
    scale_size_continuous(range = c(6, 18), name = "|Mean Delta|") +
    scale_edge_color_gradient2(low = "#3B4CC0", mid = "#EEEEEE", high = "#DC3977",
                               midpoint = 0, name = "Spearman Rho") +
    scale_edge_linetype_manual(values = c("q < 0.05" = "solid", "q \u2265 0.05" = "dashed"),
                               drop = FALSE, name = "BH-FDR") +

    guides(
      fill = guide_legend(order = 1, override.aes = list(size = 8)),
      shape = guide_legend(order = 1),
      size = guide_legend(order = 2),
      edge_color = guide_edge_colorbar(order = 3),
      edge_linetype = guide_legend(order = 5, keywidth = unit(1.6, "cm"))
    ) +

    labs(title = paste("Delta Correlation Network:", week_label),
         subtitle = paste0("Prevalence: 50% within group | Edges: Top ", TOP_N)) +
    theme_graph() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 16)
    )

  return(list(plot = p, nodes = node_info, edges = final_edges))
}

# --- 4. Execution Loop ---
times   <- list(W2 = 2, W4 = 4)
panel   <- c(W2 = "Figure5C", W4 = "Figure5D")
tabs_s6 <- list()

PANEL_DIR <- file.path(OUTPUT_DIR, "panels")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)

for(name in names(times)) {
  res <- process_and_save(times[[name]], paste(name, "vs W0"), name)
  if(!is.null(res)) {
    ggsave(file.path(RESULT_DIR, paste0(panel[[name]], ".png")), res$plot, width = 12, height = 10, dpi = 600)
    ggsave(file.path(RESULT_DIR, paste0(panel[[name]], ".svg")), res$plot, width = 12, height = 10, device = "svg")
    saveRDS(res$plot, file.path(PANEL_DIR, paste0(panel[[name]], ".rds")))
    cat("Success: PNG and SVG generated for", name, "->", panel[[name]], "\n")

    tabs_s6[[paste0("Table_S6_Nodes_", name)]] <- res$nodes
    tabs_s6[[paste0("Table_S6_Edges_", name)]] <- res$edges
  }
}

# --- 5. Table S6: one workbook, nodes + edges per timepoint ------------------
if (length(tabs_s6) > 0) {
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

  for (nm in names(tabs_s6)) add_sheet(wb, nm, tabs_s6[[nm]])
  saveWorkbook(wb, file.path(OUTPUT_DIR, "Table_S6.xlsx"), overwrite = TRUE)
  cat("Wrote", file.path(OUTPUT_DIR, "Table_S6.xlsx"), "\n")
}

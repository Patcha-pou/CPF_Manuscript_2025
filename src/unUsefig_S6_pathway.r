# ==============================================================================
# SCRIPT: Figure S6 - Functional Abundance (W0 vs W2 and W0 vs W4)
# ==============================================================================
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, purrr, ggpicrust2, tidyverse, ggprism, patchwork, funrar)

# --- 1. Set Path & Load Data ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
FILE_PATH  <- file.path(INPUT_DIR, "CPF_Downstream_analysis.xlsx")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

clinical <- read_excel(FILE_PATH, sheet = "Clinical")
ko_abundance_W2 <- read_delim(file.path(OUTPUT_DIR, "Table_S6_Abundance_W0_vs_W2.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
ko_abundance_W4 <- read_delim(file.path(OUTPUT_DIR, "Table_S6_Abundance_W0_vs_W4.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

Fiber_added_all <- filter(clinical, Treatments == "Fiber_added")

fetch_kegg_names <- function() {
  message("กำลังดึงชื่อ Pathway จาก KEGG API...")
  tryCatch({
    kegg_list <- read.delim("https://rest.kegg.jp/list/pathway/ko", header = FALSE, col.names = c("path_id", "pathway_name"))
    kegg_list$path_id <- gsub("path:", "", kegg_list$path_id)
    return(kegg_list)
  }, error = function(e) {
    message("ไม่สามารถเชื่อมต่อ KEGG ได้ จะใช้รหัส ID แทน")
    return(NULL)
  })
}
kegg_map <- fetch_kegg_names()

Fiber_0_2 <- filter(Fiber_added_all, Weeks %in% c("0", "2"))
Fiber_0_2$SampleID <- as.character(Fiber_0_2$SampleID)

df_0_2 <- ko_abundance_W2 %>% select(1, any_of(Fiber_0_2$SampleID))
colnames(df_0_2)[1] <- "#NAME"

res_0_2 <- ggpicrust2(data = df_0_2, metadata = Fiber_0_2, group = "Weeks", pathway = "KO",
                      daa_method = "LinDA", ko_to_kegg = TRUE, reference = "0", p.adjust = "none")

res_0_2$daa_results_df <- res_0_2$daa_results_df %>%
  select(-any_of(c("pathway_name", "pathway_description"))) %>%
  left_join(kegg_map, by = c("feature" = "path_id")) %>%
  mutate(pathway_name = ifelse(is.na(pathway_name), feature, pathway_name))

Fiber_0_4 <- filter(Fiber_added_all, Weeks %in% c("0", "4"))
Fiber_0_4$SampleID <- as.character(Fiber_0_4$SampleID)

df_0_4 <- ko_abundance_W4 %>% select(1, any_of(Fiber_0_4$SampleID))
colnames(df_0_4)[1] <- "#NAME"

res_0_4 <- ggpicrust2(data = df_0_4, metadata = Fiber_0_4, group = "Weeks", pathway = "KO",
                      daa_method = "LinDA", ko_to_kegg = TRUE, reference = "0", p.adjust = "none")

res_0_4$daa_results_df <- res_0_4$daa_results_df %>%
  select(-any_of(c("pathway_name", "pathway_description"))) %>%
  left_join(kegg_map, by = c("feature" = "path_id")) %>%
  mutate(pathway_name = ifelse(is.na(pathway_name), feature, pathway_name))

print("--- Result W0 vs W4 ---")
print(head(res_0_4$daa_results_df[, c("feature", "pathway_name", "p_values")]))
pacman::p_load(ggplot2, patchwork, dplyr)

create_lollipop <- function(daa_data, title_text) {
  plot_df <- daa_data %>%
    filter(p_values < 0.05) %>%
    arrange(p_values) %>%
    head(15) %>%
    mutate(direction = ifelse(log2_fold_change > 0, "Increased", "Decreased"),
           pathway_name = reorder(pathway_name, log2_fold_change))
  
  if(nrow(plot_df) == 0) {
    plot_df <- daa_data %>%
      arrange(p_values) %>%
      head(10) %>%
      mutate(direction = ifelse(log2_fold_change > 0, "Increased", "Decreased"),
             pathway_name = reorder(pathway_name, log2_fold_change))
  }
  
  ggplot(plot_df, aes(x = pathway_name, y = log2_fold_change)) +
    geom_segment(aes(x = pathway_name, xend = pathway_name, y = 0, yend = log2_fold_change), 
                 color = "grey80", linewidth = 0.6) +
    geom_point(aes(color = direction), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
    coord_flip() +
    theme_prism(base_size = 9) +
    scale_color_manual(values = c("Increased" = "#E64B35", "Decreased" = "#4DBBD5")) +
    labs(title = title_text, x = NULL, y = "Log2 Fold Change") +
    theme(legend.position = "bottom", 
          plot.title = element_text(size = 10, face = "bold"))
}

p1 <- create_lollipop(res_0_2$daa_results_df, "Fiber Added: Week 0 vs 2")

p2 <- create_lollipop(res_0_4$daa_results_df, "Fiber Added: Week 0 vs 4")

combined_plot <- p1 + p2 + 
  plot_layout(guides = "collect") + 
  plot_annotation(title = "Differential Functional Pathways (Raw p-value < 0.05)",
                  subtitle = "Comparing Baseline (W0) with W2 and W4",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) & 
  theme(legend.position = "bottom")

print(combined_plot)

ggsave(file.path(RESULT_DIR, "Figure_S6_Combined_Lollipop.png"),
       combined_plot, width = 12, height = 7, dpi = 300)
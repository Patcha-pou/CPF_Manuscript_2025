# ==============================================================================
# SCRIPT: Figure S10 - Paired Boxplot with FDR Correction
# ==============================================================================
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, ggplot2, ggprism, ggpubr,
               rstatix, svglite, scales)

# --- 0. Paths (same convention as the other scripts) ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 1. Load & Prepare Data ---
clinical <- read_excel(file.path(INPUT_DIR, "CPF_Downstream_analysis.xlsx"), sheet = "Clinical")
ko_raw   <- read_delim(file.path(INPUT_DIR, "KO_prepT.tsv"), delim = "\t",
                       show_col_types = FALSE)

Fiber_added_meta <- clinical %>%
  filter(Treatments == "Fiber_added") %>%
  select(SampleID, Subject, Weeks) %>%
  mutate(SampleID = as.character(SampleID),
         Weeks    = parse_number(as.character(Weeks)))

ko_labels <- c(
  "K00826" = "K00826: BCAT",
  "K00128" = "K00128: ALDH"
)

DETECT_MIN <- 0.1

full_pathway_kos <- names(ko_labels)

missing_kos <- setdiff(full_pathway_kos, ko_raw[[1]])
if (length(missing_kos) > 0)
  warning("KO not found in KO_prepT.tsv: ", paste(missing_kos, collapse = ", "))

ko_rel <- ko_raw %>%
  rename(ko_id = 1)

raw_combined <- ko_rel %>%
  filter(ko_id %in% full_pathway_kos) %>%
  pivot_longer(cols = -ko_id, names_to = "SampleID", values_to = "Abundance") %>%
  inner_join(Fiber_added_meta, by = "SampleID") %>%
  filter(Weeks %in% c(0, 2, 4)) %>%
  mutate(
    Weeks   = factor(as.character(Weeks), levels = c("0", "2", "4")),
    feature = factor(ko_labels[ko_id], levels = ko_labels)
  )

complete_sub <- raw_combined %>%
  distinct(Subject, Weeks) %>%
  count(Subject) %>%
  filter(n == 3) %>%
  pull(Subject)

raw_combined <- raw_combined %>%
  filter(Subject %in% complete_sub) %>%
  arrange(feature, Weeks, Subject)

detect <- raw_combined %>%
  group_by(ko_id, feature) %>%
  summarise(detect_rate = mean(Abundance > 0), .groups = "drop")

print(detect)

dropped <- detect %>% filter(detect_rate < DETECT_MIN)
if (nrow(dropped) > 0)
  message("Dropped (detected in <", DETECT_MIN * 100, "% of samples): ",
          paste(dropped$ko_id, round(dropped$detect_rate, 2), collapse = "; "))

raw_combined <- raw_combined %>%
  filter(ko_id %in% detect$ko_id[detect$detect_rate >= DETECT_MIN]) %>%
  mutate(feature = droplevels(feature))

stopifnot(nrow(raw_combined) > 0)
n_feat <- nlevels(raw_combined$feature)

cat("Paired subjects (Fiber_added, W0/W2/W4):", length(complete_sub), "\n")

stat_test <- raw_combined %>%
  group_by(feature) %>%
  rstatix::wilcox_test(Abundance ~ Weeks, paired = TRUE,
                       comparisons = list(c("0", "2"), c("0", "4"))) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "Weeks", scales = "free_y")

stat_test <- stat_test %>%
  left_join(raw_combined %>%
              group_by(feature) %>%
              summarise(.ymin = min(Abundance), .ymax = max(Abundance),
                        .groups = "drop"),
            by = "feature") %>%
  group_by(feature) %>%
  mutate(y.position = .ymax + (.ymax - .ymin) * (0.07 + 0.16 * (row_number() - 1))) %>%
  ungroup() %>%
  select(-.ymin, -.ymax)

p_paired <- ggplot(raw_combined, aes(x = Weeks, y = Abundance, fill = Weeks)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5) +
  geom_line(aes(group = Subject), color = "grey80", alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(color = Weeks), size = 2, alpha = 0.8) +
  
  stat_pvalue_manual(stat_test,
                     label = "FDR = {signif(p.adj, 2)}",
                     tip.length = 0.02,
                     hide.ns = TRUE) +
  
  facet_wrap(~ feature, scales = "free_y", ncol = n_feat) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_fill_manual(values = c("0" = "#BDC3C7", "2" = "#E64B35", "4" = "#C0392B")) +
  scale_color_manual(values = c("0" = "#7F8C8D", "2" = "#C0392B", "4" = "#96281B")) +
  theme_prism(base_size = 12) +
  labs(title = "BCAA Catabolic Potential (Paired Comparison)",
       subtitle = "P-values adjusted using Benjamini-Hochberg (FDR) method",
       x = "Study Week", y = "Relative Abundance (Log10)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, lineheight = 1.05))

# --- 4. Save ---
FIG_W <- max(5, 4.2 * n_feat)
FIG_H <- 5.5
ggsave(file.path(RESULT_DIR, "Fig_S10.png"), p_paired, width = FIG_W, height = FIG_H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S10.svg"), p_paired, width = FIG_W, height = FIG_H, device = "svg")
print(p_paired)
print(stat_test %>% select(feature, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif))

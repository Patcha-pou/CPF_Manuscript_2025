# ==============================================================================
# SCRIPT: Figure S9 - Temporal changes in BSCFA-associated taxa
# ==============================================================================
pacman::p_load(readxl, dplyr, tidyr, tibble, ggplot2, scales, svglite)

# --- 1. Path Configuration ---
BASE_DIR <- "//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025"
if (dir.exists(BASE_DIR)) setwd(BASE_DIR)

INPUT_DIR  <- "input/"
RESULT_DIR <- "result/"
FILE_PATH  <- file.path(INPUT_DIR, "CPF_Downstream_analysis.xlsx")
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

microbe_df <- read_excel(FILE_PATH, sheet = "Microbes (Rarefied)")
clinical   <- read_excel(FILE_PATH, sheet = "Clinical")

TARGET_GENERA <- c("Clostridium_sensu_stricto_1",
                   "Christensenellaceae_R-7_group",
                   "Lachnospiraceae_FCS020_group",
                   "Lachnospiraceae_UCG-010",
                   "Roseburia")
TREAT_LEVELS  <- c("Asiatica", "Fiber_added", "Placebo")

# --- 2. Genus counts -> long format -------------------------------------------
mic_genus <- microbe_df %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

missing <- setdiff(TARGET_GENERA, mic_genus$Genus)
if (length(missing) > 0)
  stop("Genus not found in 'Microbes (Rarefied)': ", paste(missing, collapse = ", "),
       "\nAvailable close matches: ",
       paste(grep(paste(sub("_.*", "", missing), collapse = "|"),
                  mic_genus$Genus, value = TRUE), collapse = ", "))

plot_df <- mic_genus %>%
  filter(Genus %in% TARGET_GENERA) %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Count") %>%
  inner_join(clinical %>%
               mutate(SampleID = as.character(SampleID)) %>%
               select(SampleID, Treatments, Weeks),
             by = "SampleID") %>%
  mutate(
    Genus      = factor(Genus, levels = TARGET_GENERA),
    Treatments = factor(Treatments, levels = TREAT_LEVELS),
    Weeks      = factor(as.character(Weeks), levels = c("0", "2", "4"))
  )

n_zero <- sum(plot_df$Count == 0)
if (n_zero > 0)
  message("Dropped ", n_zero, " zero counts (", round(100 * n_zero / nrow(plot_df), 1),
          "%) — not representable on a log10 axis")
plot_df <- plot_df %>% filter(Count > 0)

# --- 3. Plot ------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = Weeks, y = Count, fill = Genus)) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.35, colour = "grey20") +
  facet_grid(Treatments ~ Genus) +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = label_number(accuracy = 1, big.mark = "")) +
  labs(x = "Timeline", y = expression(log[10]~"(counts)")) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey85", colour = "grey40"),
    strip.text       = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

# --- 4. Save ------------------------------------------------------------------
ggsave(file.path(RESULT_DIR, "Fig_S9.png"), p, width = 13, height = 7, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S9.svg"), p, width = 13, height = 7, device = "svg")

print(p)

plot_df %>%
  group_by(Treatments, Genus, Weeks) %>%
  summarise(n = n(), median = median(Count), .groups = "drop") %>%
  as.data.frame() %>%
  print()

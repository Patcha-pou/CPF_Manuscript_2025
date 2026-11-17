# ==============================================================================
# SCRIPT: Phyloseq preparation (Filtered: Complete Subjects Only)
# ==============================================================================

# --- 1. Load Libraries ---
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(phyloseq, dplyr)

# --- 2. Set Paths ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR          <- "input/"
OUTPUT_DIR         <- "output/"
RESULT_DIR         <- "result/"
INPUT_RARE         <- file.path(INPUT_DIR, "table.from_biom_w_taxonomy.txt")
INPUT_META         <- file.path(INPUT_DIR, "METADATA.txt")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 3. Load and Process Data ---

# 3.1 Metadata: Filter for Complete Subjects (W0, W2, W4)
sample_df <- read.delim(INPUT_META, row.names = 1, comment.char = "#", check.names = FALSE)

complete_subjects <- sample_df %>%
  group_by(Subject) %>%
  filter(n_distinct(Weeks) == 3) %>%
  pull(Subject) %>% unique()

sample_df <- sample_df %>% filter(Subject %in% complete_subjects)

# 3.2 Rarefied ASV table (8,000 reads per sample; same samples as the genus table used by MaAsLin2)
otu_rare <- read.delim(INPUT_RARE, row.names = 1, check.names = FALSE)
otu_rare <- as.matrix(otu_rare[, sapply(otu_rare, is.numeric)])
missing_samples <- setdiff(rownames(sample_df), colnames(otu_rare))
if (length(missing_samples) > 0) stop("Samples missing from ASV table: ", paste(missing_samples, collapse = ", "))
otu_rare <- otu_rare[, rownames(sample_df)]

# --- 4. Create & Save Phyloseq ---
ps_rare_final <- phyloseq(otu_table(otu_rare, taxa_are_rows = TRUE), sample_data(sample_df))

saveRDS(ps_rare_final, file = file.path(OUTPUT_DIR, "ps_rare.rds"))

message("Filtered to ", length(complete_subjects), " complete subjects, ", nsamples(ps_rare_final), " samples.")

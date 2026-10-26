# ==============================================================================
# SCRIPT: Phyloseq preparation (Filtered: Complete Subjects Only)
# ==============================================================================

# --- 1. Load Libraries ---
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(phyloseq, vegan, data.table, ggplot2, tidyr, dplyr, tibble, ape)

# --- 2. Set Paths ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR          <- "input/"
OUTPUT_DIR         <- "output/"
RESULT_DIR         <- "result/"
INPUT_RARE         <- file.path(INPUT_DIR, "feature-table.txt")
INPUT_TAXA         <- file.path(INPUT_DIR, "taxonomy.tsv")
INPUT_META         <- file.path(INPUT_DIR, "METADATA.txt")
INPUT_TREE         <- file.path(INPUT_DIR, "tree.nwk")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 3. Helper Functions ---
load_otu <- function(path) {
  df <- read.delim(path, skip = 1, check.names = FALSE)
  rownames(df) <- df[,1]
  df <- df[,-1]
  df <- df[, !colnames(df) %in% "ConsensusLineage"]
  return(as.matrix(df))
}

clean_taxa <- function(tax_mat) {
  tax_mat <- gsub(".*__", "", tax_mat)
  tax_mat <- trimws(tax_mat)
  unassigned_patterns <- c("", "NA", "None", "uncultured", "unassigned", "_", "__")
  for (i in 1:nrow(tax_mat)) {
    for (j in 2:ncol(tax_mat)) {
      curr <- tax_mat[i, j]
      if (is.na(curr) || curr %in% unassigned_patterns || nchar(curr) <= 1) {
        prev_rank <- tax_mat[i, j-1]
        tax_mat[i, j] <- if(grepl("_unclassified", prev_rank)) prev_rank else paste0(prev_rank, "_unclassified")
      }
    }
  }
  return(tax_mat)
}

# --- 4. Load and Process Data ---

# 4.1 Metadata: Filter for Complete Subjects (W0, W2, W4)
sample_df <- read.delim(INPUT_META, row.names = 1, comment.char = "#", check.names = FALSE)

complete_subjects <- sample_df %>%
  group_by(Subject) %>%
  filter(n_distinct(Weeks) == 3) %>%
  pull(Subject) %>% unique()

sample_df <- sample_df %>% filter(Subject %in% complete_subjects)
META <- sample_data(sample_df)

# 4.2 OTU Table & Taxonomy
otu_rare <- load_otu(INPUT_RARE)
otu_rare <- otu_rare[, intersect(colnames(otu_rare), rownames(sample_df))] 

tax_raw  <- read.delim(INPUT_TAXA, check.names = FALSE)
tax_mat_asv <- tax_raw %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";", fill = "right") %>%
  column_to_rownames(var = "Feature ID") %>% 
  dplyr::select(-Confidence) %>% 
  as.matrix() %>% 
  clean_taxa()

phytree <- ape::read.tree(file = INPUT_TREE)

# --- 5. Create & Save Phyloseq ---
ps_rare_final <- phyloseq(otu_table(otu_rare, taxa_are_rows = TRUE), 
                          tax_table(tax_mat_asv), 
                          META, phytree)

saveRDS(ps_rare_final, file = file.path(OUTPUT_DIR, "ps_rare.rds"))

message(paste("Filtered to", length(complete_subjects), "complete subjects."))
message("Success: All samples without full timepoints (W0,W2,W4) removed.")
# ==============================================================================
# SCRIPT: MaAsLin2 — one run per omics layer x treatment arm
# ==============================================================================
library(Maaslin2)
library(missForest)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- 1. Metadata (shared by both layers) -------------------------------------
input_metadata <- read.csv(file.path(INPUT_DIR, "input_metadata.csv"), row.names = 1)

# --- 2. Feature tables -------------------------------------------------------
mic_data <- read.csv(file.path(INPUT_DIR, "input_data_rel.csv"), row.names = 1,
                     check.names = FALSE)

met_data <- read.csv(file.path(INPUT_DIR, "input_data_raw.csv"), row.names = 1,
                     check.names = FALSE)

# --- 3. Helpers --------------------------------------------------------------
to_numeric <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df
}

impute_if_needed <- function(df, seed = 123) {
  if (!anyNA(df)) return(df)
  set.seed(seed)
  missForest(df)$ximp
}

ARMS <- c("Placebo", "Asiatica", "Fiber_added")

run_layer <- function(data, layer) {
  data <- impute_if_needed(to_numeric(data))

  for (arm in ARMS) {
    ids <- rownames(input_metadata)[input_metadata$Treatments == arm]
    ids <- intersect(ids, rownames(data))

    out <- file.path(sub("/+$", "", OUTPUT_DIR), layer, paste0("output_", arm))
    if (!dir.exists(out)) dir.create(out, recursive = TRUE)

    message(layer, " / ", arm, ": n = ", length(ids), " -> ", out)

    Maaslin2(
      data.frame(data[ids, , drop = FALSE]),
      data.frame(input_metadata[ids, , drop = FALSE]),
      out,
      transform = "LOG",
      fixed_effects = c('Weeks', 'Age', 'BMI', 'Gender', 'Energy'),
      random_effects = c('Subject'),
      normalization = 'NONE',
      reference = 'Weeks,W0',
      standardize = FALSE)
  }
}

# --- 4. Run ------------------------------------------------------------------
run_layer(mic_data, "Mic")
run_layer(met_data, "Met")

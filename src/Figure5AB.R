# ==============================================================================
# SCRIPT: Figure 5A + 5B (UpSet plots, left column of the Figure 5 layout)
# ==============================================================================
library(ggplot2)
library(ComplexUpset)
library(patchwork)
library(readxl)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
TABLE_S5   <- file.path(OUTPUT_DIR, "Table_S5.xlsx")
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

P_CUT <- 0.05

FONT_BASE     <- 14
COUNT_SIZE_MM <- FONT_BASE / ggplot2::.pt
ITALIC_TAXON  <- "C. asiatica"
Y_LAB         <- "Number of correlations"

set_labeller <- function(x) {
  pattern <- gsub("\\.", "\\\\.", ITALIC_TAXON)
  out <- vapply(as.character(x), function(s) {
    if (!grepl(pattern, s)) return(sprintf("'%s'", s))
    prefix <- trimws(sub(paste0(pattern, ".*$"), "", s))
    if (nzchar(prefix)) {
      sprintf("paste('%s ', italic('%s'))", prefix, ITALIC_TAXON)
    } else {
      sprintf("italic('%s')", ITALIC_TAXON)
    }
  }, character(1), USE.NAMES = FALSE)
  parse(text = out)
}

upset_themes_md <- upset_default_themes(
  text      = element_text(size = FONT_BASE),
  axis.text = element_text(size = FONT_BASE)
)

GROUP_LABELS <- c(
  Placebo     = "Placebo",
  Asiatica    = "C. asiatica",
  Fiber_added = "Fiber-enriched C. asiatica"
)
target_sets <- unname(GROUP_LABELS)

build_membership <- function(week) {
  sheets_avail <- excel_sheets(TABLE_S5)
  
  sig <- lapply(names(GROUP_LABELS), function(g) {
    sheet <- paste0("Table_S5_Full_W", week, "_", g)
    if (!sheet %in% sheets_avail)
      stop("missing sheet '", sheet, "' in ", TABLE_S5,
           " — run FigureS7S8_TableS5.R first", call. = FALSE)
    df <- read_excel(TABLE_S5, sheet = sheet)
    df <- df[!is.na(df$p) & df$p < P_CUT, ]
    paste(df$Microbe, df$Metabolite, sep = " ~ ")
  })
  names(sig) <- names(GROUP_LABELS)
  
  pairs <- sort(unique(unlist(sig)))
  mat <- as.data.frame(
    lapply(sig, function(s) as.integer(pairs %in% s)),
    check.names = FALSE
  )
  colnames(mat) <- unname(GROUP_LABELS[names(sig)])
  rownames(mat) <- pairs
  
  message("Week ", week, ": ", nrow(mat), " unique pairs | ",
          paste(sprintf("%s=%d", colnames(mat), colSums(mat)), collapse = ", "))
  mat
}

top_data    <- build_membership(2)
bottom_data <- build_membership(4)

create_final_upset <- function(data, sets_vector) {
  upset(
    data,
    sets_vector,
    width_ratio = 0.2,
    sort_intersections = 'ascending',
    sort_sets = FALSE,
    labeller = set_labeller,
    base_annotations = list(
      'Intersection Size' = intersection_size(
        counts = TRUE,
        fill = '#333333',
        bar_number_threshold = 1,
        text = list(size = COUNT_SIZE_MM)
      ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        ylab(Y_LAB) +
        theme(
          panel.grid = element_blank(),
          axis.line.y = element_line(color = "black")
        )
    ),
    set_sizes = (
      upset_set_size() +
        geom_bar(fill = "#333333") +
        labs(x = "Set Size (Total Associations)") +
        theme(axis.line.x = element_line(color = "black"))
    ),
    stripes = upset_stripes(colors = c("grey95", "white")),
    themes = upset_themes_md
  )
}

p1 <- create_final_upset(top_data, target_sets)
p2 <- create_final_upset(bottom_data, target_sets)


ggsave(file.path(RESULT_DIR, "Figure5A.png"), p1, width = 10, height = 7, dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure5A.svg"), p1, width = 10, height = 7)

ggsave(file.path(RESULT_DIR, "Figure5B.png"), p2, width = 10, height = 7, dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure5B.svg"), p2, width = 10, height = 7)

PANEL_DIR <- file.path(OUTPUT_DIR, "panels")
if (!dir.exists(PANEL_DIR)) dir.create(PANEL_DIR, recursive = TRUE)
saveRDS(p1, file.path(PANEL_DIR, "Figure5A.rds"))
saveRDS(p2, file.path(PANEL_DIR, "Figure5B.rds"))

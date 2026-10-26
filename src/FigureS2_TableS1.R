# ==============================================================================
# SCRIPT: Figure S2 + Table S1 - Stool frequency and consistency
# ==============================================================================

# --- 1. Load Libraries ---
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lmerTest)
library(openxlsx)
library(svglite)

# --- 2. Set Paths ---
BASE_DIR <- "//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025"
if (dir.exists(BASE_DIR)) setwd(BASE_DIR)
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 3. Load & Prep Data ---
group_map <- c("Placebo" = "Placebo",
               "Asiatica_juice" = "C. asiatica",
               "Fiber_added_Asiatica_juice" = "fiber-enriched C. asiatica")

color_palette <- c("Placebo" = "#999999",
                   "C. asiatica" = "#e6c60d",
                   "fiber-enriched C. asiatica" = "#5ca08e")

italic_labels <- c(
  "Placebo" = "Placebo",
  "C. asiatica" = expression(italic("C. asiatica")),
  "fiber-enriched C. asiatica" = expression(paste("fiber-enriched ", italic("C. asiatica")))
)

VARS <- c(Frequency = "Stool frequency", Consistency = "Stool consistency (Bristol score)")

stool <- read_excel(file.path(INPUT_DIR, "inputMet_with_enterotype_stool.xlsx"), sheet = "METADATA") %>%
  transmute(ID, Subject = factor(Subject),
            Group = factor(group_map[Treatments], levels = names(color_palette)),
            Week = factor(paste0("W", Weeks), levels = c("W0", "W2", "W4")),
            Frequency, Consistency)
stopifnot(!anyNA(stool$Group))

stool_long <- stool %>%
  pivot_longer(all_of(names(VARS)), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = names(VARS)))

# --- 4. Statistics (LMM, participant as random effect) ---
tidy_coef <- function(fit, terms) {
  cf <- summary(fit)$coefficients
  cf <- cf[rownames(cf) %in% terms, , drop = FALSE]
  data.frame(term = rownames(cf), Estimate = cf[, "Estimate"], SE = cf[, "Std. Error"],
             CI95_lower = cf[, "Estimate"] - 1.96 * cf[, "Std. Error"],
             CI95_upper = cf[, "Estimate"] + 1.96 * cf[, "Std. Error"],
             df = cf[, "df"], p_value = cf[, "Pr(>|t|)"], row.names = NULL)
}

descriptive <- stool_long %>%
  group_by(Variable, Group, Week) %>%
  summarise(N = sum(!is.na(Value)), Median = median(Value, na.rm = TRUE),
            Q1 = quantile(Value, 0.25, na.rm = TRUE), Q3 = quantile(Value, 0.75, na.rm = TRUE),
            Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop")

within_group <- stool_long %>%
  group_by(Variable, Group) %>%
  group_modify(~ {
    fit <- lmer(Value ~ Week + (1 | Subject), data = .x)
    tidy_coef(fit, c("WeekW2", "WeekW4")) %>%
      mutate(Comparison = sub("Week(W\\d)", "\\1 vs W0", term), N_participants = n_distinct(.x$Subject))
  }) %>%
  ungroup() %>%
  group_by(Variable) %>% mutate(q_value = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  select(Variable, Group, Comparison, Estimate, SE, CI95_lower, CI95_upper, df, p_value, q_value, N_participants)

fit_did <- function(d, ref) {
  d$Group <- relevel(droplevels(d$Group), ref = ref)
  fit <- lmer(Value ~ Group * Week + (1 | Subject), data = d)
  int <- grep("^Group.*:Week", rownames(summary(fit)$coefficients), value = TRUE)
  tidy_coef(fit, int) %>%
    mutate(Contrast = paste(sub("^Group(.*):Week.*$", "\\1", term), "vs", ref),
           Week = sub("^.*:Week", "", term))
}

between_group <- stool_long %>%
  group_by(Variable) %>%
  group_modify(~ bind_rows(fit_did(.x, "Placebo"),
                           fit_did(.x, "C. asiatica") %>% filter(grepl("^fiber-enriched", Contrast)))) %>%
  ungroup() %>%
  group_by(Variable) %>% mutate(q_value = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  mutate(N_participants = n_distinct(stool$Subject)) %>%
  select(Variable, Contrast, Week, DiD = Estimate, SE, CI95_lower, CI95_upper, df, p_value, q_value, N_participants)

overall <- stool_long %>%
  group_by(Variable) %>%
  group_modify(~ {
    a <- anova(lmer(Value ~ Group * Week + (1 | Subject), data = .x))
    data.frame(Term = rownames(a), as.data.frame(a), row.names = NULL, check.names = FALSE)
  }) %>%
  ungroup()

# --- 5. Table S1 ---
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", halign = "center",
                  fgFill = "#EDEDED", border = "bottom")
add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, as.data.frame(df), headerStyle = hs, withFilter = TRUE)
  freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet(wb, "Descriptive", descriptive)
add_sheet(wb, "Within_group", within_group)
add_sheet(wb, "Between_group_DiD", between_group)
add_sheet(wb, "Overall_LMM", overall)
saveWorkbook(wb, file.path(OUTPUT_DIR, "Table_S1.xlsx"), overwrite = TRUE)
cat("Wrote", file.path(OUTPUT_DIR, "Table_S1.xlsx"), "\n")

# --- 6. Figure S2 ---
FONT <- "Arial"
make_panel <- function(var) {
  ggplot(filter(stool_long, Variable == var), aes(x = Week, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8, linewidth = 0.3,
                 position = position_dodge(width = 0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.4, size = 0.6, color = "black") +
    scale_fill_manual(values = color_palette, labels = italic_labels) +
    labs(x = "Week", y = VARS[[var]], fill = NULL) +
    theme_bw(base_size = 7, base_family = FONT) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_text(size = 8, face = "bold"),
          axis.text = element_text(size = 7, colour = "black"),
          legend.text = element_text(size = 7), legend.position = "bottom")
}

p_final <- (make_panel("Frequency") | make_panel("Consistency")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"),
        plot.tag = element_text(size = 9, face = "bold", family = FONT))

W <- 180 / 25.4; H <- 70 / 25.4
ggsave(file.path(RESULT_DIR, "Fig_S2.png"), p_final, width = W, height = H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S2.svg"), p_final, width = W, height = H,
       device = svglite, system_fonts = list(Arial = "Arial"))
message("Saved Fig_S2.{png,svg}")

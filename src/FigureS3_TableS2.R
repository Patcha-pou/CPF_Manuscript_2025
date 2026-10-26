# ==============================================================================
# SCRIPT: Figure S3 + Table S2 - Alpha and beta diversity
# ==============================================================================

# --- 1. Load Libraries ---
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(phyloseq, ggplot2, dplyr, tidyr, RColorBrewer, scales,
               vegan, permute, patchwork, lmerTest, ggpubr, openxlsx, svglite)

# --- 2. Set Paths ---
setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR           <- "input/"
OUTPUT_DIR          <- "output/"
RESULT_DIR          <- "result/"
INPUT_RDS_RARE      <- file.path(OUTPUT_DIR, "ps_rare.rds")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if(!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# --- 3. Load & Prep Data ---
ps_rare      <- readRDS(INPUT_RDS_RARE)

group_map <- c("Placebo" = "Placebo", 
               "Asiatica" = "C. asiatica", 
               "Fiber_added" = "fiber-enriched C. asiatica")

color_palette <- c("Placebo" = "#999999", 
                   "C. asiatica" = "#e6c60d", 
                   "fiber-enriched C. asiatica" = "#5ca08e")

tw_order <- c(
  "Placebo_W0", "Placebo_W2", "Placebo_W4",
  "Asiatica_W0", "Asiatica_W2", "Asiatica_W4",
  "Fiber_added_W0", "Fiber_added_W2", "Fiber_added_W4"
)

process_ps_data <- function(ps) {
  sd <- data.frame(sample_data(ps))
  sd$Treatments <- factor(group_map[as.character(sd$Treatments)], levels = names(color_palette))
  sd$Weeks <- factor(as.character(sd$Weeks), levels = c("W0", "W2", "W4"))
  sd$Treatments_Weeks <- factor(as.character(sd$Treatments_Weeks), levels = tw_order)
  sample_data(ps) <- sample_data(sd)
  return(ps)
}

ps_rare <- process_ps_data(ps_rare)

italic_labels <- c(
  "Placebo" = "Placebo",
  "C. asiatica" = expression(italic("C. asiatica")),
  "fiber-enriched C. asiatica" = expression(paste("fiber-enriched ", italic("C. asiatica")))
)

# --- 4. Alpha Diversity Calculation & Table S2 Export ---
alpha_meas <- estimate_richness(ps_rare, measures = c("Shannon", "Observed"))
df_alpha <- data.frame(sample_data(ps_rare))
df_alpha$Shannon <- alpha_meas$Shannon
df_alpha$Evenness <- alpha_meas$Shannon / log(alpha_meas$Observed)

# 4.1 Linear mixed models with participant as random effect (repeated measures)
df_alpha$Subject <- factor(df_alpha$Subject)
tidy_coef <- function(fit, terms) {
  cf <- summary(fit)$coefficients
  cf <- cf[rownames(cf) %in% terms, , drop = FALSE]
  data.frame(term = rownames(cf), Estimate = cf[, "Estimate"], SE = cf[, "Std. Error"],
             CI95_lower = cf[, "Estimate"] - 1.96 * cf[, "Std. Error"],
             CI95_upper = cf[, "Estimate"] + 1.96 * cf[, "Std. Error"],
             df = cf[, "df"], p_value = cf[, "Pr(>|t|)"], row.names = NULL)
}
alpha_long <- bind_rows(lapply(c("Shannon", "Evenness"), function(v)
  data.frame(Measure = v, Subject = df_alpha$Subject, Group = df_alpha$Treatments,
             Week = df_alpha$Weeks, Value = df_alpha[[v]])))

alpha_within <- alpha_long %>%
  group_by(Measure, Group) %>%
  group_modify(~ tidy_coef(lmer(Value ~ Week + (1 | Subject), data = .x), c("WeekW2", "WeekW4")) %>%
                 mutate(Comparison = sub("Week(W\\d)", "\\1 vs W0", term), N_participants = n_distinct(.x$Subject))) %>%
  ungroup() %>%
  group_by(Measure) %>% mutate(q_value = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  select(Measure, Group, Comparison, Estimate, SE, CI95_lower, CI95_upper, df, p_value, q_value, N_participants)

fit_did <- function(d, ref) {
  d$Group <- relevel(droplevels(d$Group), ref = ref)
  fit <- lmer(Value ~ Group * Week + (1 | Subject), data = d)
  int <- grep("^Group.*:Week", rownames(summary(fit)$coefficients), value = TRUE)
  tidy_coef(fit, int) %>%
    mutate(Contrast = paste(sub("^Group(.*):Week.*$", "\\1", term), "vs", ref), Week = sub("^.*:Week", "", term))
}
alpha_did <- alpha_long %>%
  group_by(Measure) %>%
  group_modify(~ bind_rows(fit_did(.x, "Placebo"),
                           fit_did(.x, "C. asiatica") %>% filter(grepl("^fiber-enriched", Contrast)))) %>%
  ungroup() %>%
  group_by(Measure) %>% mutate(q_value = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  mutate(N_participants = n_distinct(alpha_long$Subject)) %>%
  select(Measure, Contrast, Week, DiD = Estimate, SE, CI95_lower, CI95_upper, df, p_value, q_value, N_participants)

alpha_overall <- alpha_long %>%
  group_by(Measure) %>%
  group_modify(~ { a <- anova(lmer(Value ~ Group * Week + (1 | Subject), data = .x))
                   data.frame(Term = rownames(a), as.data.frame(a), row.names = NULL, check.names = FALSE) }) %>%
  ungroup()

# 4.2 Save statistics to Table_S2.xlsx
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", halign = "center",
                  fgFill = "#EDEDED", border = "bottom")

add_sheet <- function(wb, name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs, withFilter = TRUE)
  freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

add_sheet(wb, "Alpha_LMM_within", alpha_within)
add_sheet(wb, "Alpha_LMM_DiD", alpha_did)
add_sheet(wb, "Alpha_LMM_overall", alpha_overall)


p_shannon <- ggplot(df_alpha, aes(x = Treatments_Weeks, y = Shannon, color = Treatments)) +
  geom_boxplot(aes(fill = Treatments), alpha = 0.2, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.4, size = 1.2) +
  scale_color_manual(values = color_palette, labels = italic_labels) +
  scale_fill_manual(values = color_palette, labels = italic_labels) +
  theme_bw() + labs(subtitle = "Shannon Index", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), 
        legend.position = "none")

p_evenness <- ggplot(df_alpha, aes(x = Treatments_Weeks, y = Evenness, color = Treatments)) +
  geom_boxplot(aes(fill = Treatments), alpha = 0.2, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.4, size = 1.2) +
  scale_color_manual(values = color_palette, labels = italic_labels) +
  scale_fill_manual(values = color_palette, labels = italic_labels) +
  theme_bw() + labs(subtitle = "Evenness", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), 
        legend.position = "none")

dist_matrix <- phyloseq::distance(ps_rare, method = "bray")
ord_pcoa    <- ordinate(ps_rare, method = "PCoA", distance = dist_matrix)

# PERMANOVA with permutations restricted by participant: whole participants are permuted
# (tests Treatments) and time points are permuted within participant (tests Weeks and interaction).
# Restricted permutation needs equal samples per participant, so participants with all three
# time points are used.
sd_all   <- data.frame(sample_data(ps_rare))
complete <- names(which(table(sd_all$Subject) == 3))
ps_perm  <- prune_samples(as.character(sd_all$Subject) %in% complete, ps_rare)
sd_perm  <- data.frame(sample_data(ps_perm)); sd_perm$Subject <- factor(sd_perm$Subject)
sd_perm  <- sd_perm[order(sd_perm$Subject, sd_perm$Weeks), ]
dist_perm <- as.dist(as.matrix(phyloseq::distance(ps_perm, method = "bray"))[rownames(sd_perm), rownames(sd_perm)])
ctrl <- how(within = Within(type = "free"), plots = Plots(strata = sd_perm$Subject, type = "free"), nperm = 999)
set.seed(42)
perm_model <- adonis2(dist_perm ~ Treatments * Weeks, data = sd_perm, permutations = ctrl, by = NULL)
set.seed(42)
perm_terms <- adonis2(dist_perm ~ Treatments * Weeks, data = sd_perm, permutations = ctrl, by = "terms")
p_val_perm <- perm_model$`Pr(>F)`[1]

# 4.3 Add beta diversity sheets and save Table_S2.xlsx
beta_n <- data.frame(Note = c("Participants", "Samples", "Permutations"),
                     Value = c(length(complete), nrow(sd_perm), "999, restricted by participant"))
add_sheet(wb, "Beta_PERMANOVA_model", data.frame(Term = rownames(perm_model), as.data.frame(perm_model), row.names = NULL))
add_sheet(wb, "Beta_PERMANOVA_terms", data.frame(Term = rownames(perm_terms), as.data.frame(perm_terms), row.names = NULL))
add_sheet(wb, "Beta_PERMANOVA_n", beta_n)
saveWorkbook(wb, file.path(OUTPUT_DIR, "Table_S2.xlsx"), overwrite = TRUE)

p_beta <- plot_ordination(ps_rare, ord_pcoa, color = "Treatments", shape = "Weeks") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Treatments, color = Treatments), linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = color_palette, labels = italic_labels) +
  scale_shape_manual(values = c("W0" = 16, "W2" = 17, "W4" = 15)) +
  theme_bw() + labs(subtitle = "PCoA (Bray-Curtis)") +
  annotate("text", x = -Inf, y = -Inf, label = paste0("PERMANOVA p = ", round(p_val_perm, 4)), 
           hjust = -0.1, vjust = -1.1, size = 4, fontface = "italic")

# --- 6. Assembly & Save ---
p_final <- (p_shannon | p_evenness | p_beta) + 
  plot_layout(widths = c(1, 1, 1.4)) + 
  plot_annotation(tag_levels = 'A')

ggsave(file.path(RESULT_DIR, "Fig_S3.png"), plot = p_final, width = 17, height = 7, dpi = 600)
ggsave(file.path(RESULT_DIR, "Fig_S3.svg"), plot = p_final, width = 17, height = 7, device = "svg")

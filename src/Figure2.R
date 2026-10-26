# Figure 2 - diet and blood chemistry
library(tidyverse)
library(reshape2)
library(imputeMissings)
library(ggpubr)
library(patchwork)
library(svglite)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

# ==========================================
# 1. SETUP & UTILITIES
# ==========================================

group_map <- c("Control" = "Placebo", 
               "Asiatica" = "C. asiatica", 
               "Fiber added" = "fiber-enriched C. asiatica")

color_palette <- c("Placebo" = "#999999", 
                   "C. asiatica" = "#e6c60d", 
                   "fiber-enriched C. asiatica" = "#5ca08e")

FONT <- "Arial"
theme_fig <- theme_bw(base_size = 7, base_family = FONT) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "gray", color = "black"),
        strip.text = element_text(face = "bold", size = 7, margin = margin(1, 1, 1, 1)),
        axis.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7, colour = "black"),
        plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(0, 0, 1, 4, "mm")),
        plot.title.position = "plot",
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        legend.position = "bottom")

remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 3 * IQR(x, na.rm = TRUE)
  y <- x
  y[x < (qnt[1] - H) | x > (qnt[2] + H)] <- NA
  y
}

tidy_param <- function(data, cols, label) {
  df <- data[, c(1, 2, cols)]
  colnames(df) <- c("ID", "Group", "0", "2", "4")
  df_long <- melt(df, id.vars = c("ID", "Group"),
                  variable.name = "Week", value.name = "Conc.")
  df_long$Parameter <- label
  df_long$Group <- factor(group_map[df_long$Group], levels = names(color_palette))
  return(df_long)
}

# ==========================================
# 2. DIET DATA
# ==========================================

diet_refs <- data.frame(
  Parameter = c("Energy", "Protein", "Fat", "Carbohydrate", "Sugar", "Dietary fiber"),
  intercept = c(2000, 50, 62.5, 300, 50, 25),
  unit = c("Kcal", "g", "g", "g", "g", "g")
)

Food <- read.csv(file.path(INPUT_DIR, "Food.csv"))
Food[Food == 0] <- NA
Food[, 3:20] <- impute(Food[, 3:20], method = "randomForest")

diet_long <- bind_rows(
  tidy_param(Food, c(3, 9, 15), "Energy"),
  tidy_param(Food, c(4, 10, 16), "Protein"),
  tidy_param(Food, c(5, 11, 17), "Fat"),
  tidy_param(Food, c(6, 12, 18), "Carbohydrate"),
  tidy_param(Food, c(7, 13, 19), "Sugar"),
  tidy_param(Food, c(8, 14, 20), "Dietary fiber")
)

diet_long <- diet_long %>%
  group_by(Parameter) %>%
  mutate(Conc. = remove_outliers(Conc.)) %>%
  ungroup() %>%
  left_join(diet_refs %>% select(Parameter, unit, intercept), by = "Parameter") %>%
  mutate(Parameter_label = paste0(Parameter, "\n(", unit, ")"))

diet_fig <- ggplot(diet_long, aes(x = Week, y = Conc., fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.4, size = 0.3, color = "black") +
  geom_hline(aes(yintercept = intercept), color = "red", linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~Parameter_label, ncol = 6, scales = "free_y") +
  scale_fill_manual(
    values = color_palette,
    labels = c("Placebo", 
               expression(italic("C. asiatica")), 
               expression(paste("fiber-enriched ", italic("C. asiatica"))))
  ) +
  theme_fig +
  labs(title = "Diet over weeks by treatment group", y = "Amount Intake", x = "Week", fill = "Treatments")

# ==========================================
# 3. BLOOD CHEMISTRY
# ==========================================

BloodTest <- read.csv(file.path(INPUT_DIR, "BloodTest.csv"))

blood_units <- data.frame(
  Parameter = c("LDL", "TG", "FPG", "TC", "HDL"),
  ref_max = c(130, 150, 99, 200, 50),
  ref_min = c(NA, NA, 70, NA, NA),
  unit = "mg/dL"
)

blood_long <- bind_rows(
  tidy_param(BloodTest, 15:17, "LDL"),
  tidy_param(BloodTest, 12:14, "TG"),
  tidy_param(BloodTest, 3:5, "FPG"),
  tidy_param(BloodTest, 6:8, "TC"),
  tidy_param(BloodTest, 9:11, "HDL")
) %>%
  left_join(blood_units, by = "Parameter") %>%
  mutate(Parameter_label = paste0(Parameter, "\n(", unit, ")"))

blood_fig <- ggplot(blood_long, aes(x = Week, y = Conc., fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.8, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.4, size = 0.3, color = "black") +
  geom_hline(aes(yintercept = ref_max), color = "red", linetype = "dashed", linewidth = 0.4) +
  geom_hline(data = filter(blood_long, Parameter == "FPG"), 
             aes(yintercept = 70), color = "blue", linetype = "dotted", linewidth = 0.4) +
  facet_wrap(~Parameter_label, ncol = 5, scales = "free_y") +
  scale_fill_manual(
    values = color_palette,
    labels = c("Placebo", 
               expression(italic("C. asiatica")), 
               expression(paste("fiber-enriched ", italic("C. asiatica"))))
  ) +
  theme_fig +
  labs(title = "Blood chemistry over weeks by treatment group", y = "Conc.", x = "Week", fill = "Treatments")


# ==========================================
# 4. OUTPUT
# ==========================================

Figure2 <- ( blood_fig / diet_fig ) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        plot.tag = element_text(size = 9, face = "bold", family = FONT))

# Panel size as placed in Figures_EDITED pptx at 180 mm print width
W <- 180 / 25.4; H <- 101.2 / 25.4
ggsave(file.path(RESULT_DIR, "Figure2.png"), plot = Figure2, width = W, height = H, dpi = 600)
ggsave(file.path(RESULT_DIR, "Figure2.svg"), plot = Figure2, width = W, height = H,
       device = svglite, system_fonts = list(Arial = "Arial"))
# Dietary fiber intake by country, grouped by region -------------------------

library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)

setwd("//wsl.localhost/Ubuntu/home/somsri/CPF/Manuscript_2025")
INPUT_DIR  <- "input/"
OUTPUT_DIR <- "output/"
RESULT_DIR <- "result/"
if (!dir.exists(RESULT_DIR)) dir.create(RESULT_DIR, recursive = TRUE)

file   <- file.path(INPUT_DIR, "Fiber.xlsx")
outdir <- RESULT_DIR
REC  <- 25
THAI <- "#2a78d6"

# --- 1. Read raw data from the workbook ------------------------------------
raw <- read_excel(file, sheet = "Sheet1")

dat <- raw %>%
  transmute(
    country_raw = .[[1]],
    intake      = as.numeric(`Fiber intake (g/day)`)
  ) %>%
  filter(!is.na(intake))

# --- 2. Clean country names ------------------------------------------------
rename_map <- c("Thai" = "Thailand", "The Netherlands" = "Netherlands")
dat <- dat %>%
  mutate(country = str_squish(str_replace(country_raw, "\\s*\\(.*$", "")),
         country = recode(country, !!!rename_map))

region_of <- c(
  Brazil = "Latin America", Chile = "Latin America", Mexico = "Latin America",
  Thailand = "Asia", Japan = "Asia", Korea = "Asia", China = "Asia",
  Singapore = "Asia", India = "Asia",
  Sweden = "Europe", Norway = "Europe", France = "Europe", Greece = "Europe",
  `United Kingdom` = "Europe", Germany = "Europe", Italy = "Europe",
  Netherlands = "Europe", Spain = "Europe", Denmark = "Europe",
  Australia = "Oceania"
)
dat$region <- region_of[dat$country]
stopifnot(!any(is.na(dat$region)))

region_colors <- c("Latin America" = "#eda100", "Asia" = "#2a78d6",
                   "Europe" = "#1baf7a", "Oceania" = "#e87ba4")

# --- 4. Layout: regions ordered by mean; within region ordered by value ----
gap <- 1
region_levels <- dat %>% group_by(region) %>%
  summarise(m = mean(intake), .groups = "drop") %>% arrange(m) %>% pull(region)

plot_df <- data.frame(); bands <- data.frame(); y <- 0
for (r in region_levels) {
  sub <- dat %>% filter(region == r) %>% arrange(intake)
  y0 <- y
  sub$pos <- y + seq_len(nrow(sub)) - 1
  y <- y + nrow(sub)
  bands <- rbind(bands, data.frame(region = r, ymin = y0 - 0.5, ymax = y - 0.5,
                                   ymid = (y0 - 0.5 + y - 0.5) / 2))
  y <- y + gap
  plot_df <- rbind(plot_df, sub)
}
plot_df$barcol  <- ifelse(plot_df$country == "Thailand", THAI,
                          region_colors[plot_df$region])
plot_df$is_thai <- plot_df$country == "Thailand"
thai_pos <- plot_df$pos[plot_df$country == "Thailand"]

# --- 5. Plot ---------------------------------------------------------------
p <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = region),
            alpha = 0.06) +
  scale_fill_manual(values = region_colors, guide = "none") +
  geom_vline(xintercept = REC, linetype = "22", colour = "#6b6a66", linewidth = 0.5) +
  annotate("text", x = REC, y = max(plot_df$pos) + 1.1,
           label = "Recommended 25 g/day", size = 3.1, colour = "#4a4945", vjust = 0) +
  geom_rect(data = plot_df,
            aes(xmin = 0, xmax = intake, ymin = pos - 0.36, ymax = pos + 0.36),
            fill = plot_df$barcol) +
  geom_text(data = plot_df,
            aes(x = intake + 0.4, y = pos, label = sprintf("%.1f", intake),
                fontface = ifelse(is_thai, "bold", "plain")),
            hjust = 0, size = 3, colour = "#222222") +
  annotate("text", x = 11.6, y = thai_pos, label = "This study",
           hjust = 0, fontface = "italic", size = 2.9, colour = THAI) +
  geom_text(data = bands, aes(x = -11.5, y = ymid, label = region, colour = region),
            angle = 90, fontface = "bold", size = 3.3, show.legend = FALSE) +
  scale_colour_manual(values = region_colors, guide = "none") +
  scale_y_continuous(breaks = plot_df$pos, labels = plot_df$country, expand = c(0, 0.9)) +
  scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Dietary fiber intake (g/day)", y = NULL,
       title = "Dietary fiber intake falls short of recommendations across countries",
       subtitle = "Population mean intake; ordered within region by value, regions ordered by mean",
       caption = "Source: Fiber.xlsx (compiled population studies). Recommended intake 25 g/day.") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.3),
    axis.ticks         = element_blank(),
    axis.line.x        = element_line(colour = "#c3c2b7", linewidth = 0.4),
    axis.text.y        = element_text(size = 9, colour = "#333230"),
    plot.title    = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 9.5, colour = "#6b6a66"),
    plot.caption  = element_text(size = 7, colour = "#8a8984", hjust = 0),
    plot.margin   = margin(10, 14, 8, 40)
  )


ggsave(file.path(outdir, "Fig_S1.svg"), p, width = 8.4, height = 8.9)
ggsave(file.path(outdir, "Fig_S1.png"), p, width = 8.4, height = 8.9, dpi = 600)

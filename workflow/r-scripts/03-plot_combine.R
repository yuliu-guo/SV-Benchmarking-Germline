library(tidyverse)
library(here)
library(patchwork)
library(ggh4x)
library(showtext)
library(grid) # for textGrob()
library(gt)
library(svglite)

here::i_am("r-scripts/03-plot_combine.R")

main_folder <- here("pedigree-results", "truvari", "type-ignored")
fig1b_file <- here(main_folder, "plots_gt_sample.Rds")
fig1c1_file <- here(main_folder, "fig1_f1_all.Rds")
fig1c2_file <- here(main_folder, "fig1_recall_all.Rds")
fig1c3_file <- here(main_folder, "fig1_recall.Rds")
fig1c4_file <- here(main_folder, "fig1_count_split_TP.Rds")

# -------------------------
# Data + Figure A (caller metrics)
# -------------------------
source(here("r-scripts", "01-truvari-summary.R"))

data %>%
  group_by(caller) %>%
  summarize(
    f1 = mean(f1), f1_sd = sd(f1, na.rm = TRUE),
    recall = mean(recall), recall_sd = sd(recall, na.rm = TRUE),
    precision = mean(precision), precision_sd = sd(precision, na.rm = TRUE)
  )

gglayer_med_bar <- c(
  gglayer_bar_theme,
  stat_summary(aes(fill = caller), fun = "mean", geom = "col"),
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1.5, alpha = 0.5)
)

f1 <- ggplot(data, aes(x = fct_relevel(caller, f1_order), y = f1)) +
  gglayer_med_bar +
  labs(y = "F1") +
  scale_y_continuous(limits = c(0, max(data$precision)))

recall <- ggplot(data, aes(x = fct_relevel(caller, f1_order), y = recall)) +
  gglayer_med_bar +
  labs(y = "Recall") +
  scale_y_continuous(limits = c(0, max(data$precision)))

precision <- ggplot(data, aes(x = fct_relevel(caller, f1_order), y = precision)) +
  gglayer_med_bar +
  labs(y = "Precision") +
  scale_y_continuous(limits = c(0, max(data$precision)))

fig1a <- f1 / recall / precision +
  plot_layout(guides = "collect", axes = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 7, margin = margin(b = 4, l = -4), face = "bold"),
    legend.text = element_text(size = 6, hjust = -1, margin = margin(l = -4, b = -2, t = 2)),
    legend.key.spacing.x = unit(0, units = "mm"),
    legend.text.position = "top",
    legend.margin = margin(l = 5, r = -25),
    legend.box.spacing = unit(1, units = "mm")
  )

# -------------------------
# Figure B table (gt)
# -------------------------
gt_table <- readRDS(file = fig1b_file) %>% opt_vertical_padding(scale = .8)
fig1b <- gt_table %>% opt_horizontal_padding(scale = .8)

# -------------------------
# Figure B panels (SVLEN binned)
# -------------------------
if (!file.exists(fig1c1_file) | !file.exists(fig1c2_file) | !file.exists(fig1c3_file) | !file.exists(fig1c4_file)) {
  source(here("r-scripts", "r-scripts/02-vcf-investigation.R"), local = TRUE)
}

fig1c1 <- readRDS(file = fig1c1_file) + labs(y = "# of True Positives", fill = "Caller")
fig1c2 <- readRDS(file = fig1c2_file) + labs(y = "Recall", fill = "Caller")
fig1c3 <- readRDS(file = fig1c3_file) + labs(y = " Recall ", fill = "Caller")
fig1c4 <- readRDS(file = fig1c4_file) + labs(fill = "Caller")

# -------------------------
# Fonts / theme
# -------------------------
font_add("CMU Serif", regular = here("r-scripts", "cmunrm.ttf"))
showtext_auto()

theme_set(theme_bw(base_size = 11))
theme_update(text = element_text(family = "CMU Serif"))

# -------------------------
# Dimensions
# -------------------------
width <- 210
height <- 297

# -------------------------
# Layouts
# -------------------------
# B sublayout (within section B)
layout_fig1c <- c(
  area(1, 1, 1, 3), area(1, 4, 1, 6),
  area(2, 1, 2, 6),
  area(3, 1, 3, 6)
)

# Overall layout uses 6 columns.
# We'll add TWO header rows (A header, B header).
height_stretch <- 2

layout_overall <- c(
  # Row 0: Section A header (full width)
  area(1, 1, 1, 6),

  # Rows 1-? : Section A content
  area(2, 1, 1 + 3 * height_stretch, 2), # A.1 (fig1a)
  area(2, 3, 1 + 3 * height_stretch, 6), # A.2 (fig1b)

  # Next row: Section B header (full width)
  area(2 + 3 * height_stretch, 1, 2 + 3 * height_stretch, 6),

  # B.1 / B.2
  area(3 + 3 * height_stretch, 1, 3 + 4 * height_stretch - 1, 3), # B.1
  area(3 + 3 * height_stretch, 4, 3 + 4 * height_stretch - 1, 6), # B.2

  # B.3
  area(3 + 4 * height_stretch, 1, 3 + 5 * height_stretch - 1, 6), # B.3

  # B.4
  area(3 + 5 * height_stretch, 1, 3 + 6 * height_stretch - 1, 6), # B.4

  # Legend row
  area(3 + 6 * height_stretch, 1, 3 + 7 * height_stretch - 1, 6) # guide_area()
)

# -------------------------
# Section headers as grobs
# -------------------------
hdr_A <- wrap_elements(
  grid::textGrob(
    "A. Performance Metrics by Caller",
    x = 0, hjust = 0,
    gp = grid::gpar(fontface = "bold", fontsize = 12, fontfamily = "CMU Serif")
  )
)

hdr_B <- wrap_elements(
  grid::textGrob(
    "B. Performance Metrics Binned by SV Length",
    x = 0, hjust = 0,
    gp = grid::gpar(fontface = "bold", fontsize = 12, fontfamily = "CMU Serif")
  )
)

# -------------------------
# Explicit tags (robust)
# -------------------------
p_A1 <- wrap_elements(plot = fig1a, clip = FALSE) + plot_annotation(tag = "A.1")
p_A2 <- wrap_table(fig1b, panel = "body") + plot_annotation(tag = "A.2")

p_B1 <- fig1c1 + labs(tag = "B.1")
p_B2 <- fig1c2 + labs(tag = "B.2")
p_B3 <- fig1c3 + labs(tag = "B.3")
p_B4 <- fig1c4 + labs(tag = "B.4")

# -------------------------
# Combine
# -------------------------
patchwork_noC <-
  hdr_A +
    p_A1 + p_A2 +
    hdr_B +
    p_B1 + p_B2 +
    p_B3 +
    p_B4 +
    guide_area() +
    plot_layout(guides = "collect", axes = "collect", design = layout_overall) &
    theme(
      plot.tag = element_text(face = "bold"),
      plot.tag.position = c(0, 1),
      legend.position = "bottom",
      legend.byrow = TRUE,
      legend.direction = "horizontal"
    ) &
    guides(fill = guide_legend(nrow = 2, ncol = 7, byrow = TRUE))

# Keep your special-case legend behavior for fig1a
patchwork_noC[[2]] <- patchwork_noC[[2]] + theme(legend.position = "right") + plot_layout(guides = "keep")

# -------------------------
# Save
# -------------------------
ggsave(
  patchwork_noC,
  filename = "fig1_overall.pdf",
  device = cairo_pdf,
  width = width, height = height,
  units = "mm", dpi = 300
)

ggsave(
  patchwork_noC,
  filename = "fig1_overall.svg",
  device = svglite,
  width = width, height = height,
  units = "mm", dpi = 300,
  web_fonts = svglite::fonts_as_import("cmu-serif")
)

ggsave(
  patchwork_noC,
  filename = "fig1_overall.png",
  width = width, height = height,
  units = "mm", dpi = 300
)

patchwork_noC

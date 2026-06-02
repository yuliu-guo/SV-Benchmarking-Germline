#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
  library(janitor)
  library(VariantAnnotation)
  library(showtext)
  library(showtext)
  library(systemfonts)
})

here::i_am("r-scripts/04-compare_full.R")

# =============================================================================
# Configuration
# =============================================================================

data_folder <- "type-ignored"
focus_sample <- "NA12878"

# Explicit caller set (no extras)
callers <- c(
  "cnvpytor", "delly", "dysgu", "gridss", "manta",
  "octopus", "popdel", "svaba", "tardis", "tiddit"
)

linear_summary_csv <- here::here(
  "pedigree-results", "truvari", data_folder, "summary_statistics.csv"
)
graph_summary_csv <- here::here(
  "graph-results", "truvari", data_folder, "summary_statistics.csv"
)

truth_vcf_gz <- here::here("merged_hg38.svs.sort.oa.vcf.gz")

color_file <- here::here("pedigree-results", "truvari", "color_pal.Rds")
color_pal <- readRDS(color_file)

# Ordering must come from color_pal (restricted to explicit callers)
caller_levels <- intersect(names(color_pal), callers)

# 1) Name agreement checks
cat("Explicit callers:\n")
print(callers)
cat("color_pal names (raw):\n")
print(names(color_pal))

cat("Dropped by intersect (explicit not in color_pal):\n")
print(setdiff(callers, names(color_pal)))

cat("In color_pal but not in explicit callers:\n")
print(setdiff(names(color_pal), callers))


# 4) If a caller is missing: show what vcf_read found for it
# (inside your existing vcf_read it prints vcf paths; so just re-run for one missing caller)



callers <- caller_levels



save_folder <- here::here("graph-linear")
dir.create(save_folder, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Tunables (journal half-page figure)
# =============================================================================

FIG_MM_W <- 180
FIG_MM_H <- 120
EXPORT_DPI <- 600
EXPORT_SCALE <- 0.80


FONT_FAMILY <- "CMU Serif"

# Slightly larger global sizing
THEME_BASE_SIZE <- 12 # was 11
FINAL_TEXT_SIZE <- 9 # was 8
FINAL_TITLE_SIZE <- 9 # was 9
FINAL_AXIS_TITLE <- 9 # was 8
FINAL_AXIS_TEXT <- 7 # was 7
FINAL_LEGEND_TEXT <- 8 # was 7
FINAL_LEGEND_TITLE <- 9 # was 8

# Strokes (linewidth)
LW_GRID_MAJOR <- 0.18
LW_BORDER <- 0.25
LW_TICKS <- 0.25

LW_BAR_EDGE <- 0.25
LW_ERR <- 0.25
ERR_WIDTH <- 0.18

LW_UNION_EDGE <- 0.18
POINT_SIZE <- 2.2
POINT_STROKE <- 0.25

BOX_LW <- 0.25
JIT_SIZE <- 0.05
JIT_ALPHA <- 0.25

# Legend sizing (mm)
LEGEND_KEY_H_MM <- 2.4
LEGEND_KEY_W_MM <- 3.2
LEGEND_BOX_SP_MM <- 1.5
PLOT_MARGIN_MM <- 1.5

# Colors
COL_LINEAR <- "grey"
COL_GRAPH <- "#E69F00"
COL_BOTH <- "red"

# =============================================================================
# Fonts + global theme defaults (interactive)
# =============================================================================


FONT_FAMILY <- "CMU Serif"
font_add(FONT_FAMILY, regular = here::here("r-scripts", "cmunrm.ttf"))
showtext_opts(dpi = EXPORT_DPI)
showtext_auto()

theme_set(theme_bw(base_size = THEME_BASE_SIZE, base_family = FONT_FAMILY))
theme_update(text = element_text(family = FONT_FAMILY))

options(repr.plot.width = 17, repr.plot.height = 6, scipen = 999)

theme_journal_halfpage <- theme(
  text = element_text(size = FINAL_TEXT_SIZE),
  plot.title = element_text(size = FINAL_TITLE_SIZE, face = "bold", hjust = 0),
  axis.title = element_text(size = FINAL_AXIS_TITLE),
  axis.text = element_text(size = FINAL_AXIS_TEXT),
  legend.title = element_text(size = FINAL_LEGEND_TITLE),
  legend.text = element_text(size = FINAL_LEGEND_TEXT),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(linewidth = LW_GRID_MAJOR),
  panel.border = element_rect(linewidth = LW_BORDER),
  axis.ticks = element_line(linewidth = LW_TICKS),
  plot.margin = margin(PLOT_MARGIN_MM, PLOT_MARGIN_MM, PLOT_MARGIN_MM, PLOT_MARGIN_MM, "mm"),
  legend.key.height = unit(LEGEND_KEY_H_MM, "mm"),
  legend.key.width = unit(LEGEND_KEY_W_MM, "mm"),
  legend.box.spacing = unit(LEGEND_BOX_SP_MM, "mm"),
  legend.position = "right",
  legend.justification = c(0.5, 0.5),
  legend.box.just = "center"
)

BASE_THEME <- theme_bw(base_size = THEME_BASE_SIZE, base_family = FONT_FAMILY)

theme_set(BASE_THEME)

# =============================================================================
# Utilities
# =============================================================================

read_summary_long <- function(path) {
  read.csv(path) %>%
    pivot_longer(
      cols = matches("^(f1|recall|precision)_"),
      names_to = c(".value", "sample"),
      names_pattern = "^(f1|recall|precision)_(.*)$"
    )
}

summarise_mean_se <- function(df, value_col) {
  df %>%
    group_by(caller, Alignment) %>%
    summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_col]])),
      se = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(0, mean - se),
      ymax = mean + se
    )
}

make_key <- function(df) {
  df %>% mutate(key = paste(seqnames, start, end, REF, ALT, sep = "|"))
}

read_truth_sample <- function(truth_vcf_gz, sample_id) {
  vcf <- readVcf(truth_vcf_gz, genome = "hg38", param = ScanVcfParam(samples = sample_id))

  rr <- as.data.frame(rowRanges(vcf), row.names = NULL) %>%
    transmute(
      seqnames = as.character(seqnames),
      start = start,
      end = end,
      REF = as.character(REF),
      ALT = as.character(unlist(ALT))
    )

  gt <- tryCatch(geno(vcf)$GT, error = function(e) NULL)
  if (!is.null(gt) && ncol(gt) >= 1) {
    g <- gt[, 1]
    keep <- !(is.na(g) | g %in% c("0/0", "0|0", "./.", " ./.", " .|.", "0"))
    rr <- rr[keep, , drop = FALSE]
  }

  rr %>%
    make_key() %>%
    distinct(key)
}

canon_vec <- function(x) {
  if (is.list(x)) {
    vapply(x, function(z) paste(as.character(z), collapse = ","), character(1))
  } else {
    as.character(x)
  }
}


# =============================================================================
# Section 1: F1 / Recall / Precision bars (summary_statistics.csv)
# =============================================================================

pedigree <- read_summary_long(linear_summary_csv)
graph <- read_summary_long(graph_summary_csv)

stats <- bind_rows(list(Linear = pedigree, Graph = graph), .id = "Alignment") %>%
  mutate(caller = str_to_lower(as.character(caller))) %>%
  filter(caller %in% callers) %>%
  mutate(
    caller    = factor(caller, levels = caller_levels),
    Alignment = factor(Alignment, levels = c("Graph", "Linear"))
  )

pd <- position_dodge(width = 0.9)

bar_layer_bundle <- list(
  geom_col(
    position = pd,
    color = "black", linewidth = LW_BAR_EDGE, alpha = 0.9
  ),
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    position = pd,
    width = ERR_WIDTH, linewidth = LW_ERR
  ),
  scale_fill_manual(values = c(Graph = COL_GRAPH, Linear = COL_LINEAR), guide = "none")
)

f1_df <- summarise_mean_se(stats, "f1")
recall_df <- summarise_mean_se(stats, "recall")
precision_df <- summarise_mean_se(stats, "precision")

COMBINE_stats_bar <- (
  ggplot(f1_df, aes(x = caller, y = mean, fill = Alignment)) +
    bar_layer_bundle +
    labs(y = "F1", x = NULL) +
    scale_x_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "A. Performance by Caller")
) / (
  ggplot(recall_df, aes(x = caller, y = mean, fill = Alignment)) +
    bar_layer_bundle +
    labs(y = "Recall", x = NULL) +
    scale_x_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
) / (
  ggplot(precision_df, aes(x = caller, y = mean, fill = Alignment)) +
    bar_layer_bundle +
    labs(y = "Precision", x = "Caller") +
    scale_x_discrete(drop = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
) + plot_layout(axes = "collect")

COMBINE_stats_bar <- COMBINE_stats_bar & theme(legend.position = "none")
COMBINE_stats_bar
ggsave(
  plot = COMBINE_stats_bar,
  filename = here::here(save_folder, "barplot_compare_both.png"),
  width = 8, height = 16, units = "in", dpi = 150
)

# =============================================================================
# Section 2: Shared truth hits (tp-base) between Graph vs Linear by caller
# =============================================================================
right_axis_title_margin_pt <- 1 # tune this (pt)

source(here::here("r-scripts", "00-vcf_process.R"))

linear_df <- vcf_process_all_raw(
  callers = callers,
  truvari_folder = here::here("pedigree-results", "truvari"),
  data_folder = data_folder,
  samples = focus_sample,
  what = "tp-base",
  restrict = FALSE,
  debug = TRUE
) %>%
  dplyr::select(!(TruScore:Multi)) %>%
  mutate(caller = str_to_lower(as.character(caller))) %>%
  make_key() %>%
  distinct(caller, key) %>%
  mutate(linear = TRUE)

graph_df <- vcf_process_all_raw(
  callers = callers,
  truvari_folder = here::here("graph-results", "truvari"),
  data_folder = data_folder,
  samples = focus_sample,
  what = "tp-base",
  restrict = FALSE,
  debug = TRUE
) %>%
  dplyr::select(!(TruScore:Multi)) %>%
  mutate(caller = str_to_lower(as.character(caller))) %>%
  make_key() %>%
  distinct(caller, key) %>%
  mutate(graph = TRUE)


# 2) Did we actually read tp-base rows per caller?
cat("\nRows per caller in tp-base:\n")
print(linear_df %>% dplyr::count(caller, name = "n_linear"))
print(graph_df %>% dplyr::count(caller, name = "n_graph"))


truth_keys <- read_truth_sample(truth_vcf_gz, focus_sample)

truth_by_caller <- tidyr::crossing(
  caller = callers,
  key = truth_keys$key
) %>%
  mutate(caller = str_to_lower(as.character(caller))) %>%
  left_join(linear_df, by = c("caller", "key")) %>%
  left_join(graph_df, by = c("caller", "key")) %>%
  mutate(
    linear = coalesce(linear, FALSE),
    graph  = coalesce(graph, FALSE),
    caller = factor(caller, levels = caller_levels)
  )

hit_df <- truth_by_caller %>%
  mutate(any_hit = graph | linear) %>%
  filter(any_hit) %>%
  mutate(
    class = case_when(
      graph & linear ~ "Both",
      graph & !linear ~ "Graph",
      linear & !graph ~ "Linear"
    ),
    class = factor(class, levels = c("Both", "Graph", "Linear"))
  )

summary_hit <- hit_df %>%
  dplyr::count(caller, class, name = "n") %>%
  group_by(caller) %>%
  mutate(
    n_hit = sum(n),
    prop  = n / n_hit
  ) %>%
  ungroup() %>%
  mutate(caller = factor(caller, levels = caller_levels))

order_df <- hit_df %>%
  group_by(caller) %>%
  summarise(
    n_both = sum(graph & linear),
    n_hit = n(),
    jaccard = n_both / n_hit,
    .groups = "drop"
  ) %>%
  mutate(caller = factor(caller, levels = caller_levels))

COMBINE_union_plot <- ggplot(summary_hit, aes(x = caller, y = prop, fill = class)) +
  geom_col(color = "black", linewidth = LW_UNION_EDGE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  geom_point(
    data = order_df,
    aes(x = caller, y = jaccard),
    inherit.aes = FALSE,
    shape = 21, size = POINT_SIZE, stroke = POINT_STROKE,
    fill = "white", color = "black"
  ) +
  labs(
    x = NULL,
    y = "Shared among TP variants  \n(graph U linear)",
    fill = NULL,
    title = "B. Shared TP Variants between Alignments"
  ) +
  scale_fill_manual(values = c("Both" = COL_BOTH, "Graph" = COL_GRAPH, "Linear" = COL_LINEAR)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = right_axis_title_margin_pt)),
  )

# =============================================================================
# Section 3: SVLEN comparison (tp-comp + fp) — Graph vs Linear
# =============================================================================

what_calls <- c("tp-comp", "fp")



add_svlen <- function(df,
                      overwrite = FALSE,
                      pos_candidates = c("start", "POS", "pos"),
                      end_candidates = c("END", "END_1", "END_2", "END_3", "end"),
                      svlen_col = "SVLEN") {
  stopifnot(is.data.frame(df))

  to_int <- function(x) suppressWarnings(as.integer(as.character(x)))

  pos_cols <- intersect(pos_candidates, names(df))
  end_cols <- intersect(end_candidates, names(df))

  if (length(pos_cols) == 0) {
    stop("add_svlen(): couldn't find a POS column. Tried: ", paste(pos_candidates, collapse = ", "))
  }
  if (length(end_cols) == 0) {
    warning(
      "add_svlen(): couldn't find an END column (or GRanges end). Tried: ",
      paste(end_candidates, collapse = ", ")
    )
    return(df)
  }

  pos <- do.call(dplyr::coalesce, lapply(pos_cols, function(nm) to_int(df[[nm]])))
  end <- do.call(dplyr::coalesce, lapply(end_cols, function(nm) to_int(df[[nm]])))

  svlen_calc <- ifelse(!is.na(end) & !is.na(pos), end - pos, NA_integer_)

  if (!svlen_col %in% names(df)) {
    df[[svlen_col]] <- NA_integer_
  }

  fill_mask <- !is.na(svlen_calc) & (overwrite | is.na(df[[svlen_col]]))
  df[[svlen_col]][fill_mask] <- svlen_calc[fill_mask]

  # Keep your derived fields consistent with the earlier patch, if present (or create them)
  if (!"SVLEN_abs" %in% names(df)) df$SVLEN_abs <- NA_integer_
  df$SVLEN_abs[fill_mask] <- abs(df[[svlen_col]][fill_mask])

  if (!"SVLEN_all" %in% names(df)) df$SVLEN_all <- NA_character_
  df$SVLEN_all[fill_mask] <- as.character(df[[svlen_col]][fill_mask])

  if (!"SVLEN_n" %in% names(df)) df$SVLEN_n <- 0L
  df$SVLEN_n[fill_mask] <- 1L

  # Optional provenance (nice for debugging)
  if (!"SVLEN_source" %in% names(df)) df$SVLEN_source <- NA_character_
  df$SVLEN_source[fill_mask] <- "computed_END_minus_POS"

  attr(df, "add_svlen_filled_n") <- sum(fill_mask)
  df
}



linear_calls <- vcf_process_all_raw(
  callers = callers,
  truvari_folder = here::here("pedigree-results", "truvari"),
  data_folder = data_folder,
  samples = focus_sample,
  what = what_calls,
  restrict = FALSE,
  debug = FALSE
) %>%
  mutate(run = "linear", caller = str_to_lower(as.character(caller)))

linear_calls <- linear_calls %>%
  add_svlen()

graph_calls <- vcf_process_all_raw(
  callers = callers,
  truvari_folder = here::here("graph-results", "truvari"),
  data_folder = data_folder,
  samples = focus_sample,
  what = what_calls,
  restrict = FALSE,
  debug = FALSE
) %>%
  mutate(run = "graph", caller = str_to_lower(as.character(caller)))


graph_calls <- graph_calls %>%
  add_svlen()


svlen_df <- bind_rows(linear_calls, graph_calls) %>%
  filter(is.finite(SVLEN_abs), SVLEN_abs >= 50) %>%
  mutate(
    caller = factor(caller, levels = caller_levels),
    run    = factor(run, levels = c("graph", "linear"))
  )
# 3) Did SVLEN survive filtering?
cat("\nRows per caller in svlen_df:\n")
print(svlen_df %>% dplyr::count(caller, run, name = "n_svlen"))
COMBINE_svlen_plot <- ggplot(svlen_df, aes(x = caller, y = SVLEN_abs, fill = run)) +
  geom_point(
    aes(group = run, color = run),
    position = position_jitterdodge(
      dodge.width = 0.8,
      jitter.width = 0.6,
      jitter.height = 0
    ),
    size = 0.01,
    stroke = NA,
    alpha = JIT_ALPHA
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = BOX_LW) +
  scale_y_log10(labels = comma_format()) +
  labs(
    x = "Caller",
    y = "SVLEN (bp)",
    fill = NULL,
    title = "C. SV Length by Caller"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c(linear = COL_LINEAR, graph = COL_GRAPH), guide = "none") +
  scale_color_manual(values = c(linear = COL_LINEAR, graph = COL_GRAPH), guide = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = right_axis_title_margin_pt)),
  )

svlen_summary <- svlen_df %>%
  group_by(caller, run) %>%
  summarise(
    n_calls = n(),
    min_svlen = min(SVLEN_abs, na.rm = TRUE),
    median_svlen = median(SVLEN_abs, na.rm = TRUE),
    iqr_svlen = IQR(SVLEN_abs, na.rm = TRUE),
    p90_svlen = quantile(SVLEN_abs, 0.90, na.rm = TRUE),
    p95_svlen = quantile(SVLEN_abs, 0.95, na.rm = TRUE),
    max_svlen = max(SVLEN_abs, na.rm = TRUE),
  ) %>%
  arrange(match(caller, factor(caller, levels = caller_levels)), run)

cat("\n=== SVLEN summary by caller and run (tp-comp + fp) ===\n")
print(svlen_summary, n = Inf)
# =============================================================================
# Section 4: Final patchwork + export (half-page)
# =============================================================================

right_col <- (COMBINE_union_plot / COMBINE_svlen_plot) &
  plot_layout(axes = "collect")

FINAL_FIG <- (COMBINE_stats_bar | right_col) +
  plot_layout(widths = c(1.0, 1.25), guides = "collect") &
  theme_journal_halfpage &
  theme(text = element_text(family = FONT_FAMILY))

out_base <- file.path(save_folder, "fig_graph_linear_compare")

ggsave(
  filename = paste0(out_base, ".tiff"),
  plot = FINAL_FIG,
  width = FIG_MM_W, height = FIG_MM_H, units = "mm",
  dpi = EXPORT_DPI, compression = "lzw", bg = "white",
  scale = EXPORT_SCALE,
  limitsize = FALSE
)

ggsave(
  filename = paste0(out_base, ".png"),
  plot = FINAL_FIG,
  width = FIG_MM_W, height = FIG_MM_H, units = "mm",
  dpi = EXPORT_DPI, bg = "white",
  scale = EXPORT_SCALE,
  limitsize = FALSE
)

ggsave(
  filename = paste0(out_base, ".svg"),
  plot = FINAL_FIG,
  width = FIG_MM_W, height = FIG_MM_H, units = "mm",
  dpi = EXPORT_DPI, bg = "white",
  scale = EXPORT_SCALE,
  limitsize = FALSE
)

tryCatch(
  {
    ggsave(
      filename = paste0(out_base, ".pdf"),
      plot = FINAL_FIG,
      width = FIG_MM_W, height = FIG_MM_H, units = "mm",
      device = cairo_pdf, bg = "white",
      scale = EXPORT_SCALE,
      limitsize = FALSE
    )
  },
  error = function(e) {
    message("cairo_pdf failed: ", conditionMessage(e))
    ggsave(
      filename = paste0(out_base, ".pdf"),
      plot = FINAL_FIG,
      width = FIG_MM_W, height = FIG_MM_H, units = "mm",
      bg = "white",
      scale = EXPORT_SCALE,
      limitsize = FALSE
    )
  }
)


print(FINAL_FIG)

# =============================================================================
# Section 5: Console summaries
# =============================================================================

shared_stats_by_caller <- hit_df %>%
  group_by(caller) %>%
  summarise(
    n_union = n(),
    n_both = sum(graph & linear),
    n_graph_only = sum(graph & !linear),
    n_linear_only = sum(linear & !graph),
    jaccard = n_both / n_union,
    shared_of_graph = n_both / (n_both + n_graph_only),
    shared_of_linear = n_both / (n_both + n_linear_only),
    imbalance = (n_graph_only - n_linear_only) / n_union,
    .groups = "drop"
  ) %>%
  arrange(match(caller, factor(caller, levels = caller_levels)))

cat("\n=== Shared variants among union-hits (truth variants hit by >=1) ===\n")
print(shared_stats_by_caller, n = Inf)

results_truth <- truth_by_caller %>%
  group_by(caller) %>%
  summarise(
    n_truth = n(),
    n_both = sum(graph & linear),
    n_graph_only = sum(graph & !linear),
    n_linear_only = sum(linear & !graph),
    n_none = sum(!graph & !linear),
    recall_graph = (n_both + n_graph_only) / n_truth,
    recall_linear = (n_both + n_linear_only) / n_truth,
    phi = {
      x <- as.integer(graph)
      y <- as.integer(linear)
      if (length(unique(x)) < 2 || length(unique(y)) < 2) NA_real_ else cor(x, y)
    },
    .groups = "drop"
  ) %>%
  mutate(delta_recall = recall_graph - recall_linear) %>%
  arrange(match(caller, factor(caller, levels = caller_levels)))

cat("\n=== Truth-universe recall (tp-base on truth variants) ===\n")
print(results_truth, n = Inf)



stats_pretty <- stats %>%
  dplyr::select(-any_of("X")) %>%
  group_by(Alignment, caller) %>%
  summarise(
    across(c(f1, precision, recall), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Alignment,
    values_from = c(f1, precision, recall),
    names_glue = "{.value}_{Alignment}"
  ) %>%
  mutate(
    d_f1        = f1_Graph - f1_Linear,
    d_precision = precision_Graph - precision_Linear,
    d_recall    = recall_Graph - recall_Linear,
    pct_f1      = 100 * d_f1 / f1_Linear
  ) %>%
  arrange(match(caller, factor(caller, levels = caller_levels)))

cat("\n=== Mean metrics (Graph - Linear) from summary_statistics.csv ===\n")
print(stats_pretty, n = Inf)

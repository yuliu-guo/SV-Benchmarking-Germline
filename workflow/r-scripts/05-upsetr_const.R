library(tidyverse) # on server
# need openssl 1.1.1, make sure to load it
library(janitor)
library(eulerr)
library(here)
library(ComplexUpset)
library(jsonlite)
library(patchwork)
here::i_am("r-scripts/05-upsetr_const.R")
callers <- NA
parse_consistency_json <- function(svtype, filepath) {
  cons_json <- read_json(filepath)
  vcf_paths <- unlist(cons_json$vcfs)
  callers <- vapply(strsplit(vcf_paths, "/"), `[`, 4, FUN.VALUE = character(1))
  names(callers) <- vcf_paths

  group_data <- map_dfr(cons_json$detailed, function(grp) {
    bits <- strsplit(grp$group, "")[[1]]
    row <- setNames(as.logical(as.integer(bits)), callers)
    tibble(!!!row, n = grp$total)
  })

  group_data %>%
    uncount(n) %>%
    mutate(
      name = str_c("variant_", row_number()),
      SVTYPE = svtype
    ) %>%
    relocate(name)
}
focus_sample <- "NA12878"

# Set file paths
# pedigree-results/truvari-consistency/type-ignored/NA12878/DEL/consistency.json
data_folder <- here("pedigree-results", "truvari", "type-ignored", focus_sample)
files <- set_names(
  c(paste0(data_folder, "/INS", "_consistency.json"), paste0(data_folder, "/DEL", "_consistency.json"), paste0(data_folder, "/INV", "_consistency.json")),
  c("INS", "DEL", "INV")
)
files
file.exists(unlist(files))

# Combine all
all_data <- map2_dfr(names(files), files, parse_consistency_json)


# ---- Caller columns ----
upset_callers <- setdiff(colnames(all_data), c("name", "SVTYPE"))
callers <- upset_callers
upset_input <- all_data

# ---- Sanity check: Compare with JSON counts ----
validate_totals <- function(file, svtype) {
  total_data <- sum(upset_input$SVTYPE == svtype)
  cat(sprintf("%s: Parsed = %d\n", svtype, total_data))
}

walk2(files, names(files), validate_totals)
# Should print:
# INS: JSON =  XXXX | Parsed = XXXX
# DEL: JSON =  XXXX | Parsed = XXXX
# INV: JSON =  XXXX | Parsed = XXX

sv_colors <- c(
  "DEL" = "#D55E00", # Muted red
  "INS" = "#0072B2", # Muted blue
  "INV" = "#CC79A7" # Bright magenta
)
print("upset_input")
head(upset_input)

summary_file <- here("pedigree-results", "truvari", "type-ignored", "data.Rds")
color_file <- here("pedigree-results", "truvari", "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  print("Files don't exist, running 01 script....")
  source(here("r-scripts", "01-truvari-summary.R"))
}

# font_add("CMU Serif", regular = ttf) # CM Roman equivalent
# showtext_auto()


theme_set(theme_bw(base_size = 11))
# theme_update(text = element_text(family = "CMU Serif"))
# load the f1 summary data
summary_data <- readRDS(file = summary_file)
print(head(summary_data))
# get the color and orders of the callers
color_pal <- readRDS(color_file)

upset_color_pal <- as.character(color_pal)
names(upset_color_pal) <- names(color_pal)

# hack
color_metadata <- data.frame(
  set = callers,
  callers = callers
)

upset_callers <- callers[callers %in% names(upset_input)]
print(upset_callers)

# only keep colors for callers that are in your input
real_callers <- intersect(names(color_pal), upset_callers)
cp <- color_pal[real_callers]

my_queries <- map2(
  names(cp), cp,
  ~ upset_query(set = .x, fill = .y)
)

## SIMPLE INTERSECTION PLOTS

data_folder <- here("pedigree-results", "truvari", "type-ignored")


upset_input$SVTYPE <- factor(upset_input$SVTYPE, levels = c("INV", "INS", "DEL"))
sv_colors <- rev(sv_colors)

freqplot <- function(data, plot_title = "upset-freqsort",
                     min_size = 30, show_totals = FALSE,
                     inter_sort = c("cardinality", "degree"),
                     do_queries = TRUE, show_legend = TRUE,
                     top_title = "Common True Positive Variants Based on Truth Set",
                     left_title = "True Positives per Caller",
                     bottom_title = "Intersections") {
  valid_data <- data %>% dplyr::select(-where(~ is.logical(.x) && sum(.x) < min_size))
  valid_callers <- upset_callers[upset_callers %in% names(valid_data)]
  valid_color_pal <- color_pal[names(color_pal) %in% valid_callers]
  valid_callers <- names(valid_color_pal)
  queries <- NULL
  if (do_queries) {
    queries <- Map(
      function(set, col) upset_query(set = set, fill = col, only_components = c("intersections_matrix")),
      names(valid_color_pal), # the set names
      unname(valid_color_pal) # the colours (strip names)
    )
  }
  tp_labels <- NULL
  if (show_totals) {
    tp_labels <- list(
      annotate(
        geom = "text", x = Inf, y = Inf,
        label = paste("Total True Positives:", nrow(data)),
        vjust = 1.1, hjust = 1.1,
        size = 5, family = "CMU Modern"
      ),
      annotate(
        geom = "text", x = Inf, y = Inf, vjust = 1.5, hjust = 1.3,
        label = paste(
          "DEL ", nrow(data %>% filter(SVTYPE == "DEL")),
          "\nINS ", nrow(data %>% filter(SVTYPE == "INS")),
          "\nINV ", nrow(data %>% filter(SVTYPE == "INV")), " "
        ),
        size = 4, family = "CMU Modern"
      )
    )
  }
  guide_value <- "none"
  if (show_legend) {
    guide_value <- "legend"
  }

  freqsort <- upset(data, rev(valid_callers),
    name = bottom_title,
    width_ratio = 0.1, min_size = min_size,

    # build the upper (bar) plot
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = SVTYPE),
        text = list(family = "CMU Serif", size = 3, color = "black"),
        bar_number_threshold = 1,
      ) +
        ylab(top_title) +
        scale_fill_manual(values = sv_colors, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        tp_labels +
        theme(axis.title = element_text(family = "CMU Serif", size = 14))
    ),
    queries = queries,

    # build the left (count) plot
    set_sizes = (
      upset_set_size(
        geom = geom_bar(
          aes(fill = SVTYPE, x = group),
          width = 0.8
        ),
        position = "left"
      ) + ylab(left_title) +
        labs(fill = "SV Type") +
        scale_y_reverse(expand = expansion(mult = c(0, 0))) +
        scale_fill_manual(values = sv_colors, guide = guide_value)),

    # color the stripes
    stripes = adjustcolor(unname(rev(valid_color_pal)), alpha.f = 0.2),

    # color the matrix!
    matrix = intersection_matrix(
      geom = geom_point(
        aes(alpha = value),
        shape = "circle filled", size = 3
      )
    ) + xlab(bottom_title) +
      scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1), guide = "none"),

    # remove extra grid lines and modify some spacing
    themes = upset_modify_themes(
      list(
        "intersections_matrix" = theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(0, "in"),
          text = element_text(family = "CMU Serif", size = 20)
        ),
        "overall_sizes" = theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(x = 0, units = "in"),
        )
      )
    ),
    guides = "collect",
    sort_sets = FALSE, sort_intersections_by = inter_sort,
  )

  # row 1
  freqsort[[1]] <- freqsort[[1]] & theme(plot.margin = margin(b = 0, r = -10, l = 5))
  freqsort[[2]] <- freqsort[[2]] & theme(plot.margin = margin(b = 0, l = -10))
  # row 2
  freqsort[[3]] <- freqsort[[3]] & theme(plot.margin = margin(t = 0, r = -10, b = 0, l = 5))
  freqsort[[4]] <- freqsort[[4]] & theme(plot.margin = margin(t = 0, r = 0, b = 0, l = -10))

  freqsort
  ggsave(plot = freqsort, width = 16, height = 9, here(data_folder, paste0(plot_title, ".pdf")))
  ggsave(plot = freqsort, width = 16, height = 9, here(data_folder, paste0(plot_title, ".png")))
  return(freqsort)
}
freqsort <- freqplot(upset_input, show_totals = TRUE, min_size = 25)
freqsort
upset_inv <- upset_input %>% filter(SVTYPE == "INV")
inversions_upset <- freqplot(upset_inv,
  plot_title = "inv_only",
  min_size = 1, show_legend = FALSE,
  show_totals = FALSE, top_title = "", left_title = "", bottom_title = ""
)
inversions_upset

all_data <- upset_input %>%
  filter(SVTYPE == "INV") %>%
  pivot_longer(cols = upset_callers, names_to = "caller", names_repair = "minimal") %>%
  filter(value) %>%
  mutate(caller = factor(caller, names(color_pal)))
inv_counts <- all_data %>%
  group_by(caller) %>%
  summarize(n = n()) %>%
  mutate(outline = case_when(
    caller == "svaba" ~ "TRUE",
    caller == "tardis" ~ "TRUE",
    TRUE ~ "FALSE"
  ))

inversion_plot <- ggplot(inv_counts, aes(x = caller, y = n, color = outline)) +
  geom_col(fill = sv_colors["INV"]) +
  scale_color_manual(values = c("TRUE" = "#000", "FALSE" = "#fff"), guide = "none") +
  geom_text(aes(label = n), vjust = 1.5) +
  labs(y = "Number of INVs Called", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(family = "CMU Serif", angle = 45, vjust = 1, hjust = 1))
inversion_plot
#ggsave(width = 16, height = 9, here(data_folder, paste0("inversion_plot", ".pdf")))
ggsave(width = 16, height = 9, here(data_folder, paste0("inversion_plot", ".png")))

# (freqsort | inversions_upset) + plot_layout(guides = "collect")

wrap_elements(freqsort) + inset_element(
  ggplotify::as.ggplot(inversions_upset) + theme(
    plot.background =
      element_rect(fill = "white", color = "black", linewidth = 2)
  ),
  left = 0.5, bottom = 0.46, top = 0.905, right = 0.9, align_to = "full"
)
ggsave(width = 16, height = 9, here(data_folder, paste0("inset_inv", ".pdf")))
ggsave(width = 16, height = 9, here(data_folder, paste0("inset_inv", ".png")))


wrap_elements(freqsort) + inset_element(
  ggplotify::as.ggplot(inversion_plot) + theme(
    plot.background =
      element_rect(fill = "white", color = "black"),
    plot.margin = margin(b = -2, t = 1)
  ),
  left = 0.5, bottom = 0.5,
  top = 0.85, right = 0.78,
  align_to = "full"
)
ggsave(width = 16, height = 9, here(data_folder, paste0("inset_inv_smol", ".pdf")))
ggsave(width = 16, height = 9, here(data_folder, paste0("inset_inv_smol", ".png")))


# print out values
count_unique <- function(col_name, df = NULL, svtype = NULL) {
  if (is.null(df)) {
    df <- upset_input
  }
  df %>%
    filter(
      # only keep rows where this column is TRUE
      .data[[col_name]],
      # only one TRUE among all boolean columns
      rowSums(across(where(is.logical))) == 1,
      # optionally filter SVTYPE
      if (!is.null(svtype)) SVTYPE == svtype else TRUE
    ) %>%
    nrow()
}

print(count_unique("dragen"))
print(count_unique("dysgu"))
print(count_unique("octopus"))
print(count_unique("svaba"))
print(count_unique("tardis"))

print(count_unique("dragen") / (unique(upset_input$name) %>% length()))

print(count_unique("dragen", svtype = "INS"))
print(count_unique("dragen", svtype = "INS") / count_unique("dragen"))

print(count_unique("dysgu", svtype = "INS"))
print(count_unique("dysgu", svtype = "DEL"))
print(count_unique("dysgu", svtype = "INV"))

print(count_unique("tardis", svtype = "INS"))
print(count_unique("tardis", svtype = "DEL"))
print(count_unique("tardis", svtype = "INV"))

print(count_unique("delly"))
print(count_unique("delly"), svtype = "INS")

print(count_unique("gridss"), svtype = "INS")


print(count_unique("lumpy"))
print(count_unique("lumpy"), svtype = "INS")

print(count_unique("smoove"))
print(count_unique("wham"))
print(count_unique("wham"), svtype = "INS")

count_one <- function(col_name, df = NULL, svtype = NULL) {
  if (is.null(df)) {
    df <- upset_input
  }
  df %>%
    filter(
      # only keep rows where this column is TRUE
      .data[[col_name]],
      # optionally filter SVTYPE
      if (!is.null(svtype)) SVTYPE == svtype else TRUE
    ) %>%
    nrow()
}

count_double <- function(col_name_1, col_name_2, df = NULL, svtype = NULL) {
  if (is.null(df)) {
    df <- upset_input
  }
  df %>%
    filter(
      # only keep rows where this column is TRUE
      .data[[col_name_1]] & .data[[col_name_2]],
      # optionally filter SVTYPE
      if (!is.null(svtype)) SVTYPE == svtype else TRUE
    ) %>%
    nrow()
}

report_match <- function(var_list, title = NULL) {
  var_1 <- var_list[1]
  var_2 <- var_list[2]
  print(paste(var_1, count_one(var_1)))
  print(paste(var_2, count_one(var_2)))

  print(count_double(var_1, var_2))

  print(count_double(var_1, var_2) / count_one(var_1))
  print(count_double(var_1, var_2) / count_one(var_2))

  mat <- upset_input %>% select(var_list)
  fit <- euler(mat, shape = "ellipse")
  plot(fit,
    quantities = TRUE, type = c("counts", "percent"),
    edges = color_pal[var_list],
    fills = adjustcolor(color_pal[var_list], alpha.f = 0.5),
    main = title
  )
  return(fit)
}
report_match(c("smoove", "lumpy"), title = "Smoove and Lumpy TP Overlaps")
report_match(c("dragen", "manta"), title = "Dragen and Manta TP Overlaps")

report_match(names(color_pal)[1:3], title = "Top 3 Callers")
report_match(names(color_pal)[1:5], title = "Top 5 Callers")

# note - these commands take a long time and output a mess 
report_match(names(color_pal)[1:10], title = "Top 10 Callers")
#report_match(names(color_pal), title = "Abandon all hope ye who enter here")


report_match(c("dragen", "dysgu", "octopus"), title = "Top 3 Callers with the Most 'Unique' Calls")


count_double("dragen", "manta") / count_one("dragen")
count_double("dragen", "manta") / count_one("manta")

mat <- upset_input %>% select(c("dragen", "manta"))
fit <- euler(mat, shape = "ellipse")
plot(fit, quantities = TRUE)

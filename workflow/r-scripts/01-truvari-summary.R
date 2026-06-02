library(jsonlite) # on server
library(tidyverse) # on server
library(ggrepel) # on server
library(forcats) # on server
library(ggforce)
library(VariantAnnotation) 
library(here)

# take in our inputs
library("optparse")
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = "pedigree-results",
    help = "directory with results [default = %default", metavar = "character"
  ),
  make_option(c("-n", "--names"),
    type = "character", default = "NA12878,NA12879,NA12881,NA12882",
    help = "comma separated names [default= %default]", metavar = "character"
  ),
  make_option(c("-f", "--files"), type = "character", default = "summary.json", help = "summary file to use", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
sample_names <- strsplit(opt$names, ",")[[1]]
print(sample_names)
# setup our file structure
here::i_am("workflow/r-scripts/01-truvari-summary.R")
truvari_folder <- here(opt$input, "truvari")


theme_set(theme_bw(base_size = 11))
# theme_update(text=element_text(size=10,family="Arial"))
options(repr.plot.width = 17, repr.plot.height = 6)

# available color palettes
# colors <- palettes_d_names |> filter(novelty == FALSE, length >= 16,  type == "qualitative") |> arrange(type, length)
# for(i in 1:nrow(colors)){
# color <- colors[i,]
# color_name <- str_c(color$package, "::", color$palette)
# print(str_c(color_name, " - ", color$type, " n = ", color$length))
# print(paletteer_d(palette=color_name, n=16,type = ))
# }

color_pal_set <- "ggthemes::Tableau_20" # set to your fav color name

# 1 Summary Statistics
## 1.0 Getting Summary Files
dbl_vars <- c(
  "precision", "recall",
  "f1"
)

data_folders <- c("type-ignored")
print(paste("Data folders are: ", data_folders))
all_runs <- data.frame()

get_name <- function(x) {
  return(str_replace_all(str_replace(x, "_", " with "), "-", " "))
}
# Define our graph defaults so the graphs all look the same
gglayer_bar_theme <- list(
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)),
  xlab("Callers")
) # this will have our colors added to it later
gglayer_sep_bar <- c(geom_col(aes(fill = caller)), facet_grid(. ~ sample))
gglayer_med_bar <- c(
  gglayer_bar_theme,
  stat_summary(aes(fill = caller), fun = "median", geom = "col"),
  geom_point(aes(shape = sample))
)

gglayer_scatter <- list(
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "lightgrey"
  ),
  geom_point(aes(shape = sample)),
  geom_text_repel(
    max.overlaps = 10, # skip overplotted labels (sad)
    size = 3, aes(
      label = if_else( # only label one of the samples
        sample == levels(sample)[2],
        caller,
        ""
      )
    ), show.legend = FALSE
  ), labs(x = "Precison", y = "Recall"),
  theme(legend.text = element_text(size = 10))
)

for (data_folder_name in data_folders) {
  tryCatch({
    folder_name <- get_name(data_folder_name)
    data_folder <- here(truvari_folder, data_folder_name)
    print(paste("DOING ANALYSIS ON TRUVARI RUN:", data_folder, " ", folder_name))
    callers <-
      list.dirs(
        path = data_folder,
        full.names = FALSE,
        recursive = FALSE
      )
    files <-
      list.files(
        path = data_folder,
        recursive = TRUE,
        pattern = paste0("^", opt$files),
        full.names = FALSE
      )
    files <- files[!grepl("phab_bench", files)] # get a list of all the summary files
    files <- files[!grepl("temp", files)] # get a list of all the summary files
    print(files)
    temp <- data.frame(t(sapply(here(data_folder, files), fromJSON))) %>%
      rownames_to_column(var = "file")
    print("Loaded raw input...")
    print(temp)
    print("Processing raw input....")
    data <- temp %>%
      dplyr::select(-"gt_matrix") %>% # get rid of genotyping nonsense
      dplyr::select(-contains("gt")) 
    head(data) 
    data <- data %>% # get rid of genotyping nonsense
        unnest(cols = dplyr::where(is.list)) %>% 
        mutate(across(all_of(dbl_vars), as.double)) %>%
        # mutate(across(!dbl_vars, as.integer)) %>% # convert everything
      mutate(file = str_remove(file, data_folder)) %>%
      mutate(file = str_remove(file, "/summary.json")) %>% # caller names from files
      separate(file, c(NA, "caller", "sample"), sep = "/")
    print(data)
    print("Formatting raw input!")
    data <- data %>%
      mutate(sample = str_remove(sample, "/")) %>%
      mutate(sample = case_when(
        str_count(sample, sample_names[1]) > 0 ~ sample_names[1],
        str_count(sample, sample_names[2]) > 0 ~ sample_names[2],
        str_count(sample, sample_names[3]) > 0 ~ sample_names[3],
        TRUE ~ sample
      )) %>% # map the sample names to shorter versions
      mutate(sample = as_factor(sample)) %>%
      mutate(caller = as_factor(str_to_lower(caller)))
    print(data)
    data <- data %>%
      mutate(f1 = ifelse(is.na(f1), 0, f1)) %>% # map NA f1s to 0
      distinct()
    print(levels(data$sample))


    # Define our levels and get the color palette
    if (data_folder_name == data_folders[1]) {
      f1_order <- levels((data |> mutate(caller = fct_reorder(caller, -f1)))$caller)
      print("Order of callers:")
      print(f1_order)
      pr_order <- levels((data |> mutate(caller = fct_reorder(caller, -precision)))$caller)
      re_order <- levels((data |> mutate(caller = fct_reorder(caller, -recall)))$caller)
      #color_pal <- paletteer_d(color_pal_set, n = length(f1_order)) # original command
      # no paletteer workaround 
      color_pal <- c('#4E79A7FF','#A0CBE8FF','#F28E2BFF','#FFBE7DFF','#59A14FFF','#8CD17DFF','#B6992DFF','#F1CE63FF','#499894FF','#86BCB6FF','#E15759FF','#FF9D9AFF','#79706EFF','#BAB0ACFF','#D37295FF','#FABFD2FF','#B07AA1FF','#D4A6C8FF','#9D7660FF','#D7B5A6FF')
      color_pal <- color_pal[1:length(f1_order)]
      names(color_pal) <- f1_order
      color_string <- c(as.character(color_pal))
      df <- data.frame(callers = names(color_pal), colors = color_string)
      write.csv(df, file = here(truvari_folder, "color_pal.csv"))
      saveRDS(color_pal, file = here(truvari_folder, "color_pal.Rds"))
      gglayer_bar_theme <- c(gglayer_bar_theme, scale_fill_manual(values = color_pal, guide = "none"))
      gglayer_sep_bar <- c(gglayer_sep_bar, gglayer_bar_theme)
      gglayer_med_bar <- c(gglayer_med_bar, gglayer_bar_theme)

      gglayer_scatter <- c(gglayer_scatter, scale_color_manual(values = color_pal))
    }
    data <- data |>
      mutate(caller = fct_relevel(caller, f1_order)) |>
      arrange(caller) # reorder by f1

    print(paste(
      "Saving summary data to ",
      here(data_folder, "summary_statistics.csv")
    ))
    to_print <- data %>%
      dplyr::select(sample, caller, f1, precision, recall) %>%
      pivot_wider(
        names_from = sample,
        values_from = c(f1, precision, recall),
        names_sort = TRUE
      )
    print(to_print)
    write.csv(to_print,
      file = here(data_folder, "summary_statistics.csv")
    )
    library(gt)
    gt(to_print) %>% tab_footnote(
      footnote = md("Truth set restricted only to deletions"),
      locations = cells_body(
        columns = caller,
        rows = caller == "popdel"
      )
    )

    library(gt)
    library(stringr)

    # create formatting for table
    to_print_aug <- to_print %>%
      mutate(
        f1_mean        = rowMeans(dplyr::select(., starts_with("f1_")), na.rm = TRUE),
        precision_mean = rowMeans(dplyr::select(., starts_with("precision_")), na.rm = TRUE),
        recall_mean    = rowMeans(dplyr::select(., starts_with("recall_")), na.rm = TRUE)
      )
    sample_names_aug <- c(sample_names, "mean")

    # make basic table
    gt_sample <- purrr::reduce(
      sample_names_aug,
      .init = gt(to_print_aug),
      .f    = ~ .x %>% tab_spanner(label = .y, columns = ends_with(.y))
    )
    gt_sample <- gt_sample %>%
      cols_label_with(
        fn = ~ gsub(paste0("_(", paste(sample_names_aug, collapse = "|"), ")$"), "", .)
      ) %>%
      cols_label_with(fn = ~ gsub("precision", "P", .)) %>%
      cols_label_with(fn = ~ gsub("recall", "R", .)) %>%
      cols_label_with(fn = ~ gsub("f1", "F1", .)) %>%
      fmt_number(decimals = 2) %>% 
      data_color(columns=f1_NA12879:recall_mean, palette = c("#E7F4E9FF","#C7E5C9FF", "#A5D6A6FF", "#80C684FF", "#66BA6AFF",
                                                                         "#4CAE50FF", "#439F46FF", "#388D3BFF", "#2D7D32FF","#1A5E1FFF")) # ggsci::green_material
    # add better styling
    last_recall_sample <- paste0("recall_", tail(sample_names, 1))
    gt_sample <- gt_sample %>%
      tab_style(
        style = list(
          cell_fill(color = "grey90"),
          cell_text(weight = "bold")
        ),
        locations = list(
          cells_column_labels(columns = ends_with("_mean")),
          cells_column_spanners(spanners = "mean")
        )
      ) %>%
      tab_style(
        style = cell_borders(sides = "r", color = "white", weight = px(2)),
        locations = list(
          cells_body(columns = starts_with("recall")),
          cells_column_labels(columns = starts_with("recall"))
        )
      ) %>%
      tab_style(
        style = cell_borders(sides = "l", color = "white", weight = px(3)),
        locations = list(
          cells_body(columns = "f1_mean"),
          cells_column_labels(columns = "f1_mean"),
          cells_column_spanners(spanners = "mean")
        )
      ) %>%
      opt_table_font(size = 10)

    # add some extra text
    gt_sample %>% tab_footnote(footnote = "Assessed on a deletion-only subset of the truth set", locations = cells_body(columns = caller, rows = caller == "popdel"))
    gt_sample
    saveRDS(gt_sample, file = here(data_folder, "plots_gt_sample.Rds"))
    #gtsave(gt_sample, filename = here(data_folder, "table_by_sample.png"))
    #gtsave(gt_sample, filename = here(data_folder, "table_by_sample.pdf"))

    sample_split <- gt_sample %>% gt_split(col_slice_at = starts_with("recall"))
    saveRDS(data, file = here(data_folder, "data.Rds"))
    print("Starting plotting...")
    print(data)

    print("Making precision recall plot")
    pr_plot <-
      ggplot(
        data,
        aes(
          x = precision,
          y = recall,
          color = caller,
          label = caller
        )
      ) +
      gglayer_scatter +
      labs(title = str_c("Truvari run: ", folder_name))

    pr_plot
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_recall.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_recall.pdf")
    )
    print("Making precision recall square plot")
    pr_plot +
      coord_cartesian(
        xlim = c(0, max(data$precision, data$recall)),
        ylim = c(0, max(data$precision, data$recall)),
      )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_recall_square.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_recall_square.pdf")
    )

    ## 3.0 F-scores
    print("Making F1 bar plots....")
    print(data)
    ggplot(
      data,
      aes(x = caller, y = f1)
    ) +
      gglayer_sep_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "f1_barplot_sep.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "f1_barplot_sep.pdf")
    )

    ggplot(data, aes(x = caller, y = f1)) +
      gglayer_med_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "f1_barplot_median.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "f1_barplot_median.pdf")
    )

    print("Making precision bar plots....")
    ggplot(
      data,
      aes(x = fct_relevel(caller, pr_order), y = precision)
    ) +
      gglayer_sep_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_barplot_sep.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_barplot_sep.pdf")
    )
    ggplot(data, aes(x = fct_relevel(caller, pr_order), y = precision)) +
      gglayer_med_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_barplot_median.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "precision_barplot_median.pdf")
    )

    print("Making recall bar plots....")
    ggplot(
      data,
      aes(x = fct_relevel(caller, re_order), y = recall)
    ) +
      gglayer_sep_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "recall_barplot_sep.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "recall_barplot_sep.pdf")
    )
    ggplot(data, aes(x = fct_relevel(caller, re_order), y = recall)) +
      gglayer_med_bar +
      labs(title = str_c("Truvari run: ", folder_name))
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "recall_barplot_median.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "recall_barplot_median.pdf")
    )

    # make a combined plot
    ggplot(data %>% pivot_longer(cols = c(f1, precision, recall)), aes(x = caller, y = value)) +
      gglayer_sep_bar +
      facet_grid(sample ~ name)
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "combined_barplot_sep.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "combined_barplot_sep.pdf")
    )

    combined_barplot_median <- ggplot(data %>% pivot_longer(cols = c(f1, precision, recall)), aes(x = caller, y = value)) +
      gglayer_med_bar +
      facet_row(~name)
    combined_barplot_median
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "combined_barplot_median.png")
    )
    ggplot2::ggsave(
      width = 17, height = 6,
      file = here(data_folder, "combined_barplot_median.pdf")
    )
    combined_barplot_median <- ggplot(data %>% pivot_longer(cols = c(f1, precision, recall)), aes(x = caller, y = value)) +
      gglayer_med_bar +
      facet_col(~name)
    all_runs <- rbind(all_runs, data |>
      mutate(run = folder_name))
  })
}

# add facets to the plots
gglayer_sep_bar_facet <- c(
  gglayer_sep_bar,
  facet_grid(run ~ sample)
)
gglayer_med_bar_facet <- c(gglayer_med_bar, facet_wrap(run ~ .))
print("DOING COMPARISON ANALYSIS BETWEEN TRUVARI RUNS")
print("Comparison scatter plots...")
pr_plot <-
  ggplot(
    all_runs,
    aes(
      x = precision,
      y = recall,
      color = caller
    )
  ) +
  gglayer_scatter +
  facet_wrap(run ~ .)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = pr_plot, file = here(truvari_folder, "precision_recall.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = pr_plot, file = here(truvari_folder, "precision_recall.pdf")
)
pr_plot <- pr_plot +
  coord_cartesian(
    xlim = c(0, max(all_runs$precision, all_runs$recall)),
    ylim = c(0, max(all_runs$precision, all_runs$recall)),
  )
ggplot2::ggsave(
  width = 17, height = 6,
  plot = pr_plot, file = here(truvari_folder, "precision_recall_square.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = pr_plot, file = here(truvari_folder, "precision_recall_square.pdf")
)

pr_plot$layers <- pr_plot$layers[1:2] # remove text repel
pr_plot + geom_mark_ellipse(aes(fill = caller), alpha = 0.5, expand = unit(3, "mm")) +
  scale_fill_manual(values = color_pal) +
  geom_text_repel(
    max.overlaps = 20,
    box.padding = 1,
    min.segment.length = 1,
    show.legend = FALSE,
    nudge_y = -0.0015, # adjustment so in default, tiddit or octopus label isnt inside delly
    aes(
      label = if_else(
        sample == levels(sample)[1],
        caller,
        ""
      )
    )
  )
ggplot2::ggsave(width = 17, height = 6, file = here(truvari_folder, "precision_recall_square_area.png"))

ggplot2::ggsave(width = 17, height = 6, file = here(truvari_folder, "precision_recall_square_area.pdf"))
## 3.0 F-scores



print("Comparison bar plots....")
f1_barplot_sep <- ggplot(
  all_runs,
  aes(x = caller, y = f1)
) + gglayer_sep_bar_facet
ggplot2::ggsave(
  width = 17, height = 6,
  plot = f1_barplot_sep, file = here(truvari_folder, "f1_barplot_sep.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = f1_barplot_sep, file = here(truvari_folder, "f1_barplot_sep.pdf")
)

f1_barplot_median <- ggplot(all_runs, aes(x = caller, y = f1)) +
  labs(x = "Callers", y = "F1") +
  gglayer_med_bar_facet
f1_barplot_median
ggplot2::ggsave(
  width = 17, height = 6,
  plot = f1_barplot_median, file = here(truvari_folder, "f1_barplot_median.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = f1_barplot_median, file = here(truvari_folder, "f1_barplot_median.pdf")
)

precision_barplot_sep <- ggplot(
  all_runs,
  aes(x = fct_relevel(caller, pr_order), y = precision)
) + gglayer_sep_bar_facet
ggplot2::ggsave(
  width = 17, height = 6,
  plot = precision_barplot_sep, file = here(truvari_folder, "precision_barplot_sep.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = precision_barplot_sep, file = here(truvari_folder, "precision_barplot_sep.pdf")
)

precision_barplot_median <- ggplot(all_runs, aes(x = fct_relevel(caller, pr_order), y = precision)) +
  gglayer_med_bar_facet
precision_barplot_median
ggplot2::ggsave(
  width = 17, height = 6,
  plot = precision_barplot_median, file = here(truvari_folder, "precision_barplot_median.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = precision_barplot_median, file = here(truvari_folder, "precision_barplot_median.pdf")
)

recall_barplot_sep <- ggplot(
  all_runs,
  aes(x = fct_relevel(caller, re_order), y = recall)
) + gglayer_sep_bar_facet
ggplot2::ggsave(
  width = 17, height = 6,
  plot = recall_barplot_sep, file = here(truvari_folder, "recall_barplot_sep.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = recall_barplot_sep, file = here(truvari_folder, "recall_barplot_sep.pdf")
)

recall_barplot_median <- ggplot(all_runs, aes(x = fct_relevel(caller, re_order), y = recall)) + gglayer_med_bar_facet
recall_barplot_median
ggplot2::ggsave(
  width = 17, height = 6,
  plot = recall_barplot_median, file = here(truvari_folder, "recall_barplot_median.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  plot = recall_barplot_median, file = here(truvari_folder, "recall_barplot_median.pdf")
)

## 4.0 Combined Plots
wider_data <- pivot_longer(all_runs, cols = c(f1, precision, recall), names_to = "stat", values_to = "value")
ggplot(wider_data, aes(x = caller, y = value)) +
  facet_grid(stat ~ run, scales = "free_y", switch = "y") +
  stat_summary(aes(fill = caller), fun = "median", geom = "col") +
  geom_point(aes(shape = sample)) +
  gglayer_bar_theme +
  labs(x = "Callers", y = " ")
ggplot2::ggsave(
  width = 17, height = 6,
  file = here(truvari_folder, "stat_summary.png")
)
ggplot2::ggsave(
  width = 17, height = 6,
  file = here(truvari_folder, "stat_summary.pdf")
)

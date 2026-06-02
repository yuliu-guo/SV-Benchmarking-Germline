#!/usr/bin/env Rscript
library(tidyverse) # on server
# need openssl 1.1.1, make sure to load it
library(StructuralVariantAnnotation) # on server
library(janitor)
library(ggrepel)
library(here)
library("optparse")
library(patchwork)

library(dplyr)
library(tidyr)
library(GenomicRanges)
library(Gviz)
library(data.table)



# options(error=traceback)
option_list <- list(
  make_option(c("-i", "--input"),
    type = "character", default = "pedigree-results",
    help = "directory with results [default = %default", metavar = "character"
  ),
  make_option(c("-r", "--rds"),
    type = "character", default = "FALSE",
    help = "if you want to reuse or name [default = %default", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

here::i_am("r-scripts/02-vcf-investigation.R")
truvari_folder <- here(opt$input, "truvari")

# load vcf processing functions
source(here("r-scripts", "00-vcf_process.R"))

args <- commandArgs(trailingOnly = TRUE)
# useRDS <- opt$rds
useRDS <- TRUE
options(repr.plot.width = 16, repr.plot.height = 9)
library(showtext)
font_add("CMU Serif", regular = here("r-scripts", "cmunrm.ttf")) # CM Roman equivalent
showtext_auto()


theme_set(theme_bw(base_size = 11))
theme_update(text = element_text(family = "CMU Serif"))


# Loading Data ------------------------------------------------------------



data_folder <-
  "type-ignored/"
truvari_name <- str_replace(str_remove(data_folder, "/"), "-", " ")
folders <-
  list.dirs(
    path = here(truvari_folder, data_folder),
    full.names = FALSE,
    recursive = TRUE
  )
folders <- folders[!grepl("phab_bench", folders)] # get a list of all the summary files
folders <- folders[!grepl("temp", folders)] # get a list of all the summary files
print(folders)
callers <- unique(str_split_i(str_remove(folders, data_folder), "/", 1))[-1]
samples <- unique(str_split_i(str_remove(folders, data_folder), "/", 2))[-1]
print(samples)
print(callers)
tp_type <- "tp-base"
what <- c(tp_type, "fn")
sample_name <- samples[1]

# quit(save="ask")

summary_file <- here(truvari_folder, data_folder, "data.Rds")
color_file <- here(truvari_folder, "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  source(here("r-scripts", "01-truvari-report.R"))
}
# load the f1 summary data
summary_data <- readRDS(file = summary_file)
# get the color and orders of the callers
callers <- levels(summary_data$caller)
samples <- levels(summary_data$sample)
color_pal <- readRDS(color_file)
print(color_pal)

samples <- c("NA12878")

caller_i <- "manta"
sample_i <- sample_name
overall_data <- data.frame()
rds_location <- here(truvari_folder, data_folder, tp_type, "vcf_data.Rds")
if (useRDS == "TRUE" & file.exists(rds_location)) { # use that hard won vcf data :D
  print("LOADING OLD RDS DATA.....")
  overall_data <- readRDS(file = rds_location)
  print(head(overall_data, 1))
} else {
  for (caller_i in callers) {
    for (sample_i in samples) {
      print(caller_i)
      vcfs <- vcf_read(sample_i, caller_i,
        truvari_folder = truvari_folder,
        data_folder = data_folder,
        include = str_c(what, collapse = "|")
      )
      print(str_c("VCFs loaded for sample ", sample_i, " and caller ", caller_i))
      if (length(vcfs) == length(what)) { # make sure we're loading all 3

        for (i in seq(1, length(vcfs))) {
          print(str_c("Processing ", names(vcfs)[i]))
          data <- vcf_process(vcfs[[i]], status = names(vcfs)[i]) |>
            mutate(sample = sample_i) |>
            mutate(caller = caller_i)

          if (length(overall_data) > 1) {
            data <- data |> dplyr::select(any_of(names(overall_data)))
            print(head(data, 1))
          }
          overall_data <- bind_rows(overall_data, data) %>% remove_empty()
        }
      }
    }
  }

  # save that hard won vcf data :)
  print("Saving RDS data...")
  print(head(overall_data, 1))
  saveRDS(overall_data, file = rds_location)
}



# based on hue_pal()2, then grey
group.colors <- c(
  "True Positive" = "#00BFC4",
  "False Positive" = "#F8766D",
  "False Negative" = "grey"
)
names(overall_data)
overall_data <- dplyr::rename(overall_data, c("type" = "SVTYPE"))


if (length(unique(overall_data$caller) != length(names(color_pal)))) {
  print(str_c("We are missing ", unique(overall_data$caller)[!(unique(overall_data$caller) %in% names(color_pal))], " in the summary data"))
  print(str_c("We are missing ", names(color_pal)[!(names(color_pal) %in% unique(overall_data$caller))], " in the vcfs"))
}

overall_data <- overall_data %>%
  mutate(sample = as_factor(overall_data$sample)) %>%
  mutate(caller = factor(str_to_lower(caller), levels = names(color_pal), )) %>%
  filter(!is.na(caller)) # remove anything not mapped properly
print(levels(overall_data$sample))

# remove extra typing
overall_data <- overall_data |>
  mutate(type = case_when(caller == "tardis" ~ str_split_i(type, ":", 1), type == "DUP:TANDEM" ~ "DUP", TRUE ~ type))

### 0 Sanity Check
print("Running sanity checks for read in data....")
sanity <- overall_data %>%
  group_by(caller, sample, status) |>
  dplyr::count() |>
  pivot_wider(
    names_from = status,
    values_from = n,
    values_fill = 0,
  ) |>
  clean_names() |>
  full_join(
    summary_data |>
      dplyr::select(caller, sample,
        tp_comp = TP.comp,
        fn = FN, fp = FP, f1, precision, recall
      ),
    by = c("caller", "sample"), suffix = c(".vcf", ".json")
  ) |>
  mutate(across(fn.vcf:f1, ~ ifelse(is.na(.x), 0, .x)))

# print("Removing any caller sample set found here:")
# sanity[!(sanity$tp_comp.json == sanity$tp_comp.vcf), ]
# sanity[!(sanity$fn.json == sanity$fn.vcf), ]
# sanity[!(sanity$fp.json == sanity$fp.vcf), ]


overall_data <- overall_data |>
  #  anti_join(sanity[!(sanity$fn.json == sanity$fn.vcf), ], by = c("caller", "sample")) |>
  #  anti_join(sanity[!(sanity$fp.json == sanity$fp.vcf), ], by = c("caller", "sample")) |>
  #  anti_join(sanity[!(sanity$tp_comp.json == sanity$tp_comp.vcf), ], by = c("caller", "sample")) |>
  mutate(status = case_when(status == "fn" ~ "False Negative", status == "fp" ~ "False Positive", status == tp_type ~ "True Positive"))

head(overall_data)


# Type Counts -------------------------------------------------------------


type_counts <- overall_data %>%
  group_by(sample, caller, status, type) %>%
  dplyr::count() |>
  mutate(caller = factor(caller, levels = unique(overall_data$caller))) |>
  mutate(n_label = case_when(n > 1000 ~ NA, TRUE ~ n))

print(levels(type_counts$caller))

sample_names <- unique(overall_data$sample)
for (sample_name in sample_names) {
  type_count_graph <- ggplot(
    type_counts %>% filter(sample == sample_name),
    aes(x = reorder(type, -n), y = n, fill = status)
  ) +
    geom_col() +
    scale_fill_manual(values = group.colors) +
    scale_color_manual(values = group.colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title = str_c("Types called and missed, with ", truvari_name, " on sample ", sample_name), x = "Type", y = "Counts")


  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_raw.png"))

  type_count_graph + facet_wrap(vars(caller))
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts.png"))

  type_count_graph + facet_wrap(vars(caller), scales = "free_x") +
    geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freex_label.png"))


  type_count_graph + facet_wrap(vars(caller), scales = "free_x")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freexy.png"))

  type_count_graph + facet_wrap(vars(caller), scales = "free_x") + geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_freexy_label.png"))
  print(sample_names)
  type_count_graph <- ggplot(
    type_counts %>% filter(sample == sample_name) %>%
      filter(status != "False Negative"),
    aes(x = reorder(type, -n), y = n, fill = status)
  ) +
    geom_col() +
    scale_fill_manual(values = group.colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_wrap(vars(caller), scales = "free_x") +
    labs(title = str_c("Types called with ", truvari_name, " on sample ", sample_name), x = "Type", y = "Counts")
  type_count_graph
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts.png"))

  type_count_graph +
    scale_color_manual(values = group.colors) +
    geom_label_repel(max.overlaps = 50, nudge_y = 100, aes(label = n_label, color = status), fill = "white", direction = "y")
  ggplot2::ggsave(create.dir = TRUE, width = 16, height = 9, here(truvari_folder, data_folder, tp_type, sample_name, "type_counts_label.png"))
}
overall_data <- overall_data %>% dplyr::select(rowname, seqnames, start, end, width, FILTER, SVLEN, SVLEN_1, type, sample, caller, status)



## SANITY CHECK: BASIC VALUES MATCH?
calculated_summary <- overall_data %>%
  group_by(caller, status) %>%
  dplyr::count() %>%
  pivot_wider(names_from = status, values_from = n) %>%
  mutate(recall = `True Positive` / (`True Positive` + `False Negative`))
calculated_summary
summary_data %>%
  dplyr::filter(sample == sample_name) %>%
  dplyr::select(caller, FN, TP.base, recall)

## NOW LETS DO IT
quantile(abs(overall_data$SVLEN))
# bin the SVLEN
overall_data <- overall_data %>% mutate(
  SVLEN =
    factor(case_when(
      abs(SVLEN) < 100 ~ "50-100",
      abs(SVLEN) < 500 ~ "100-500",
      abs(SVLEN) < 1000 ~ "500-1000",
      TRUE ~ "1000+"
    ), levels = c("50-100", "100-500", "500-1000", "1000+"))
)

save_folder <- here(truvari_folder, data_folder)

binned_data <- overall_data %>%
  group_by(SVLEN, caller, status) %>%
  dplyr::count() %>%
  pivot_wider(names_from = status, values_from = n) %>%
  mutate(recall = `True Positive` / (`True Positive` + `False Negative`))


f1_all <- ggplot(binned_data, aes(x = SVLEN, y =  `True Positive`, fill = caller)) +
  geom_col(position = "dodge") +
  labs(x = "SVLEN", y = "# of True Positives") +
  scale_fill_manual(values = color_pal)
saveRDS(f1_all, file = here(save_folder, "fig1_f1_all.Rds"))
recall_all <- ggplot(binned_data, aes(x = SVLEN, y = recall, fill = caller)) +
  geom_col(position = "dodge") +
  labs(x = "SVLEN", y = "Recall") +
  scale_fill_manual(values = color_pal)
saveRDS(recall_all, file = here(save_folder, "fig1_recall_all.Rds"))

binned_data_type <- overall_data %>%
  group_by(SVLEN, caller, type, status) %>%
  dplyr::count() %>%
  pivot_wider(names_from = status, values_from = n) %>%
  mutate(recall = `True Positive` / (`True Positive` + `False Negative`)) %>%
  mutate(SVTYPE = factor(type, levels = c("DEL", "INS", "INV"))) %>%
  filter(!is.na(SVTYPE))

recall_split <- ggplot(binned_data_type, aes(x = SVLEN, y = recall, fill = caller)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color_pal) +
  facet_wrap(~SVTYPE, scales = "fixed")
recall_split
saveRDS(recall_split, file = here(save_folder, "fig1_recall.Rds"))

library(ggh4x)
count_split_TP <- ggplot(binned_data_type, aes(x = SVLEN, y = `True Positive`, fill = caller)) +
  geom_col(position = "dodge") +
  labs(y = "# of True Positives") +
  scale_fill_manual(values = color_pal) +
  facet_wrap(~SVTYPE, scales = "free_y", axes = "all_y") +
  scale_y_continuous(limits = c(0, 2035)) +
  scale_y_facet(PANEL == 2, labels = NULL, limits = c(0, 2035)) +
  scale_y_facet(SVTYPE == "INV", limits = c(0, 20.35))
saveRDS(count_split_TP, file = here(save_folder, "fig1_count_split_TP.Rds"))


library(patchwork)
recall_split / count_split_TP + plot_layout(guides = "collect", axes = "collect")

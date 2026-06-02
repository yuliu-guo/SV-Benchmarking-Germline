library(tidyverse) # on server
# need openssl 1.1.1, make sure to load it
library(StructuralVariantAnnotation) # on server
library(janitor)
library(eulerr)
library(here)
library(ComplexUpset)


here::i_am("r-scripts/upsetr.R")
source(here("r-scripts", "00-vcf_process.R"))
options(repr.plot.width = 16, repr.plot.height = 9)

theme_set(theme_bw())

truvari_run <- "type-ignored"
truvari_folder <- here("pedigree-results", "truvari")
data_folder <- here("pedigree-results", "truvari", truvari_run)
folders <-
  list.dirs(
    path = data_folder,
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
sample <- samples[1]

# quit(save="ask")

summary_file <- here(data_folder, "data.Rds")
if (!file.exists(summary_file)) {
  source(here("r-scripts", "01-truvari-report.R"))
}
data <- readRDS(file = summary_file)

matches <- data.frame()
for (caller in callers) {
  tryCatch(
    {
      print(caller)
      vcf <- vcf_read(sample, caller, include = "tp-base", data_folder = truvari_run)[[1]]
      tp_base <- vcf_process(vcf, sample = "NA12879", caller = caller, status = "tp-base") %>% mutate(id = str_c(
        seqnames, "_",
        start, "-", end,
        "_", SVTYPE
      ))
      tp_base <- tp_base %>% mutate(caller = caller)
      matches <- bind_rows(matches, tp_base)
    },
    error = function(e) {
      message("An Error Occurred")
      print(e)
    }
  )
}

head(matches)
summary_file <- here(data_folder, "data.Rds")
color_file <- here("pedigree-results", "truvari", "color_pal.Rds")
if (!file.exists(summary_file) | !file.exists(color_file)) { # make sure these files exist
  source(here("r-scripts", "01-truvari-summary.R"))
}
# load the f1 summary data
summary_data <- readRDS(file = summary_file)
# get the color and orders of the callers
color_pal <- readRDS(color_file)


matches <- matches %>% dplyr::select(c(id, caller, QUAL, PctSizeSimilarity, PctSeqSimilarity, SVTYPE, width, TruScore, PctRecOverlap))
upset_input <- data.frame(matches %>% mutate(value = TRUE) %>%
  pivot_wider(
    id_cols = id, names_from = caller,
    values_from = value, values_fill = FALSE, values_fn = function(x) any(x),
    unused_fn = c(
      QUAL = mean, PctSizeSimilarity = mean,
      PctSeqSimilarity = mean,
      SVTYPE = unique, width = unique,
      NumClusterSVs = mean, TruScore = mean,
      PctRecOverlap = mean
    )
  )) %>%
  dplyr::rename(name = id) # id caues a duplicate error
head(upset_input)
upset_color_pal <- as.character(color_pal)
names(upset_color_pal) <- names(color_pal)

# hack
color_metadata <- data.frame(
  set = callers,
  callers = callers
)

upset_callers <- callers[callers %in% names(upset_input)]

# only keep colors for callers that are in your input
real_callers <- intersect(names(color_pal), upset_callers)
cp <- color_pal[real_callers]

my_queries <- map2(
  names(cp), cp,
  ~ upset_query(set = .x, fill = .y)
)


sv_colors <- c(
  "DEL" = "#D55E00", # Muted red
  "INS" = "#0072B2", # Muted blue
  "INV" = "#FF0000" # Bright red
)
upset_input <- upset_input %>%
  mutate(SVTYPE = factor(SVTYPE, levels = c("INV", "DEL", "INS")))
## SIMPLE INTERSECTION PLOTS

upset_general <- function(sort = c("degree", "cardinality"), min_size = 30,
                          save_name = NA) {
  upset(upset_input, upset_callers,
    name = "callers",
    width_ratio = 0.1, min_size = min_size,
    base_annotations = list(
      "Intersection size" = intersection_size(
        counts = TRUE,
        mapping = aes(fill = SVTYPE)
      ) + ylab("Common True Positive Variants Based on Truth Set") + scale_fill_manual(values = sv_colors)
    ),
    guides = "collect",
    sort_intersections_by = sort,
    matrix = (intersection_matrix(geom = geom_point(shape = 21, size = 3, color = "#ffffff00"))
    + scale_color_manual(
        values = color_pal,
        guide = guide_legend(override.aes = list(shape = "circle"))
      )),
    queries = my_queries,
    set_sizes = (
      upset_set_size(
        geom = geom_bar(
          aes(fill = SVTYPE, x = group),
          width = 0.8
        ),
        position = "left"
      ) + ylab("Overall True Positives\nper Caller") + scale_fill_manual(guide = "none", values = sv_colors))
  )
  if (is.na(save_name)) {
    save_name <- str_c("upset-", sort[[1]], "sort_callercolored_cutoff", min_size, ".png")
  }
  ggsave(width = 16, height = 9, here(data_folder, save_name))
}


### colored by caller
upset_general()

upset_general(min_size = 5)

upset_general(sort = c("cardinality", "degree"))


### TYPE SPECIFIC
upset(upset_input %>% filter(SVTYPE == "INS"), upset_callers,
  name = "callers",
  width_ratio = 0.1, min_size = 30,
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE) + ylab("Common True Positive Insertions Based on Truth Set")
  ),
  set_sizes = (upset_set_size() + ylab("True Positives Insertions\nper Caller")),
  guides = "collect",
  sort_intersections_by = c("degree", "cardinality")
)
ggsave(width = 16, height = 9, here(data_folder, "upset-degreesort-insertions.png"))

upset(upset_input %>% filter(SVTYPE == "DEL"), upset_callers,
  name = "callers",
  width_ratio = 0.1, min_size = 30,
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE) + ylab("Common True Positive Deletions Based on Truth Set")
  ),
  set_sizes = (upset_set_size() + ylab("True Positive Deletions\nper Caller")),
  guides = "collect",
  sort_intersections_by = c("degree", "cardinality")
)
ggsave(width = 16, height = 9, here(data_folder, "upset-degreesort-deletions.png"))


### NO COLORS
upset(upset_input, upset_callers,
  name = "callers",
  width_ratio = 0.1, min_size = 30,
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = SVTYPE)
    ) + ylab("Common True Positive Variants Based on Truth Set") + scale_fill_manual(values = sv_colors)
  ),
  set_sizes = (
    upset_set_size(
      geom = geom_bar(
        aes(fill = SVTYPE, x = group),
        width = 0.8
      ),
      position = "left"
    ) + ylab("Overall True Positives\nper Caller") + scale_fill_manual(values = sv_colors, guide = "none")),
  guides = "collect",
  sort_intersections_by = c("degree", "cardinality")
)
ggsave(width = 16, height = 9, here(data_folder, "upset-degreesort.png"))

upset(upset_input, upset_callers,
  name = "callers",
  width_ratio = 0.1, min_size = 30,
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = SVTYPE)
    ) + scale_fill_manual(values = sv_colors) + ylab("Common True Positive Variants Based on Truth Set")
  ),
  set_sizes = (
    upset_set_size(
      geom = geom_bar(
        aes(fill = SVTYPE, x = group),
        width = 0.8
      ),
      position = "left"
    ) + ylab("Overall True Positives per Caller") + scale_fill_manual(values = sv_colors, guide = "none")),
  guides = "collect",
  sort_intersections_by = c("cardinality", "degree"),
)
ggsave(width = 16, height = 9, here(data_folder, "upset-freqsort.png"))

library(VariantAnnotation)
library(tidyverse)
library(here)

theme_set(theme_bw())
here::i_am("workflow/r-scripts/00-vcf_process.R")
options(repr.plot.width = 8, repr.plot.height = 4.5)

vcf_read <- function(sample = NA, caller = NA, data_folder = "type-ignored",
                     include = "tp-comp|fp|fn", restrict = TRUE, param = ScanVcfParam(), phab = FALSE, debug = FALSE) {
  tryCatch(
    {
      vcfs <- NULL
      folder_name <- here(truvari_folder, data_folder)

      vcf_files <-
        list.files(
          path = folder_name,
          recursive = TRUE,
          pattern = str_c("*.vcf.gz$"),
          full.names = TRUE
        ) # get a list of all the summary files

      if (!is.na(caller) & !is.na(sample)) {
        vcf_files <- vcf_files[grep(str_c(caller, "/", sample), vcf_files)]
      } else if (is.na(caller)) {
        vcf_files <- vcf_files[grep(sample, vcf_files)]
      } else if (is.na(sample)) {
        vcf_files <- vcf_files[grep(caller, vcf_files)]
      }
      vcf_files <- vcf_files[grep(include, vcf_files, invert = FALSE)]
      vcf_files <- vcf_files[grep("phab", vcf_files, invert = !phab)]
      print(paste("vcf files for", sample, caller))
      print(vcf_files)
      vcfs <- lapply(vcf_files, readVcf, genome = "hg38", param = param)
      names(vcfs) <- str_remove(str_split_i(str_remove(vcf_files, folder_name), pattern = "/", 4), ".vcf.gz")

      for (vcf in vcfs) {
        if (restrict) {
          # some callers (octopus)adont label with a st
          if ("SVTYPE" %in% names(vcf@info)) {
            vcf@info <- vcf@info[c("SVTYPE", "PctSizeSimilarity", "PctSeqSimilarity")]
          } else {
            vcf@info <- vcf@info[c("PctSizeSimilarity", "PctSeqSimilarity")]
            vcf@info$SVTYPE <- NA
          }
        }
        if (debug) {
          print("INFO names")
          print(names(vcf@info))
        }
        seqlevels(vcf) <-
          str_remove(seqlevels(vcf), "chr") # change the chr names for plotting
        # try(end(vcf) <- vcf@info$END)
        # properly add the ends
      }
      return(vcfs)
    },
    error = function(e) {
      message("An Error Occurred")
      print(e)
    }
  )
  return(NA)
}

vcf_process <- function(vcf, sample = NA, caller = NA, status = NA, restrict = FALSE, debug = FALSE, sample_names = NA) {
  if (debug) {
    print(paste("Processing vcf for", sample, caller, status))
  }
  test <- as.data.frame(rowRanges(vcf), optional = TRUE, row.names = NULL)
  type <- as.data.frame(vcf@info) %>% unnest_wider(col = where(is.list), names_sep = "_", simplify = TRUE)
  if (restrict) {
    type <- type %>% dplyr::select(any_of(c(
      "SVTYPE", "PctSizeSimilarity",
      "TruScore", "MatchId", "SVLEN"
    )))
    if (debug) {
      print("Restricted")
    }
  }
  if (debug) {
    print(names(data))
  }
  data <- bind_cols(test, type)

  data <- data %>%
    rownames_to_column() %>%
    mutate(sample = sample) %>%
    mutate(caller = caller) %>%
    mutate(status = status) %>%
    dplyr::select(-any_of("SC")) # i dont know what this is, but its problematic
  return(data)
}
vcf_process_all_raw <- function(truvari_folder = "pedigree-results/truvari", data_folder = "type-ignored",
                                callers = NA, samples = NA,
                                useRDS = FALSE, RDSname = "vcf_data",
                                param = ScanVcfParam(), debug = TRUE,
                                what = c("tp-comp", "fp", "fn"), restrict = TRUE) {
  truvari_name <- str_replace(str_remove(data_folder, "/"), "-", " ")
  folders <-
    list.dirs(
      path = here::here(truvari_folder, data_folder),
      full.names = FALSE,
      recursive = TRUE
    )
  print(folders)
  if (length(callers) == 1 && is.na(callers)) {
    callers <- unique(str_split_i(str_remove(folders, data_folder), "/", 1))[-1]
  }
  if (length(samples) == 1 && is.na(samples)) {
    samples <- unique(str_split_i(str_remove(folders, data_folder), "/", 2))[-1]
  }
  print(samples)
  print(callers)
  overall_data <- data.frame()
  RDSfile <- here::here(truvari_folder, str_c(data_folder, "/", RDSname, ".Rds"))
  if (useRDS & file.exists(RDSfile)) { # use that hard won vcf data :D
    print("LOADING OLD RDS DATA.....")
    overall_data <- readRDS(file = RDSfile)
  } else {
    for (caller_i in callers) {
      for (sample_i in samples) {
        print(caller_i)
        vcfs <- vcf_read(sample_i, caller_i,
          debug = TRUE,
          data_folder = data_folder, param = param,
          include = str_c(what, collapse = "|"), restrict = restrict
        )
        print("VCFs loaded")
        if (length(vcfs) == length(what)) { # make sure we're loading all 3

          for (i in seq(1, length(vcfs))) {
            print(str_c("Processing ", names(vcfs)[i]))
            data <- vcf_process(vcfs[[i]],
              restrict = restrict, debug = TRUE,
              sample = sample_i, caller = caller_i, status = names(vcfs)[i]
            )
            if (length(overall_data) > 1) {
              data <- data |> dplyr::select(any_of(names(overall_data)))
            }
            overall_data <- bind_rows(overall_data, data) 
            # adapted from janitor::remove_empty 
            overall_data <- overall_data[rowSums(is.na(overall_data)) != ncol(overall_data),,drop=FALSE]
            overall_data <- overall_data[,colSums(is.na(overall_data)) != nrow(overall_data),drop=FALSE]
          }
        }
      }
    }
    # save that hard won vcf data :)
    try(saveRDS(overall_data, file = RDSfile))
    return(overall_data)
  }


  vcf_process_all <- function(truvari_folder = "pedigree-results/truvari", data_folder = "type-ignored",
                              callers = NA, samples = NA,
                              useRDS = FALSE, RDSname = "vcf_data",
                              param = ScanVcfParam(),
                              what = c("tp-comp", "fp", "fn"), restrict = TRUE) {
    overall_data <- vcf_process_all_raw(truvari_folder, data_folder, callers, samples, useRDS, RDSname, param, what, restrict = restrict)
    overall_data <- rename(overall_data, c("SVTYPE" = "type"))

    print("Loading color pal....")
    # get the color and orders of the callers
    color_pal <- readRDS(here(truvari_folder, "color_pal.Rds"))
    print(color_pal)

    if (length(unique(overall_data$caller) != length(names(color_pal)))) {
      print(str_c("We are missing ", unique(overall_data$caller)[!(unique(overall_data$caller) %in% names(color_pal))], " in the summary data"))
      print(str_c("We are missing ", names(color_pal)[!(names(color_pal) %in% unique(overall_data$caller))], " in the vcfs"))
    }

    overall_data <- overall_data %>%
      mutate(sample = as_factor(sample)) %>%
      mutate(caller = factor(str_to_lower(caller), levels = names(color_pal), )) %>%
      filter(!is.na(caller)) %>% # remove anything not mapped properly
      mutate(status = case_when(status == "fn" ~ "False Negative", status == "fp" ~ "False Positive", status == "tp-comp" ~ "True Positive"))
    print(levels(overall_data$sample))
    sample_names <- c("200915", "200921", "NA24385")
    focus_sample <- sample_names[1]
    levels(overall_data$sample) <- sample_names

    # remove extra typing
    overall_data <- overall_data |>
      mutate(type = case_when(caller == "tardis" ~ str_split_i(type, ":", 1), type == "DUP:TANDEM" ~ "DUP", TRUE ~ type))


    # save that hard won vcf data :)
    saveRDS(overall_data, file = RDSfile)
  }

  return(overall_data)
}

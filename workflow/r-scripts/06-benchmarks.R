library(here)
library(tidyverse)
library(tidyr)
library(purrr)
library(lubridate)
library(ggpubr)

theme_set(theme_minimal())
options(scipen = 999)
# read in all the data
here::i_am("r-scripts/03-benchmarks.R")
data_folder <- here("pedigree-results", "benchmarks")
files <- list.files(path = data_folder, full.names = TRUE, pattern = "*.tsv", recursive = TRUE)
print(files)
data <- files |>
  map(read_delim, col_types = "d_ddddddd", delim = "\t") |>
  reduce(bind_rows)
data$file <- str_remove(str_remove(files, data_folder), ".tsv")

# get the callers and samples
data <- data |>
  mutate(file = str_remove(file, data_folder)) |>
  separate(file, c(NA, "caller", "sample"), sep = "/")
data <- data |> mutate(sample = as.factor(sample))

# some callers are spread across multiple benchmarks
data <- data |> separate(caller, c("caller", "subprocess"), sep = "-")

print(data)
# sum up the subprocesses
data <- data |>
  group_by(caller, sample) |>
  summarise(s = sum(s), cpu_time = sum(cpu_time), max_rss = max(max_rss), max_vms = max(max_vms), max_uss = max(max_uss), max_pss = max(max_pss), io_in = sum(io_in), io_out = sum(io_out), mean_load = mean(mean_load))

# pretty format seconds
data <- data |> mutate(time = lubridate::seconds_to_period(s), cpu_time_pretty = lubridate::seconds_to_period(cpu_time))

write.csv(data, file = here(data_folder, "caller-summary.csv"))

print(data)
print(levels(data$sample))

# grab the colors
color_pal <- readRDS(here("pedigree-results", "truvari", "color_pal.Rds"))
# and subset only to valid
print(unique(data$caller))
print(names(color_pal))
data <- data %>%
  filter(caller %in% names(color_pal)) %>%
  mutate(caller = factor(caller, levels = names(color_pal)))
# and plot!
# Running time (s)
runtime <- ggplot(data, aes(fill = caller, x = caller, y = s)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "Running Time (s)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "max_sec.png"))
# CPU Time (s)
runtime <- ggplot(data, aes(fill = caller, x = caller, y = cpu_time)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "CPU Time (s)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "cpu_time.png"))
# CPU Time (s)
runtime <- ggplot(data, aes(fill = caller, x = caller, y = cpu_time_pretty)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "CPU Time") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "cpu_time_pretty.png"))


run_to_cpu <- ggplot(data, aes(color = caller, x = s, y = cpu_time, shape = sample)) +
  geom_point() +
  scale_color_manual(values = color_pal) +
  labs(x = " Run time (s)", y = "CPU time (s)")
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "run_to_cpu.png"))

cpu_to_mem <- ggplot(data, aes(color = caller, x = cpu_time, y = max_rss, shape = sample)) +
  geom_point() +
  scale_color_manual(values = color_pal) +
  labs(x = "CPU time (s)", y = "RAM allocated (Mb)")
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "cpu_to_rss.png"))

# RSS - RAM allocated during allocation, including preloaded libraries
rss <- ggplot(data, aes(fill = caller, x = caller, y = max_rss)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "RAM allocated (Mb)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "max_rss.png"))

ggplot(data, aes(fill = caller, x = caller, y = max_vms)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "VMS") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "max_vms.png"))

# Unique Set Size libraries and pages allocated only to this process
ggplot(data, aes(fill = caller, x = caller, max_uss)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "USS (Mb)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "max_uss.png"))

# Mean load
load <- ggplot(data, aes(fill = caller, x = caller, y = mean_load)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "Mean Load") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "mead_load.png"))

# IO in
ioin <- ggplot(data, aes(fill = caller, x = caller, y = io_in)) +
  stat_summary(fun = "mean", geom = "col") +
  scale_fill_manual(values = color_pal, guide = "none") +
  geom_point(aes(shape = sample)) +
  labs(x = "Callers", y = "Max IO In") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggplot2::ggsave(width = 16, height = 9, file = here(data_folder, "max_ioin.png"))

figure <- ggarrange(runtime, rss, load, ioin, ncol = 2, nrow = 2, common.legend = TRUE)
figure
ggexport(figure, filename = here(data_folder, "summary.png"), width = 1600, height = 900)

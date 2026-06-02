# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(tidyr)


# Read the datasets
union_data <- read_csv("exhaustive-union-result.csv") # can be found in 10.5281/zenodo.20410726
consensus_data <- read_csv("exhaustive-consensus-result.csv") # can be found in 10.5281/zenodo.20410726

# Add support indicator
union_data$support <- 1
consensus_data$support <- 2

# Merge the datasets
combined_data <- rbind(union_data, consensus_data)

# Convert support to factor with labels
combined_data$support_label <- factor(combined_data$support, 
                                      levels = c(1, 2), 
                                      labels = c("Union", "Concensus"))

# Round F1 scores to 2 decimal places
combined_data$f1_rounded <- round(combined_data$f1, 3)


top_10_consensus <- consensus_data %>%
  arrange(desc(f1)) %>%
  head(10)

# Round F1 scores to 2 decimal places
top_10_consensus$f1_rounded <- round(top_10_consensus$f1, 4)
top_10_consensus$combo <- gsub("_vote2$", "", top_10_consensus$combo)


# Create precision-recall scatter plot
consensus_pr_plot <- ggplot(top_10_consensus, aes(x = precision, y = recall, color = f1)) +
  geom_point(size = 4, alpha = 0.8, shape = 17) +
  geom_text_repel(aes(label = paste0(combo, "\n(F1: ", f1_rounded, ")")), 
                  size = 4.5, 
                  max.overlaps = 25,
                  box.padding = 0.5,
                  point.padding = 0.3) +
  labs(title = "top 10 concensus combinations by F1",
       x = "Precision",
       y = "Recall",
       color = "F1 Score") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_gradient(low = "#FFB6C1", high = "#E31A1C", 
                       name = "F1 Score")


print(consensus_pr_plot)


top_10_union <- union_data  %>%
  arrange(desc(f1)) %>%
  head(10)

# Round F1 scores to 2 decimal places
top_10_union$f1_rounded <- round(top_10_union$f1, 4)


# Create precision-recall scatter plot
consensus_pr_plot_union <- ggplot(top_10_union, aes(x = precision, y = recall, color = f1)) +
  geom_point(size = 4, alpha = 0.8, shape = 17) +
  geom_text_repel(aes(label = paste0(combo, "\n(F1: ", f1_rounded, ")")), 
                  size = 4.5, 
                  max.overlaps = 25,
                  box.padding = 0.5,
                  point.padding = 0.3) +
  labs(title = "top 10 union combinations by F1",
       x = "Precision",
       y = "Recall",
       color = "F1 Score") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_gradient(low = "#ADD8E6", high = "#1F78B4", 
                       name = "F1 Score")


print(consensus_pr_plot_union)


# ──────────────────────────────────────────────────────
INDIVIDUAL_CALLERS <- c("Dysgu", "Manta", "Octopus", "DELLY", "SvABA")

UNION_COMBOS <- c(
  "Manta+Dysgu+Octopus+Smoove+LUMPY+PopDel+Tardis"
)

CONSENSUS_COMBOS <- c(
  "Manta+Dysgu+Octopus+DELLY+Smoove+SvABA+LUMPY+PopDel+Tardis_vote2"
)
# ──────────────────────────────────────────────────────────────────────────────

individual_rows <- union_data %>%
  filter(combo %in% INDIVIDUAL_CALLERS) %>%
  mutate(group = "Individual", type = combo)

union_rows <- union_data %>%
  filter(combo %in% UNION_COMBOS) %>%
  mutate(group = "Union", type = paste0("Union: ", combo))

consensus_rows <- consensus_data %>%
  filter(combo %in% CONSENSUS_COMBOS) %>%
  mutate(group = "Consensus", type = paste0("Consensus: ", gsub("_vote2$", "", combo)))

plot_df <- bind_rows(individual_rows, union_rows, consensus_rows) %>%
  select(type, group, precision, recall, f1) %>%
  pivot_longer(cols = c(recall, precision, f1),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("recall", "precision", "f1"),
                    labels = c("Recall", "Precision", "F1")),
    group  = factor(group, levels = c("Individual", "Union", "Consensus")),
    type   = factor(type, levels = c(
      unique(individual_rows$type),
      unique(union_rows$type),
      unique(consensus_rows$type)
    ))
  )

metric_colours <- c(
  "Recall"    = "#E41A1C",
  "Precision" = "#377EB8",
  "F1"        = "goldenrod3"
)

# single separator between individuals and combos
n_individuals <- length(unique(individual_rows$type))

comparison_plot <- ggplot(plot_df, aes(x = type, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.8) +
  geom_vline(xintercept = n_individuals + 0.5,
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  scale_fill_manual(values = metric_colours) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Individual Callers vs Combinations",
    x     = NULL,
    y     = "score",
    fill  = "metric"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 14),
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 9),
    axis.title.y    = element_text(size = 12),
    legend.position = "top",
    legend.title    = element_text(size = 11)
  )

print(comparison_plot)



# Create precision-recall scatter plot
consensus_all <- ggplot(consensus_data, aes(x = precision, y = recall, color = f1)) +
  geom_point(size = 0.7, alpha = 0.8, shape = 17) +
  labs(title = "consensus combinations",
       x = "Precision",
       y = "Recall",
       color = "F1 Score") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_gradient(low = "#FFB6C1", high = "#E31A1C", 
                       name = "F1 Score")


ggsave("consensus_all.png", plot = consensus_all, 
       width = 8, height = 6, units = "in", 
       dpi = 300, bg = "white")

print(consensus_all)


# Create precision-recall scatter plot
union_all <- ggplot(union_data, aes(x = precision, y = recall, color = f1)) +
  geom_point(size = 1, alpha = 0.8, shape = 17) +
  labs(x = "Precision",
       y = "Recall",
       color = "F1 Score") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  scale_color_gradient(low = "#ADD8E6", high = "#1F78B4", 
                       name = "F1 Score")


ggsave("union_all.png", plot = union_all, 
       width = 8, height = 6, units = "in", 
       dpi = 300, bg = "white")



print(union_all)
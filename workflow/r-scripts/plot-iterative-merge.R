library(tidyverse)
library(ggplot2)
library(ggrepel)

iterations <- read_excel("iterative-merge-result.xlsx") #can be found in 10.5281/zenodo.20410726
data <- iterations %>% na.omit() 

highlight_combinations <- c(
  "Dysgu", 
  "Dysgu+Octopus", 
  "Dysgu+Octopus+Manta", 
  "Dysgu+Octopus+Manta+Tardis",
  "Dysgu+Octopus+Manta+Tardis+SvABA",
  "Dysgu+Octopus+Manta+Tardis+SvABA+CNVpytor",
  "Dysgu+Octopus+Manta+Tardis+SvABA+CNVpytor+Tiddit",
  "Dysgu+Octopus+Manta+Tardis+SvABA+CNVpytor+Tiddit+PopDel"
)

# Add a column to mark highlighted combinations
data$is_highlighted <- data$combinations %in% highlight_combinations

# Add formatted labels for recall and f1
data$recall_label <- paste0(data$combinations, " (", sprintf("%.3f", data$recall), ")")
data$f1_label <- paste0(data$combinations, " (", sprintf("%.3f", data$f1), ")")

# Create plot styling function to maintain consistent look between plots
create_styled_plot <- function(data, metric_name, metric_column, label_column, metric_color) {
  # Create the plot
  plot <- ggplot() +
    # Add grey points for all non-highlighted combinations
    geom_point(data = subset(data, !is_highlighted), 
               aes(x = iteration, y = !!sym(metric_column)),
               color = "darkgray", alpha = 0.7, size = 2) +
    
    # Add grey labels for non-highlighted combinations - with adjusted parameters for spacing
    geom_text_repel(
      data = subset(data, !is_highlighted),
      aes(x = iteration, y = !!sym(metric_column), label = !!sym(label_column)),
      size = 3.0,
      color = "darkgray",
      segment.color = "lightgray",
      segment.size = 0.2,
      box.padding = 0.5,       
      point.padding = 0.3,     
      force = 12,              
      max.overlaps = 8,        
      nudge_x = 0.2,           
      direction = "both",      
      hjust = 0                
    ) +
    
    # Add colored points for highlighted combinations
    geom_point(data = subset(data, is_highlighted), 
               aes(x = iteration, y = !!sym(metric_column)),
               color = metric_color, size = 3.5) +
    
    # Add lines connecting the highlighted combinations
    geom_line(data = subset(data, is_highlighted),
              aes(x = iteration, y = !!sym(metric_column), group = 1),
              color = metric_color, linewidth = 0.8) +
    
    # Add colored labels for highlighted combinations - with adjusted parameters for spacing
    geom_text_repel(
      data = subset(data, is_highlighted),
      aes(x = iteration, y = !!sym(metric_column), label = !!sym(label_column)),
      color = metric_color,
      fontface = "bold",
      size = 4.0,
      segment.color = metric_color,
      segment.size = 0.4,
      box.padding = 0.8,       
      point.padding = 0.5,     
      force = 20,              
      nudge_x = 0.3,           
      direction = "y",         
      min.segment.length = 0   
    ) +
    
    # Adjust scales and labels - wider x-axis to accommodate spread-out labels
    scale_x_continuous(breaks = 1:7, limits = c(0.7, 7.5)) +  
    scale_y_continuous(limits = c(0.0, 0.5)) +
    
    labs(title = paste(metric_name, "of Selected Tool Combinations Across Iterations"),
         x = "Iteration (Number of Combined Tools)",
         y = metric_name) +
    
    theme_bw() +
    theme( 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(10, 30, 10, 10)  
    )
  
  return(plot)
}

# Create individual plots with the metric-specific labels
recall_plot <- create_styled_plot(data, "Recall", "recall", "recall_label", "#E41A1C")
f1_plot <- create_styled_plot(data, "F1 Score", "f1", "f1_label", "goldenrod3")

print(f1_plot)

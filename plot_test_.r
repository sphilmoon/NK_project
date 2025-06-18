# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

# Read CSV using base R
data <- read.csv("asc_mock_data.csv", stringsAsFactors = FALSE)

# Set time_point as factor with proper levels for ordering
data$time_point <- factor(data$time_point, levels = c(
  "0M", "1M", "1M+3D", "1M+6D",
  "3M", "3M+3D", "3M+6D",
  "6M", "6M+3D", "6M+7D",
  "12M", "12M+3D", "12M+6D"
))

# Shape and fill definitions
strain_shapes <- c("O1M" = 21, "Asia1" = 24)
strain_fills  <- c("O1M" = "#999999", "Asia1" = "#555555")

# Plot function per group
plot_group <- function(group_label, arrow_x = NULL, arrow_label = NULL, arrow_color = "black") {
  group_data <- data %>% filter(group == group_label)

  p <- ggplot(group_data, aes(x = time_point, y = asc_count, fill = strain, shape = strain)) +
    geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = asc_count - sd, ymax = asc_count + sd),
                  position = position_dodge(width = 0.8), width = 0.25) +
    geom_point(position = position_dodge(width = 0.8), size = 3, color = "black") +
    scale_shape_manual(values = strain_shapes) +
    scale_fill_manual(values = strain_fills) +
    ylim(0, 150) +
    labs(
      y = if (group_label == "Group 1") "Serotype-specific ASC per million cells" else NULL,
      x = if (group_label == "Group 3") "Months (and Days) post prime vaccination" else NULL
    ) +
    ggtitle(group_label) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 1, face = "bold"),
      legend.position = "none"
    )

  # Add optional arrow and label
  if (!is.null(arrow_x) && !is.null(arrow_label)) {
    p <- p +
      annotate("segment", x = arrow_x, xend = arrow_x, y = 150, yend = 90,
               arrow = arrow(length = unit(0.25, "cm")), color = arrow_color, size = 1.2) +
      annotate("text", x = arrow_x, y = 155, label = arrow_label, angle = 45, hjust = 0)
  }

  return(p)
}

# Create plots with annotations
p1 <- plot_group("Group 1", arrow_x = 1, arrow_label = "Prime", arrow_color = "magenta") +
      annotate("segment", x = 2, xend = 2, y = 150, yend = 90,
               arrow = arrow(length = unit(0.25, "cm")), color = "magenta", size = 1.2) +
      annotate("text", x = 2, y = 155, label = "Boost", angle = 45, hjust = 0)

p2 <- plot_group("Group 2", arrow_x = 5, arrow_label = "Boost", arrow_color = "blue")
p3 <- plot_group("Group 3", arrow_x = 9, arrow_label = "Boost", arrow_color = "darkgreen")

# Combine plots vertically
final_plot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))

# Optional: Add legend using cowplot
legend_plot <- ggplot(data, aes(x = time_point, y = asc_count, fill = strain, shape = strain)) +
  geom_point() +
  scale_shape_manual(values = strain_shapes) +
  scale_fill_manual(values = strain_fills) +
  theme_void() +
  theme(legend.position = "bottom")

legend <- cowplot::get_legend(legend_plot)

# Save the figure
ggsave("ASC_combined_group_plot.pdf", plot = final_plot, width = 12, height = 12, dpi = 600)

# # Print to viewer
# print(final_plot)
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)

# Optional: Define your input CSV file
# Your CSV should have columns: time_point, strain, days_post_immunization, asc_per_million
data <- read_csv("your_input_data.csv")

# Preview expected structure
# head(data)
#   time_point      strain       days_post_immunization     asc_per_million
#   "Day 7"         "Strain A"    7                        105
#   "Day 7"         "Strain B"    7                        98
#   "Day 14"        "Strain A"    14                       130
#   ...

# Ensure 'time_point' is treated as a factor with desired order
data$time_point <- factor(data$time_point, levels = c("Day 7", "Day 14", "Day 21"))

# Set consistent color palette for strains
strain_colors <- c("Strain A" = "#1f77b4", "Strain B" = "#ff7f0e")

# Create a plotting function
plot_group <- function(time) {
  ggplot(filter(data, time_point == time), aes(x = factor(days_post_immunization),
                                               y = asc_per_million,
                                               fill = strain)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    scale_fill_manual(values = strain_colors) +
    ylim(0, 150) +
    labs(
      title = paste("Immunization Time Point:", time),
      x = "Days Post Immunization",
      y = "ASC per Million"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10)
    )
}

# Generate individual plots
plot_day7 <- plot_group("Day 7")
plot_day14 <- plot_group("Day 14")
plot_day21 <- plot_group("Day 21")

# Combine with patchwork and add shared legend
final_plot <- (plot_day7 / plot_day14 / plot_day21) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

# Save as high-resolution PDF
ggsave("ASC_combined_figure.pdf", plot = final_plot, width = 10, height = 12, dpi = 600)
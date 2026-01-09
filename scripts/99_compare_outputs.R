# scripts/99_compare_outputs.R
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork) # For combining plots

# 1. Load Data
my_data <- read_csv("data/selected_environmental_variables.csv")
paper_data <- read_csv("data/paper_guy_selected_environmental_variables.csv") # Assuming you have this

# 2. Create Plot for YOUR Data (PC1)
p1 <- ggplot(my_data, aes(x=x, y=y, color=PC1)) +
  geom_point(size=0.1) +
  scale_color_viridis_c(option = "magma") +
  coord_fixed() +
  labs(title = "My Data: PC1", subtitle = "PCA Axis 1 (Climate)") +
  theme_minimal()

# 3. Create Plot for YOUR Data (Rugosity)
p2 <- ggplot(my_data, aes(x=x, y=y, color=rugosity)) +
  geom_point(size=0.1) +
  scale_color_viridis_c(option = "cividis", limits = c(0, 50)) + # Cap at 50 to see contrast
  coord_fixed() +
  labs(title = "My Data: Rugosity", subtitle = "Static Terrain") +
  theme_minimal()

# 4. Create Plot for PAPER Data (SST Mean - proxy for PC1)
p3 <- ggplot(paper_data, aes(x=x, y=y, color=sstmean)) +
  geom_point(size=0.1) +
  scale_color_viridis_c(option = "plasma") +
  coord_fixed() +
  labs(title = "Paper Guy: SST Mean", subtitle = "Raw Variable") +
  theme_minimal()

# 5. Combine and Save
final_plot <- p1 / p3
ggsave("outputs/comparison_check.png", final_plot, width=12, height=10)

print(final_plot)
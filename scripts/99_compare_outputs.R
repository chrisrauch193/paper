# scripts/99_compare_outputs.R
# ------------------------------------------------------------------------------
# COMPARISON PLOTTER (Smart Mode)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, readr, patchwork, viridis)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs", "figures")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# 1. Load Data
cat("Loading datasets...\n")
my_data <- read_csv(file.path(DATA_DIR, "selected_environmental_variables.csv"), show_col_types = FALSE)

# Check if Paper Guy data exists (for side-by-side)
paper_path <- file.path(DATA_DIR, "paper_guy_selected_environmental_variables.csv")
has_paper_data <- file.exists(paper_path)
if(has_paper_data) {
  paper_data <- read_csv(paper_path, show_col_types = FALSE)
  # Ensure 0-360 shift for comparison if he isn't already
  if(min(paper_data$x, na.rm=T) < 0) paper_data$x <- paper_data$x + 360
}

# 2. Determine Plotting Logic (Auto-Detect Mode)
if ("PC1" %in% names(my_data)) {
  # --- EXPANSION MODE (PCA) ---
  cat("Detected EXPANSION MODE (Plotting PC1)...\n")
  plot_var <- "PC1"
  plot_label <- "My Data (PC1 - Expansion)"
  
} else {
  # --- REPLICATION MODE (Raw Vars) ---
  cat("Detected REPLICATION MODE (Plotting Raw Variables)...\n")
  # Find the temperature variable (usually thetao or sst) to compare with his SST
  temp_col <- grep("thetao|sst", names(my_data), value = TRUE)[1]
  
  if(is.na(temp_col)) temp_col <- names(my_data)[3] # Fallback to 3rd column if no temp found
  
  plot_var <- temp_col
  plot_label <- paste("My Data (", temp_col, " - Replication)")
}

# 3. Create YOUR Plot
p1 <- ggplot(my_data, aes(x=x, y=y, color=.data[[plot_var]])) +
  geom_point(size=0.1) +
  scale_color_viridis_c(option = "plasma") +
  coord_fixed(xlim=c(30, 260), ylim=c(-50, 50)) + # Standard View
  theme_minimal() +
  labs(title = plot_label, subtitle = paste(nrow(my_data), "points"))

# 4. Create PAPER GUY Plot (if available)
if(has_paper_data) {
  # He usually names it 'sstmean' or similar. Check his names.
  his_cols <- names(paper_data)
  his_var <- grep("sst|temp|thetao", his_cols, value=TRUE, ignore.case=TRUE)[1]
  
  if(is.na(his_var)) his_var <- his_cols[3] # Fallback
  
  p2 <- ggplot(paper_data, aes(x=x, y=y, color=.data[[his_var]])) +
    geom_point(size=0.1) +
    scale_color_viridis_c(option = "plasma") +
    coord_fixed(xlim=c(30, 260), ylim=c(-50, 50)) +
    theme_minimal() +
    labs(title = paste("Paper Guy Data (", his_var, ")"), 
         subtitle = paste(nrow(paper_data), "points"))
  
  # Combine
  final_plot <- p1 / p2
} else {
  final_plot <- p1
}

# 5. Save
out_file <- file.path(OUTPUT_DIR, "Env_Comparison_Check.png")
ggsave(out_file, final_plot, width=10, height=12, bg="white")
cat("Saved comparison map to:", out_file, "\n")

print(final_plot)
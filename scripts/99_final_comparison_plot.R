# scripts/99_final_comparison_plot.R
# ------------------------------------------------------------------------------
# FINAL COMPARISON: REPLICATION VS PAPER GUY (Direct CSV Plotting)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, readr, ggrepel, viridis, patchwork, stringr)

# --- CONFIGURATION ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs", "figures")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# INPUT FILES (Polygon Definitions)
USER_FILE <- file.path(DATA_DIR, "meow_ecos_df.csv")       # Your file
PAPER_FILE <- file.path(DATA_DIR, "paper_guy_meow_ecos_df.csv") # His file

# --- 1. DATA LOADING & PREP FUNCTION ---
prepare_data <- function(filepath, label_tag) {
  
  if(!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Load CSV
  df <- read_csv(filepath, show_col_types = FALSE)
  
  # Ensure standard columns (Handle minor naming diffs if any)
  # Expected: long, lat, PROVINCE, REALM, id (for polygon grouping)
  colnames(df) <- tolower(colnames(df)) # standardize to lowercase first
  
  # Rename back to standard
  df <- df %>% 
    rename(long = long, lat = lat, PROVINCE = province, REALM = realm, group = id)
  
  # 1. Shift Longitude to 0-360 (if not already)
  # If min long is negative, shift it.
  if(min(df$long, na.rm=TRUE) < 0) {
    df <- df %>% mutate(long = ifelse(long < 0, long + 360, long))
  }
  
  # 2. Calculate Centroids (For Labeling)
  # We group by Province to get one label per province
  centroids <- df %>% 
    group_by(PROVINCE, REALM) %>% 
    summarise(
      X_cent = mean(long), 
      Y_cent = mean(lat), 
      .groups = "drop"
    ) %>% 
    arrange(X_cent) %>% # Sort West to East
    mutate(Map_Index = row_number()) # Assign 1..N
  
  # 3. Clean Names for Key
  centroids$Clean_Name <- str_trunc(centroids$PROVINCE, 22, "right")
  
  # 4. Join Index back to main dataframe
  df_final <- df %>% 
    left_join(centroids %>% dplyr::select(PROVINCE, Map_Index), by = "PROVINCE")
  
  return(list(polygons = df_final, labels = centroids, tag = label_tag))
}

# --- 2. LOAD DATASETS ---
cat("Loading User Data...\n")
user_data <- prepare_data(USER_FILE, "My Replication")

cat("Loading Paper Guy Data...\n")
paper_data <- prepare_data(PAPER_FILE, "Paper Guy Original")

# --- 3. PLOTTING FUNCTION ---
create_map <- function(data_list) {
  
  polys <- data_list$polygons
  labs  <- data_list$labels
  title <- data_list$tag
  
  # Standard Colors
  realm_cols <- c(
    "Central Indo-Pacific" = "#E6AB02", 
    "Eastern Indo-Pacific" = "#A6761D",
    "Western Indo-Pacific" = "#66A61E", 
    "Temperate Northern Pacific" = "#7570B3",
    "Temperate Australasia" = "#E7298A", 
    "Temperate Southern Africa" = "#D95F02"
  )
  
  # MAP PLOT
  p_map <- ggplot() +
    # Background World (0-360 Pacific View)
    borders("world2", colour="gray90", fill="gray95", size=0.1) +
    
    # The Regions
    geom_polygon(data = polys, aes(x = long, y = lat, group = group, fill = REALM),
                 color = "white", size = 0.1, alpha = 0.9) +
    
    # The Labels (1, 2, 3...)
    geom_text_repel(data = labs, aes(x = X_cent, y = Y_cent, label = Map_Index),
                    size = 3, fontface = "bold",
                    min.segment.length = 0, # Always draw line
                    segment.color = "gray40", segment.size = 0.3,
                    box.padding = 0.4,
                    bg.color = "white", bg.r = 0.15) +
    
    scale_fill_manual(values = realm_cols) +
    
    # *** EXACT FRAMING ***
    coord_fixed(xlim = c(30, 240), ylim = c(-40, 40), expand = FALSE) +
    
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.background = element_rect(fill = "#E0F3F8", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank()
    )
  
  # KEY PLOT (Table below map)
  # Split into 4 columns
  n_rows <- ceiling(nrow(labs) / 4)
  labs$col_group <- rep(1:4, each = n_rows, length.out = nrow(labs))
  
  p_key <- ggplot(labs, aes(x = 0, y = -Map_Index, label = paste0(Map_Index, ". ", Clean_Name))) +
    geom_text(aes(color = REALM), hjust = 0, size = 2.8, fontface = "bold") +
    scale_color_manual(values = realm_cols) +
    facet_wrap(~col_group, scales = "free") +
    theme_void() +
    theme(legend.position = "none", strip.text = element_blank()) +
    scale_x_continuous(expand = c(0, 0.1))
  
  # Combine Map (Top) and Key (Bottom)
  final_plot <- p_map / p_key + plot_layout(heights = c(4, 1))
  
  return(final_plot)
}

# --- 4. GENERATE PLOTS ---

cat("Generating Plot 1 (User)...\n")
plot1 <- create_map(user_data)
ggsave(file.path(OUTPUT_DIR, "Comparison_1_User.png"), plot1, width = 12, height = 10, dpi = 300)

cat("Generating Plot 2 (Paper Guy)...\n")
plot2 <- create_map(paper_data)
ggsave(file.path(OUTPUT_DIR, "Comparison_2_PaperGuy.png"), plot2, width = 12, height = 10, dpi = 300)

# --- 5. SUMMARY COMPARISON PRINT ---

cat("\n======================================================\n")
cat("                COMPARISON SUMMARY\n")
cat("======================================================\n")

u_provs <- user_data$labels %>% arrange(PROVINCE) %>% pull(PROVINCE)
p_provs <- paper_data$labels %>% arrange(PROVINCE) %>% pull(PROVINCE)

cat("\n[COUNTS]\n")
cat("User Provinces:      ", length(u_provs), "\n")
cat("Paper Guy Provinces: ", length(p_provs), "\n")

cat("\n[DIFFERENCES]\n")
extra_in_user <- setdiff(u_provs, p_provs)
missing_in_user <- setdiff(p_provs, u_provs)

if(length(extra_in_user) > 0) {
  cat("Extra in User (Remove these for exact replication):\n")
  print(extra_in_user)
} else {
  cat("Extra in User: NONE\n")
}

if(length(missing_in_user) > 0) {
  cat("Missing in User (Add these):\n")
  print(missing_in_user)
} else {
  cat("Missing in User: NONE\n")
}

cat("\n======================================================\n")
cat("Done. Check 'outputs/figures' for Comparison_1_User.png and Comparison_2_PaperGuy.png\n")

# Display in RStudio
print(plot1)
print(plot2)
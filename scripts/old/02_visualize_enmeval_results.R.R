# scripts/new_pipeline/02_visualize_enmeval_results.R
#-------------------------------------------------------------------------------
# 1. Visualize NEW ENMeval Pipeline Results
# 2. Extract Metrics (AUC/CBI) and Summary
#-------------------------------------------------------------------------------

# --- Setup ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, ggplot2, sf, rnaturalearth, dplyr, readr, tidyr, viridis)

# Define paths
base_dir <- "/home/bi-server-kyoto/a0236995/paper/" # Your new project root
rds_dir  <- file.path(base_dir, "Rdata")
fig_dir  <- file.path(base_dir, "figures_new_pipeline")
res_dir  <- file.path(base_dir, "results_new_pipeline")

if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if(!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# Load World Map
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- USER CUSTOM COLOR SCALE ---
custom_colors <- c("#FEF9C3", "#FDBA74", "#D946EF", "#4F46E5")

#-------------------------------------------------------------------------------
# PART 1: VISUALIZATION FUNCTION (Adapted for SpatRaster & Unwrapping)
#-------------------------------------------------------------------------------
plot_and_save <- function(model_obj, model_type, output_folder) {
  
  if(is.null(model_obj$maps)) {
    warning(paste("No maps found for", model_type)); return(NULL)
  }
  
  # --- UNWRAP LOGIC ---
  # Check if it's a packed object (from terra::wrap) and unwrap it
  maps <- model_obj$maps
  if (inherits(maps, "PackedSpatRaster")) {
    maps <- terra::unwrap(maps)
  }
  
  species_names <- names(maps)
  
  cat(paste("Plotting", length(species_names), "maps for", model_type, "...\n"))
  
  for (sp in species_names) {
    # Convert SpatRaster layer to dataframe
    # terra::as.data.frame with xy=TRUE
    df <- as.data.frame(maps[[sp]], xy = TRUE)
    colnames(df) <- c("x", "y", "suitability")
    
    # Plot
    p <- ggplot() +
      geom_sf(data = world, fill = "grey80", color = "white", size = 0.1) +
      geom_raster(data = df, aes(x = x, y = y, fill = suitability)) +
      
      # Custom Colors
      scale_fill_gradientn(colors = custom_colors, limits = c(0, 1), name = "Probability") +
      
      coord_sf(xlim = c(30, 240), ylim = c(-40, 40), expand = FALSE) +
      labs(title = paste(model_type, "-", gsub("_", " ", sp)),
           subtitle = "New ENMeval Pipeline",
           x = NULL, y = NULL) +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white", color = NA),
            legend.position = "bottom")
    
    # Save
    filename <- file.path(output_folder, paste0(model_type, "_", sp, ".png"))
    ggsave(filename, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
  }
}

#-------------------------------------------------------------------------------
# PART 2: METRICS EXTRACTION FUNCTION
#-------------------------------------------------------------------------------
extract_metrics <- function(model_obj, model_type) {
  
  # In our new pipeline, we explicitly stored a dataframe called 'eval'
  if(is.null(model_obj$eval)) {
    warning(paste("No evaluation dataframe found for", model_type)); return(NULL)
  }
  
  df <- model_obj$eval
  df$ModelType <- model_type
  
  # Reorder columns
  df <- df %>% 
    dplyr::select(ModelType, Species, everything())
  
  return(df)
}

#-------------------------------------------------------------------------------
# EXECUTION ROUTINE
#-------------------------------------------------------------------------------

all_metrics <- list()

# 1. Host Anemones
if(file.exists(file.path(rds_dir, "anemENMs.RDS"))) {
  cat("\n--- Processing Anemone ENMs ---\n")
  res <- readRDS(file.path(rds_dir, "anemENMs.RDS"))
  plot_and_save(res, "Phase1_Anemone", fig_dir)
  all_metrics[["Phase1"]] <- extract_metrics(res, "Phase1_Anemone")
}

# 2. Clownfish Environment Only
if(file.exists(file.path(rds_dir, "amphENMs.RDS"))) {
  cat("\n--- Processing Clownfish Env Models ---\n")
  res <- readRDS(file.path(rds_dir, "amphENMs.RDS"))
  plot_and_save(res, "Phase2_Clownfish_Env", fig_dir)
  all_metrics[["Phase2"]] <- extract_metrics(res, "Phase2_Clownfish_Env")
}

# # 3. Clownfish Biotic Only (Host Only)
# if(file.exists(file.path(rds_dir, "amphBioticOnly.RDS"))) {
#   cat("\n--- Processing Clownfish Biotic Only Models ---\n")
#   res <- readRDS(file.path(rds_dir, "amphBioticOnly.RDS"))
#   plot_and_save(res, "Phase3_Clownfish_HostOnly", fig_dir)
#   all_metrics[["Phase3"]] <- extract_metrics(res, "Phase3_Clownfish_HostOnly")
# }
# 
# # 4. Clownfish Combined
# if(file.exists(file.path(rds_dir, "amphEBMs.RDS"))) {
#   cat("\n--- Processing Clownfish Combined Models ---\n")
#   res <- readRDS(file.path(rds_dir, "amphEBMs.RDS"))
#   plot_and_save(res, "Phase4_Clownfish_Combined", fig_dir)
#   all_metrics[["Phase4"]] <- extract_metrics(res, "Phase4_Clownfish_Combined")
# }

#-------------------------------------------------------------------------------
# PART 3: SAVE METRICS SUMMARY
#-------------------------------------------------------------------------------
if(length(all_metrics) > 0) {
  cat("\n--- Summarizing Metrics ---\n")
  
  final_stats <- bind_rows(all_metrics) %>%
    arrange(ModelType, Species)
  
  # Save CSV
  csv_path <- file.path(res_dir, "new_pipeline_metrics.csv")
  write_csv(final_stats, csv_path)
  
  cat("Success! Metrics saved to:", csv_path, "\n")
  print(head(final_stats))
  
} else {
  cat("No metrics extracted.\n")
}
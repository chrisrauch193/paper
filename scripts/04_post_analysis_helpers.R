# scripts/04_post_analysis_helpers.R
# ------------------------------------------------------------------------------
# POST-ANALYSIS HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# --- 1. CONFIG & THEMES ---
get_garcia_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30", margin = margin(b = 10)),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      axis.text = element_text(size = 10, color = "grey40")
    )
}

# --- 2. ROBUST FILE SAVING ---
save_and_print <- function(plot_obj, filename, output_dir, width = 10, height = 8) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  full_path <- file.path(output_dir, filename)
  ggsave(full_path, plot = plot_obj, width = width, height = height, dpi = 300, bg = "white")
  return(plot_obj)
}

# --- 3. SPECIES METADATA ---
get_specialist_status <- function(species_name) {
  # Define your specialists here based on your thesis
  specialists <- c("Amphiprion_frenatus", "Amphiprion_biaculeatus", "Amphiprion_nigripes", "Amphiprion_latezonatus") 
  if (species_name %in% specialists) return("Specialist")
  return("Generalist")
}

# --- 4. ROBUST STACKING (The Ant Paper Method) ---
stack_richness <- function(species_list, pred_dir, pattern_suffix = "_mean.tif") {
  
  valid_rasters <- list()
  
  for (sp in species_list) {
    # Clean name for file matching
    sp_clean <- gsub(" ", "_", sp)
    f_path <- file.path(pred_dir, paste0(sp_clean, pattern_suffix))
    
    if (file.exists(f_path)) {
      r <- terra::rast(f_path)
      valid_rasters[[sp]] <- r
    } else {
      warning(paste("Missing prediction for:", sp, "at", f_path))
    }
  }
  
  if (length(valid_rasters) == 0) return(NULL)
  
  # Stack and Sum
  r_stack <- terra::rast(valid_rasters)
  richness <- sum(r_stack, na.rm = TRUE)
  names(richness) <- "Richness"
  return(richness)
}

# --- 5. DELTA CALCULATION ---
calculate_delta <- function(curr, fut) {
  # Align geometry if slightly off due to different projections/saving
  if (!terra::compareGeom(curr, fut, stopOnError = FALSE)) {
    fut <- terra::resample(fut, curr, method = "bilinear")
  }
  delta <- fut - curr
  names(delta) <- "Delta_Suitability"
  return(delta)
}
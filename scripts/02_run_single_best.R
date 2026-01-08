# scripts/02_run_single_best.R
# Runs the full pipeline ONCE to generate:
# 1. Best Model Objects (for future projection)
# 2. Prediction Rasters (for plotting maps)
# 3. Host Maps (needed by the 03_run_iterations looper)

# --- Setup ---
if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, stringr, ENMeval)

# Load shared functions
source("scripts/01_model_functions.R")

# --- Config ---
OUT_DIR_MOD <- "outputs/models/single_run"
OUT_DIR_PRED <- "outputs/predictions/current"
dir.create(OUT_DIR_MOD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_PRED, recursive = TRUE, showWarnings = FALSE)

# We also keep the 'Rdata' folder for backward compatibility with your looper logic
dir.create("Rdata", showWarnings = FALSE) 

# --- Load Data ---
cat("Loading Data...\n")
env_df <- read_csv("data/selected_environmental_variables.csv", show_col_types = FALSE)
env_rast <- terra::rast(env_df, type = "xyz", crs = "EPSG:4326")

amph_occ <- read_csv("data/amph_occ_env_final_dataset.csv", show_col_types = FALSE)
anem_occ <- read_csv("data/anem_occ_env_final_dataset.csv", show_col_types = FALSE)

int_mat <- read.csv("data/interaction_matrix.csv", row.names = 1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

#===============================================================================
# PHASE 1: Host Anemones (Env Only)
# This is critical because Phase 3 & 4 depend on these maps.
#===============================================================================
cat("\n--- Running PHASE 1: Host Anemones ---\n")
anem_species <- unique(anem_occ$species)
anem_maps <- terra::rast()
anem_eval <- data.frame()

for (sp in anem_species) {
  cat("Processing:", sp, "...\n")
  res <- run_species_model(sp, anem_occ, env_rast, "env", verbose = FALSE)
  
  if (!is.null(res)) {
    # Generate Prediction Raster
    # Note: run_species_model returns the model object, we must predict here to save TIF
    best_mod <- res$model@models[[as.numeric(rownames(res$stats))]]
    pred_map <- enm.maxnet@predict(best_mod, env_rast, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
    names(pred_map) <- gsub(" ", "_", sp)
    
    # Add to stack and save
    anem_maps <- c(anem_maps, pred_map)
    anem_eval <- bind_rows(anem_eval, res$stats %>% mutate(species = sp))
    
    # Save individual raster
    terra::writeRaster(pred_map, file.path(OUT_DIR_PRED, paste0("anem_", names(pred_map), ".tif")), overwrite = TRUE)
  }
}

# SAVE PHASE 1 RESULTS (Critical for Looper)
# We save wrapped maps to avoid pointer dead ends
anem_results_packed <- list(maps = terra::wrap(anem_maps), eval = anem_eval)
saveRDS(anem_results_packed, file.path(OUT_DIR_MOD, "anemENMs.rds"))
# Save copy to Rdata/ for the looper script to find easily
saveRDS(anem_results_packed, "Rdata/anemENMs.RDS") 
cat("Phase 1 Complete. Host maps saved.\n")

#===============================================================================
# PHASE 2: Fish (Env Only)
#===============================================================================
cat("\n--- Running PHASE 2: Fish (Env Only) ---\n")
fish_species <- unique(amph_occ$species)
fish_env_maps <- terra::rast()
fish_env_eval <- data.frame()

for (sp in fish_species) {
  cat("Processing:", sp, "...\n")
  res <- run_species_model(sp, amph_occ, env_rast, "env", verbose = FALSE)
  
  if (!is.null(res)) {
    best_mod <- res$model@models[[as.numeric(rownames(res$stats))]]
    pred_map <- enm.maxnet@predict(best_mod, env_rast, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
    names(pred_map) <- gsub(" ", "_", sp)
    
    fish_env_maps <- c(fish_env_maps, pred_map)
    fish_env_eval <- bind_rows(fish_env_eval, res$stats %>% mutate(species = sp))
    
    terra::writeRaster(pred_map, file.path(OUT_DIR_PRED, paste0("fish_env_", names(pred_map), ".tif")), overwrite = TRUE)
  }
}

saveRDS(list(maps = terra::wrap(fish_env_maps), eval = fish_env_eval), file.path(OUT_DIR_MOD, "amphENMs.rds"))
cat("Phase 2 Complete.\n")

#===============================================================================
# PHASE 3: Fish (Biotic Only)
#===============================================================================
cat("\n--- Running PHASE 3: Fish (Biotic Only) ---\n")
fish_bio_maps <- terra::rast()
fish_bio_eval <- data.frame()

# Use the UNWRAPPED maps from Phase 1
# (anem_maps is already a live SpatRaster in this session, no need to unwrap)
# Ensure names match matrix
names(anem_maps) <- gsub(" ", "_", names(anem_maps))
names(anem_maps) <- gsub("\\.", "_", names(anem_maps))

for (sp in fish_species) {
  cat("Processing:", sp, "...\n")
  
  # Helper function creates the layer
  bio_layer <- get_biotic_layer(sp, int_mat, anem_maps)
  
  if (!is.null(bio_layer)) {
    res <- run_species_model(sp, amph_occ, env_rast, "biotic_only", host_map_stack = bio_layer, verbose = FALSE)
    
    if (!is.null(res)) {
      best_mod <- res$model@models[[as.numeric(rownames(res$stats))]]
      pred_map <- enm.maxnet@predict(best_mod, bio_layer, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
      names(pred_map) <- gsub(" ", "_", sp)
      
      fish_bio_maps <- c(fish_bio_maps, pred_map)
      fish_bio_eval <- bind_rows(fish_bio_eval, res$stats %>% mutate(species = sp))
      
      terra::writeRaster(pred_map, file.path(OUT_DIR_PRED, paste0("fish_bio_", names(pred_map), ".tif")), overwrite = TRUE)
    }
  } else {
    cat("  Skipping (No host map found)\n")
  }
}

saveRDS(list(maps = terra::wrap(fish_bio_maps), eval = fish_bio_eval), file.path(OUT_DIR_MOD, "amphBioticOnly.rds"))
cat("Phase 3 Complete.\n")

#===============================================================================
# PHASE 4: Fish (Combined)
#===============================================================================
cat("\n--- Running PHASE 4: Fish (Combined) ---\n")
fish_cmb_maps <- terra::rast()
fish_cmb_eval <- data.frame()

for (sp in fish_species) {
  cat("Processing:", sp, "...\n")
  
  bio_layer <- get_biotic_layer(sp, int_mat, anem_maps)
  
  if (!is.null(bio_layer)) {
    # run_species_model handles the stacking internally if type="combined"
    res <- run_species_model(sp, amph_occ, env_rast, "combined", host_map_stack = bio_layer, verbose = FALSE)
    
    if (!is.null(res)) {
      # Re-stack for prediction
      full_stack <- c(env_rast, bio_layer)
      best_mod <- res$model@models[[as.numeric(rownames(res$stats))]]
      pred_map <- enm.maxnet@predict(best_mod, full_stack, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
      names(pred_map) <- gsub(" ", "_", sp)
      
      fish_cmb_maps <- c(fish_cmb_maps, pred_map)
      fish_cmb_eval <- bind_rows(fish_cmb_eval, res$stats %>% mutate(species = sp))
      
      terra::writeRaster(pred_map, file.path(OUT_DIR_PRED, paste0("fish_combined_", names(pred_map), ".tif")), overwrite = TRUE)
    }
  } else {
    cat("  Skipping (No host map found)\n")
  }
}

saveRDS(list(maps = terra::wrap(fish_cmb_maps), eval = fish_cmb_eval), file.path(OUT_DIR_MOD, "amphEBMs.rds"))
# Copy for legacy scripts if needed
saveRDS(list(maps = terra::wrap(fish_cmb_maps), eval = fish_cmb_eval), "Rdata/amphEBMs.RDS")

cat("\n--- SINGLE RUN COMPLETE ---\n")
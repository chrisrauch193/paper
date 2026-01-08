# scripts/new_pipeline/01_enmeval_models.R
#-------------------------------------------------------------------------------
# REPRODUCTION OF JIMENEZ PIPELINE USING ENMEVAL
# 1. Host SDMs (Base)
# 2. Fish SDMs (Base Env)
# 3. Fish SDMs (Biotic / Combined)
#-------------------------------------------------------------------------------

# --- 1. Setup ---
if(!require("pacman")) install.packages("pacman")
pacman::p_load(ENMeval, terra, dplyr, readr, sf, tools, stringr)

# Paths
setwd("/home/bi-server-kyoto/a0236995/paper/") # Your project root
dir.create("./Rdata", showWarnings = FALSE)
dir.create("./results", showWarnings = FALSE)

# --- 2. Load Data ---
cat("Loading Data...\n")

# Environment (Convert CSV -> RasterStack)
# Assumes the CSV is a regular grid
env_df <- read_csv("./data/selected_environmental_variables.csv", show_col_types = FALSE)
env_rast <- terra::rast(env_df, type = "xyz", crs = "EPSG:4326")

# Occurrences
anem_occ <- read_csv("./data/anem_occ_env_final_dataset.csv", show_col_types = FALSE)
amph_occ <- read_csv("./data/amph_occ_env_final_dataset.csv", show_col_types = FALSE)

# Interaction Matrix
int_mat <- read.csv("./data/interaction_matrix.csv", row.names = 1)

# --- 3. Define ENMeval Wrapper Function (FIXED ARGUMENTS) ---
run_enmeval_batch <- function(species_list, occ_df, env_stack, bg_n = 10000) {
  
  results_list <- list()
  maps_stack <- terra::rast()
  eval_df <- data.frame()
  
  for (sp in species_list) {
    cat("\nProcessing:", sp, "...\n")
    
    # 1. Filter Occs
    sp_occ <- occ_df %>% filter(species == sp) %>% dplyr::select(x, y)
    
    # Ensure coordinates are numeric
    sp_occ$x <- as.numeric(sp_occ$x)
    sp_occ$y <- as.numeric(sp_occ$y)
    
    if (nrow(sp_occ) < 5) {
      cat("  Skipping: Too few occurrences (<5)\n")
      next
    }
    
    # 2. Background Points
    # Sample from env_stack (returns dataframe with cell values, we just need coords)
    # Using 'xy=TRUE' returns x, y columns.
    bg_sample <- terra::spatSample(env_stack, size = bg_n, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
    bg <- bg_sample[, c("x", "y")] # Ensure just coords
    
    # 3. Run ENMeval
    tryCatch({
      # FIXED ARGUMENT NAMES HERE: occs, envs, bg
      mod <- ENMevaluate(occs = sp_occ, 
                         envs = env_stack, 
                         bg = bg, 
                         algorithm = 'maxnet', 
                         partitions = 'block', 
                         tune.args = list(fc = c("L", "LQ", "H"), rm = 1:2),
                         parallel = TRUE, numCores = 4, quiet = TRUE)
      
      # 4. Select Best Model (Lowest AICc)
      best_res <- mod@results %>% filter(delta.AICc == 0) %>% slice(1)
      if(nrow(best_res) == 0) best_res <- mod@results %>% arrange(desc(auc.val.avg)) %>% slice(1)
      
      # 5. Predict
      # ENMeval 2.0 stores predictions for the best model automatically? No, we must generate.
      # But first, get the best model object
      best_model_obj <- mod@models[[as.numeric(rownames(best_res))]]
      
      # Generate prediction raster
      # Use enm.maxnet@predict directly or terra::predict
      # Note: 'envs' must be a RasterStack/SpatRaster
      pred_map <- enm.maxnet@predict(best_model_obj, env_stack, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
      names(pred_map) <- sp
      
      # 6. Store Results
      maps_stack <- c(maps_stack, pred_map)
      
      stats <- data.frame(
        Species = sp,
        AUC_mean = best_res$auc.val.avg,
        CBI_mean = best_res$cbi.val.avg,
        Settings = paste0("fc:", best_res$fc, "_rm:", best_res$rm)
      )
      eval_df <- bind_rows(eval_df, stats)
      
      results_list[[sp]] <- list(model = mod, best_settings = best_res)
      
      cat("  Success. AUC:", round(best_res$auc.val.avg, 3), "\n")
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
    })
  }
  
  return(list(maps = terra::wrap(maps_stack), eval = eval_df, details = results_list))
}

#===============================================================================
# PHASE 1: Sea Anemones (anemENM)
#===============================================================================
cat("\n--- Running PHASE 1: Host Anemones ---\n")
anem_species <- unique(anem_occ$species)

anemENM <- run_enmeval_batch(anem_species, anem_occ, env_rast)

# Save
saveRDS(anemENM, "./Rdata/anemENMs.RDS")
cat("Anemone models saved.\n")

#===============================================================================
# PHASE 2: Clownfish Env Only (amphENM)
#===============================================================================
cat("\n--- Running PHASE 2: Clownfish (Env Only) ---\n")
amph_species <- unique(amph_occ$species)

amphENM <- run_enmeval_batch(amph_species, amph_occ, env_rast)

# Save
saveRDS(amphENM, "./Rdata/amphENMs.RDS")
cat("Clownfish Env models saved.\n")

#===============================================================================
# PHASE 3: Clownfish HOST ONLY (Biotic Only)
# Hypothesis: "Is the fish distribution driven PURELY by host availability?"
# Input: ONLY the Max Suitability of Hosts. No Temp, No Salinity.
#===============================================================================
cat("\n--- Running PHASE 3: Clownfish (Host Only) ---\n")

# --- STEP A: PREPARE HOST MAPS ---
# We must unwrap the Phase 1 maps to use them for calculation
cat("Preparing Host Maps...\n")
if (inherits(anemENM$maps, "PackedSpatRaster")) {
  host_model_stack <- terra::unwrap(anemENM$maps)
} else {
  host_model_stack <- anemENM$maps
}

# Force map names to match the interaction matrix format (Underscores)
names(host_model_stack) <- gsub(" ", "_", names(host_model_stack))
names(host_model_stack) <- gsub("\\.", "_", names(host_model_stack))
cat("Host Maps Ready:", length(names(host_model_stack)), "species available.\n")

# Initialize Containers
biotic_maps <- terra::rast()
biotic_eval <- data.frame()
amph_species <- unique(amph_occ$species)

# --- STEP B: LOOP THROUGH FISH ---
for (fish in amph_species) {
  cat("\nProcessing Host-Only:", fish, "...\n")
  
  # 1. Identify Hosts
  # Clean fish name to match matrix rownames (Underscores)
  fish_clean <- gsub(" ", "_", fish)
  
  # Check existence in matrix
  if (!fish_clean %in% rownames(int_mat)) { 
    cat("  Skipping: Fish name '", fish_clean, "' not found in interaction matrix.\n", sep="")
    next 
  }
  
  # Get hosts with ANY association (> 0)
  hosts <- names(int_mat)[int_mat[fish_clean, ] > 0]
  
  # Keep only hosts we actually modeled in Phase 1
  valid_hosts <- hosts[hosts %in% names(host_model_stack)]
  
  if (length(valid_hosts) == 0) { 
    cat("  Skipping: No modeled hosts found. (Required:", paste(hosts, collapse=", "), ")\n")
    next 
  }
  
  cat("  Hosts found:", paste(valid_hosts, collapse=", "), "\n")
  
  # 2. Create Biotic Layer (Max Suitability)
  # Logic: "If the best host is 80% suitable here, the site is 80% suitable biotically."
  if (length(valid_hosts) == 1) {
    biotic_layer <- host_model_stack[[valid_hosts]]
  } else {
    sub_stack <- host_model_stack[[valid_hosts]]
    biotic_layer <- terra::app(sub_stack, fun = max, na.rm = TRUE)
  }
  names(biotic_layer) <- "host_suitability"
  
  # 3. Stack (ONLY BIOTIC)
  env_input <- biotic_layer 
  
  # 4. Prepare Occurrences & Background
  sp_occ <- amph_occ %>% filter(species == fish) %>% dplyr::select(x, y)
  sp_occ$x <- as.numeric(sp_occ$x); sp_occ$y <- as.numeric(sp_occ$y)
  
  if (nrow(sp_occ) < 5) {
    cat("  Skipping: Too few occurrences (<5)\n")
    next
  }
  
  bg <- terra::spatSample(env_input, size = 10000, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
  bg <- bg[, c("x", "y")]
  
  # 5. Run ENMeval
  tryCatch({
    mod <- ENMevaluate(occs = sp_occ, envs = env_input, bg = bg, 
                       algorithm = 'maxnet', partitions = 'block', 
                       tune.args = list(fc = c("L", "LQ"), rm = 1:2), # Simpler tuning for single variable
                       parallel = TRUE, numCores = 4, quiet = TRUE)
    
    # Select Best Model
    best_res <- mod@results %>% filter(delta.AICc == 0) %>% slice(1)
    if(nrow(best_res) == 0) best_res <- mod@results %>% arrange(desc(auc.val.avg)) %>% slice(1)
    best_model <- mod@models[[as.numeric(rownames(best_res))]]
    
    # Predict
    pred_map <- enm.maxnet@predict(best_model, env_input, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
    names(pred_map) <- fish
    biotic_maps <- c(biotic_maps, pred_map)
    
    # Save Stats
    stats <- data.frame(
      Species = fish, 
      AUC_mean = best_res$auc.val.avg, 
      CBI_mean = best_res$cbi.val.avg, 
      Settings = paste0("fc:", best_res$fc, "_rm:", best_res$rm)
    )
    biotic_eval <- bind_rows(biotic_eval, stats)
    cat("  Success. AUC:", round(best_res$auc.val.avg, 3), "\n")
    
  }, error = function(e) { cat("  ERROR:", e$message, "\n") })
}

# Save Phase 3 (Wrapped)
amphBioticOnly <- list(maps = terra::wrap(biotic_maps), eval = biotic_eval)
saveRDS(amphBioticOnly, "./Rdata/amphBioticOnly.RDS")
cat("Host-Only models saved.\n")


#===============================================================================
# PHASE 4: Clownfish COMBINED (Env + Host)
# Hypothesis: "Realized Niche = Abiotic Niche + Host Availability"
# Input: Environment Stack + Host Suitability Layer
#===============================================================================
cat("\n--- Running PHASE 4: Clownfish Combined (Biotic + Env) ---\n")

combined_maps <- terra::rast()
combined_eval <- data.frame()

# We reuse 'host_model_stack' from Phase 3, so we don't need to unwrap again.

for (fish in amph_species) {
  cat("\nProcessing Combined:", fish, "...\n")
  
  # 1. Identify Hosts (Same logic as above)
  fish_clean <- gsub(" ", "_", fish)
  
  if (!fish_clean %in% rownames(int_mat)) { next }
  
  hosts <- names(int_mat)[int_mat[fish_clean, ] > 0]
  valid_hosts <- hosts[hosts %in% names(host_model_stack)]
  
  if (length(valid_hosts) == 0) { next }
  
  # 2. Create Biotic Layer
  if (length(valid_hosts) == 1) {
    biotic_layer <- host_model_stack[[valid_hosts]]
  } else {
    sub_stack <- host_model_stack[[valid_hosts]]
    biotic_layer <- terra::app(sub_stack, fun = max, na.rm = TRUE)
  }
  names(biotic_layer) <- "host_suitability"
  
  # 3. Stack (COMBINED: Env + Biotic)
  # env_rast is Global, biotic_layer is Global. They stack perfectly.
  env_input <- c(env_rast, biotic_layer)
  
  # 4. Prepare Occurrences & Background
  sp_occ <- amph_occ %>% filter(species == fish) %>% dplyr::select(x, y)
  sp_occ$x <- as.numeric(sp_occ$x); sp_occ$y <- as.numeric(sp_occ$y)
  
  if (nrow(sp_occ) < 5) { next }
  
  bg <- terra::spatSample(env_input, size = 10000, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
  bg <- bg[, c("x", "y")]
  
  # 5. Run ENMeval
  tryCatch({
    mod <- ENMevaluate(occs = sp_occ, envs = env_input, bg = bg, 
                       algorithm = 'maxnet', partitions = 'block', 
                       tune.args = list(fc = c("L", "LQ", "H"), rm = 1:2),
                       parallel = TRUE, numCores = 4, quiet = TRUE)
    
    # Select Best & Predict
    best_res <- mod@results %>% filter(delta.AICc == 0) %>% slice(1)
    if(nrow(best_res) == 0) best_res <- mod@results %>% arrange(desc(auc.val.avg)) %>% slice(1)
    best_model <- mod@models[[as.numeric(rownames(best_res))]]
    
    pred_map <- enm.maxnet@predict(best_model, env_input, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
    names(pred_map) <- fish
    combined_maps <- c(combined_maps, pred_map)
    
    # Save Stats
    stats <- data.frame(
      Species = fish, 
      AUC_mean = best_res$auc.val.avg, 
      CBI_mean = best_res$cbi.val.avg, 
      Settings = paste0("fc:", best_res$fc, "_rm:", best_res$rm)
    )
    combined_eval <- bind_rows(combined_eval, stats)
    cat("  Success. AUC:", round(best_res$auc.val.avg, 3), "\n")
    
  }, error = function(e) { cat("  ERROR:", e$message, "\n") })
}

# Save Phase 4 (Wrapped)
amphEBM <- list(maps = terra::wrap(combined_maps), eval = combined_eval)
saveRDS(amphEBM, "./Rdata/amphEBMs.RDS")
cat("Combined models saved.\n")

cat("\n--- PIPELINE COMPLETE ---\n")
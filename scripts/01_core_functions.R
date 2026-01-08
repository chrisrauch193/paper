# scripts/01_core_functions.R
# Core logic for SDM Pipeline (Test & Production)

library(terra)
library(ENMeval)
library(dplyr)
library(readr)

# --- 1. Worker wrapper for training one iteration ---
train_model_worker <- function(species, occ_df, env_stack, seed, verbose=FALSE) {
  # Set seed for reproducibility of background points/folds
  set.seed(seed)
  
  sp_occ <- occ_df %>% filter(species == !!species) %>% dplyr::select(x, y)
  
  # Safety check
  if(nrow(sp_occ) < 5) return(NULL)
  
  # Sample background (10k)
  bg <- terra::spatSample(env_stack, size = 10000, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
  bg <- bg[, c("x", "y")]
  
  tryCatch({
    # ENMeval settings (Block partition for rigorous testing)
    mod <- ENMevaluate(occs = sp_occ, envs = env_stack, bg = bg, 
                       algorithm = 'maxnet', partitions = 'block', 
                       tune.args = list(fc = c("L", "LQ", "H"), rm = 1:2),
                       parallel = FALSE, quiet = !verbose)
    
    # Extract Best Model (Lowest AICc)
    best_res <- mod@results %>% filter(delta.AICc == 0) %>% slice(1)
    if(nrow(best_res) == 0) best_res <- mod@results %>% arrange(desc(auc.val.avg)) %>% slice(1)
    
    # Return compact object to save memory
    return(list(
      model_obj = mod@models[[as.numeric(rownames(best_res))]], # The maxnet object
      stats = best_res,                                         # The metrics
      best_settings = best_res$tune.args
    ))
    
  }, error = function(e) {
    return(list(error = e$message))
  })
}

# --- 2. Helper to construct the Biotic Layer ---
get_biotic_layer <- function(fish_name, int_mat, host_map_stack) {
  fish_clean <- gsub(" ", "_", fish_name)
  
  if (!fish_clean %in% rownames(int_mat)) return(NULL)
  
  # Find associated hosts with interaction > 0
  hosts <- names(int_mat)[int_mat[fish_clean, ] > 0]
  
  # Intersect with hosts we actually have maps for
  valid_hosts <- hosts[hosts %in% names(host_map_stack)]
  
  if (length(valid_hosts) == 0) return(NULL)
  
  # Create Biotic Layer (Max Suitability logic)
  subset_stack <- host_map_stack[[valid_hosts]]
  
  if (terra::nlyr(subset_stack) > 1) {
    # If multiple hosts, take the MAX value per pixel
    bio_layer <- terra::app(subset_stack, max, na.rm = TRUE)
  } else {
    bio_layer <- subset_stack
  }
  
  names(bio_layer) <- "host_suitability"
  return(bio_layer)
}

# --- 3. Logging Helper ---
write_log <- function(path, msg) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n"), 
      file = path, append = TRUE)
}
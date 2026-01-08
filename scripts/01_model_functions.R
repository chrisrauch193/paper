# scripts/01_model_functions.R
# Core functions for ENMeval pipeline
# Loaded by 02 and 03 scripts

library(ENMeval)
library(terra)
library(dplyr)

# --- 1. Universal Wrapper for One Species ---
run_species_model <- function(species_name, occ_df, env_stack, 
                              model_type = c("env", "biotic_only", "combined"), 
                              seed = NULL, bg_n = 10000, 
                              host_map_stack = NULL, verbose = FALSE) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # Filter Occurrences
  sp_occ <- occ_df %>% 
    filter(species == species_name) %>% 
    dplyr::select(x, y) %>%
    mutate(x = as.numeric(x), y = as.numeric(y))
  
  if (nrow(sp_occ) < 5) return(NULL) # Skip rare species
  
  # Define Environment based on Model Type
  if (model_type == "env") {
    mod_env <- env_stack
  } else if (model_type == "biotic_only") {
    if (is.null(host_map_stack)) stop("Host stack missing for biotic model")
    mod_env <- host_map_stack # Expecting pre-calculated max suitability layer
  } else if (model_type == "combined") {
    if (is.null(host_map_stack)) stop("Host stack missing for combined model")
    mod_env <- c(env_stack, host_map_stack)
  }
  
  # Sample Background
  bg <- terra::spatSample(mod_env, size = bg_n, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
  bg <- bg[, c("x", "y")]
  
  # Run ENMeval
  # Note: using 'random' partition for iterations is often safer for statistical 
  # comparison of AUC unless you strictly control block assignment across folds.
  # But 'block' is better for spatial independence. Let's stick to 'block'.
  
  tryCatch({
    mod <- ENMevaluate(occs = sp_occ, envs = mod_env, bg = bg, 
                       algorithm = 'maxnet', partitions = 'block', 
                       tune.args = list(fc = c("L", "LQ", "H"), rm = 1:2),
                       parallel = FALSE, # We parallelize the LOOP, not the internal tuning
                       quiet = !verbose)
    
    # Select Best Model (AICc)
    best_res <- mod@results %>% filter(delta.AICc == 0) %>% slice(1)
    if(nrow(best_res) == 0) best_res <- mod@results %>% arrange(desc(auc.val.avg)) %>% slice(1)
    
    return(list(model = mod, stats = best_res, settings = best_res$tune.args))
    
  }, error = function(e) {
    if(verbose) cat("Error in", species_name, ":", e$message, "\n")
    return(NULL)
  })
}

# --- 2. Helper to Prepare Biotic Layer ---
get_biotic_layer <- function(fish_name, int_mat, host_maps) {
  fish_clean <- gsub(" ", "_", fish_name)
  if (!fish_clean %in% rownames(int_mat)) return(NULL)
  
  hosts <- names(int_mat)[int_mat[fish_clean, ] > 0]
  valid_hosts <- hosts[hosts %in% names(host_maps)]
  
  if (length(valid_hosts) == 0) return(NULL)
  
  if (length(valid_hosts) == 1) {
    layer <- host_maps[[valid_hosts]]
  } else {
    sub_stack <- host_maps[[valid_hosts]]
    layer <- terra::app(sub_stack, fun = max, na.rm = TRUE)
  }
  names(layer) <- "host_suitability"
  return(layer)
}
# scripts/02_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL "JOURNAL STANDARD" SDM PIPELINE RUNNER
# ------------------------------------------------------------------------------

# --- 1. CONFIGURATION ---------------------------------------------------------
RUN_MODE <- "TEST" 

STATIC_VARS <- c("rugosity", "bathymetry", "slope", "aspect") 

if (RUN_MODE == "FINAL") {
  N_CORES     <- 30 
  N_HOST_BOOT <- 10 
  N_FISH_BOOT <- 40 
  OUTPUT_ROOT <- "outputs/final_run"
} else {
  N_CORES     <- 30
  N_HOST_BOOT <- 2
  N_FISH_BOOT <- 2
  OUTPUT_ROOT <- "outputs/test_run"
}

# --- 2. SETUP -----------------------------------------------------------------
if(!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, terra, dplyr, readr, ENMeval, maxnet, dismo)

source("scripts/01_core_functions.R")

# Directories
DIRS <- c("models", 
          "predictions/current/hosts", 
          "predictions/current/fish_env", 
          "predictions/current/fish_combined",
          "logs") # NEW: Specific log directory

for(d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive=TRUE, showWarnings=FALSE)

MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
file.create(MASTER_LOG)
write_log(MASTER_LOG, paste("--- PIPELINE STARTED IN", RUN_MODE, "MODE ---"))

# --- 3. DATA LOADING ----------------------------------------------------------

# A. Load Environmental Data
ENV_PATH <- "data/selected_environmental_variables.csv"
cat("Loading Current Environmental Data...\n")

env_df <- read_csv(ENV_PATH, show_col_types=FALSE)
current_env <- terra::rast(env_df, type="xyz", crs="EPSG:4326")

packed_env <- terra::wrap(current_env)

# B. Load Future Scenarios
FUT_DIR   <- "data/env/future_pca" 
FUT_FILES <- list.files(FUT_DIR, full=TRUE, pattern="\\.tif$")
scenario_names <- tools::file_path_sans_ext(basename(FUT_FILES))

for(scen in scenario_names) {
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"), recursive=TRUE, showWarnings=FALSE)
}

# C. Load Occurrences (Force DataFrame)
amph_occ <- as.data.frame(read_csv("data/amph_occ_env_final_dataset.csv", show_col_types=FALSE))
anem_occ <- as.data.frame(read_csv("data/anem_occ_env_final_dataset.csv", show_col_types=FALSE))
int_mat  <- read.csv("data/interaction_matrix.csv", row.names=1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

# --- 4. CLUSTER INIT ----------------------------------------------------------
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

clusterExport(cl, c("get_best_params", "fit_bootstrap_worker", "get_biotic_layer", "write_log", 
                    "packed_env", "FUT_FILES", "scenario_names", "STATIC_VARS",
                    "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG", 
                    "N_HOST_BOOT", "N_FISH_BOOT"))

clusterEvalQ(cl, { 
  library(terra); library(dplyr); library(maxnet); library(ENMeval); library(readr)
  
  prepare_future_stack <- function(fut_file, curr_stack, static_names) {
    fut_rast <- terra::rast(fut_file)
    wanted_static <- intersect(names(curr_stack), static_names)
    if(length(wanted_static) > 0) {
      missing_static <- setdiff(wanted_static, names(fut_rast))
      if(length(missing_static) > 0) {
        static_layers <- curr_stack[[missing_static]]
        if(!terra::compareGeom(fut_rast, static_layers, stopOnError=FALSE)) {
          static_layers <- terra::resample(static_layers, fut_rast, method="near")
        }
        return(c(fut_rast, static_layers))
      }
    }
    return(fut_rast)
  }
})

# ==============================================================================
# PHASE 1: HOSTS
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 1: HOSTS STARTED ---")
cat("--- PHASE 1: HOSTS ---\n")

anem_species <- unique(anem_occ$species)

host_results <- foreach(sp = anem_species, .packages=c('terra','dplyr','maxnet')) %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  write_log(MASTER_LOG, paste("START Host:", sp))
  
  # Define Species Log File
  sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
  file.create(sp_log)
  
  env_stack <- terra::unwrap(packed_env)
  
  tryCatch({
    sp_dat <- anem_occ %>% dplyr::filter(species == sp)
    if(nrow(sp_dat) < 5) stop("Not enough occurrences")
    
    future_stacks <- list()
    for(i in seq_along(FUT_FILES)) {
      scen <- scenario_names[i]
      f_stack <- prepare_future_stack(FUT_FILES[i], env_stack, STATIC_VARS)
      vars_needed <- names(env_stack)
      if(all(vars_needed %in% names(f_stack))) {
        future_stacks[[scen]] <- f_stack[[vars_needed]]
      }
    }
    
    params <- get_best_params(sp_dat, env_stack)
    if(is.null(params)) stop("Tuning failed")
    
    # PASS LOG FILE HERE
    results <- fit_bootstrap_worker(sp_dat, env_stack, future_stacks, params, n_boot=N_HOST_BOOT, debug_log=sp_log)
    
    terra::writeRaster(results$current, 
                       file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, ".tif")), 
                       overwrite=TRUE)
    
    for(scen in names(results$future)) {
      terra::writeRaster(results$future[[scen]], 
                         file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, ".tif")), 
                         overwrite=TRUE)
    }
    
    write_log(MASTER_LOG, paste("FINISH Host:", sp))
    return(paste("SUCCESS:", sp))
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Host:", sp, "->", e$message))
    return(paste("FAILED:", sp, "-", e$message))
  })
}

print(host_results)

# ==============================================================================
# PHASE 2: FISH
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 2: FISH STARTED ---")
cat("--- PHASE 2: FISH ---\n")

host_files <- list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=TRUE, pattern="\\.tif$")
if(length(host_files) == 0) stop("CRITICAL: Phase 1 failed (No Host Maps).")

packed_hosts_curr <- terra::wrap(terra::rast(host_files))
clusterExport(cl, "packed_hosts_curr")

fish_species <- unique(amph_occ$species)

fish_results <- foreach(sp = fish_species, .packages=c('terra','dplyr','maxnet')) %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  write_log(MASTER_LOG, paste("START Fish:", sp))
  
  # Species Log
  sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
  file.create(sp_log)
  
  env_stack_curr  <- terra::unwrap(packed_env)
  host_stack_curr <- terra::unwrap(packed_hosts_curr)
  
  tryCatch({
    sp_dat <- amph_occ %>% dplyr::filter(species == sp)
    if(nrow(sp_dat) < 5) stop("Not enough occurrences")
    
    # A. Current
    biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat)
    if(is.null(biotic_curr)) stop("No matching hosts")
    
    full_stack_curr <- c(env_stack_curr, biotic_curr)
    names(full_stack_curr) <- c(names(env_stack_curr), "biotic_suitability")
    
    # B. Futures
    future_stacks <- list()
    for(i in seq_along(FUT_FILES)) {
      scen <- scenario_names[i]
      f_env <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
      f_env <- f_env[[names(env_stack_curr)]] 
      
      host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
      host_files_fut <- list.files(host_dir_fut, full.names=TRUE, pattern="\\.tif$")
      
      if(length(host_files_fut) > 0) {
        host_stack_fut <- terra::rast(host_files_fut)
        biotic_fut <- get_biotic_layer(sp, host_stack_fut, int_mat)
        if(!is.null(biotic_fut)) {
          f_comb <- c(f_env, biotic_fut)
          names(f_comb) <- c(names(f_env), "biotic_suitability")
          future_stacks[[scen]] <- f_comb
        }
      }
    }
    
    params <- get_best_params(sp_dat, full_stack_curr)
    if(is.null(params)) stop("Tuning failed")
    
    # PASS LOG FILE HERE
    results <- fit_bootstrap_worker(sp_dat, full_stack_curr, future_stacks, params, n_boot=N_FISH_BOOT, debug_log=sp_log)
    
    terra::writeRaster(results$current, 
                       file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, ".tif")), 
                       overwrite=TRUE)
    
    for(scen in names(results$future)) {
      terra::writeRaster(results$future[[scen]], 
                         file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, ".tif")), 
                         overwrite=TRUE)
    }
    
    write_log(MASTER_LOG, paste("FINISH Fish:", sp))
    return(paste("SUCCESS:", sp))
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish:", sp, "->", e$message))
    return(paste("FAILED:", sp, "-", e$message))
  })
}

print(fish_results)
stopCluster(cl)
write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
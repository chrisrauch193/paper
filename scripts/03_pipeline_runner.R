# scripts/03_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL SDM PIPELINE: EXECUTION
# ------------------------------------------------------------------------------

# 1. SETUP
rm(list=ls())
gc()
if(!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, terra, dplyr, readr, ENMeval, maxnet, dismo)

# 2. LOAD MODULES
source("scripts/00_config.R")
source("scripts/01_functions_core.R")
source("scripts/02_functions_model.R")

# 3. DIRECTORIES
DIRS <- c("models_stats", "models_tuning", "models", # New 'models' folder for iteration RDS
          "predictions/current/hosts", "predictions/current/fish_combined", "predictions/current/fish_env_only", "logs")
for(d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive=TRUE, showWarnings=FALSE)

# 4. AUTO-WIPE
if(WIPE_PREDICTIONS) {
  try(unlink(file.path(OUTPUT_ROOT, "models_stats", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "models", "*")), silent=T) # Wipe iteration stats
  try(unlink(file.path(OUTPUT_ROOT, "models_tuning", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/hosts", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/fish_combined", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", "*")), silent=T)
}

# 5. LOGGING
MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
if(file.exists(MASTER_LOG)) unlink(MASTER_LOG)
file.create(MASTER_LOG)
write_log(MASTER_LOG, paste("--- PIPELINE START | ID:", RUN_ID, "| BG:", BG_SAMPLING_METHOD, "---"))

# 6. DATA
cat("Loading Data...\n")
if(!file.exists(ENV_PATH)) stop("Env Raster missing.")
current_env <- terra::rast(ENV_PATH)
packed_env  <- terra::wrap(current_env) 

FUT_FILES <- list.files(FUT_DIR, full=TRUE, pattern="\\.tif$")
scenario_names <- tools::file_path_sans_ext(basename(FUT_FILES))
for(scen in scenario_names) {
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"), recursive=T, showWarnings=F)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"), recursive=T, showWarnings=F)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only"), recursive=T, showWarnings=F)
}

amph_occ <- as.data.frame(read_csv("data/amph_occ_env_final_dataset.csv", show_col_types=FALSE))
anem_occ <- as.data.frame(read_csv("data/anem_occ_env_final_dataset.csv", show_col_types=FALSE))
int_mat  <- read.csv("data/interaction_matrix.csv", row.names=1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

# 7. CLUSTER
try(stopCluster(cl), silent=TRUE)
cl <- makeCluster(N_CORES)
registerDoParallel(cl)
clusterEvalQ(cl, { rm(list=ls()); gc() })

# !!! CRITICAL FIX: EXPLICIT EXPORT !!!
clusterExport(cl, c("packed_env", "FUT_FILES", "scenario_names", "STATIC_VARS",
                    "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG", 
                    "N_HOST_BOOT", "N_FISH_BOOT", "TUNE_ARGS", 
                    "USE_SPATIAL_THINNING", "USE_SPATIAL_TUNING", "BG_SAMPLING_METHOD",
                    # Explicitly list functions to ensure they cross over
                    "write_log", "thin_occurrences", "get_bias_corrected_background", 
                    "get_random_background", "get_biotic_layer", "get_best_params", "fit_bootstrap_worker"),
              envir = environment()) 

clusterEvalQ(cl, { 
  library(terra); library(dplyr); library(maxnet); library(ENMeval); library(readr); library(dismo)
  
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
write_log(MASTER_LOG, "--- PHASE 1: HOSTS ---")
if(!is.null(TARGET_HOSTS)) anem_run_list <- intersect(unique(anem_occ$species), TARGET_HOSTS) else anem_run_list <- unique(anem_occ$species)

foreach(sp = anem_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  # Stats file checks
  stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_Host_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_Host_params.csv"))
  pred_file  <- file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, ".tif"))
  
  if(file.exists(stats_file) && file.exists(pred_file)) return(NULL)
  
  write_log(MASTER_LOG, paste("START Host:", sp))
  sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
  env_stack <- terra::unwrap(get("packed_env"))
  
  tryCatch({
    sp_dat <- anem_occ %>% dplyr::filter(species == sp)
    if(nrow(sp_dat) < 5) stop("Not enough occurrences")
    
    if(USE_SPATIAL_THINNING) sp_dat <- thin_occurrences(sp_dat, env_stack)
    
    if(BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
      bg_coords <- get_bias_corrected_background(sp_dat, env_stack, method=BG_SAMPLING_METHOD)
    } else {
      bg_coords <- get_random_background(sp_dat, env_stack)
    }
    
    future_stacks <- list()
    for(i in seq_along(FUT_FILES)) {
      scen <- scenario_names[i]
      f_stack <- prepare_future_stack(FUT_FILES[i], env_stack, STATIC_VARS)
      vars_needed <- names(env_stack)
      if(all(vars_needed %in% names(f_stack))) future_stacks[[scen]] <- f_stack[[vars_needed]]
    }
    
    if(file.exists(tune_file)) {
      write_log(MASTER_LOG, paste("  > Skipped Tuning:", sp))
      params <- read_csv(tune_file, show_col_types=FALSE)
    } else {
      write_log(MASTER_LOG, paste("  > Tuning:", sp))
      params <- get_best_params(sp_dat, env_stack, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if(is.null(params)) stop("Tuning failed")
      params$fc <- tolower(params$fc)
      write_csv(as.data.frame(params), tune_file)
    }
    
    write_log(MASTER_LOG, paste("  > Running Model:", sp))
    # Note: Passed output_dir for saving RDS iterations
    results <- fit_bootstrap_worker(sp_dat, env_stack, future_stacks, bg_coords, params, n_boot=N_HOST_BOOT, 
                                    sp_name=sp_clean, model_type="Host", output_dir=OUTPUT_ROOT, debug_log=sp_log)
    
    write_csv(results$stats, stats_file)
    terra::writeRaster(results$current, pred_file, overwrite=TRUE)
    for(scen in names(results$future)) terra::writeRaster(results$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, ".tif")), overwrite=TRUE)
    
    write_log(MASTER_LOG, paste("FINISH Host:", sp))
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Host:", sp, "->", e$message))
  })
}

# ==============================================================================
# PHASE 2: FISH (ENV ONLY + COMBINED)
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 2: FISH ---")
host_files <- list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=TRUE, pattern="\\.tif$")
if(length(host_files) == 0) stop("CRITICAL: Phase 1 failed.")

hosts_curr <- terra::rast(host_files)
names(hosts_curr) <- tools::file_path_sans_ext(basename(host_files)) 
packed_hosts_curr <- terra::wrap(hosts_curr)
clusterExport(cl, "packed_hosts_curr", envir = environment())

if(!is.null(TARGET_FISH)) fish_run_list <- intersect(unique(amph_occ$species), TARGET_FISH) else fish_run_list <- unique(amph_occ$species)

foreach(sp = fish_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
  env_stack_curr  <- terra::unwrap(get("packed_env"))
  host_stack_curr <- terra::unwrap(get("packed_hosts_curr"))
  
  tryCatch({
    sp_dat <- amph_occ %>% dplyr::filter(species == sp)
    if(nrow(sp_dat) < 5) stop("Not enough occurrences")
    if(USE_SPATIAL_THINNING) sp_dat <- thin_occurrences(sp_dat, env_stack_curr)
    
    # ---------------------------------------------------------
    # A. ENV ONLY MODEL
    # ---------------------------------------------------------
    write_log(MASTER_LOG, paste("START Fish EnvOnly:", sp))
    
    if(BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
      bg_coords <- get_bias_corrected_background(sp_dat, env_stack_curr, method=BG_SAMPLING_METHOD)
    } else {
      bg_coords <- get_random_background(sp_dat, env_stack_curr)
    }
    
    # Prepare Futures (Env Only)
    future_stacks_env <- list()
    for(i in seq_along(FUT_FILES)) {
      scen <- scenario_names[i]
      f_stack <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
      vars_needed <- names(env_stack_curr)
      if(all(vars_needed %in% names(f_stack))) future_stacks_env[[scen]] <- f_stack[[vars_needed]]
    }
    
    # Tuning (Env Only)
    tune_file_env <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishEnvOnly_params.csv"))
    if(file.exists(tune_file_env)) {
      params_env <- read_csv(tune_file_env, show_col_types=FALSE)
    } else {
      params_env <- get_best_params(sp_dat, env_stack_curr, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if(!is.null(params_env)) {
        params_env$fc <- tolower(params_env$fc)
        write_csv(as.data.frame(params_env), tune_file_env)
      }
    }
    
    if(!is.null(params_env)) {
      results_env <- fit_bootstrap_worker(sp_dat, env_stack_curr, future_stacks_env, bg_coords, params_env, n_boot=N_FISH_BOOT, 
                                          sp_name=sp_clean, model_type="EnvOnly", output_dir=OUTPUT_ROOT, debug_log=sp_log)
      
      write_csv(results_env$stats, file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishEnvOnly_stats.csv")))
      terra::writeRaster(results_env$current, file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", paste0(sp_clean, ".tif")), overwrite=TRUE)
      for(scen in names(results_env$future)) terra::writeRaster(results_env$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only", paste0(sp_clean, ".tif")), overwrite=TRUE)
    }
    
    # ---------------------------------------------------------
    # B. COMBINED MODEL (Env + Host)
    # ---------------------------------------------------------
    write_log(MASTER_LOG, paste("START Fish Combined:", sp))
    
    biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat, debug_path=MASTER_LOG)
    if(!is.null(biotic_curr)) {
      full_stack_curr <- c(env_stack_curr, biotic_curr)
      names(full_stack_curr) <- c(names(env_stack_curr), "biotic_suitability")
      
      # Re-generate Background for full stack (technically same coords, but safer to re-extract)
      # Note: We reuse 'bg_coords' geometry, but extraction happens inside worker
      
      # Futures (Combined)
      future_stacks_comb <- list()
      for(i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        f_env <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
        f_env <- f_env[[names(env_stack_curr)]] 
        host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
        host_files_fut <- list.files(host_dir_fut, full.names=TRUE, pattern="\\.tif$")
        
        if(length(host_files_fut) > 0) {
          host_stack_fut <- terra::rast(host_files_fut)
          names(host_stack_fut) <- tools::file_path_sans_ext(basename(host_files_fut))
          biotic_fut <- get_biotic_layer(sp, host_stack_fut, int_mat) 
          if(!is.null(biotic_fut)) {
            f_comb <- c(f_env, biotic_fut)
            names(f_comb) <- c(names(f_env), "biotic_suitability")
            future_stacks_comb[[scen]] <- f_comb
          }
        }
      }
      
      # Tuning (Combined)
      tune_file_comb <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishCombined_params.csv"))
      if(file.exists(tune_file_comb)) {
        params_comb <- read_csv(tune_file_comb, show_col_types=FALSE)
      } else {
        params_comb <- get_best_params(sp_dat, full_stack_curr, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
        if(!is.null(params_comb)) {
          params_comb$fc <- tolower(params_comb$fc)
          write_csv(as.data.frame(params_comb), tune_file_comb)
        }
      }
      
      if(!is.null(params_comb)) {
        results_comb <- fit_bootstrap_worker(sp_dat, full_stack_curr, future_stacks_comb, bg_coords, params_comb, n_boot=N_FISH_BOOT, 
                                             sp_name=sp_clean, model_type="Combined", output_dir=OUTPUT_ROOT, debug_log=sp_log)
        
        write_csv(results_comb$stats, file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishCombined_stats.csv")))
        terra::writeRaster(results_comb$current, file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, ".tif")), overwrite=TRUE)
        for(scen in names(results_comb$future)) terra::writeRaster(results_comb$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, ".tif")), overwrite=TRUE)
      }
    }
    
    write_log(MASTER_LOG, paste("FINISH Fish:", sp))
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish:", sp, "->", e$message))
  })
}

stopCluster(cl)
write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
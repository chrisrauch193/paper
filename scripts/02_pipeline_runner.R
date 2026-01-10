# scripts/02_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL SDM PIPELINE: CLEAN CLUSTER START
# ------------------------------------------------------------------------------

pipeline_main <- function() {
  
  # --- 1. CONFIGURATION ---
  RUN_MODE <- "TEST" 
  
  TARGET_HOSTS <- c("Entacmaea_quadricolor", "Heteractis_magnifica")
  TARGET_FISH  <- c("Amphiprion_clarkii", "Amphiprion_frenatus")
  
  STATIC_VARS <- c("rugosity", "bathymetry", "slope", "aspect") 
  
  if (RUN_MODE == "FINAL") {
    N_CORES     <- 30 
    N_HOST_BOOT <- 10 
    N_FISH_BOOT <- 40 
    OUTPUT_ROOT <- "outputs/final_run"
  } else {
    N_CORES     <- 4
    N_HOST_BOOT <- 2
    N_FISH_BOOT <- 2
    OUTPUT_ROOT <- "outputs/test_run"
  }
  
  # --- 2. SETUP ---
  if(!require("pacman")) install.packages("pacman")
  pacman::p_load(foreach, doParallel, terra, dplyr, readr, ENMeval, maxnet, dismo)
  
  # Source locally for the main thread
  source("scripts/01_core_functions.R")
  
  DIRS <- c("models_stats", "predictions/current/hosts", "predictions/current/fish_combined", "predictions/current/fish_env_only", "logs")
  for(d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive=TRUE, showWarnings=FALSE)
  
  MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
  if(!file.exists(MASTER_LOG)) write_log(MASTER_LOG, paste("--- PIPELINE INIT:", RUN_MODE, "---"))
  
  # --- 3. DATA LOADING ---
  ENV_PATH <- "data/final_env_stack.tif"
  if(!file.exists(ENV_PATH)) stop("Env Raster missing.")
  
  cat("Loading Environmental Raster...\n")
  current_env <- terra::rast(ENV_PATH)
  packed_env <- terra::wrap(current_env) 
  
  # Load Futures
  FUT_DIR   <- "data/env/future_pca" 
  FUT_FILES <- list.files(FUT_DIR, full=TRUE, pattern="\\.tif$")
  scenario_names <- tools::file_path_sans_ext(basename(FUT_FILES))
  for(scen in scenario_names) {
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"), recursive=T, showWarnings=F)
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"), recursive=T, showWarnings=F)
  }
  
  # Load Occurrences
  amph_occ <- as.data.frame(read_csv("data/amph_occ_env_final_dataset.csv", show_col_types=FALSE))
  anem_occ <- as.data.frame(read_csv("data/anem_occ_env_final_dataset.csv", show_col_types=FALSE))
  int_mat  <- read.csv("data/interaction_matrix.csv", row.names=1)
  colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
  rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))
  
  # --- 4. CLUSTER INIT ---
  try(stopCluster(cl), silent=TRUE)  # optional; fine to keep
  cl <- makeCluster(N_CORES)
  
  # IMPORTANT: on.exit MUST be inside a function (or local({...})) to work as intended
  on.exit({
    try(stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  registerDoParallel(cl)
  
  # 1. Clean Workers
  clusterEvalQ(cl, { rm(list=ls()); gc() })
  
  # 2. Export Data
  clusterExport(cl, c("packed_env", "FUT_FILES", "scenario_names", "STATIC_VARS",
                      "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG", 
                      "N_HOST_BOOT", "N_FISH_BOOT"))
  
  # 3. Load Libraries & Source Functions Fresh
  clusterEvalQ(cl, { 
    library(terra); library(dplyr); library(maxnet); library(ENMeval); library(readr); library(dismo)
    source("scripts/01_core_functions.R") # Loads the Anonymization logic
    
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
  write_log(MASTER_LOG, "--- PHASE 1: HOSTS CHECK ---")
  
  all_anem_species <- unique(anem_occ$species)
  if(!is.null(TARGET_HOSTS)) anem_run_list <- intersect(all_anem_species, TARGET_HOSTS) else anem_run_list <- all_anem_species
  
  host_results <- foreach(sp = anem_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
    sp_clean <- gsub(" ", "_", sp)
    stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_Host_stats.csv"))
    
    if(file.exists(stats_file)) return(paste("SKIPPED:", sp))
    
    write_log(MASTER_LOG, paste("START Host:", sp))
    sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
    env_stack <- terra::unwrap(get("packed_env"))
    
    tryCatch({
      sp_dat <- anem_occ %>% dplyr::filter(species == sp)
      if(nrow(sp_dat) < 5) stop("Not enough occurrences")
      
      write_log(MASTER_LOG, paste("  > Tuning parameters for:", sp))
      
      future_stacks <- list()
      for(i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        f_stack <- prepare_future_stack(FUT_FILES[i], env_stack, STATIC_VARS)
        vars_needed <- names(env_stack)
        if(all(vars_needed %in% names(f_stack))) future_stacks[[scen]] <- f_stack[[vars_needed]]
      }
      
      params <- get_best_params(sp_dat, env_stack)
      if(is.null(params)) stop("Tuning failed")
      params$fc <- tolower(params$fc)
      
      write_log(MASTER_LOG, paste("  > Tuning Complete:", sp, "|", params$fc, params$rm))
      
      results <- fit_bootstrap_worker(sp_dat, env_stack, future_stacks, params, n_boot=N_HOST_BOOT, debug_log=sp_log)
      
      write_csv(results$stats, stats_file)
      terra::writeRaster(results$current, file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, ".tif")), overwrite=TRUE)
      for(scen in names(results$future)) {
        terra::writeRaster(results$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, ".tif")), overwrite=TRUE)
      }
      
      write_log(MASTER_LOG, paste("FINISH Host:", sp))
      return(paste("SUCCESS:", sp))
      
    }, error = function(e) {
      write_log(MASTER_LOG, paste("ERROR Host:", sp, "->", e$message))
      return(paste("FAILED:", sp, "-", e$message))
    })
  }
  
  # ==============================================================================
  # PHASE 2: FISH
  # ==============================================================================
  write_log(MASTER_LOG, "--- PHASE 2: FISH CHECK ---")
  
  host_files <- list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=TRUE, pattern="\\.tif$")
  if(length(host_files) == 0) stop("CRITICAL: Phase 1 failed (No Host Maps).")
  
  packed_hosts_curr <- terra::wrap(terra::rast(host_files))
  clusterExport(cl, "packed_hosts_curr")
  
  all_fish_species <- unique(amph_occ$species)
  if(!is.null(TARGET_FISH)) fish_run_list <- intersect(all_fish_species, TARGET_FISH) else fish_run_list <- all_fish_species
  
  fish_results <- foreach(sp = fish_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
    sp_clean <- gsub(" ", "_", sp)
    stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishCombined_stats.csv"))
    
    if(file.exists(stats_file)) return(paste("SKIPPED:", sp))
    
    write_log(MASTER_LOG, paste("START Fish:", sp))
    sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
    
    env_stack_curr  <- terra::unwrap(get("packed_env"))
    host_stack_curr <- terra::unwrap(get("packed_hosts_curr"))
    
    tryCatch({
      sp_dat <- amph_occ %>% dplyr::filter(species == sp)
      if(nrow(sp_dat) < 5) stop("Not enough occurrences")
      
      biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat)
      if(is.null(biotic_curr)) stop("No matching hosts found")
      
      full_stack_curr <- c(env_stack_curr, biotic_curr)
      names(full_stack_curr) <- c(names(env_stack_curr), "biotic_suitability")
      
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
      
      write_log(MASTER_LOG, paste("  > Tuning parameters for:", sp))
      params <- get_best_params(sp_dat, full_stack_curr)
      if(is.null(params)) stop("Tuning failed")
      
      write_log(MASTER_LOG, paste("  > Tuning Complete:", sp, "|", params$fc, params$rm))
      
      results <- fit_bootstrap_worker(sp_dat, full_stack_curr, future_stacks, params, n_boot=N_FISH_BOOT, debug_log=sp_log)
      
      write_csv(results$stats, stats_file)
      terra::writeRaster(results$current, file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, ".tif")), overwrite=TRUE)
      for(scen in names(results$future)) {
        terra::writeRaster(results$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, ".tif")), overwrite=TRUE)
      }
      
      write_log(MASTER_LOG, paste("FINISH Fish:", sp))
      return(paste("SUCCESS:", sp))
      
    }, error = function(e) {
      write_log(MASTER_LOG, paste("ERROR Fish:", sp, "->", e$message))
      return(paste("FAILED:", sp, "-", e$message))
    })
  }
  
  print(host_results)
  print(fish_results)
  stopCluster(cl)
  write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
}

pipeline_main()

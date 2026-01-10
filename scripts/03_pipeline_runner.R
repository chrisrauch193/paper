# scripts/03_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL SDM PIPELINE: GOLD STANDARD (3-Stage Fish + Future Biotic)
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
DIRS <- c("models_stats", "models_tuning", "models", 
          "predictions/current/hosts", 
          "predictions/current/fish_env_only", 
          "predictions/current/fish_host_only", 
          "predictions/current/fish_combined", "logs")
for(d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive=TRUE, showWarnings=FALSE)

# 4. AUTO-WIPE
if(WIPE_PREDICTIONS) {
  try(unlink(file.path(OUTPUT_ROOT, "models_stats", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "models", "*")), silent=T) 
  try(unlink(file.path(OUTPUT_ROOT, "models_tuning", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/hosts", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/fish_host_only", "*")), silent=T)
  try(unlink(file.path(OUTPUT_ROOT, "predictions/current/fish_combined", "*")), silent=T)
  # Note: Future folders are wiped/recreated in the loop below
}

# 5. LOGGING
MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
if(file.exists(MASTER_LOG)) unlink(MASTER_LOG)
file.create(MASTER_LOG)
write_log(MASTER_LOG, paste("--- PIPELINE START | ID:", RUN_ID, "| BG:", BG_SAMPLING_METHOD, "---"))

# 6. DATA LOADING
cat("Loading Data...\n")
if(!file.exists(ENV_PATH)) stop("Env Raster missing.")
current_env <- terra::rast(ENV_PATH)
packed_env  <- terra::wrap(current_env) 

FUT_FILES <- list.files(FUT_DIR, full=TRUE, pattern="\\.tif$")
if(length(FUT_FILES) == 0) warning("No future layers found in ", FUT_DIR)

scenario_names <- tools::file_path_sans_ext(basename(FUT_FILES))
# Create future directories for all scenarios (2050, 2100, etc.)
for(scen in scenario_names) {
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"), recursive=T, showWarnings=F)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only"), recursive=T, showWarnings=F)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only"), recursive=T, showWarnings=F)
  dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"), recursive=T, showWarnings=F)
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

clusterExport(cl, c("packed_env", "FUT_FILES", "scenario_names", "STATIC_VARS",
                    "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG", 
                    "N_HOST_BOOT", "N_FISH_BOOT", "TUNE_ARGS", 
                    "USE_SPATIAL_THINNING", "USE_SPATIAL_TUNING", "BG_SAMPLING_METHOD",
                    "write_log", "thin_occurrences", "get_bias_corrected_background", "get_random_background",
                    "get_biotic_layer", "get_best_params", "fit_bootstrap_worker"),
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
  stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_Host_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_Host_params.csv"))
  pred_file  <- file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, "_mean.tif"))
  
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
    
    results <- fit_bootstrap_worker(sp_dat, env_stack, future_stacks, bg_coords, params, n_boot=N_HOST_BOOT, 
                                    sp_name=sp_clean, model_type="Host", output_dir=OUTPUT_ROOT, debug_log=sp_log)
    
    write_csv(results$stats, stats_file)
    terra::writeRaster(results$current$mean, file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
    terra::writeRaster(results$current$sd,   file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
    
    for(scen in names(results$future)) {
      terra::writeRaster(results$future[[scen]]$mean, file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
      terra::writeRaster(results$future[[scen]]$sd,   file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
    }
    
    write_log(MASTER_LOG, paste("FINISH Host:", sp))
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Host:", sp, "->", e$message))
  })
}

# ==============================================================================
# PHASE 2: FISH (Models A, B, C)
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 2: FISH ---")
host_files <- list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=TRUE, pattern="_mean\\.tif$")
if(length(host_files) == 0) stop("CRITICAL: Phase 1 failed.")

hosts_curr <- terra::rast(host_files)
names(hosts_curr) <- tools::file_path_sans_ext(basename(host_files)) %>% gsub("_mean", "", .)
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
    # A. ENV ONLY MODEL (PC1-4 + Rugosity ONLY)
    # ---------------------------------------------------------
    write_log(MASTER_LOG, paste("START Fish EnvOnly:", sp))
    
    # Define variables: PCs + Rugosity. Drop extraneous vars (e.g. Slope/Bathy) if unused.
    # We grab names from the first Future file to know which PCs exist.
    ref_fut <- terra::rast(FUT_FILES[1])
    vars_to_keep_env <- c(names(ref_fut), STATIC_VARS)
    env_stack_model_a <- env_stack_curr[[intersect(names(env_stack_curr), vars_to_keep_env)]]
    
    # Background
    if(BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
      bg_coords <- get_bias_corrected_background(sp_dat, env_stack_model_a, method=BG_SAMPLING_METHOD)
    } else {
      bg_coords <- get_random_background(sp_dat, env_stack_model_a)
    }
    
    # Futures
    future_stacks_env <- list()
    for(i in seq_along(FUT_FILES)) {
      scen <- scenario_names[i]
      f_stack <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
      if(all(names(env_stack_model_a) %in% names(f_stack))) {
        future_stacks_env[[scen]] <- f_stack[[names(env_stack_model_a)]]
      }
    }
    
    # Tune & Run
    tune_file_env <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishEnvOnly_params.csv"))
    if(file.exists(tune_file_env)) {
      params_env <- read_csv(tune_file_env, show_col_types=FALSE)
    } else {
      params_env <- get_best_params(sp_dat, env_stack_model_a, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if(!is.null(params_env)) {
        params_env$fc <- tolower(params_env$fc)
        write_csv(as.data.frame(params_env), tune_file_env)
      }
    }
    
    if(!is.null(params_env)) {
      res_env <- fit_bootstrap_worker(sp_dat, env_stack_model_a, future_stacks_env, bg_coords, params_env, n_boot=N_FISH_BOOT, 
                                      sp_name=sp_clean, model_type="EnvOnly", output_dir=OUTPUT_ROOT, debug_log=sp_log)
      
      write_csv(res_env$stats, file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishEnvOnly_stats.csv")))
      terra::writeRaster(res_env$current$mean, file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
      terra::writeRaster(res_env$current$sd,   file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
      for(scen in names(res_env$future)) {
        terra::writeRaster(res_env$future[[scen]]$mean, file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
        terra::writeRaster(res_env$future[[scen]]$sd,   file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
      }
    }
    
    # ---------------------------------------------------------
    # PREP HOST DATA (Current)
    # ---------------------------------------------------------
    biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat, debug_path=MASTER_LOG)
    
    if(!is.null(biotic_curr)) {
      
      # ---------------------------------------------------------
      # B. HOST ONLY MODEL (Host + Rugosity)
      # ---------------------------------------------------------
      write_log(MASTER_LOG, paste("START Fish HostOnly:", sp))
      
      if("rugosity" %in% names(env_stack_curr)) {
        rug_layer <- env_stack_curr[["rugosity"]]
        host_only_stack <- c(rug_layer, biotic_curr)
        names(host_only_stack) <- c("rugosity", "biotic_suitability")
        
        # Futures for Host Only (Static Rugosity + Dynamic Future Host)
        future_stacks_hostonly <- list()
        for(i in seq_along(FUT_FILES)) {
          scen <- scenario_names[i]
          # 1. Get Static Rugosity
          f_env_raw <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
          f_rug <- f_env_raw[["rugosity"]]
          
          # 2. Get Future Biotic
          host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
          host_files_fut <- list.files(host_dir_fut, full.names=TRUE, pattern="_mean\\.tif$")
          if(length(host_files_fut) > 0) {
            h_stack_f <- terra::rast(host_files_fut)
            names(h_stack_f) <- tools::file_path_sans_ext(basename(host_files_fut)) %>% gsub("_mean", "", .)
            biotic_fut <- get_biotic_layer(sp, h_stack_f, int_mat)
            if(!is.null(biotic_fut)) {
              f_comb <- c(f_rug, biotic_fut)
              names(f_comb) <- c("rugosity", "biotic_suitability")
              future_stacks_hostonly[[scen]] <- f_comb
            }
          }
        }
        
        tune_file_host <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishHostOnly_params.csv"))
        if(file.exists(tune_file_host)) {
          params_host <- read_csv(tune_file_host, show_col_types=FALSE)
        } else {
          params_host <- get_best_params(sp_dat, host_only_stack, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
          if(!is.null(params_host)) {
            params_host$fc <- tolower(params_host$fc)
            write_csv(as.data.frame(params_host), tune_file_host)
          }
        }
        
        if(!is.null(params_host)) {
          res_ho <- fit_bootstrap_worker(sp_dat, host_only_stack, future_stacks_hostonly, bg_coords, params_host, n_boot=N_FISH_BOOT, 
                                         sp_name=sp_clean, model_type="HostOnly", output_dir=OUTPUT_ROOT, debug_log=sp_log)
          
          write_csv(res_ho$stats, file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishHostOnly_stats.csv")))
          terra::writeRaster(res_ho$current$mean, file.path(OUTPUT_ROOT, "predictions/current/fish_host_only", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
          terra::writeRaster(res_ho$current$sd,   file.path(OUTPUT_ROOT, "predictions/current/fish_host_only", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
          for(scen in names(res_ho$future)) {
            terra::writeRaster(res_ho$future[[scen]]$mean, file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
            terra::writeRaster(res_ho$future[[scen]]$sd,   file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
          }
        }
      }
      
      # ---------------------------------------------------------
      # C. COMBINED MODEL (Env + Host)
      # ---------------------------------------------------------
      write_log(MASTER_LOG, paste("START Fish Combined:", sp))
      
      # Stack: Env Only (PCs+Rug) + Biotic
      full_stack_curr <- c(env_stack_model_a, biotic_curr)
      names(full_stack_curr) <- c(names(env_stack_model_a), "biotic_suitability")
      
      # Futures (Env + Future Host)
      future_stacks_comb <- list()
      for(i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        # 1. Get Env Part
        f_stack_raw <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
        if(all(names(env_stack_model_a) %in% names(f_stack_raw))) {
          f_env <- f_stack_raw[[names(env_stack_model_a)]]
          
          # 2. Get Biotic Part (Future)
          host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
          host_files_fut <- list.files(host_dir_fut, full.names=TRUE, pattern="_mean\\.tif$")
          if(length(host_files_fut) > 0) {
            h_stack_f <- terra::rast(host_files_fut)
            names(h_stack_f) <- tools::file_path_sans_ext(basename(host_files_fut)) %>% gsub("_mean", "", .)
            biotic_fut <- get_biotic_layer(sp, h_stack_f, int_mat) 
            if(!is.null(biotic_fut)) {
              f_comb <- c(f_env, biotic_fut)
              names(f_comb) <- c(names(f_env), "biotic_suitability")
              future_stacks_comb[[scen]] <- f_comb
            }
          }
        }
      }
      
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
        res_comb <- fit_bootstrap_worker(sp_dat, full_stack_curr, future_stacks_comb, bg_coords, params_comb, n_boot=N_FISH_BOOT, 
                                         sp_name=sp_clean, model_type="Combined", output_dir=OUTPUT_ROOT, debug_log=sp_log)
        
        write_csv(res_comb$stats, file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishCombined_stats.csv")))
        terra::writeRaster(res_comb$current$mean, file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
        terra::writeRaster(res_comb$current$sd,   file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
        for(scen in names(res_comb$future)) {
          terra::writeRaster(res_comb$future[[scen]]$mean, file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, "_mean.tif")), overwrite=TRUE)
          terra::writeRaster(res_comb$future[[scen]]$sd,   file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, "_sd.tif")), overwrite=TRUE)
        }
      }
    }
    
    write_log(MASTER_LOG, paste("FINISH Fish:", sp))
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish:", sp, "->", e$message))
  })
}

stopCluster(cl)
write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
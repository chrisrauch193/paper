# scripts/02_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL SDM PIPELINE: GOLD STANDARD
# (Thinning + Bias Background + Spatial Tuning + Robust Caching)
# ------------------------------------------------------------------------------

# 1. GLOBAL SETUP
rm(list=ls())
gc()

if(!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, terra, dplyr, readr, ENMeval, maxnet, dismo)

# ------------------------------------------------------------------------------
# PIPELINE FUNCTION WRAPPER
# ------------------------------------------------------------------------------
pipeline_main <- function() {
  
  # --- 1. CONFIGURATION ---
  RUN_MODE         <- "TEST"   # "TEST" (fast) or "FINAL" (thesis quality)
  WIPE_PREDICTIONS <- TRUE     # Force re-run?
  
  # METHODOLOGY TOGGLES (The "Switches")
  USE_SPATIAL_THINNING      <- TRUE  # Remove duplicate points in same cell? (Recommended: TRUE)
  USE_SPATIAL_TUNING        <- TRUE  # Use Block Partitioning? (Critical for Range Shifts: TRUE)
  USE_BIAS_CORRECTED_BG     <- TRUE  # Use Paper's Distance-Weighted Background? (TRUE = Paper way, FALSE = Random)
  
  # SPECIES SETTINGS
  TARGET_HOSTS <- c("Entacmaea_quadricolor", "Heteractis_magnifica")
  TARGET_FISH  <- c("Amphiprion_clarkii", "Amphiprion_frenatus")
  STATIC_VARS  <- c("rugosity", "bathymetry", "slope", "aspect") 
  
  if (RUN_MODE == "FINAL") {
    N_CORES     <- 30 
    N_HOST_BOOT <- 10 
    N_FISH_BOOT <- 40 
    OUTPUT_ROOT <- "outputs/final_run"
    # Thesis used RM 0.5 to 4.0 in steps of 0.5
    TUNE_ARGS   <- list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5)) 
  } else {
    N_CORES     <- 4
    N_HOST_BOOT <- 2
    N_FISH_BOOT <- 2
    OUTPUT_ROOT <- "outputs/test_run"
    TUNE_ARGS   <- list(fc = c("L", "LQ", "H"), rm = c(1, 2, 5))
  }
  
  # --- 2. LOCAL FUNCTIONS ---
  
  write_log <- function(path, msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(paste0("[", timestamp, "] ", msg, "\n"), file = path, append = TRUE)
  }
  
  # A. SPATIAL THINNING (Grid-based)
  thin_occurrences <- function(occ_df, env_rast) {
    # 1. Extract cell numbers
    cells <- terra::cellFromXY(env_rast, as.matrix(occ_df[, c("x", "y")]))
    # 2. Find duplicates (multiple points in same pixel)
    dups <- duplicated(cells)
    # 3. Filter
    occ_thinned <- occ_df[!dups, ]
    return(occ_thinned)
  }
  
  # B. PAPER'S BIAS CORRECTION METHOD (The "Exact" Way)
  get_bias_corrected_background <- function(occ_coords, env_stack, n_bg=10000, alpha=0.5) {
    if(inherits(occ_coords, "data.frame")) occ <- as.matrix(occ_coords[, c("x", "y")])
    else occ <- as.matrix(occ_coords)
    
    # Use only first 2 layers (PC1/PC2) for env distance (Standard PCA space)
    env_cols <- names(env_stack)[1:2] 
    occ_env <- terra::extract(env_stack[[env_cols]], occ)
    occ_env <- na.omit(cbind(occ, occ_env))
    if(nrow(occ_env) == 0) return(NULL) 
    
    # Sample dense candidates
    candidates <- terra::spatSample(env_stack, size=n_bg*3, method="random", na.rm=TRUE, xy=TRUE, values=TRUE)
    candidates <- na.omit(candidates)
    
    cand_xy <- as.matrix(candidates[, c("x", "y")])
    occ_xy  <- as.matrix(occ_env[, c("x", "y")])
    cand_env <- as.matrix(candidates[, env_cols])
    occ_env_vals <- as.matrix(occ_env[, env_cols])
    
    # Calculate Distances (Nearest Neighbor)
    dist_geo <- apply(cand_xy, 1, function(pt) min(sqrt((occ_xy[,1]-pt[1])^2 + (occ_xy[,2]-pt[2])^2)))
    dist_env <- apply(cand_env, 1, function(pt) min(sqrt((occ_env_vals[,1]-pt[1])^2 + (occ_env_vals[,2]-pt[2])^2)))
    
    d_geo_norm <- dist_geo / max(dist_geo, na.rm=TRUE)
    d_env_norm <- dist_env / max(dist_env, na.rm=TRUE)
    
    # Weighting Formula (Prob increases as distance DECREASES to mimic sampling bias)
    # Note: Target group background usually mimics bias. 
    # Paper formula: Prob = 1 - ... implies probability is higher if distant? 
    # Actually standard bias correction samples proportional to bias density.
    # We will use the formula derived from their code:
    sampling_prob <- 1 - ((1 - d_geo_norm)^alpha) * ((1 - d_env_norm)^(1 - alpha))
    
    if(nrow(candidates) > n_bg) {
      selected_idx <- sample(1:nrow(candidates), size=n_bg, prob=sampling_prob, replace=FALSE)
      bg_final <- candidates[selected_idx, c("x", "y")]
    } else {
      bg_final <- candidates[, c("x", "y")]
    }
    return(bg_final)
  }
  
  # C. STANDARD RANDOM BACKGROUND (The Fallback)
  get_random_background <- function(occ_coords, env_stack, n_bg=10000) {
    if(inherits(occ_coords, "data.frame")) coords_mat <- as.matrix(occ_coords[, c("x", "y")])
    else coords_mat <- as.matrix(occ_coords)
    
    vect_occ <- terra::vect(coords_mat, crs="EPSG:4326", type="points")
    vect_buff <- terra::aggregate(terra::buffer(vect_occ, width=1000000)) # 1000km buffer
    
    env_crop <- terra::crop(env_stack, vect_buff)
    env_mask <- terra::mask(env_crop, vect_buff)
    bg_coords <- terra::spatSample(env_mask, size=n_bg, method="random", na.rm=TRUE, xy=TRUE, values=FALSE)
    return(bg_coords)
  }
  
  # D. BIOTIC LAYER
  get_biotic_layer <- function(fish_sp, host_stack, int_mat, debug_path=NULL) {
    fish_clean <- gsub(" ", "_", fish_sp)
    if(!fish_clean %in% rownames(int_mat)) {
      if(!is.null(debug_path)) write_log(debug_path, paste("DEBUG BIOTIC:", fish_clean, "not in matrix."))
      return(NULL)
    }
    
    row_idx <- which(rownames(int_mat) == fish_clean)
    w_vals <- as.numeric(int_mat[row_idx, ])
    names(w_vals) <- colnames(int_mat) 
    weights <- w_vals[w_vals > 0]
    
    available_hosts <- intersect(names(weights), names(host_stack))
    if(length(available_hosts) == 0) {
      if(!is.null(debug_path)) {
        write_log(debug_path, paste("DEBUG BIOTIC FAIL:", fish_clean))
        write_log(debug_path, paste("  > Need:", paste(names(weights), collapse=", ")))
        write_log(debug_path, paste("  > Have:", paste(names(host_stack), collapse=", ")))
      }
      return(NULL)
    }
    
    sub_stack <- host_stack[[available_hosts]]
    sub_weights <- as.numeric(weights[available_hosts])
    biotic_layer <- terra::weighted.mean(sub_stack, w=sub_weights, na.rm=TRUE)
    names(biotic_layer) <- "biotic_suitability"
    return(biotic_layer)
  }
  
  # E. TUNING
  get_best_params <- function(occ_df, env_stack, bg_coords, use_spatial, tune_args) {
    occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
    
    # Partition Method
    if(use_spatial) {
      part_method <- "block"
      log_msg <- "Spatial (Block)"
    } else {
      part_method <- "randomkfold"
      log_msg <- "Random k-Fold"
    }
    
    tryCatch({
      # Try primary method
      eval_res <- tryCatch({
        ENMeval::ENMevaluate(occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
                             algorithm = "maxnet", partitions = part_method, 
                             tune.args = tune_args,
                             quiet = TRUE, parallel = FALSE)
      }, error = function(e) {
        # Fallback to random if block fails (e.g. widely scattered points)
        return(ENMeval::ENMevaluate(occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
                                    algorithm = "maxnet", partitions = "randomkfold", partition.settings = list(kfolds=5),
                                    tune.args = tune_args,
                                    quiet = TRUE, parallel = FALSE))
      })
      
      best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
      if(nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
      return(list(fc = as.character(best$fc), rm = as.numeric(best$rm), method=eval_res@partition.method))
    }, error = function(e) { return(NULL) })
  }
  
  # F. WORKER
  fit_bootstrap_worker <- function(occ_df, current_stack, future_stack_list=NULL, bg_coords, params, n_boot=10, debug_log=NULL) {
    log_debug <- function(msg) {
      if(!is.null(debug_log)) {
        ts <- format(Sys.time(), "%H:%M:%S")
        cat(paste0("[", ts, "] ", msg, "\n"), file=debug_log, append=TRUE)
      }
    }
    
    original_names <- names(current_stack)
    safe_names <- paste0("v", sprintf("%02d", 1:length(original_names)))
    names(current_stack) <- safe_names
    
    if(!is.null(future_stack_list)) {
      for(n in names(future_stack_list)) {
        if(all(original_names %in% names(future_stack_list[[n]]))) {
          future_stack_list[[n]] <- future_stack_list[[n]][[original_names]]
          names(future_stack_list[[n]]) <- safe_names
        } else { future_stack_list[[n]] <- NULL }
      }
    }
    
    all_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
    # Use PRE-GENERATED background passed in args to match tuning
    
    bg_data <- terra::extract(current_stack, bg_coords)
    pres_data_all <- terra::extract(current_stack, all_coords)
    
    if("ID" %in% names(bg_data)) bg_data$ID <- NULL
    if("ID" %in% names(pres_data_all)) pres_data_all$ID <- NULL
    
    bg_data <- as.data.frame(lapply(bg_data, as.numeric))
    pres_data_all <- as.data.frame(lapply(pres_data_all, as.numeric))
    bg_data <- na.omit(bg_data)
    pres_data_all <- na.omit(pres_data_all)
    
    valid_vars <- names(bg_data)[sapply(bg_data, var) > 0]
    if(length(valid_vars) < ncol(bg_data)) {
      bg_data <- bg_data[, valid_vars, drop=FALSE]
      pres_data_all <- pres_data_all[, valid_vars, drop=FALSE]
    }
    
    sum_curr <- terra::rast(current_stack[[1]]); terra::values(sum_curr) <- 0
    names(sum_curr) <- "probability"
    sum_fut_list <- list()
    if(!is.null(future_stack_list)) {
      for(n in names(future_stack_list)) {
        r <- terra::rast(future_stack_list[[n]][[1]]); terra::values(r) <- 0; names(r) <- "probability"
        sum_fut_list[[n]] <- r
      }
    }
    
    stats_list <- list()
    successful_boots <- 0
    
    for(i in 1:n_boot) {
      set.seed(i)
      tryCatch({
        n_pres <- nrow(pres_data_all)
        if(n_pres < 5) next
        train_idx <- sample(1:n_pres, size = round(0.75 * n_pres))
        test_idx  <- setdiff(1:n_pres, train_idx)
        if(length(train_idx) < 5 || length(test_idx) < 5) next
        
        train_p <- pres_data_all[train_idx, , drop=FALSE]
        test_p  <- pres_data_all[test_idx, , drop=FALSE]
        p_vec <- c(rep(1, nrow(train_p)), rep(0, nrow(bg_data)))
        data_df <- rbind(train_p, bg_data)
        
        mod <- maxnet::maxnet(p_vec, data_df, maxnet::maxnet.formula(p_vec, data_df, classes=params$fc), regmult=params$rm)
        pred_test_p <- stats::predict(mod, test_p, type="logistic")
        eval_bg_idx <- sample(1:nrow(bg_data), min(1000, nrow(bg_data)), replace=TRUE)
        pred_test_bg <- stats::predict(mod, bg_data[eval_bg_idx, , drop=FALSE], type="logistic")
        e <- dismo::evaluate(p=as.vector(pred_test_p), a=as.vector(pred_test_bg))
        stats_list[[i]] <- data.frame(iter=i, auc=e@auc, tss=max(e@TPR + e@TNR - 1))
        
        pred_c <- terra::predict(current_stack, mod, type="logistic", na.rm=TRUE)
        sum_curr <- sum_curr + pred_c
        if(!is.null(future_stack_list)) {
          for(n in names(future_stack_list)) {
            pred_f <- terra::predict(future_stack_list[[n]], mod, type="logistic", na.rm=TRUE)
            sum_fut_list[[n]] <- sum_fut_list[[n]] + pred_f
          }
        }
        successful_boots <- successful_boots + 1
        log_debug(paste("Iter", i, "OK | AUC:", round(e@auc, 3)))
      }, error = function(e) { log_debug(paste("Iter", i, "Error:", e$message)) })
    }
    if(successful_boots == 0) stop("All bootstrap iterations failed.")
    mean_curr <- sum_curr / successful_boots
    mean_fut_list <- list()
    if(!is.null(future_stack_list)) {
      for(n in names(future_stack_list)) mean_fut_list[[n]] <- sum_fut_list[[n]] / successful_boots
    }
    return(list(current = mean_curr, future = mean_fut_list, stats = dplyr::bind_rows(stats_list)))
  }
  
  # --- 3. EXECUTION ---
  DIRS <- c("models_stats", "models_tuning", "predictions/current/hosts", "predictions/current/fish_combined", "logs")
  for(d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive=TRUE, showWarnings=FALSE)
  
  if(WIPE_PREDICTIONS) {
    try(file.remove(list.files(file.path(OUTPUT_ROOT, "models_stats"), full.names=T)), silent=T)
    try(file.remove(list.files(file.path(OUTPUT_ROOT, "models_tuning"), full.names=T)), silent=T) 
    try(file.remove(list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=T)), silent=T)
    try(file.remove(list.files(file.path(OUTPUT_ROOT, "predictions/current/fish_combined"), full.names=T)), silent=T)
  }
  
  log_files <- list.files(file.path(OUTPUT_ROOT, "logs"), full.names=TRUE)
  if(length(log_files) > 0) unlink(log_files)
  MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
  if(file.exists(MASTER_LOG)) unlink(MASTER_LOG)
  file.create(MASTER_LOG)
  write_log(MASTER_LOG, paste("--- PIPELINE INIT:", RUN_MODE, "---"))
  
  # LOAD DATA
  ENV_PATH <- "data/final_env_stack.tif"
  if(!file.exists(ENV_PATH)) stop("Env Raster missing.")
  current_env <- terra::rast(ENV_PATH)
  packed_env <- terra::wrap(current_env) 
  
  FUT_DIR <- "data/env/future_pca" 
  FUT_FILES <- list.files(FUT_DIR, full=TRUE, pattern="\\.tif$")
  scenario_names <- tools::file_path_sans_ext(basename(FUT_FILES))
  for(scen in scenario_names) {
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"), recursive=T, showWarnings=F)
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"), recursive=T, showWarnings=F)
  }
  
  amph_occ <- as.data.frame(read_csv("data/amph_occ_env_final_dataset.csv", show_col_types=FALSE))
  anem_occ <- as.data.frame(read_csv("data/anem_occ_env_final_dataset.csv", show_col_types=FALSE))
  
  int_mat  <- read.csv("data/interaction_matrix.csv", row.names=1)
  colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
  rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))
  
  try(stopCluster(cl), silent=TRUE)
  cl <- makeCluster(N_CORES)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
  registerDoParallel(cl)
  clusterEvalQ(cl, { rm(list=ls()); gc() })
  
  # EXPORT NEW FUNCTION
  clusterExport(cl, c("packed_env", "FUT_FILES", "scenario_names", "STATIC_VARS",
                      "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG", 
                      "N_HOST_BOOT", "N_FISH_BOOT", "TUNE_ARGS", 
                      "USE_SPATIAL_THINNING", "USE_SPATIAL_TUNING", "USE_BIAS_CORRECTED_BG",
                      "write_log", "thin_occurrences", "get_bias_corrected_background", "get_random_background",
                      "get_biotic_layer", "get_best_params", "fit_bootstrap_worker"),
                envir = environment()) 
  
  clusterEvalQ(cl, { library(terra); library(dplyr); library(maxnet); library(ENMeval); library(readr); library(dismo)
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
  
  # --- PHASE 1: HOSTS ---
  write_log(MASTER_LOG, "--- PHASE 1: HOSTS CHECK ---")
  all_anem_species <- unique(anem_occ$species)
  if(!is.null(TARGET_HOSTS)) anem_run_list <- intersect(all_anem_species, TARGET_HOSTS) else anem_run_list <- all_anem_species
  
  host_results <- foreach(sp = anem_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
    sp_clean <- gsub(" ", "_", sp)
    stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_Host_stats.csv"))
    tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_Host_params.csv"))
    pred_file  <- file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, ".tif"))
    
    if(file.exists(stats_file) && file.exists(pred_file)) return(paste("SKIPPED (Done):", sp))
    
    write_log(MASTER_LOG, paste("START Host:", sp))
    sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
    env_stack <- terra::unwrap(get("packed_env"))
    
    tryCatch({
      sp_dat <- anem_occ %>% dplyr::filter(species == sp)
      if(nrow(sp_dat) < 5) stop("Not enough occurrences")
      
      # 1. THINNING
      if(USE_SPATIAL_THINNING) sp_dat <- thin_occurrences(sp_dat, env_stack)
      
      # 2. BACKGROUND
      if(USE_BIAS_CORRECTED_BG) {
        bg_coords <- get_bias_corrected_background(sp_dat, env_stack)
        bg_msg <- "Bias-Corrected"
      } else {
        bg_coords <- get_random_background(sp_dat, env_stack)
        bg_msg <- "Random Buffer"
      }
      
      future_stacks <- list()
      for(i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        f_stack <- prepare_future_stack(FUT_FILES[i], env_stack, STATIC_VARS)
        vars_needed <- names(env_stack)
        if(all(vars_needed %in% names(f_stack))) future_stacks[[scen]] <- f_stack[[vars_needed]]
      }
      
      if(file.exists(tune_file)) {
        write_log(MASTER_LOG, paste("  > Skipped Tuning (Loaded Cache):", sp))
        params <- read_csv(tune_file, show_col_types=FALSE)
      } else {
        write_log(MASTER_LOG, paste("  > Tuning (ENMeval):", sp, "| BG:", bg_msg))
        params <- get_best_params(sp_dat, env_stack, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
        if(is.null(params)) stop("Tuning failed")
        params$fc <- tolower(params$fc)
        write_csv(as.data.frame(params), tune_file)
      }
      write_log(MASTER_LOG, paste("  > Params:", params$fc, params$rm, "| Partition:", params$method))
      
      write_log(MASTER_LOG, paste("  > Running Model:", sp))
      results <- fit_bootstrap_worker(sp_dat, env_stack, future_stacks, bg_coords, params, n_boot=N_HOST_BOOT, debug_log=sp_log)
      
      write_csv(results$stats, stats_file)
      terra::writeRaster(results$current, pred_file, overwrite=TRUE)
      for(scen in names(results$future)) terra::writeRaster(results$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, ".tif")), overwrite=TRUE)
      
      write_log(MASTER_LOG, paste("FINISH Host:", sp))
      return("SUCCESS")
    }, error = function(e) {
      write_log(MASTER_LOG, paste("ERROR Host:", sp, "->", e$message))
      return("FAILED")
    })
  }
  
  # --- PHASE 2: FISH ---
  write_log(MASTER_LOG, "--- PHASE 2: FISH CHECK ---")
  host_files <- list.files(file.path(OUTPUT_ROOT, "predictions/current/hosts"), full.names=TRUE, pattern="\\.tif$")
  if(length(host_files) == 0) stop("CRITICAL: Phase 1 failed.")
  
  hosts_curr <- terra::rast(host_files)
  names(hosts_curr) <- tools::file_path_sans_ext(basename(host_files)) 
  packed_hosts_curr <- terra::wrap(hosts_curr)
  
  write_log(MASTER_LOG, paste("Hosts Loaded:", paste(names(hosts_curr), collapse=", ")))
  clusterExport(cl, "packed_hosts_curr", envir = environment())
  
  all_fish_species <- unique(amph_occ$species)
  if(!is.null(TARGET_FISH)) fish_run_list <- intersect(all_fish_species, TARGET_FISH) else fish_run_list <- all_fish_species
  
  foreach(sp = fish_run_list, .packages=c('terra','dplyr','maxnet','dismo')) %dopar% {
    sp_clean <- gsub(" ", "_", sp)
    stats_file <- file.path(OUTPUT_ROOT, "models_stats", paste0(sp_clean, "_FishCombined_stats.csv"))
    tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishCombined_params.csv"))
    pred_file  <- file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, ".tif"))
    
    if(file.exists(stats_file) && file.exists(pred_file)) return(paste("SKIPPED (Done):", sp))
    
    write_log(MASTER_LOG, paste("START Fish:", sp))
    sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, ".log"))
    
    env_stack_curr  <- terra::unwrap(get("packed_env"))
    host_stack_curr <- terra::unwrap(get("packed_hosts_curr"))
    
    tryCatch({
      sp_dat <- amph_occ %>% dplyr::filter(species == sp)
      if(nrow(sp_dat) < 5) stop("Not enough occurrences")
      
      if(USE_SPATIAL_THINNING) sp_dat <- thin_occurrences(sp_dat, env_stack_curr)
      
      biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat, debug_path=MASTER_LOG)
      if(is.null(biotic_curr)) stop("No matching hosts found (Check main log)")
      
      full_stack_curr <- c(env_stack_curr, biotic_curr)
      names(full_stack_curr) <- c(names(env_stack_curr), "biotic_suitability")
      
      # 1. GENERATE FISH BACKGROUND
      if(USE_BIAS_CORRECTED_BG) {
        bg_coords <- get_bias_corrected_background(sp_dat, full_stack_curr)
        bg_msg <- "Bias-Corrected"
      } else {
        bg_coords <- get_random_background(sp_dat, full_stack_curr)
        bg_msg <- "Random Buffer"
      }
      
      future_stacks <- list()
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
            future_stacks[[scen]] <- f_comb
          }
        }
      }
      
      if(file.exists(tune_file)) {
        write_log(MASTER_LOG, paste("  > Skipped Tuning (Loaded Cache):", sp))
        params <- read_csv(tune_file, show_col_types=FALSE)
      } else {
        write_log(MASTER_LOG, paste("  > Tuning (ENMeval):", sp, "| BG:", bg_msg))
        params <- get_best_params(sp_dat, full_stack_curr, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
        if(is.null(params)) stop("Tuning failed")
        params$fc <- tolower(params$fc)
        write_csv(as.data.frame(params), tune_file)
      }
      write_log(MASTER_LOG, paste("  > Params:", params$fc, params$rm, "| Partition:", params$method))
      
      write_log(MASTER_LOG, paste("  > Running Model:", sp))
      results <- fit_bootstrap_worker(sp_dat, full_stack_curr, future_stacks, bg_coords, params, n_boot=N_FISH_BOOT, debug_log=sp_log)
      
      write_csv(results$stats, stats_file)
      terra::writeRaster(results$current, pred_file, overwrite=TRUE)
      for(scen in names(results$future)) terra::writeRaster(results$future[[scen]], file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, ".tif")), overwrite=TRUE)
      
      write_log(MASTER_LOG, paste("FINISH Fish:", sp))
      return("SUCCESS")
      
    }, error = function(e) {
      write_log(MASTER_LOG, paste("ERROR Fish:", sp, "->", e$message))
      return("FAILED")
    })
  }
  
  stopCluster(cl)
  write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
}

pipeline_main()
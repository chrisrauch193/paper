# scripts/01_core_functions.R
# ------------------------------------------------------------------------------
# CORE FUNCTIONS: MATCHING DEBUGGER LOGIC
# ------------------------------------------------------------------------------

library(terra)
library(dplyr)
library(maxnet)
library(ENMeval)
library(dismo)

# 1. SIMPLE LOGGER
write_log <- function(path, msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] ", msg, "\n"), file = path, append = TRUE)
}

# 2. BUFFERED BACKGROUND
get_buffered_background <- function(occ_coords, env_stack, n_bg=10000) {
  if(inherits(occ_coords, "data.frame")) {
    coords_mat <- as.matrix(occ_coords[, c("x", "y")])
  } else {
    coords_mat <- as.matrix(occ_coords)
  }
  
  vect_occ <- terra::vect(coords_mat, crs="EPSG:4326", type="points")
  vect_buff <- terra::buffer(vect_occ, width=1000000) 
  vect_buff <- terra::aggregate(vect_buff)
  
  env_crop <- terra::crop(env_stack, vect_buff)
  env_mask <- terra::mask(env_crop, vect_buff)
  
  bg_coords <- terra::spatSample(env_mask, size=n_bg, method="random", 
                                 na.rm=TRUE, xy=TRUE, values=FALSE)
  return(bg_coords)
}

# 3. PARAMETER TUNING
get_best_params <- function(occ_df, env_stack) {
  occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  
  pt_vals <- terra::extract(env_stack, occ_coords)
  if("ID" %in% names(pt_vals)) pt_vals$ID <- NULL
  valid_rows <- complete.cases(pt_vals)
  occ_coords <- occ_coords[valid_rows, ]
  
  if(nrow(occ_coords) < 5) return(NULL) 
  
  bg_coords <- get_buffered_background(occ_coords, env_stack)
  if(nrow(bg_coords) < 100) return(NULL)
  
  tryCatch({
    eval_res <- ENMeval::ENMevaluate(
      occs = occ_coords, bg = bg_coords, envs = env_stack,
      algorithm = "maxnet", partitions = "randomkfold",
      partition.settings = list(kfolds = 5),
      tune.args = list(fc = c("L", "LQ", "H"), rm = c(1, 2, 5)),
      quiet = TRUE, parallel = FALSE
    )
    
    best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
    if(nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
    return(list(fc = as.character(best$fc), rm = as.numeric(best$rm)))
  }, error = function(e) { return(NULL) })
}

# 4. BOOTSTRAP WORKER (The "Nuclear Debug" Logic)
fit_bootstrap_worker <- function(occ_df, current_stack, future_stack_list=NULL, params, n_boot=10, debug_log=NULL) {
  
  log_debug <- function(msg) {
    if(!is.null(debug_log)) {
      ts <- format(Sys.time(), "%H:%M:%S")
      cat(paste0("[", ts, "] ", msg, "\n"), file=debug_log, append=TRUE)
    }
  }
  
  log_debug("--- WORKER STARTED ---")
  
  # --- STEP 1: ANONYMIZE PREDICTORS (Matches Debugger) ---
  original_names <- names(current_stack)
  safe_names <- paste0("v", sprintf("%02d", 1:length(original_names)))
  names(current_stack) <- safe_names
  log_debug(paste("Anonymized Vars:", paste(safe_names, collapse=", ")))
  
  # Handle Futures
  if(!is.null(future_stack_list)) {
    for(n in names(future_stack_list)) {
      if(all(original_names %in% names(future_stack_list[[n]]))) {
        future_stack_list[[n]] <- future_stack_list[[n]][[original_names]]
        names(future_stack_list[[n]]) <- safe_names
      } else {
        future_stack_list[[n]] <- NULL 
      }
    }
  }
  
  # --- STEP 2: DATA EXTRACTION & CLEANING ---
  all_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  bg_coords  <- get_buffered_background(all_coords, current_stack)
  
  bg_data <- terra::extract(current_stack, bg_coords)
  pres_data_all <- terra::extract(current_stack, all_coords)
  
  # Clean ID
  if("ID" %in% names(bg_data)) bg_data$ID <- NULL
  if("ID" %in% names(pres_data_all)) pres_data_all$ID <- NULL
  
  # Force Numeric (Matches Debugger)
  bg_data <- as.data.frame(lapply(bg_data, as.numeric))
  pres_data_all <- as.data.frame(lapply(pres_data_all, as.numeric))
  
  # Remove NAs
  bg_data <- na.omit(bg_data)
  pres_data_all <- na.omit(pres_data_all)
  
  log_debug(paste("Data Ready. P:", nrow(pres_data_all), "BG:", nrow(bg_data)))
  
  # Init Output
  sum_curr <- terra::rast(current_stack[[1]]); terra::values(sum_curr) <- 0
  names(sum_curr) <- "probability"
  
  sum_fut_list <- list()
  if(!is.null(future_stack_list)) {
    for(n in names(future_stack_list)) {
      r <- terra::rast(future_stack_list[[n]][[1]]); terra::values(r) <- 0
      names(r) <- "probability"
      sum_fut_list[[n]] <- r
    }
  }
  
  stats_list <- list()
  successful_boots <- 0
  
  log_debug(paste("predict() comes from:", paste(find("predict"), collapse=" -> ")))
  log_debug(paste("predict binding is:", paste0(deparse(utils::getAnywhere("predict")$where), collapse=" | ")))
  
  # --- STEP 3: LOOP ---
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
      
      # FIT MAXNET
      mod <- maxnet::maxnet(
        p_vec, data_df,
        maxnet::maxnet.formula(p_vec, data_df, classes=params$fc),
        regmult=params$rm
      )
      
      # EVALUATE
      pred_test_p <- stats::predict(mod, test_p, type="logistic")
      
      eval_bg_idx <- sample(1:nrow(bg_data), min(1000, nrow(bg_data)), replace=TRUE)
      pred_test_bg <- stats::predict(mod, bg_data[eval_bg_idx, , drop=FALSE], type="logistic")
      
      e <- dismo::evaluate(p=as.vector(pred_test_p), a=as.vector(pred_test_bg))
      
      stats_list[[i]] <- data.frame(
        iteration = i,
        auc = e@auc,
        tss = max(e@TPR + e@TNR - 1),
        n_train = length(train_idx),
        n_test = length(test_idx)
      )
      
      # PROJECT
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

# 5. BIOTIC LAYER (Same as before)
get_biotic_layer <- function(fish_sp, host_stack, int_mat) {
  fish_clean <- gsub(" ", "_", fish_sp)
  if(!fish_clean %in% rownames(int_mat)) return(NULL)
  weights <- int_mat[fish_clean, ]
  weights <- weights[weights > 0]
  available_hosts <- intersect(names(weights), names(host_stack))
  if(length(available_hosts) == 0) return(NULL)
  sub_stack <- host_stack[[available_hosts]]
  sub_weights <- as.numeric(weights[available_hosts])
  biotic_layer <- terra::weighted.mean(sub_stack, w=sub_weights, na.rm=TRUE)
  names(biotic_layer) <- "biotic_suitability"
  return(biotic_layer)
}
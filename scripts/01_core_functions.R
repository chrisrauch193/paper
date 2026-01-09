# scripts/01_core_functions.R
# ------------------------------------------------------------------------------
# CORE FUNCTIONS (DEFENSIBLE: BUFFERED BACKGROUND)
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

# --- HELPER: GENERATE BUFFERED BACKGROUND ---
# This creates a "donut" or buffer around points to sample background from.
# Prevents inflation of AUC by restricting training to accessible areas.
get_buffered_background <- function(occ_coords, env_stack, buffer_dist_km=1000, n_bg=10000) {
  
  # Create a SpatVector from coordinates
  vect_occ <- terra::vect(occ_coords, geom=c("x", "y"), crs="EPSG:4326")
  
  # Buffer (convert km to m for terra if using lon/lat, but approx is 1 deg ~ 111km)
  # For safety/speed on WGS84, we approximate: 1000km ~ 10 degrees.
  # Adjust 'width' as needed. 10 degrees is a generous dispersal limit.
  vect_buff <- terra::buffer(vect_occ, width=1000000) # 1000 km
  
  # Aggregate to a single shape to avoid overlapping sampling
  vect_buff <- terra::aggregate(vect_buff)
  
  # Sample points ONLY within this buffer
  bg_coords <- terra::spatSample(env_stack, size=n_bg, method="random", 
                                 x=vect_buff, na.rm=TRUE, xy=TRUE, values=FALSE)
  
  return(bg_coords)
}

# 2. PARAMETER TUNING
get_best_params <- function(occ_df, env_stack) {
  if(inherits(occ_df, "SpatRaster")) stop("occ_df input is a Raster!")
  
  occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  
  # Clean points falling in NA
  pt_vals <- terra::extract(env_stack, occ_coords)
  if("ID" %in% names(pt_vals)) pt_vals$ID <- NULL
  
  valid_rows <- complete.cases(pt_vals)
  occ_coords <- occ_coords[valid_rows, ]
  
  if(nrow(occ_coords) < 5) return(NULL) 
  
  # --- DEFENSIBLE FIX: BUFFERED BACKGROUND ---
  bg_coords <- get_buffered_background(occ_coords, env_stack)
  
  # If buffer fails (too small), fall back to global, but warn
  if(nrow(bg_coords) < 100) {
    bg_coords <- terra::spatSample(env_stack, size=10000, method="random", 
                                   na.rm=TRUE, xy=TRUE, values=FALSE)
  }
  
  tryCatch({
    eval_res <- ENMeval::ENMevaluate(
      occs = occ_coords,
      bg = bg_coords,
      envs = env_stack,
      algorithm = "maxnet",
      partitions = "randomkfold",
      partition.settings = list(kfolds = 5),
      tune.args = list(fc = c("L", "LQ", "H"), rm = c(1, 2, 5)),
      quiet = TRUE,
      parallel = FALSE
    )
    
    best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
    if(nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
    
    return(list(fc = as.character(best$fc), rm = as.numeric(best$rm)))
    
  }, error = function(e) { return(NULL) })
}

# 3. BOOTSTRAP WORKER
fit_bootstrap_worker <- function(occ_df, current_stack, future_stack_list=NULL, params, n_boot=10, debug_log=NULL) {
  
  log_debug <- function(msg) {
    if(!is.null(debug_log)) cat(paste0("[BOOTSTRAP] ", msg, "\n"), file=debug_log, append=TRUE)
  }
  log_debug(paste("Init: FC=", params$fc, " RM=", params$rm))
  
  sum_curr <- terra::rast(current_stack, nlyrs=1)
  terra::values(sum_curr) <- 0
  
  sum_fut_list <- list()
  if(!is.null(future_stack_list)) {
    for(n in names(future_stack_list)) {
      r <- terra::rast(future_stack_list[[n]], nlyrs=1)
      terra::values(r) <- 0
      sum_fut_list[[n]] <- r
    }
  }
  
  all_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  n_total <- nrow(all_coords)
  
  # --- DEFENSIBLE FIX: BUFFERED BACKGROUND FOR TRAINING ---
  bg_coords <- get_buffered_background(all_coords, current_stack)
  bg_data   <- terra::extract(current_stack, bg_coords)
  if("ID" %in% names(bg_data)) bg_data$ID <- NULL 
  
  successful_boots <- 0
  
  for(i in 1:n_boot) {
    tryCatch({
      # Resample
      train_idx <- sample(1:n_total, n_total, replace=TRUE)
      train_coords <- all_coords[train_idx, ]
      
      # Extract Presence
      pres_data <- terra::extract(current_stack, train_coords)
      if("ID" %in% names(pres_data)) pres_data$ID <- NULL 
      
      # Combine
      p_vec <- c(rep(1, nrow(pres_data)), rep(0, nrow(bg_data)))
      data_df <- rbind(pres_data, bg_data)
      
      # Clean NAs
      valid_rows <- complete.cases(data_df)
      data_df <- data_df[valid_rows, , drop=FALSE]
      p_vec <- p_vec[valid_rows]
      
      # Nuclear Option: Ensure Numeric
      data_df <- as.data.frame(lapply(data_df, as.numeric))
      p_vec   <- as.numeric(p_vec)
      
      if(sum(p_vec == 1) < 5) stop("Too few presence points.")
      
      # Fit
      mod <- maxnet::maxnet(p_vec, data_df, 
                            maxnet.formula(p_vec, data_df, classes=params$fc), 
                            regmult=params$rm)
      
      # Predict (PROJECT TO FULL EXTENT)
      # We trained on the buffer, but we project to the WHOLE current stack.
      pred_c <- terra::predict(current_stack, mod, type="logistic", na.rm=TRUE)
      sum_curr <- sum_curr + pred_c
      
      if(!is.null(future_stack_list)) {
        for(n in names(future_stack_list)) {
          pred_f <- terra::predict(future_stack_list[[n]], mod, type="logistic", na.rm=TRUE)
          sum_fut_list[[n]] <- sum_fut_list[[n]] + pred_f
        }
      }
      
      successful_boots <- successful_boots + 1
      
    }, error = function(e) {
      log_debug(paste("Iter", i, "FAILED ->", e$message))
    })
  }
  
  if(successful_boots == 0) {
    log_debug("CRITICAL: All bootstrap iterations failed.")
    stop("All bootstrap iterations failed.")
  }
  
  mean_curr <- sum_curr / successful_boots
  mean_fut_list <- list()
  if(!is.null(future_stack_list)) {
    for(n in names(future_stack_list)) {
      mean_fut_list[[n]] <- sum_fut_list[[n]] / successful_boots
    }
  }
  
  return(list(current = mean_curr, future = mean_fut_list))
}

get_biotic_layer <- function(fish_sp, host_stack, int_mat) {
  fish_clean <- gsub(" ", "_", fish_sp)
  if(!fish_clean %in% rownames(int_mat)) return(NULL)
  weights <- int_mat[fish_clean, ]
  available_hosts <- intersect(names(weights), names(host_stack))
  if(length(available_hosts) == 0) return(NULL)
  sub_stack <- host_stack[[available_hosts]]
  sub_weights <- as.numeric(weights[available_hosts])
  biotic_layer <- terra::weighted.mean(sub_stack, w=sub_weights, na.rm=TRUE)
  names(biotic_layer) <- "biotic_suitability"
  return(biotic_layer)
}
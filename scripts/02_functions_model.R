# scripts/02_functions_model.R
# ------------------------------------------------------------------------------
# MODELING WRAPPERS
# ------------------------------------------------------------------------------

get_best_params <- function(occ_df, env_stack, bg_coords, use_spatial, tune_args) {
  occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  part_method <- if(use_spatial) "block" else "randomkfold"
  
  tryCatch({
    eval_res <- tryCatch({
      ENMeval::ENMevaluate(occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
                           algorithm = "maxnet", partitions = part_method, 
                           tune.args = tune_args, quiet = TRUE, parallel = FALSE)
    }, error = function(e) {
      return(ENMeval::ENMevaluate(occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
                                  algorithm = "maxnet", partitions = "randomkfold", partition.settings = list(kfolds=5),
                                  tune.args = tune_args, quiet = TRUE, parallel = FALSE))
    })
    best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
    if(nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
    return(list(fc = as.character(best$fc), rm = as.numeric(best$rm), method=eval_res@partition.method))
  }, error = function(e) { return(NULL) })
}

# UPDATED: Saves RDS per iteration for Post-Analysis
fit_bootstrap_worker <- function(occ_df, current_stack, future_stack_list=NULL, bg_coords, params, n_boot=10, 
                                 sp_name, model_type, output_dir, debug_log=NULL) {
  
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
  
  # LOOP
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
      
      # SAVE ITERATION STATS (For LMM Analysis)
      iter_stats <- data.frame(species=sp_name, model=model_type, iter=i, auc=e@auc, tss=max(e@TPR + e@TNR - 1))
      stats_list[[i]] <- iter_stats
      
      # Save individual model object if needed? (Usually too big, stick to stats)
      saveRDS(iter_stats, file.path(output_dir, "models", paste0(sp_name, "_", model_type, "_iter", i, ".rds")))
      
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
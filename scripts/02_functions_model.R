# scripts/02_functions_model.R
# ------------------------------------------------------------------------------
# MODELING WRAPPERS (Tuning + Resumable Bootstrap + Ensemble Mean/SD)
# ------------------------------------------------------------------------------

get_best_params <- function(occ_df, env_stack, bg_coords, use_spatial, tune_args, seed=42) {
  occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  
  # Reproducible partitions
  set.seed(seed)
  
  part_method <- if(use_spatial) "block" else "randomkfold"
  
  tryCatch({
    eval_res <- tryCatch({
      ENMeval::ENMevaluate(
        occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
        algorithm = "maxnet", partitions = part_method,
        tune.args = tune_args, quiet = TRUE, parallel = FALSE
      )
    }, error = function(e) {
      ENMeval::ENMevaluate(
        occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
        algorithm = "maxnet", partitions = "randomkfold", partition.settings = list(kfolds=5),
        tune.args = tune_args, quiet = TRUE, parallel = FALSE
      )
    })
    
    best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
    if(nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
    
    list(fc = as.character(best$fc), rm = as.numeric(best$rm), method = eval_res@partition.method)
  }, error = function(e) NULL)
}

# ------------------------------------------------------------------------------
# RESUMABLE BOOTSTRAP WORKER
# - Skips iterations whose "iter_###.ok" exists
# - Allows increasing n_boot later: run 2 now, change to 10 later, runs 3..10
# - Saves model objects per iter, per-iter stats, and sum/sum_sq rasters
#
# Robustness fixes:
# - block-safe sum raster init (terra::init + mask), no huge values<-0
# - defensive hasValues() checks after predict() to fail loudly if terra temp/mem breaks
# ------------------------------------------------------------------------------

fit_bootstrap_worker <- function(occ_df, current_stack, future_stack_list=NULL, bg_coords, params, n_boot=10,
                                 sp_name, model_type, output_dir, debug_log=NULL, stage_dir=NULL) {
  
  log_debug <- function(msg) {
    if(!is.null(debug_log)) {
      ts <- format(Sys.time(), "%H:%M:%S")
      cat(paste0("[", ts, "] ", msg, "\n"), file=debug_log, append=TRUE)
    }
  }
  
  # ---- stage cache dir (per species + model type) ----
  if(is.null(stage_dir)) {
    stage_dir <- file.path(output_dir, "stage_cache", model_type, sp_name)
  }
  dir.create(stage_dir, recursive=TRUE, showWarnings=FALSE)
  
  # ---- sanitize params (works whether params is list or tibble) ----
  fc <- tolower(as.character(params$fc[1]))
  rm <- as.numeric(params$rm[1])
  
  # ---- helpers: atomic-ish writes (temp then rename) ----
  atomic_write_lines <- function(text, path) {
    tmp <- paste0(path, ".tmp")
    writeLines(text, tmp)
    if(file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  
  atomic_save_rds <- function(obj, path) {
    tmp <- paste0(path, ".tmp")
    saveRDS(obj, tmp)
    if(file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  
  atomic_write_raster <- function(r, path) {
    tmp <- paste0(path, ".tmp.tif")
    terra::writeRaster(r, tmp, overwrite=TRUE)
    if(file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  
  # ---- anonymize vars ----
  original_names <- names(current_stack)
  safe_names <- paste0("v", sprintf("%02d", 1:length(original_names)))
  names(current_stack) <- safe_names
  
  if(!is.null(future_stack_list)) {
    for(n in names(future_stack_list)) {
      if(all(original_names %in% names(future_stack_list[[n]]))) {
        future_stack_list[[n]] <- future_stack_list[[n]][[original_names]]
        names(future_stack_list[[n]]) <- safe_names
      } else {
        future_stack_list[[n]] <- NULL
      }
    }
    # drop NULL scenarios
    future_stack_list <- future_stack_list[!vapply(future_stack_list, is.null, logical(1))]
    if(length(future_stack_list) == 0) future_stack_list <- NULL
  }
  
  # ---- extract presence/background env values ----
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
  
  # ---- file paths ----
  sum_curr_path    <- file.path(stage_dir, "sum_current.tif")
  sum_sq_curr_path <- file.path(stage_dir, "sum_sq_current.tif")
  progress_path    <- file.path(stage_dir, "progress.rds")
  
  # future sum paths (sanitized scenario names)
  scen_key <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
  sum_fut_paths <- list()
  sum_sq_fut_paths <- list()
  if(!is.null(future_stack_list)) {
    for(sc in names(future_stack_list)) {
      k <- scen_key(sc)
      sum_fut_paths[[sc]]    <- file.path(stage_dir, paste0("sum_future_", k, ".tif"))
      sum_sq_fut_paths[[sc]] <- file.path(stage_dir, paste0("sum_sq_future_", k, ".tif"))
    }
  }
  
  iter_ok_path <- function(i) file.path(stage_dir, sprintf("iter_%03d.ok", i))
  
  model_path <- function(i) file.path(output_dir, "models", paste0(sp_name, "_", model_type, "_iter", sprintf("%03d", i), "_model.rds"))
  iter_stats_path <- function(i) file.path(output_dir, "models_stats", paste0(sp_name, "_", model_type, "_iter", sprintf("%03d", i), ".csv"))
  
  dir.create(file.path(output_dir, "models"),       recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(output_dir, "models_stats"), recursive=TRUE, showWarnings=FALSE)
  
  # ---- initialize or load progress + sums ----
  progress <- if(file.exists(progress_path)) {
    tryCatch(readRDS(progress_path), error=function(e) list(sum_completed=integer(0)))
  } else {
    list(sum_completed=integer(0))
  }
  
  # infer completed iters from ok markers (more robust than relying on RDS)
  ok_files <- list.files(stage_dir, pattern="^iter_[0-9]{3}\\.ok$", full.names=FALSE)
  ok_iters <- integer(0)
  if(length(ok_files) > 0) {
    ok_iters <- as.integer(sub("^iter_([0-9]{3})\\.ok$", "\\1", ok_files))
    ok_iters <- ok_iters[!is.na(ok_iters)]
  }
  ok_iters <- sort(unique(ok_iters[ok_iters >= 1 & ok_iters <= n_boot]))
  
  sum_completed <- sort(unique(as.integer(progress$sum_completed)))
  sum_completed <- sum_completed[sum_completed >= 1 & sum_completed <= n_boot]
  
  # block-safe zero raster preserving NA mask
  zero_like <- function(r) {
    z <- terra::init(terra::rast(r), fun = 0)
    z <- terra::mask(z, r)  # keep NA outside study area
    z
  }
  
  # load sums if they exist; otherwise rebuild from scratch if ok_iters exist
  need_rebuild <- FALSE
  if(file.exists(sum_curr_path) && file.exists(sum_sq_curr_path)) {
    sum_curr    <- terra::rast(sum_curr_path)
    sum_sq_curr <- terra::rast(sum_sq_curr_path)
  } else {
    sum_curr    <- zero_like(current_stack[[1]])
    sum_sq_curr <- zero_like(current_stack[[1]])
    need_rebuild <- length(ok_iters) > 0
  }
  
  sum_fut_list <- list()
  sum_sq_fut_list <- list()
  if(!is.null(future_stack_list)) {
    for(sc in names(future_stack_list)) {
      if(file.exists(sum_fut_paths[[sc]]) && file.exists(sum_sq_fut_paths[[sc]])) {
        sum_fut_list[[sc]]    <- terra::rast(sum_fut_paths[[sc]])
        sum_sq_fut_list[[sc]] <- terra::rast(sum_sq_fut_paths[[sc]])
      } else {
        sum_fut_list[[sc]]    <- zero_like(future_stack_list[[sc]][[1]])
        sum_sq_fut_list[[sc]] <- zero_like(future_stack_list[[sc]][[1]])
        if(length(ok_iters) > 0) need_rebuild <- TRUE
      }
    }
  }
  
  # rebuild sums if needed (rare: after crash before sums were written)
  if(need_rebuild && length(ok_iters) > 0) {
    log_debug(paste("Rebuilding sums from existing OK models:", paste(ok_iters, collapse=",")))
    sum_curr    <- zero_like(current_stack[[1]])
    sum_sq_curr <- zero_like(current_stack[[1]])
    if(!is.null(future_stack_list)) {
      for(sc in names(future_stack_list)) {
        sum_fut_list[[sc]]    <- zero_like(future_stack_list[[sc]][[1]])
        sum_sq_fut_list[[sc]] <- zero_like(future_stack_list[[sc]][[1]])
      }
    }
    sum_completed <- integer(0)
    
    for(i in ok_iters) {
      mp <- model_path(i)
      if(!file.exists(mp)) next
      mod <- tryCatch(readRDS(mp), error=function(e) NULL)
      if(is.null(mod)) next
      
      pred_c <- terra::predict(current_stack, mod, type="logistic", na.rm=TRUE)
      if(!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values during rebuild (temp/mem).")
      
      sum_curr    <- sum_curr + pred_c
      sum_sq_curr <- sum_sq_curr + (pred_c^2)
      
      if(!is.null(future_stack_list)) {
        for(sc in names(future_stack_list)) {
          pred_f <- terra::predict(future_stack_list[[sc]], mod, type="logistic", na.rm=TRUE)
          if(!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values (future=", sc, ") during rebuild."))
          sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
          sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
        }
      }
      
      sum_completed <- sort(unique(c(sum_completed, i)))
    }
    
    atomic_write_raster(sum_curr,    sum_curr_path)
    atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
    if(!is.null(future_stack_list)) {
      for(sc in names(future_stack_list)) {
        atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
        atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
      }
    }
    atomic_save_rds(list(sum_completed=sum_completed), progress_path)
  }
  
  # ---- run missing iterations ----
  for(i in seq_len(n_boot)) {
    
    # already completed and summed → skip
    if(file.exists(iter_ok_path(i)) && (i %in% sum_completed)) {
      next
    }
    
    # already OK but not in sums → incorporate by re-predicting from saved model
    if(file.exists(iter_ok_path(i)) && !(i %in% sum_completed)) {
      mp <- model_path(i)
      if(file.exists(mp)) {
        mod <- tryCatch(readRDS(mp), error=function(e) NULL)
        if(!is.null(mod)) {
          pred_c <- terra::predict(current_stack, mod, type="logistic", na.rm=TRUE)
          if(!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values while summing missing OK iter.")
          
          sum_curr    <- sum_curr + pred_c
          sum_sq_curr <- sum_sq_curr + (pred_c^2)
          
          if(!is.null(future_stack_list)) {
            for(sc in names(future_stack_list)) {
              pred_f <- terra::predict(future_stack_list[[sc]], mod, type="logistic", na.rm=TRUE)
              if(!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values (future=", sc, ") while summing missing OK iter."))
              sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
              sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
            }
          }
          
          sum_completed <- sort(unique(c(sum_completed, i)))
          atomic_write_raster(sum_curr,    sum_curr_path)
          atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
          if(!is.null(future_stack_list)) {
            for(sc in names(future_stack_list)) {
              atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
              atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
            }
          }
          atomic_save_rds(list(sum_completed=sum_completed), progress_path)
          log_debug(paste(model_type, "| Iter", i, "was OK but not summed → summed now"))
          next
        }
      }
      # fallthrough: if model missing, rerun
    }
    
    # run iteration i
    set.seed(i)
    
    tryCatch({
      n_pres <- nrow(pres_data_all)
      if(n_pres < 5) stop("Too few presences after NA omission")
      
      train_idx <- sample(1:n_pres, size = round(0.75 * n_pres))
      test_idx  <- setdiff(1:n_pres, train_idx)
      if(length(train_idx) < 5 || length(test_idx) < 5) stop("Train/test split too small")
      
      train_p <- pres_data_all[train_idx, , drop=FALSE]
      test_p  <- pres_data_all[test_idx, , drop=FALSE]
      
      p_vec <- c(rep(1, nrow(train_p)), rep(0, nrow(bg_data)))
      data_df <- rbind(train_p, bg_data)
      
      mod <- maxnet::maxnet(
        p_vec, data_df,
        maxnet::maxnet.formula(p_vec, data_df, classes=fc),
        regmult=rm
      )
      attr(mod, "orig_names") <- original_names
      attr(mod, "safe_names") <- safe_names
      attr(mod, "fc") <- fc
      attr(mod, "rm") <- rm
      attr(mod, "iter_seed") <- i
      
      pred_test_p <- stats::predict(mod, test_p, type="logistic")
      eval_bg_idx <- sample(1:nrow(bg_data), min(1000, nrow(bg_data)), replace=TRUE)
      pred_test_bg <- stats::predict(mod, bg_data[eval_bg_idx, , drop=FALSE], type="logistic")
      e <- dismo::evaluate(p=as.vector(pred_test_p), a=as.vector(pred_test_bg))
      
      iter_stats <- data.frame(
        species=sp_name, model=model_type, iter=i,
        fc=fc, rm=rm,
        n_pres=n_pres, n_bg=nrow(bg_data),
        auc=e@auc, tss=max(e@TPR + e@TNR - 1),
        seed=i
      )
      
      # PROJECTIONS
      pred_c <- terra::predict(current_stack, mod, type="logistic", na.rm=TRUE)
      if(!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values (current).")
      
      sum_curr    <- sum_curr + pred_c
      sum_sq_curr <- sum_sq_curr + (pred_c^2)
      
      if(!is.null(future_stack_list)) {
        for(sc in names(future_stack_list)) {
          pred_f <- terra::predict(future_stack_list[[sc]], mod, type="logistic", na.rm=TRUE)
          if(!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values (future=", sc, ")."))
          sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
          sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
        }
      }
      
      # SAVE MODEL + ITER STATS
      atomic_save_rds(mod, model_path(i))
      readr::write_csv(iter_stats, iter_stats_path(i))
      
      # MARK OK + UPDATE SUM CHECKPOINTS
      atomic_write_lines("OK", iter_ok_path(i))
      sum_completed <- sort(unique(c(sum_completed, i)))
      
      atomic_write_raster(sum_curr,    sum_curr_path)
      atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
      if(!is.null(future_stack_list)) {
        for(sc in names(future_stack_list)) {
          atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
          atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
        }
      }
      atomic_save_rds(list(sum_completed=sum_completed), progress_path)
      
      log_debug(paste(model_type, "| Iter", i, "OK | AUC:", round(e@auc, 3), "| TSS:", round(iter_stats$tss, 3)))
    }, error=function(e) {
      log_debug(paste(model_type, "| Iter", i, "ERROR:", e$message))
    })
  }
  
  # ---- finalize mean/sd from sums ----
  N <- length(sum_completed)
  if(N == 0) stop("All bootstrap iterations failed (no completed iters).")
  
  calc_mean_sd <- function(sum_r, sum_sq_r, N) {
    mean_r <- sum_r / N
    if(N <= 1) {
      sd_r <- sum_r * 0
    } else {
      var_r <- (sum_sq_r - (sum_r^2)/N) / (N - 1)
      var_r <- terra::clamp(var_r, lower=0)
      sd_r  <- sqrt(var_r)
    }
    names(mean_r) <- "mean_prob"
    names(sd_r)   <- "sd_prob"
    list(mean=mean_r, sd=sd_r)
  }
  
  res_curr <- calc_mean_sd(sum_curr, sum_sq_curr, N)
  
  res_fut_list <- list()
  if(!is.null(future_stack_list)) {
    for(sc in names(future_stack_list)) {
      res_fut_list[[sc]] <- calc_mean_sd(sum_fut_list[[sc]], sum_sq_fut_list[[sc]], N)
    }
  }
  
  # ---- collect stats from per-iter files ----
  stats_rows <- list()
  for(i in sum_completed) {
    p <- iter_stats_path(i)
    if(file.exists(p)) {
      stats_rows[[length(stats_rows)+1]] <- tryCatch(readr::read_csv(p, show_col_types=FALSE), error=function(e) NULL)
    }
  }
  stats_df <- dplyr::bind_rows(stats_rows)
  
  list(current=res_curr, future=res_fut_list, stats=stats_df, completed=N, target=n_boot)
}

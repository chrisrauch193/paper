# scripts/02_functions_model.R
# ------------------------------------------------------------------------------
# MODELING WRAPPERS (Tuning + Resumable Bootstrap + Ensemble Mean/SD)
# ------------------------------------------------------------------------------

get_best_params <- function(occ_df, env_stack, bg_coords, use_spatial, tune_args, seed = 42) {
  occ_coords <- occ_df %>% dplyr::select(x, y) %>% as.matrix()
  set.seed(seed)
  
  part_method <- if (isTRUE(use_spatial)) "block" else "randomkfold"
  
  eval_res <- tryCatch({
    ENMeval::ENMevaluate(
      occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
      algorithm = "maxnet", partitions = part_method,
      tune.args = tune_args, quiet = TRUE, parallel = FALSE
    )
  }, error = function(e) {
    if (isTRUE(exists("STRICT_BLOCKCV", inherits = TRUE)) &&
        isTRUE(get("STRICT_BLOCKCV", inherits = TRUE)) &&
        part_method == "block") {
      stop("ENMevaluate block partitions failed and STRICT_BLOCKCV=TRUE: ", e$message)
    }
    
    ENMeval::ENMevaluate(
      occs = occ_coords, bg = as.matrix(bg_coords), envs = env_stack,
      algorithm = "maxnet", partitions = "randomkfold",
      partition.settings = list(kfolds = 5),
      tune.args = tune_args, quiet = TRUE, parallel = FALSE
    )
  })
  
  best <- eval_res@results %>% dplyr::filter(delta.AICc == 0) %>% dplyr::slice(1)
  if (nrow(best) == 0) best <- eval_res@results %>% dplyr::arrange(desc(auc.val.avg)) %>% dplyr::slice(1)
  
  list(fc = as.character(best$fc), rm = as.numeric(best$rm), method = eval_res@partition.method)
}

# ------------------------------------------------------------------------------
# RESUMABLE BOOTSTRAP WORKER (now: spatial block holdout evaluation)
# ------------------------------------------------------------------------------

fit_bootstrap_worker <- function(
    occ_df, current_stack, future_stack_list = NULL, bg_coords, params, n_boot = 10,
    sp_name, model_type, output_dir, debug_log = NULL, stage_dir = NULL
) {
  
  log_debug <- function(msg) {
    if (!is.null(debug_log)) {
      ts <- format(Sys.time(), "%H:%M:%S")
      cat(paste0("[", ts, "] ", msg, "\n"), file = debug_log, append = TRUE)
    }
  }
  
  if (is.null(stage_dir)) stage_dir <- file.path(output_dir, "stage_cache", model_type, sp_name)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  
  fc <- tolower(as.character(params$fc[1]))
  rm <- as.numeric(params$rm[1])
  
  atomic_write_lines <- function(text, path) {
    tmp <- paste0(path, ".tmp")
    writeLines(text, tmp)
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  atomic_save_rds <- function(obj, path) {
    tmp <- paste0(path, ".tmp")
    saveRDS(obj, tmp)
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  atomic_write_raster <- function(r, path) {
    tmp <- paste0(path, ".tmp.tif")
    terra::writeRaster(r, tmp, overwrite = TRUE)
    if (file.exists(path)) unlink(path)
    file.rename(tmp, path)
  }
  
  # ---- anonymize vars (but keep provenance in stats) ----
  original_names <- names(current_stack)
  safe_names <- paste0("v", sprintf("%02d", seq_along(original_names)))
  names(current_stack) <- safe_names
  
  if (!is.null(future_stack_list)) {
    for (n in names(future_stack_list)) {
      if (all(original_names %in% names(future_stack_list[[n]]))) {
        future_stack_list[[n]] <- future_stack_list[[n]][[original_names]]
        names(future_stack_list[[n]]) <- safe_names
      } else {
        future_stack_list[[n]] <- NULL
      }
    }
    future_stack_list <- future_stack_list[!vapply(future_stack_list, is.null, logical(1))]
    if (length(future_stack_list) == 0) future_stack_list <- NULL
  }
  
  # ---- extract env values and KEEP coords (needed for spatial folds) ----
  pres_xy <- as.matrix(occ_df %>% dplyr::select(x, y))
  bg_xy   <- if (inherits(bg_coords, "data.frame")) as.matrix(bg_coords[, c("x", "y")]) else as.matrix(bg_coords)
  
  pres_ext <- terra::extract(current_stack, pres_xy)
  bg_ext   <- terra::extract(current_stack, bg_xy)
  
  if ("ID" %in% names(pres_ext)) pres_ext$ID <- NULL
  if ("ID" %in% names(bg_ext))   bg_ext$ID   <- NULL
  
  pres_df <- cbind(x = pres_xy[, 1], y = pres_xy[, 2], pres_ext)
  bg_df   <- cbind(x = bg_xy[, 1],   y = bg_xy[, 2],   bg_ext)
  
  pres_df <- pres_df[stats::complete.cases(pres_df), , drop = FALSE]
  bg_df   <- bg_df[stats::complete.cases(bg_df), , drop = FALSE]
  
  if (nrow(pres_df) < 10) stop("Too few presences after NA omission")
  if (nrow(bg_df)   < 500) stop("Too few background points after NA omission")
  
  pred_cols <- safe_names
  keep_vars <- pred_cols[sapply(as.data.frame(bg_df[, pred_cols, drop = FALSE]), stats::var) > 0]
  
  if (length(keep_vars) < 1) stop("No valid predictors after variance filtering.")
  
  pres_xy2 <- as.matrix(pres_df[, c("x", "y")])
  bg_xy2   <- as.matrix(bg_df[, c("x", "y")])
  
  pres_data_all <- as.data.frame(lapply(pres_df[, keep_vars, drop = FALSE], as.numeric))
  bg_data_all   <- as.data.frame(lapply(bg_df[, keep_vars, drop = FALSE], as.numeric))
  
  current_stack <- current_stack[[keep_vars]]
  if (!is.null(future_stack_list)) {
    for (sc in names(future_stack_list)) future_stack_list[[sc]] <- future_stack_list[[sc]][[keep_vars]]
  }
  
  orig_used <- original_names[match(keep_vars, safe_names)]
  
  # ---- paths ----
  sum_curr_path    <- file.path(stage_dir, "sum_current.tif")
  sum_sq_curr_path <- file.path(stage_dir, "sum_sq_current.tif")
  progress_path    <- file.path(stage_dir, "progress.rds")
  
  scen_key <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
  sum_fut_paths <- list()
  sum_sq_fut_paths <- list()
  if (!is.null(future_stack_list)) {
    for (sc in names(future_stack_list)) {
      k <- scen_key(sc)
      sum_fut_paths[[sc]]    <- file.path(stage_dir, paste0("sum_future_", k, ".tif"))
      sum_sq_fut_paths[[sc]] <- file.path(stage_dir, paste0("sum_sq_future_", k, ".tif"))
    }
  }
  
  iter_ok_path <- function(i) file.path(stage_dir, sprintf("iter_%03d.ok", i))
  model_path <- function(i) file.path(output_dir, "models", paste0(sp_name, "_", model_type, "_iter", sprintf("%03d", i), "_model.rds"))
  iter_stats_path <- function(i) file.path(output_dir, "models_stats", paste0(sp_name, "_", model_type, "_iter", sprintf("%03d", i), ".csv"))
  iter_permimp_path <- function(i) file.path(output_dir, "models_stats", paste0(sp_name, "_", model_type, "_iter", sprintf("%03d", i), "_permimp.csv"))
  
  dir.create(file.path(output_dir, "models"),       recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "models_stats"), recursive = TRUE, showWarnings = FALSE)
  
  # ---- init/load progress + sums ----
  progress <- if (file.exists(progress_path)) {
    tryCatch(readRDS(progress_path), error = function(e) list(sum_completed = integer(0)))
  } else list(sum_completed = integer(0))
  
  ok_files <- list.files(stage_dir, pattern = "^iter_[0-9]{3}\\.ok$", full.names = FALSE)
  ok_iters <- integer(0)
  if (length(ok_files) > 0) {
    ok_iters <- as.integer(sub("^iter_([0-9]{3})\\.ok$", "\\1", ok_files))
    ok_iters <- ok_iters[!is.na(ok_iters)]
  }
  ok_iters <- sort(unique(ok_iters[ok_iters >= 1 & ok_iters <= n_boot]))
  
  sum_completed <- sort(unique(as.integer(progress$sum_completed)))
  sum_completed <- sum_completed[sum_completed >= 1 & sum_completed <= n_boot]
  
  zero_like <- function(r) {
    z <- terra::init(terra::rast(r), fun = 0)
    terra::mask(z, r)
  }
  
  need_rebuild <- FALSE
  if (file.exists(sum_curr_path) && file.exists(sum_sq_curr_path)) {
    sum_curr    <- terra::rast(sum_curr_path)
    sum_sq_curr <- terra::rast(sum_sq_curr_path)
  } else {
    sum_curr    <- zero_like(current_stack[[1]])
    sum_sq_curr <- zero_like(current_stack[[1]])
    need_rebuild <- length(ok_iters) > 0
  }
  
  sum_fut_list <- list()
  sum_sq_fut_list <- list()
  if (!is.null(future_stack_list)) {
    for (sc in names(future_stack_list)) {
      if (file.exists(sum_fut_paths[[sc]]) && file.exists(sum_sq_fut_paths[[sc]])) {
        sum_fut_list[[sc]]    <- terra::rast(sum_fut_paths[[sc]])
        sum_sq_fut_list[[sc]] <- terra::rast(sum_sq_fut_paths[[sc]])
      } else {
        sum_fut_list[[sc]]    <- zero_like(future_stack_list[[sc]][[1]])
        sum_sq_fut_list[[sc]] <- zero_like(future_stack_list[[sc]][[1]])
        if (length(ok_iters) > 0) need_rebuild <- TRUE
      }
    }
  }
  
  if (need_rebuild && length(ok_iters) > 0) {
    log_debug(paste("Rebuilding sums from existing OK models:", paste(ok_iters, collapse = ",")))
    sum_curr    <- zero_like(current_stack[[1]])
    sum_sq_curr <- zero_like(current_stack[[1]])
    if (!is.null(future_stack_list)) {
      for (sc in names(future_stack_list)) {
        sum_fut_list[[sc]]    <- zero_like(future_stack_list[[sc]][[1]])
        sum_sq_fut_list[[sc]] <- zero_like(future_stack_list[[sc]][[1]])
      }
    }
    sum_completed <- integer(0)
    
    for (i in ok_iters) {
      mp <- model_path(i)
      if (!file.exists(mp)) next
      mod <- tryCatch(readRDS(mp), error = function(e) NULL)
      if (is.null(mod)) next
      
      pred_c <- terra::predict(current_stack, mod, type = "logistic", na.rm = TRUE)
      if (!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values during rebuild (current).")
      
      sum_curr    <- sum_curr + pred_c
      sum_sq_curr <- sum_sq_curr + (pred_c^2)
      
      if (!is.null(future_stack_list)) {
        for (sc in names(future_stack_list)) {
          pred_f <- terra::predict(future_stack_list[[sc]], mod, type = "logistic", na.rm = TRUE)
          if (!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values during rebuild (future=", sc, ")."))
          sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
          sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
        }
      }
      
      sum_completed <- sort(unique(c(sum_completed, i)))
    }
    
    atomic_write_raster(sum_curr,    sum_curr_path)
    atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
    if (!is.null(future_stack_list)) {
      for (sc in names(future_stack_list)) {
        atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
        atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
      }
    }
    atomic_save_rds(list(sum_completed = sum_completed), progress_path)
  }
  
  # ---- spatial folds (after NA/var filtering) ----
  cv_method <- if (exists("CV_METHOD", inherits = TRUE)) get("CV_METHOD", inherits = TRUE) else "block_quadrant"
  cv_k      <- if (exists("CV_FOLDS",  inherits = TRUE)) get("CV_FOLDS",  inherits = TRUE) else 4L
  
  folds <- assign_spatial_folds(pres_xy2, bg_xy2, method = cv_method, k = cv_k)
  pres_fold <- folds$pres_fold
  bg_fold   <- folds$bg_fold
  
  min_p <- if (exists("MIN_TEST_PRES", inherits = TRUE)) get("MIN_TEST_PRES", inherits = TRUE) else 5L
  min_a <- if (exists("MIN_TEST_BG",   inherits = TRUE)) get("MIN_TEST_BG",   inherits = TRUE) else 200L
  
  tab_p <- table(pres_fold)
  tab_a <- table(bg_fold)
  
  eligible <- intersect(
    as.integer(names(tab_p[tab_p >= min_p])),
    as.integer(names(tab_a[tab_a >= min_a]))
  )
  eligible <- sort(unique(eligible))
  if (length(eligible) == 0) {
    msg <- "No eligible spatial folds (too few presences/bg in every fold)."
    if (isTRUE(exists("STRICT_BLOCKCV", inherits = TRUE)) && isTRUE(get("STRICT_BLOCKCV", inherits = TRUE))) stop(msg)
    log_debug(paste("WARN:", msg, "Falling back to random eval split."))
  }
  
  # ---- run iterations ----
  for (i in seq_len(n_boot)) {
    
    if (file.exists(iter_ok_path(i)) && (i %in% sum_completed)) next
    
    if (file.exists(iter_ok_path(i)) && !(i %in% sum_completed)) {
      mp <- model_path(i)
      if (file.exists(mp)) {
        mod <- tryCatch(readRDS(mp), error = function(e) NULL)
        if (!is.null(mod)) {
          pred_c <- terra::predict(current_stack, mod, type = "logistic", na.rm = TRUE)
          if (!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values while summing missing OK iter (current).")
          sum_curr    <- sum_curr + pred_c
          sum_sq_curr <- sum_sq_curr + (pred_c^2)
          
          if (!is.null(future_stack_list)) {
            for (sc in names(future_stack_list)) {
              pred_f <- terra::predict(future_stack_list[[sc]], mod, type = "logistic", na.rm = TRUE)
              if (!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values while summing missing OK iter (future=", sc, ")."))
              sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
              sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
            }
          }
          
          sum_completed <- sort(unique(c(sum_completed, i)))
          atomic_write_raster(sum_curr,    sum_curr_path)
          atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
          if (!is.null(future_stack_list)) {
            for (sc in names(future_stack_list)) {
              atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
              atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
            }
          }
          atomic_save_rds(list(sum_completed = sum_completed), progress_path)
          log_debug(paste(model_type, "| Iter", i, "OK but not summed → summed now"))
          next
        }
      }
    }
    
    set.seed(i)
    
    tryCatch({
      
      cv_assign <- if (exists("CV_ASSIGNMENT", inherits = TRUE)) get("CV_ASSIGNMENT", inherits = TRUE) else "rotate"
      test_fold <- NA_integer_
      
      if (length(eligible) > 0 && tolower(cv_method) != "random") {
        if (tolower(cv_assign) == "random") test_fold <- sample(eligible, 1)
        else test_fold <- eligible[((i - 1L) %% length(eligible)) + 1L]
      }
      
      if (!is.na(test_fold)) {
        train_p_idx_all <- which(pres_fold != test_fold)
        test_p_idx      <- which(pres_fold == test_fold)
        train_a_idx_all <- which(bg_fold   != test_fold)
        test_a_idx      <- which(bg_fold   == test_fold)
      } else {
        n_pres <- nrow(pres_data_all)
        train_p_idx_all <- sample(seq_len(n_pres), size = round(0.75 * n_pres))
        test_p_idx      <- setdiff(seq_len(n_pres), train_p_idx_all)
        
        n_bg <- nrow(bg_data_all)
        train_a_idx_all <- sample(seq_len(n_bg), size = round(0.75 * n_bg))
        test_a_idx      <- setdiff(seq_len(n_bg), train_a_idx_all)
        test_fold <- NA_integer_
      }
      
      if (length(train_p_idx_all) < 5 || length(test_p_idx) < min_p) stop("Train/test pres split too small")
      if (length(train_a_idx_all) < 200 || length(test_a_idx) < min_a) stop("Train/test bg split too small")
      
      boot_pres <- isTRUE(exists("BOOTSTRAP_PRES_WITH_REPLACEMENT", inherits = TRUE)) &&
        isTRUE(get("BOOTSTRAP_PRES_WITH_REPLACEMENT", inherits = TRUE))
      boot_bg   <- isTRUE(exists("BOOTSTRAP_BG_WITH_REPLACEMENT", inherits = TRUE)) &&
        isTRUE(get("BOOTSTRAP_BG_WITH_REPLACEMENT", inherits = TRUE))
      
      train_p_idx <- if (boot_pres) sample(train_p_idx_all, length(train_p_idx_all), replace = TRUE) else train_p_idx_all
      train_a_idx <- if (boot_bg)   sample(train_a_idx_all, length(train_a_idx_all), replace = TRUE) else train_a_idx_all
      
      train_p <- pres_data_all[train_p_idx, , drop = FALSE]
      train_a <- bg_data_all[train_a_idx,   , drop = FALSE]
      
      test_p  <- pres_data_all[test_p_idx, , drop = FALSE]
      test_a0 <- bg_data_all[test_a_idx,   , drop = FALSE]
      
      eval_bg_max <- if (exists("EVAL_BG_MAX", inherits = TRUE)) get("EVAL_BG_MAX", inherits = TRUE) else 5000L
      if (nrow(test_a0) > eval_bg_max) {
        test_a <- test_a0[sample(seq_len(nrow(test_a0)), eval_bg_max, replace = FALSE), , drop = FALSE]
      } else test_a <- test_a0
      
      p_vec <- c(rep(1, nrow(train_p)), rep(0, nrow(train_a)))
      data_df <- rbind(train_p, train_a)
      
      mod <- maxnet::maxnet(
        p_vec, data_df,
        maxnet::maxnet.formula(p_vec, data_df, classes = fc),
        regmult = rm
      )
      
      attr(mod, "orig_names") <- orig_used
      attr(mod, "safe_names") <- keep_vars
      attr(mod, "fc") <- fc
      attr(mod, "rm") <- rm
      attr(mod, "iter_seed") <- i
      attr(mod, "cv_method") <- folds$method
      attr(mod, "cv_k")      <- folds$k
      attr(mod, "test_fold") <- test_fold
      
      pred_test_p  <- stats::predict(mod, test_p, type = "logistic")
      pred_test_bg <- stats::predict(mod, test_a, type = "logistic")
      
      e <- dismo::evaluate(p = as.vector(pred_test_p), a = as.vector(pred_test_bg))
      
      # threshold @ maxTSS
      thr_maxTSS <- NA_real_
      tss_val <- NA_real_
      tss_vec <- e@TPR + e@TNR - 1
      imax <- which.max(tss_vec)
      if (length(imax) == 1 && is.finite(tss_vec[imax])) {
        tss_val <- tss_vec[imax]
        thr_maxTSS <- e@t[imax]
      } else {
        tss_val <- max(tss_vec, na.rm = TRUE)
      }
      
      # confusion + calibration-ish
      eps <- 1e-15
      y <- c(rep(1, length(pred_test_p)), rep(0, length(pred_test_bg)))
      p <- c(as.numeric(pred_test_p), as.numeric(pred_test_bg))
      thr <- if (is.finite(thr_maxTSS)) thr_maxTSS else 0.5
      yhat <- as.integer(p >= thr)
      
      tp <- sum(y == 1 & yhat == 1)
      fn <- sum(y == 1 & yhat == 0)
      tn <- sum(y == 0 & yhat == 0)
      fp <- sum(y == 0 & yhat == 1)
      
      sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      bal_acc <- if (is.finite(sens) && is.finite(spec)) (sens + spec) / 2 else NA_real_
      
      brier <- mean((y - p)^2)
      logloss <- -mean(y * log(pmax(p, eps)) + (1 - y) * log(pmax(1 - p, eps)))
      
      boyce_bins <- if (exists("BOYCE_N_BINS", inherits = TRUE)) get("BOYCE_N_BINS", inherits = TRUE) else 20L
      boyce <- boyce_index_values(pred_test_p, pred_test_bg, n_bins = boyce_bins)
      
      iter_stats <- data.frame(
        species = sp_name,
        model   = model_type,
        iter    = i,
        fc      = fc,
        rm      = rm,
        
        cv_method = folds$method,
        cv_k      = folds$k,
        test_fold = test_fold,
        
        n_pres_total = nrow(pres_data_all),
        n_bg_total   = nrow(bg_data_all),
        n_pres_train = nrow(train_p),
        n_bg_train   = nrow(train_a),
        n_pres_test  = nrow(test_p),
        n_bg_test    = nrow(test_a),
        
        auc = e@auc,
        tss = as.numeric(tss_val),
        threshold_maxTSS = as.numeric(thr_maxTSS),
        
        sensitivity = as.numeric(sens),
        specificity = as.numeric(spec),
        balanced_accuracy = as.numeric(bal_acc),
        
        tp = as.integer(tp), fp = as.integer(fp), tn = as.integer(tn), fn = as.integer(fn),
        
        brier = as.numeric(brier),
        logloss = as.numeric(logloss),
        boyce = as.numeric(boyce),
        
        predictors = paste(orig_used, collapse = ";"),
        seed = i
      )
      
      if (isTRUE(exists("DO_PERM_IMPORTANCE", inherits = TRUE)) && isTRUE(get("DO_PERM_IMPORTANCE", inherits = TRUE))) {
        test_df_all <- rbind(test_p, test_a)
        base_auc <- e@auc
        pimps <- lapply(colnames(test_df_all), function(v) {
          tmp <- test_df_all
          tmp[[v]] <- sample(tmp[[v]])
          pp <- stats::predict(mod, tmp[seq_len(nrow(test_p)), , drop = FALSE], type = "logistic")
          aa <- stats::predict(mod, tmp[(nrow(test_p) + 1):nrow(tmp), , drop = FALSE], type = "logistic")
          ee <- dismo::evaluate(p = as.vector(pp), a = as.vector(aa))
          data.frame(var_safe = v, var_orig = orig_used[match(v, keep_vars)], auc_drop = base_auc - ee@auc)
        })
        pimps <- do.call(rbind, pimps)
        readr::write_csv(pimps, iter_permimp_path(i))
      }
      
      # Projections for ensemble
      pred_c <- terra::predict(current_stack, mod, type = "logistic", na.rm = TRUE)
      if (!terra::hasValues(pred_c)) stop("terra::predict returned raster with no values (current).")
      
      sum_curr    <- sum_curr + pred_c
      sum_sq_curr <- sum_sq_curr + (pred_c^2)
      
      if (!is.null(future_stack_list)) {
        for (sc in names(future_stack_list)) {
          pred_f <- terra::predict(future_stack_list[[sc]], mod, type = "logistic", na.rm = TRUE)
          if (!terra::hasValues(pred_f)) stop(paste0("terra::predict returned raster with no values (future=", sc, ")."))
          sum_fut_list[[sc]]    <- sum_fut_list[[sc]] + pred_f
          sum_sq_fut_list[[sc]] <- sum_sq_fut_list[[sc]] + (pred_f^2)
        }
      }
      
      atomic_save_rds(mod, model_path(i))
      readr::write_csv(iter_stats, iter_stats_path(i))
      
      atomic_write_lines("OK", iter_ok_path(i))
      sum_completed <- sort(unique(c(sum_completed, i)))
      
      atomic_write_raster(sum_curr,    sum_curr_path)
      atomic_write_raster(sum_sq_curr, sum_sq_curr_path)
      if (!is.null(future_stack_list)) {
        for (sc in names(future_stack_list)) {
          atomic_write_raster(sum_fut_list[[sc]],    sum_fut_paths[[sc]])
          atomic_write_raster(sum_sq_fut_list[[sc]], sum_sq_fut_paths[[sc]])
        }
      }
      atomic_save_rds(list(sum_completed = sum_completed), progress_path)
      
      log_debug(paste(model_type, "| Iter", i, "OK | fold:", test_fold, "| AUC:", round(e@auc, 3), "| TSS:", round(iter_stats$tss, 3)))
      
    }, error = function(e) {
      log_debug(paste(model_type, "| Iter", i, "ERROR:", e$message))
    })
  }
  
  N <- length(sum_completed)
  if (N == 0) stop("All bootstrap iterations failed (no completed iters).")
  
  calc_mean_sd <- function(sum_r, sum_sq_r, N) {
    mean_r <- sum_r / N
    if (N <= 1) {
      sd_r <- sum_r * 0
    } else {
      var_r <- (sum_sq_r - (sum_r^2) / N) / (N - 1)
      var_r <- terra::clamp(var_r, lower = 0)
      sd_r  <- sqrt(var_r)
    }
    names(mean_r) <- "mean_prob"
    names(sd_r)   <- "sd_prob"
    list(mean = mean_r, sd = sd_r)
  }
  
  res_curr <- calc_mean_sd(sum_curr, sum_sq_curr, N)
  
  res_fut_list <- list()
  if (!is.null(future_stack_list)) {
    for (sc in names(future_stack_list)) {
      res_fut_list[[sc]] <- calc_mean_sd(sum_fut_list[[sc]], sum_sq_fut_list[[sc]], N)
    }
  }
  
  stats_rows <- list()
  for (i in sum_completed) {
    p <- iter_stats_path(i)
    if (file.exists(p)) {
      stats_rows[[length(stats_rows) + 1]] <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
    }
  }
  stats_df <- dplyr::bind_rows(stats_rows)
  
  list(current = res_curr, future = res_fut_list, stats = stats_df, completed = N, target = n_boot)
}

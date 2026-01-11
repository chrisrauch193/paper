# scripts/03_pipeline_runner.R
# ------------------------------------------------------------------------------
# FINAL SDM PIPELINE: GOLD STANDARD (3-Stage Fish + Future Biotic)
# ------------------------------------------------------------------------------

rm(list = ls())
gc()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  foreach, doParallel,
  terra, dplyr, readr, ENMeval, maxnet, dismo
)

source("scripts/00_config.R")
source("scripts/01_functions_core.R")
source("scripts/02_functions_model.R")

THREAD_ENV <- c(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS    = "1",
  GDAL_NUM_THREADS       = "1"
)
do.call(Sys.setenv, as.list(THREAD_ENV))

set.seed(42)

# 3. DIRECTORIES ---------------------------------------------------------------
OUTPUT_ROOT <- file.path("outputs", RUN_ID)

DIRS <- c(
  "models_stats", "models_tuning", "models",
  "predictions/current/hosts",
  "predictions/current/fish_env_only",
  "predictions/current/fish_host_only",
  "predictions/current/fish_combined",
  "predictions/future",
  "logs",
  "status/hosts",
  "status/fish_env_only",
  "status/fish_host_only",
  "status/fish_combined",
  "tmp_terra",
  "occurrences/hosts",
  "occurrences/fish",
  "stage_cache/Host",
  "stage_cache/EnvOnly",
  "stage_cache/HostOnly",
  "stage_cache/Combined"
)
for (d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive = TRUE, showWarnings = FALSE)

MASTER_TD <- file.path(OUTPUT_ROOT, "tmp_terra", "master")
dir.create(MASTER_TD, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(
  tempdir  = MASTER_TD,
  todisk   = TRUE,
  memfrac  = 0.25,
  progress = 0
)

# 4. AUTO-WIPE -----------------------------------------------------------------
if (isTRUE(WIPE_PREDICTIONS)) {
  suppressWarnings({
    try(unlink(file.path(OUTPUT_ROOT, "models_stats"), recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "models"),       recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "models_tuning"),recursive = TRUE, force = TRUE), silent = TRUE)
    
    try(unlink(file.path(OUTPUT_ROOT, "predictions/current"), recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "predictions/future"),  recursive = TRUE, force = TRUE), silent = TRUE)
    
    try(unlink(file.path(OUTPUT_ROOT, "status"),      recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "tmp_terra"),   recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "stage_cache"), recursive = TRUE, force = TRUE), silent = TRUE)
    try(unlink(file.path(OUTPUT_ROOT, "occurrences"), recursive = TRUE, force = TRUE), silent = TRUE)
    
    for (d in DIRS) dir.create(file.path(OUTPUT_ROOT, d), recursive = TRUE, showWarnings = FALSE)
    
    dir.create(MASTER_TD, recursive = TRUE, showWarnings = FALSE)
    terra::terraOptions(
      tempdir  = MASTER_TD,
      todisk   = TRUE,
      memfrac  = 0.25,
      progress = 0
    )
  })
}

# 5. LOGGING -------------------------------------------------------------------
MASTER_LOG <- file.path(OUTPUT_ROOT, "pipeline_log.txt")
if (!file.exists(MASTER_LOG)) file.create(MASTER_LOG)
write_log(MASTER_LOG, paste0("=== PIPELINE START | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                             " | ID: ", RUN_ID, " | BG: ", BG_SAMPLING_METHOD, " ==="))

# 6. DATA LOADING --------------------------------------------------------------
cat("Loading Data...\n")

if (!file.exists(ENV_PATH)) stop("Env Raster missing at: ", ENV_PATH)
current_env_full <- terra::rast(ENV_PATH)

# Enforce predictor rule globally (PC1-4 + rugosity ONLY; PC5 excluded)
missing_preds <- setdiff(ENV_PREDICTORS, names(current_env_full))
if (length(missing_preds) > 0) stop("Missing required predictors in ENV_PATH: ", paste(missing_preds, collapse = ", "))

current_env <- current_env_full[[ENV_PREDICTORS]]
packed_env  <- terra::wrap(current_env)

write_log(MASTER_LOG, paste("Current env predictors:", paste(names(current_env), collapse = ", ")))

FUT_FILES <- list.files(FUT_DIR, full.names = TRUE, pattern = "\\.tif$")
HAS_FUT   <- length(FUT_FILES) > 0
if (!HAS_FUT) warning("No future layers found in: ", FUT_DIR)
scenario_names <- if (HAS_FUT) tools::file_path_sans_ext(basename(FUT_FILES)) else character(0)

if (HAS_FUT) {
  for (scen in scenario_names) {
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts"),          recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only"),  recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined"),  recursive = TRUE, showWarnings = FALSE)
  }
}

amph_occ <- as.data.frame(readr::read_csv("data/amph_occ_env_final_dataset.csv", show_col_types = FALSE))
anem_occ <- as.data.frame(readr::read_csv("data/anem_occ_env_final_dataset.csv", show_col_types = FALSE))

int_mat  <- utils::read.csv("data/interaction_matrix.csv", row.names = 1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

filter_species <- function(df, min_n = 10) {
  counts <- dplyr::count(df, species)
  valid  <- dplyr::filter(counts, n >= min_n) |> dplyr::pull(species)
  intersect(unique(df$species), valid)
}

valid_anemones <- filter_species(anem_occ, 10)
valid_fish     <- filter_species(amph_occ, 10)

anem_run_list <- if (!is.null(TARGET_HOSTS)) intersect(valid_anemones, TARGET_HOSTS) else valid_anemones
fish_run_list <- if (!is.null(TARGET_FISH))  intersect(valid_fish, TARGET_FISH)      else valid_fish

cat("Final Target Hosts:", length(anem_run_list), "\n")
cat("Final Target Fish:", length(fish_run_list), "\n")
write_log(MASTER_LOG, paste("Targets | Hosts:", length(anem_run_list), "| Fish:", length(fish_run_list)))

# Keep compatibility with your older naming
ENV_MODEL_VARS <- ENV_PREDICTORS

# 7. CLUSTER HELPERS -----------------------------------------------------------
prepare_future_stack <- function(fut_file, curr_stack, static_names) {
  fut_rast <- terra::rast(fut_file)
  
  wanted_static <- intersect(names(curr_stack), static_names)
  if (length(wanted_static) > 0) {
    missing_static <- setdiff(wanted_static, names(fut_rast))
    if (length(missing_static) > 0) {
      static_layers <- curr_stack[[missing_static]]
      if (!terra::compareGeom(fut_rast, static_layers, stopOnError = FALSE)) {
        static_layers <- terra::resample(static_layers, fut_rast, method = "near")
      }
      return(c(fut_rast, static_layers))
    }
  }
  fut_rast
}

start_cluster <- function(n, tag, output_root, thread_env) {
  outfile <- base::file.path(output_root, "logs", paste0("cluster_", tag, ".out"))
  base::dir.create(base::dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  cl <- NULL
  ok <- FALSE
  
  on.exit({
    if (!ok && !is.null(cl)) {
      base::message("Cluster setup failed; stopping cluster.")
      base::try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
  
  cl <- parallel::makeCluster(n, type = "PSOCK", outfile = outfile)
  doParallel::registerDoParallel(cl)
  
  parallel::clusterCall(
    cl,
    function(env, out_root) {
      do.call(Sys.setenv, as.list(env))
      
      if (requireNamespace("terra", quietly = TRUE)) {
        pid <- Sys.getpid()
        td  <- file.path(out_root, "tmp_terra", paste0("worker_", pid))
        dir.create(td, recursive = TRUE, showWarnings = FALSE)
        terra::terraOptions(
          tempdir  = td,
          todisk   = TRUE,
          memfrac  = 0.25,
          progress = 0
        )
      }
      
      Sys.getpid()
    },
    thread_env, output_root
  )
  
  parallel::clusterEvalQ(cl, {
    invisible(lapply(
      c("terra", "dplyr", "readr", "maxnet", "dismo", "ENMeval"),
      requireNamespace, quietly = TRUE
    ))
    gc()
    NULL
  })
  
  parallel::clusterCall(cl, base::Sys.getpid)
  
  ok <- TRUE
  cl
}

export_common <- function(cl, extra = character(0), envir = parent.frame()) {
  test <- base::try(parallel::clusterCall(cl, base::Sys.getpid), silent = TRUE)
  if (inherits(test, "try-error")) {
    base::stop(
      "Cluster is not alive when export_common() ran.\n",
      "Original error: ", base::conditionMessage(attr(test, "condition"))
    )
  }
  
  base_vars <- c(
    # data / paths
    "packed_env", "FUT_FILES", "scenario_names", "HAS_FUT",
    "STATIC_VARS", "ENV_MODEL_VARS", "ENV_PREDICTORS", "ENV_PCS",
    "amph_occ", "anem_occ", "int_mat", "OUTPUT_ROOT", "MASTER_LOG",
    
    # params
    "N_HOST_BOOT", "N_FISH_BOOT", "TUNE_ARGS",
    "USE_SPATIAL_THINNING", "USE_SPATIAL_TUNING", "BG_SAMPLING_METHOD",
    
    # CV / eval settings
    "CV_METHOD", "CV_FOLDS", "CV_ASSIGNMENT", "STRICT_BLOCKCV",
    "MIN_TEST_PRES", "MIN_TEST_BG", "EVAL_BG_MAX",
    "BOOTSTRAP_PRES_WITH_REPLACEMENT", "BOOTSTRAP_BG_WITH_REPLACEMENT",
    "DO_PERM_IMPORTANCE", "BOYCE_N_BINS",
    
    # BG settings
    "BG_N_BG", "BG_CAND_MULT", "BG_ALPHA", "BG_ENV_VARS", "BG_GEO_METRIC", "BG_BUFFER_M",
    
    # funcs
    "write_log", "thin_occurrences",
    "get_bias_corrected_background", "get_random_background",
    "get_biotic_layer",
    "get_best_params", "fit_bootstrap_worker",
    "prepare_future_stack",
    "read_status_file", "write_status_ok", "write_status_progress",
    
    # helpers for blockCV + boyce + haversine
    ".haversine_m", "make_block_quadrant_folds", "assign_spatial_folds", "boyce_index_values"
  )
  
  parallel::clusterExport(cl, unique(c(base_vars, extra)), envir = envir)
  invisible(TRUE)
}

FOREACH_OPTS <- list(preschedule = FALSE)

N_CORES_HOST <- min(N_CORES, max(1L, length(anem_run_list)))
N_CORES_FISH <- min(N_CORES, max(1L, length(fish_run_list)), 10L)

# ==============================================================================
# PHASE 1: HOSTS
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 1: HOSTS ---")

cl_hosts <- start_cluster(N_CORES_HOST, "hosts", OUTPUT_ROOT, THREAD_ENV)
export_common(cl_hosts)

host_results <- foreach::foreach(
  sp = anem_run_list,
  .packages = c("terra", "dplyr", "maxnet", "dismo", "readr", "ENMeval"),
  .errorhandling = "pass",
  .options.snow = FOREACH_OPTS
) %dopar% {
  
  sp_clean <- gsub(" ", "_", sp)
  
  done_file <- file.path(OUTPUT_ROOT, "status", "hosts", paste0(sp_clean, ".done"))
  st <- read_status_file(done_file)
  
  stats_file <- file.path(OUTPUT_ROOT, "models_stats",  paste0(sp_clean, "_Host_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_Host_params.csv"))
  pred_mean  <- file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, "_mean.tif"))
  pred_sd    <- file.path(OUTPUT_ROOT, "predictions/current/hosts", paste0(sp_clean, "_sd.tif"))
  
  if (st$status == "OK" && !is.na(st$n) && st$n >= N_HOST_BOOT &&
      file.exists(pred_mean) && file.exists(pred_sd) && file.exists(stats_file)) {
    write_log(MASTER_LOG, paste("SKIP Host:", sp_clean, "| already OK n=", st$n))
    return(list(species = sp_clean, stage = "Host", status = "SKIP"))
  }
  
  write_log(MASTER_LOG, paste("START Host:", sp_clean, "| target iters:", N_HOST_BOOT))
  sp_log <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, "_Host.log"))
  
  env_stack <- terra::unwrap(packed_env)  # already PC1-4 + rugosity
  
  out <- tryCatch({
    
    occ_out <- file.path(OUTPUT_ROOT, "occurrences", "hosts", paste0(sp_clean, "_occ_used.csv"))
    if (file.exists(occ_out)) {
      sp_dat <- readr::read_csv(occ_out, show_col_types = FALSE)
    } else {
      sp_dat <- dplyr::filter(anem_occ, species == sp)
      if (nrow(sp_dat) < 5) stop("Not enough occurrences")
      if (isTRUE(USE_SPATIAL_THINNING)) sp_dat <- thin_occurrences(sp_dat, env_stack)
      readr::write_csv(sp_dat, occ_out)
    }
    
    # Background
    if (BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
      bg_coords <- get_bias_corrected_background(
        sp_dat, env_stack,
        n_bg = BG_N_BG, alpha = BG_ALPHA, method = BG_SAMPLING_METHOD,
        env_cols = BG_ENV_VARS, cand_mult = BG_CAND_MULT, geo_metric = BG_GEO_METRIC,
        seed = 42
      )
    } else {
      bg_coords <- get_random_background(sp_dat, env_stack, n_bg = BG_N_BG, buffer_m = BG_BUFFER_M)
    }
    if (is.null(bg_coords) || nrow(bg_coords) < 50) stop("Background sampling failed/too small")
    
    # Futures (aligned to current env vars)
    future_stacks <- NULL
    if (isTRUE(HAS_FUT)) {
      future_stacks <- list()
      for (i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        f_stack <- prepare_future_stack(FUT_FILES[i], env_stack, STATIC_VARS)
        vars_needed <- names(env_stack)
        if (all(vars_needed %in% names(f_stack))) {
          future_stacks[[scen]] <- f_stack[[vars_needed]]
        }
      }
      if (length(future_stacks) == 0) future_stacks <- NULL
    }
    
    # Tuning (cached)
    if (file.exists(tune_file)) {
      params <- readr::read_csv(tune_file, show_col_types = FALSE)
      write_log(MASTER_LOG, paste("  > Host tuning cached:", sp_clean))
    } else {
      write_log(MASTER_LOG, paste("  > Tuning Host:", sp_clean))
      params <- get_best_params(sp_dat, env_stack, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if (is.null(params)) stop("Tuning failed")
      params$fc <- tolower(params$fc)
      readr::write_csv(as.data.frame(params), tune_file)
    }
    
    results <- fit_bootstrap_worker(
      occ_df = sp_dat,
      current_stack = env_stack,
      future_stack_list = future_stacks,
      bg_coords = bg_coords,
      params = params,
      n_boot = N_HOST_BOOT,
      sp_name = sp_clean,
      model_type = "Host",
      output_dir = OUTPUT_ROOT,
      debug_log = sp_log
    )
    
    readr::write_csv(results$stats, stats_file)
    terra::writeRaster(results$current$mean, pred_mean, overwrite = TRUE)
    terra::writeRaster(results$current$sd,   pred_sd,   overwrite = TRUE)
    
    if (!is.null(results$future) && length(results$future) > 0) {
      for (scen in names(results$future)) {
        terra::writeRaster(
          results$future[[scen]]$mean,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, "_mean.tif")),
          overwrite = TRUE
        )
        terra::writeRaster(
          results$future[[scen]]$sd,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts", paste0(sp_clean, "_sd.tif")),
          overwrite = TRUE
        )
      }
    }
    
    if (!is.null(results$completed) && results$completed >= N_HOST_BOOT) {
      write_status_ok(done_file, N_HOST_BOOT)
      write_log(MASTER_LOG, paste("FINISH Host:", sp_clean, "| OK n=", N_HOST_BOOT))
      list(species = sp_clean, stage = "Host", status = "OK", n = N_HOST_BOOT)
    } else {
      write_status_progress(done_file, results$completed, N_HOST_BOOT)
      write_log(MASTER_LOG, paste("FINISH Host:", sp_clean, "| PROGRESS n=", results$completed, "/", N_HOST_BOOT))
      list(species = sp_clean, stage = "Host", status = "PROGRESS", n = results$completed)
    }
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Host:", sp_clean, "->", e$message))
    list(species = sp_clean, stage = "Host", status = "ERROR", msg = e$message)
  })
  
  out
}

base::try(parallel::stopCluster(cl_hosts), silent = TRUE)

# ==============================================================================
# PREP HOST STACK (Current) for fish phases
# ==============================================================================
host_files <- list.files(
  file.path(OUTPUT_ROOT, "predictions/current/hosts"),
  full.names = TRUE, pattern = "_mean\\.tif$"
)
if (length(host_files) == 0) stop("CRITICAL: Phase 1 failed (no host rasters found).")

hosts_curr <- terra::rast(host_files)
names(hosts_curr) <- gsub("_mean$", "", tools::file_path_sans_ext(basename(host_files)))
packed_hosts_curr <- terra::wrap(hosts_curr)

# ==============================================================================
# PHASE 2: FISH — 3 PASSES
# ==============================================================================
write_log(MASTER_LOG, "--- PHASE 2: FISH (3 passes) ---")

cl_fish <- start_cluster(N_CORES_FISH, "fish", OUTPUT_ROOT, THREAD_ENV)
export_common(cl_fish, extra = c("packed_hosts_curr"))

# ------------------------------------------------------------------------------
# Pass A: ENV ONLY (PC1-4 + rugosity)  [MODEL A]
# ------------------------------------------------------------------------------
write_log(MASTER_LOG, "--- PHASE 2A: FISH ENV ONLY ---")

fish_env_results <- foreach::foreach(
  sp = fish_run_list,
  .packages = c("terra", "dplyr", "maxnet", "dismo", "readr", "ENMeval"),
  .errorhandling = "pass",
  .options.snow = FOREACH_OPTS
) %dopar% {
  
  sp_clean <- gsub(" ", "_", sp)
  
  done_file <- file.path(OUTPUT_ROOT, "status", "fish_env_only", paste0(sp_clean, ".done"))
  st <- read_status_file(done_file)
  
  sp_log    <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, "_EnvOnly.log"))
  
  stats_file <- file.path(OUTPUT_ROOT, "models_stats",  paste0(sp_clean, "_FishEnvOnly_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishEnvOnly_params.csv"))
  bg_file    <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_bg.rds"))
  
  out_mean <- file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", paste0(sp_clean, "_mean.tif"))
  out_sd   <- file.path(OUTPUT_ROOT, "predictions/current/fish_env_only", paste0(sp_clean, "_sd.tif"))
  
  if (st$status == "OK" && !is.na(st$n) && st$n >= N_FISH_BOOT &&
      file.exists(out_mean) && file.exists(out_sd) && file.exists(stats_file)) {
    write_log(MASTER_LOG, paste("SKIP Fish EnvOnly:", sp_clean, "| already OK n=", st$n))
    return(list(species = sp_clean, stage = "FishEnvOnly", status = "SKIP"))
  }
  
  write_log(MASTER_LOG, paste("START Fish EnvOnly:", sp_clean, "| target iters:", N_FISH_BOOT))
  env_stack_curr <- terra::unwrap(packed_env)
  
  out <- tryCatch({
    
    occ_out <- file.path(OUTPUT_ROOT, "occurrences", "fish", paste0(sp_clean, "_occ_used.csv"))
    if (file.exists(occ_out)) {
      sp_dat <- readr::read_csv(occ_out, show_col_types = FALSE)
    } else {
      sp_dat <- dplyr::filter(amph_occ, species == sp)
      if (nrow(sp_dat) < 5) stop("Not enough occurrences")
      if (isTRUE(USE_SPATIAL_THINNING)) sp_dat <- thin_occurrences(sp_dat, env_stack_curr)
      readr::write_csv(sp_dat, occ_out)
    }
    
    env_stack_model_a <- env_stack_curr[[ENV_PREDICTORS]]
    if (!all(c(ENV_PCS, "rugosity") %in% names(env_stack_model_a))) stop("EnvOnly stack missing required predictors.")
    if ("PC5" %in% names(env_stack_model_a)) stop("BUG: PC5 leaked into EnvOnly predictors.")
    
    # Background (cached; built from EnvOnly stack)
    if (file.exists(bg_file)) {
      bg_coords <- readRDS(bg_file)
    } else {
      if (BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
        bg_coords <- get_bias_corrected_background(
          sp_dat, env_stack_model_a,
          n_bg = BG_N_BG, alpha = BG_ALPHA, method = BG_SAMPLING_METHOD,
          env_cols = BG_ENV_VARS, cand_mult = BG_CAND_MULT, geo_metric = BG_GEO_METRIC,
          seed = 42
        )
      } else {
        bg_coords <- get_random_background(sp_dat, env_stack_model_a, n_bg = BG_N_BG, buffer_m = BG_BUFFER_M)
      }
      if (is.null(bg_coords) || nrow(bg_coords) < 50) stop("Background sampling failed/too small")
      saveRDS(bg_coords, bg_file)
    }
    
    # Futures (subset to EnvOnly predictors)
    future_stacks_env <- NULL
    if (isTRUE(HAS_FUT)) {
      future_stacks_env <- list()
      for (i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        f_stack <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
        if (all(names(env_stack_model_a) %in% names(f_stack))) {
          future_stacks_env[[scen]] <- f_stack[[names(env_stack_model_a)]]
        }
      }
      if (length(future_stacks_env) == 0) future_stacks_env <- NULL
    }
    
    # Tuning (cached)
    if (file.exists(tune_file)) {
      params_env <- readr::read_csv(tune_file, show_col_types = FALSE)
      write_log(MASTER_LOG, paste("  > EnvOnly tuning cached:", sp_clean))
    } else {
      params_env <- get_best_params(sp_dat, env_stack_model_a, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if (!is.null(params_env)) {
        params_env$fc <- tolower(params_env$fc)
        readr::write_csv(as.data.frame(params_env), tune_file)
      }
    }
    if (is.null(params_env)) stop("Tuning failed for EnvOnly")
    
    res_env <- fit_bootstrap_worker(
      occ_df = sp_dat,
      current_stack = env_stack_model_a,
      future_stack_list = future_stacks_env,
      bg_coords = bg_coords,
      params = params_env,
      n_boot = N_FISH_BOOT,
      sp_name = sp_clean,
      model_type = "EnvOnly",
      output_dir = OUTPUT_ROOT,
      debug_log = sp_log
    )
    
    readr::write_csv(res_env$stats, stats_file)
    terra::writeRaster(res_env$current$mean, out_mean, overwrite = TRUE)
    terra::writeRaster(res_env$current$sd,   out_sd,   overwrite = TRUE)
    
    if (!is.null(res_env$future) && length(res_env$future) > 0) {
      for (scen in names(res_env$future)) {
        terra::writeRaster(
          res_env$future[[scen]]$mean,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only", paste0(sp_clean, "_mean.tif")),
          overwrite = TRUE
        )
        terra::writeRaster(
          res_env$future[[scen]]$sd,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_env_only", paste0(sp_clean, "_sd.tif")),
          overwrite = TRUE
        )
      }
    }
    
    if (!is.null(res_env$completed) && res_env$completed >= N_FISH_BOOT) {
      write_status_ok(done_file, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish EnvOnly:", sp_clean, "| OK n=", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishEnvOnly", status = "OK", n = N_FISH_BOOT)
    } else {
      write_status_progress(done_file, res_env$completed, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish EnvOnly:", sp_clean, "| PROGRESS n=", res_env$completed, "/", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishEnvOnly", status = "PROGRESS", n = res_env$completed)
    }
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish EnvOnly:", sp_clean, "->", e$message))
    list(species = sp_clean, stage = "FishEnvOnly", status = "ERROR", msg = e$message)
  })
  
  out
}

# ------------------------------------------------------------------------------
# Pass B: HOST ONLY (rugosity + biotic)  [MODEL B]
# ------------------------------------------------------------------------------
write_log(MASTER_LOG, "--- PHASE 2B: FISH HOST ONLY ---")

fish_host_results <- foreach::foreach(
  sp = fish_run_list,
  .packages = c("terra", "dplyr", "maxnet", "dismo", "readr", "ENMeval"),
  .errorhandling = "pass",
  .options.snow = FOREACH_OPTS
) %dopar% {
  
  sp_clean <- gsub(" ", "_", sp)
  
  done_file <- file.path(OUTPUT_ROOT, "status", "fish_host_only", paste0(sp_clean, ".done"))
  st <- read_status_file(done_file)
  
  if (st$status == "SKIP_NO_BIOTIC") {
    write_log(MASTER_LOG, paste("SKIP Fish HostOnly:", sp_clean, "| status SKIP_NO_BIOTIC"))
    return(list(species = sp_clean, stage = "FishHostOnly", status = "SKIP_NO_BIOTIC"))
  }
  
  sp_log    <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, "_HostOnly.log"))
  
  stats_file <- file.path(OUTPUT_ROOT, "models_stats",  paste0(sp_clean, "_FishHostOnly_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishHostOnly_params.csv"))
  bg_file    <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_bg.rds"))
  
  out_mean <- file.path(OUTPUT_ROOT, "predictions/current/fish_host_only", paste0(sp_clean, "_mean.tif"))
  out_sd   <- file.path(OUTPUT_ROOT, "predictions/current/fish_host_only", paste0(sp_clean, "_sd.tif"))
  
  if (st$status == "OK" && !is.na(st$n) && st$n >= N_FISH_BOOT &&
      file.exists(out_mean) && file.exists(out_sd) && file.exists(stats_file)) {
    write_log(MASTER_LOG, paste("SKIP Fish HostOnly:", sp_clean, "| already OK n=", st$n))
    return(list(species = sp_clean, stage = "FishHostOnly", status = "SKIP"))
  }
  
  write_log(MASTER_LOG, paste("START Fish HostOnly:", sp_clean, "| target iters:", N_FISH_BOOT))
  
  env_stack_curr  <- terra::unwrap(packed_env)
  host_stack_curr <- terra::unwrap(packed_hosts_curr)
  
  out <- tryCatch({
    
    occ_out <- file.path(OUTPUT_ROOT, "occurrences", "fish", paste0(sp_clean, "_occ_used.csv"))
    if (file.exists(occ_out)) {
      sp_dat <- readr::read_csv(occ_out, show_col_types = FALSE)
    } else {
      sp_dat <- dplyr::filter(amph_occ, species == sp)
      if (nrow(sp_dat) < 5) stop("Not enough occurrences")
      if (isTRUE(USE_SPATIAL_THINNING)) sp_dat <- thin_occurrences(sp_dat, env_stack_curr)
      readr::write_csv(sp_dat, occ_out)
    }
    
    if (!("rugosity" %in% names(env_stack_curr))) stop("Missing rugosity layer for HostOnly")
    rug_layer <- env_stack_curr[["rugosity"]]
    
    biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat, debug_path = MASTER_LOG)
    if (is.null(biotic_curr)) {
      write_log(MASTER_LOG, paste("SKIP Fish HostOnly:", sp_clean, "-> no biotic layer"))
      writeLines("SKIP_NO_BIOTIC", done_file)
      return(list(species = sp_clean, stage = "FishHostOnly", status = "SKIP_NO_BIOTIC"))
    }
    
    host_only_stack <- c(rug_layer, biotic_curr)
    names(host_only_stack) <- c("rugosity", "biotic_suitability")
    if ("PC5" %in% names(host_only_stack)) stop("BUG: PC5 leaked into HostOnly predictors.")
    
    # Background (cached; built from EnvOnly env stack)
    if (file.exists(bg_file)) {
      bg_coords <- readRDS(bg_file)
    } else {
      env_stack_model_a <- env_stack_curr[[ENV_PREDICTORS]]
      if (BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
        bg_coords <- get_bias_corrected_background(
          sp_dat, env_stack_model_a,
          n_bg = BG_N_BG, alpha = BG_ALPHA, method = BG_SAMPLING_METHOD,
          env_cols = BG_ENV_VARS, cand_mult = BG_CAND_MULT, geo_metric = BG_GEO_METRIC,
          seed = 42
        )
      } else {
        bg_coords <- get_random_background(sp_dat, env_stack_model_a, n_bg = BG_N_BG, buffer_m = BG_BUFFER_M)
      }
      if (is.null(bg_coords) || nrow(bg_coords) < 50) stop("Background sampling failed/too small")
      saveRDS(bg_coords, bg_file)
    }
    
    # Futures: rugosity + future biotic from future host predictions
    future_stacks_hostonly <- NULL
    if (isTRUE(HAS_FUT)) {
      future_stacks_hostonly <- list()
      
      for (i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        
        f_env_raw <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
        if (!("rugosity" %in% names(f_env_raw))) next
        f_rug <- f_env_raw[["rugosity"]]
        
        host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
        host_files_fut <- list.files(host_dir_fut, full.names = TRUE, pattern = "_mean\\.tif$")
        if (length(host_files_fut) == 0) next
        
        h_stack_f <- terra::rast(host_files_fut)
        names(h_stack_f) <- gsub("_mean$", "", tools::file_path_sans_ext(basename(host_files_fut)))
        
        biotic_fut <- get_biotic_layer(sp, h_stack_f, int_mat)
        if (is.null(biotic_fut)) next
        
        f_comb <- c(f_rug, biotic_fut)
        names(f_comb) <- c("rugosity", "biotic_suitability")
        future_stacks_hostonly[[scen]] <- f_comb
      }
      
      if (length(future_stacks_hostonly) == 0) future_stacks_hostonly <- NULL
    }
    
    # Tuning (cached)
    if (file.exists(tune_file)) {
      params_host <- readr::read_csv(tune_file, show_col_types = FALSE)
      write_log(MASTER_LOG, paste("  > HostOnly tuning cached:", sp_clean))
    } else {
      params_host <- get_best_params(sp_dat, host_only_stack, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if (!is.null(params_host)) {
        params_host$fc <- tolower(params_host$fc)
        readr::write_csv(as.data.frame(params_host), tune_file)
      }
    }
    if (is.null(params_host)) stop("Tuning failed for HostOnly")
    
    res_ho <- fit_bootstrap_worker(
      occ_df = sp_dat,
      current_stack = host_only_stack,
      future_stack_list = future_stacks_hostonly,
      bg_coords = bg_coords,
      params = params_host,
      n_boot = N_FISH_BOOT,
      sp_name = sp_clean,
      model_type = "HostOnly",
      output_dir = OUTPUT_ROOT,
      debug_log = sp_log
    )
    
    readr::write_csv(res_ho$stats, stats_file)
    terra::writeRaster(res_ho$current$mean, out_mean, overwrite = TRUE)
    terra::writeRaster(res_ho$current$sd,   out_sd,   overwrite = TRUE)
    
    if (!is.null(res_ho$future) && length(res_ho$future) > 0) {
      for (scen in names(res_ho$future)) {
        terra::writeRaster(
          res_ho$future[[scen]]$mean,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only", paste0(sp_clean, "_mean.tif")),
          overwrite = TRUE
        )
        terra::writeRaster(
          res_ho$future[[scen]]$sd,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_host_only", paste0(sp_clean, "_sd.tif")),
          overwrite = TRUE
        )
      }
    }
    
    if (!is.null(res_ho$completed) && res_ho$completed >= N_FISH_BOOT) {
      write_status_ok(done_file, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish HostOnly:", sp_clean, "| OK n=", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishHostOnly", status = "OK", n = N_FISH_BOOT)
    } else {
      write_status_progress(done_file, res_ho$completed, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish HostOnly:", sp_clean, "| PROGRESS n=", res_ho$completed, "/", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishHostOnly", status = "PROGRESS", n = res_ho$completed)
    }
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish HostOnly:", sp_clean, "->", e$message))
    list(species = sp_clean, stage = "FishHostOnly", status = "ERROR", msg = e$message)
  })
  
  out
}

# ------------------------------------------------------------------------------
# Pass C: COMBINED (PC1-4 + rugosity + biotic)  [MODEL C]
# ------------------------------------------------------------------------------
write_log(MASTER_LOG, "--- PHASE 2C: FISH COMBINED ---")

fish_comb_results <- foreach::foreach(
  sp = fish_run_list,
  .packages = c("terra", "dplyr", "maxnet", "dismo", "readr", "ENMeval"),
  .errorhandling = "pass",
  .options.snow = FOREACH_OPTS
) %dopar% {
  
  sp_clean <- gsub(" ", "_", sp)
  
  done_file <- file.path(OUTPUT_ROOT, "status", "fish_combined", paste0(sp_clean, ".done"))
  st <- read_status_file(done_file)
  
  if (st$status == "SKIP_NO_BIOTIC") {
    write_log(MASTER_LOG, paste("SKIP Fish Combined:", sp_clean, "| status SKIP_NO_BIOTIC"))
    return(list(species = sp_clean, stage = "FishCombined", status = "SKIP_NO_BIOTIC"))
  }
  
  sp_log    <- file.path(OUTPUT_ROOT, "logs", paste0(sp_clean, "_Combined.log"))
  
  stats_file <- file.path(OUTPUT_ROOT, "models_stats",  paste0(sp_clean, "_FishCombined_stats.csv"))
  tune_file  <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_FishCombined_params.csv"))
  bg_file    <- file.path(OUTPUT_ROOT, "models_tuning", paste0(sp_clean, "_bg.rds"))
  
  out_mean <- file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, "_mean.tif"))
  out_sd   <- file.path(OUTPUT_ROOT, "predictions/current/fish_combined", paste0(sp_clean, "_sd.tif"))
  
  if (st$status == "OK" && !is.na(st$n) && st$n >= N_FISH_BOOT &&
      file.exists(out_mean) && file.exists(out_sd) && file.exists(stats_file)) {
    write_log(MASTER_LOG, paste("SKIP Fish Combined:", sp_clean, "| already OK n=", st$n))
    return(list(species = sp_clean, stage = "FishCombined", status = "SKIP"))
  }
  
  write_log(MASTER_LOG, paste("START Fish Combined:", sp_clean, "| target iters:", N_FISH_BOOT))
  
  env_stack_curr  <- terra::unwrap(packed_env)
  host_stack_curr <- terra::unwrap(packed_hosts_curr)
  
  out <- tryCatch({
    
    occ_out <- file.path(OUTPUT_ROOT, "occurrences", "fish", paste0(sp_clean, "_occ_used.csv"))
    if (file.exists(occ_out)) {
      sp_dat <- readr::read_csv(occ_out, show_col_types = FALSE)
    } else {
      sp_dat <- dplyr::filter(amph_occ, species == sp)
      if (nrow(sp_dat) < 5) stop("Not enough occurrences")
      if (isTRUE(USE_SPATIAL_THINNING)) sp_dat <- thin_occurrences(sp_dat, env_stack_curr)
      readr::write_csv(sp_dat, occ_out)
    }
    
    env_stack_model_a <- env_stack_curr[[ENV_PREDICTORS]]
    if ("PC5" %in% names(env_stack_model_a)) stop("BUG: PC5 leaked into Combined env predictors.")
    
    biotic_curr <- get_biotic_layer(sp, host_stack_curr, int_mat, debug_path = MASTER_LOG)
    if (is.null(biotic_curr)) {
      write_log(MASTER_LOG, paste("SKIP Fish Combined:", sp_clean, "-> no biotic layer"))
      writeLines("SKIP_NO_BIOTIC", done_file)
      return(list(species = sp_clean, stage = "FishCombined", status = "SKIP_NO_BIOTIC"))
    }
    
    full_stack_curr <- c(env_stack_model_a, biotic_curr)
    names(full_stack_curr) <- c(names(env_stack_model_a), "biotic_suitability")
    
    # Background (cached; built from EnvOnly env stack)
    if (file.exists(bg_file)) {
      bg_coords <- readRDS(bg_file)
    } else {
      if (BG_SAMPLING_METHOD %in% c("paper_exact", "nearest_neighbor")) {
        bg_coords <- get_bias_corrected_background(
          sp_dat, env_stack_model_a,
          n_bg = BG_N_BG, alpha = BG_ALPHA, method = BG_SAMPLING_METHOD,
          env_cols = BG_ENV_VARS, cand_mult = BG_CAND_MULT, geo_metric = BG_GEO_METRIC,
          seed = 42
        )
      } else {
        bg_coords <- get_random_background(sp_dat, env_stack_model_a, n_bg = BG_N_BG, buffer_m = BG_BUFFER_M)
      }
      if (is.null(bg_coords) || nrow(bg_coords) < 50) stop("Background sampling failed/too small")
      saveRDS(bg_coords, bg_file)
    }
    
    # Futures: env + future biotic from future host predictions
    future_stacks_comb <- NULL
    if (isTRUE(HAS_FUT)) {
      future_stacks_comb <- list()
      
      for (i in seq_along(FUT_FILES)) {
        scen <- scenario_names[i]
        
        f_stack_raw <- prepare_future_stack(FUT_FILES[i], env_stack_curr, STATIC_VARS)
        if (!all(names(env_stack_model_a) %in% names(f_stack_raw))) next
        f_env <- f_stack_raw[[names(env_stack_model_a)]]
        
        host_dir_fut <- file.path(OUTPUT_ROOT, "predictions/future", scen, "hosts")
        host_files_fut <- list.files(host_dir_fut, full.names = TRUE, pattern = "_mean\\.tif$")
        if (length(host_files_fut) == 0) next
        
        h_stack_f <- terra::rast(host_files_fut)
        names(h_stack_f) <- gsub("_mean$", "", tools::file_path_sans_ext(basename(host_files_fut)))
        
        biotic_fut <- get_biotic_layer(sp, h_stack_f, int_mat)
        if (is.null(biotic_fut)) next
        
        f_comb <- c(f_env, biotic_fut)
        names(f_comb) <- c(names(f_env), "biotic_suitability")
        future_stacks_comb[[scen]] <- f_comb
      }
      
      if (length(future_stacks_comb) == 0) future_stacks_comb <- NULL
    }
    
    # Tuning (cached)
    if (file.exists(tune_file)) {
      params_comb <- readr::read_csv(tune_file, show_col_types = FALSE)
      write_log(MASTER_LOG, paste("  > Combined tuning cached:", sp_clean))
    } else {
      params_comb <- get_best_params(sp_dat, full_stack_curr, bg_coords, USE_SPATIAL_TUNING, TUNE_ARGS)
      if (!is.null(params_comb)) {
        params_comb$fc <- tolower(params_comb$fc)
        readr::write_csv(as.data.frame(params_comb), tune_file)
      }
    }
    if (is.null(params_comb)) stop("Tuning failed for Combined")
    
    res_comb <- fit_bootstrap_worker(
      occ_df = sp_dat,
      current_stack = full_stack_curr,
      future_stack_list = future_stacks_comb,
      bg_coords = bg_coords,
      params = params_comb,
      n_boot = N_FISH_BOOT,
      sp_name = sp_clean,
      model_type = "Combined",
      output_dir = OUTPUT_ROOT,
      debug_log = sp_log
    )
    
    readr::write_csv(res_comb$stats, stats_file)
    terra::writeRaster(res_comb$current$mean, out_mean, overwrite = TRUE)
    terra::writeRaster(res_comb$current$sd,   out_sd,   overwrite = TRUE)
    
    if (!is.null(res_comb$future) && length(res_comb$future) > 0) {
      for (scen in names(res_comb$future)) {
        terra::writeRaster(
          res_comb$future[[scen]]$mean,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, "_mean.tif")),
          overwrite = TRUE
        )
        terra::writeRaster(
          res_comb$future[[scen]]$sd,
          file.path(OUTPUT_ROOT, "predictions/future", scen, "fish_combined", paste0(sp_clean, "_sd.tif")),
          overwrite = TRUE
        )
      }
    }
    
    if (!is.null(res_comb$completed) && res_comb$completed >= N_FISH_BOOT) {
      write_status_ok(done_file, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish Combined:", sp_clean, "| OK n=", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishCombined", status = "OK", n = N_FISH_BOOT)
    } else {
      write_status_progress(done_file, res_comb$completed, N_FISH_BOOT)
      write_log(MASTER_LOG, paste("FINISH Fish Combined:", sp_clean, "| PROGRESS n=", res_comb$completed, "/", N_FISH_BOOT))
      list(species = sp_clean, stage = "FishCombined", status = "PROGRESS", n = res_comb$completed)
    }
    
  }, error = function(e) {
    write_log(MASTER_LOG, paste("ERROR Fish Combined:", sp_clean, "->", e$message))
    list(species = sp_clean, stage = "FishCombined", status = "ERROR", msg = e$message)
  })
  
  out
}

base::try(parallel::stopCluster(cl_fish), silent = TRUE)
write_log(MASTER_LOG, paste0("=== PIPELINE FINISHED | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ==="))
cat("Pipeline finished. See log:", MASTER_LOG, "\n")

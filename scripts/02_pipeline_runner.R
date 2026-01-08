# scripts/02_pipeline_runner.R
# ------------------------------------------------------------------------------
# MASTER SDM PIPELINE RUNNER
# 1. HOSTS: Run N iterations -> Ensemble (Mean) -> Project to Future.
# 2. FISH:  Run N iterations -> Ensemble (Mean) -> Project to Future.
# ------------------------------------------------------------------------------

# --- 1. CONFIGURATION ---------------------------------------------------------
RUN_MODE <- "TEST" # "TEST" (2 iters) or "FINAL" (10 host / 50 fish)

if (RUN_MODE == "TEST") {
  N_CORES_USE   <- 30
  N_ITER_HOST   <- 2
  N_ITER_FISH   <- 2
  OUTPUT_ROOT   <- "outputs/test_run"
} else if (RUN_MODE == "FINAL") {
  N_CORES_USE   <- 30
  N_ITER_HOST   <- 10   # Stable predictor
  N_ITER_FISH   <- 50   # Statistical significance
  OUTPUT_ROOT   <- "outputs/final_run"
}

# --- 2. SETUP -----------------------------------------------------------------
if(!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, terra, dplyr, readr, stringr, ENMeval, tools)

source("scripts/01_core_functions.R")

# Create Directories
DIR_MODELS <- file.path(OUTPUT_ROOT, "models")
DIR_PREDS  <- file.path(OUTPUT_ROOT, "predictions")
DIR_LOGS   <- file.path(OUTPUT_ROOT, "logs")
DIR_DETAILS <- file.path(DIR_LOGS, "details")

# Ensure subdirectories exist
dir.create(file.path(DIR_PREDS, "current", "hosts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_PREDS, "current", "fish_combined"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DIR_PREDS, "current", "fish_env"), recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MODELS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_DETAILS, recursive = TRUE, showWarnings = FALSE)

MASTER_LOG <- file.path(DIR_LOGS, "pipeline_progress.txt")
cat(paste0("--- START: ", RUN_MODE, " MODE ---\n"), file = MASTER_LOG, append = FALSE)

# --- 3. DATA LOADING ----------------------------------------------------------
cat("Loading Data...\n")
ENV_PATH <- "data/selected_environmental_variables.csv"

# Future Env
FUTURE_DIR <- "data/env/future"
FUTURE_SCENARIO_PATHS <- list()
if (dir.exists(FUTURE_DIR)) {
  fut_files <- list.files(FUTURE_DIR, pattern = "\\.tif$", full.names = TRUE)
  for (f in fut_files) {
    scen_name <- tools::file_path_sans_ext(basename(f))
    FUTURE_SCENARIO_PATHS[[scen_name]] <- f
    dir.create(file.path(DIR_PREDS, scen_name, "hosts"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(DIR_PREDS, scen_name, "fish_combined"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(DIR_PREDS, scen_name, "fish_env"), recursive = TRUE, showWarnings = FALSE)
  }
}

amph_occ <- read_csv("data/amph_occ_env_final_dataset.csv", show_col_types = FALSE)
anem_occ <- read_csv("data/anem_occ_env_final_dataset.csv", show_col_types = FALSE)
int_mat  <- read.csv("data/interaction_matrix.csv", row.names = 1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

# --- 4. CLUSTER INIT ----------------------------------------------------------
cl <- makeCluster(N_CORES_USE)
registerDoParallel(cl)

clusterExport(cl, c("train_model_worker", "get_biotic_layer", "write_log",
                    "ENV_PATH", "FUTURE_SCENARIO_PATHS", 
                    "amph_occ", "anem_occ", "int_mat",
                    "DIR_MODELS", "DIR_PREDS", "DIR_DETAILS", "MASTER_LOG"))

clusterEvalQ(cl, {
  library(terra); library(ENMeval); library(dplyr); library(readr)
  
  get_env_raster <- function() {
    df <- read_csv(ENV_PATH, show_col_types = FALSE)
    terra::rast(df, type = "xyz", crs = "EPSG:4326")
  }
  
  save_tif <- function(r, folder, species, scenario) {
    sp_clean <- gsub(" ", "_", species)
    path <- file.path(DIR_PREDS, scenario, folder, paste0(sp_clean, ".tif"))
    terra::writeRaster(r, path, overwrite=TRUE)
    return(path)
  }
  
  run_with_log <- function(species, task_name, expr) {
    log_file <- file.path(DIR_DETAILS, paste0(gsub(" ", "_", species), "_", task_name, ".txt"))
    con <- file(log_file, open="wt")
    sink(con, type="output"); sink(con, type="message")
    cat(paste("---", task_name, ":", species, "---\n"))
    res <- tryCatch(expr, error = function(e) { cat("\nERROR:", e$message); NULL })
    sink(type="message"); sink(type="output"); close(con)
    return(res)
  }
})

# ==============================================================================
# PHASE 1: HOST SDMs (Ensemble)
# ==============================================================================
cat("\n--- PHASE 1: HOSTS ---\n")
anem_species <- unique(anem_occ$species)
clusterExport(cl, "N_ITER_HOST")

host_results <- foreach(sp = anem_species, .packages = c('terra', 'dplyr'), .errorhandling="pass") %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  write_log(MASTER_LOG, paste("START Host:", sp))
  
  run_with_log(sp, "Phase1_Host", {
    env_curr <- get_env_raster()
    stack_curr <- terra::rast()
    list_fut <- list()
    for(s in names(FUTURE_SCENARIO_PATHS)) list_fut[[s]] <- terra::rast()
    
    for(i in 1:N_ITER_HOST) {
      out <- train_model_worker(sp, anem_occ, env_curr, seed = i)
      if(!is.null(out) && is.null(out$error)) {
        pred <- enm.maxnet@predict(out$model_obj, env_curr, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
        stack_curr <- c(stack_curr, pred)
        
        for(scen in names(FUTURE_SCENARIO_PATHS)) {
          env_f <- terra::rast(FUTURE_SCENARIO_PATHS[[scen]])
          names(env_f) <- names(env_curr)
          pred_f <- enm.maxnet@predict(out$model_obj, env_f, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
          list_fut[[scen]] <- c(list_fut[[scen]], pred_f)
        }
      }
    }
    
    if(terra::nlyr(stack_curr) > 0) {
      save_tif(terra::app(stack_curr, mean, na.rm=TRUE), "hosts", sp, "current")
      for(scen in names(list_fut)) {
        if(terra::nlyr(list_fut[[scen]]) > 0) {
          save_tif(terra::app(list_fut[[scen]], mean, na.rm=TRUE), "hosts", sp, scen)
        }
      }
      return(TRUE)
    } else { stop("No models converged") }
  })
  write_log(MASTER_LOG, paste("DONE Host:", sp))
}

# --- REASSEMBLE HOSTS ---
cat("Reloading Host Ensembles...\n")
host_map_stack <- terra::rast()
for(sp in anem_species) {
  sp_clean <- gsub(" ", "_", sp)
  f <- file.path(DIR_PREDS, "current", "hosts", paste0(sp_clean, ".tif"))
  if(file.exists(f)) host_map_stack <- c(host_map_stack, terra::rast(f))
}
host_map_packed <- terra::wrap(host_map_stack)
clusterExport(cl, "host_map_packed")

# ==============================================================================
# PHASE 2: FISH SDMs (Stats & Ensembles)
# ==============================================================================
cat("\n--- PHASE 2: FISH ---\n")
fish_species <- unique(amph_occ$species)
clusterExport(cl, "N_ITER_FISH")

fish_results <- foreach(sp = fish_species, .packages = c('terra', 'dplyr'), .errorhandling="pass") %dopar% {
  sp_clean <- gsub(" ", "_", sp)
  write_log(MASTER_LOG, paste("START Fish:", sp))
  
  run_with_log(sp, "Phase2_Fish", {
    env_curr <- get_env_raster()
    host_maps <- terra::unwrap(host_map_packed)
    bio_layer <- get_biotic_layer(sp, int_mat, host_maps)
    
    # Stacks for Ensembling
    stack_cmb_curr <- terra::rast()
    stack_env_curr <- terra::rast()
    list_cmb_fut <- list(); for(s in names(FUTURE_SCENARIO_PATHS)) list_cmb_fut[[s]] <- terra::rast()
    list_env_fut <- list(); for(s in names(FUTURE_SCENARIO_PATHS)) list_env_fut[[s]] <- terra::rast()
    
    for(i in 1:N_ITER_FISH) {
      
      # 1. ENV ONLY
      out_env <- train_model_worker(sp, amph_occ, env_curr, seed = i)
      if(!is.null(out_env) && is.null(out_env$error)) {
        saveRDS(out_env$stats, file.path(DIR_MODELS, paste0(sp_clean, "_EnvOnly_iter_", i, ".rds")))
        
        # Accumulate Current
        pred <- enm.maxnet@predict(out_env$model_obj, env_curr, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
        stack_env_curr <- c(stack_env_curr, pred)
        
        # Accumulate Future
        for(scen in names(FUTURE_SCENARIO_PATHS)) {
          env_f <- terra::rast(FUTURE_SCENARIO_PATHS[[scen]])
          names(env_f) <- names(env_curr)
          pred_f <- enm.maxnet@predict(out_env$model_obj, env_f, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
          list_env_fut[[scen]] <- c(list_env_fut[[scen]], pred_f)
        }
      }
      
      # 2. COMBINED
      if(!is.null(bio_layer)) {
        combined_stack <- c(env_curr, bio_layer)
        out_cmb <- train_model_worker(sp, amph_occ, combined_stack, seed = i)
        
        if(!is.null(out_cmb) && is.null(out_cmb$error)) {
          saveRDS(out_cmb$stats, file.path(DIR_MODELS, paste0(sp_clean, "_Combined_iter_", i, ".rds")))
          
          # Accumulate Current
          pred <- enm.maxnet@predict(out_cmb$model_obj, combined_stack, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
          stack_cmb_curr <- c(stack_cmb_curr, pred)
          
          # Accumulate Future (Requires Future Hosts)
          for(scen in names(FUTURE_SCENARIO_PATHS)) {
            env_f <- terra::rast(FUTURE_SCENARIO_PATHS[[scen]])
            names(env_f) <- names(env_curr)
            
            # Re-build Biotic Layer using Future Host Ensembles
            hosts <- names(int_mat)[int_mat[sp_clean, ] > 0]
            stack_fut_hosts <- terra::rast()
            for(h in hosts) {
              f <- file.path(DIR_PREDS, scen, "hosts", paste0(h, ".tif"))
              if(file.exists(f)) stack_fut_hosts <- c(stack_fut_hosts, terra::rast(f))
            }
            
            if(terra::nlyr(stack_fut_hosts) > 0) {
              if(terra::nlyr(stack_fut_hosts) == 1) bio_f <- stack_fut_hosts else bio_f <- terra::app(stack_fut_hosts, max, na.rm=TRUE)
              names(bio_f) <- "host_suitability"
              combined_f <- c(env_f, bio_f)
              
              pred_f <- enm.maxnet@predict(out_cmb$model_obj, combined_f, other.settings = list(pred.type = "cloglog", doClamp = TRUE))
              list_cmb_fut[[scen]] <- c(list_cmb_fut[[scen]], pred_f)
            }
          }
        }
      }
    }
    
    # Save Final Ensembles (Means)
    if(terra::nlyr(stack_env_curr) > 0) save_tif(terra::app(stack_env_curr, mean, na.rm=TRUE), "fish_env", sp, "current")
    for(s in names(list_env_fut)) if(terra::nlyr(list_env_fut[[s]]) > 0) save_tif(terra::app(list_env_fut[[s]], mean, na.rm=TRUE), "fish_env", sp, s)
    
    if(terra::nlyr(stack_cmb_curr) > 0) save_tif(terra::app(stack_cmb_curr, mean, na.rm=TRUE), "fish_combined", sp, "current")
    for(s in names(list_cmb_fut)) if(terra::nlyr(list_cmb_fut[[s]]) > 0) save_tif(terra::app(list_cmb_fut[[s]], mean, na.rm=TRUE), "fish_combined", sp, s)
  })
  write_log(MASTER_LOG, paste("DONE Fish:", sp))
}

stopCluster(cl)
write_log(MASTER_LOG, "--- PIPELINE FINISHED ---")
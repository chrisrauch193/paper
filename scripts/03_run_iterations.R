# scripts/03_run_iterations.R
# Runs 50 iterations of models for statistical comparison
# Uses parallel processing with terra-safe wrapping

# --- Setup ---
if(!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, terra, dplyr, readr, stringr, ENMeval)

# Define functions locally first so they exist in global env
source("scripts/01_model_functions.R")

# --- Config ---
N_CORES <- 20 
N_ITERATIONS <- 50
OUT_DIR <- "outputs/models/iterations"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
cat("Loading Data...\n")
env_df <- read_csv("data/selected_environmental_variables.csv", show_col_types = FALSE)
env_rast <- terra::rast(env_df, type = "xyz", crs = "EPSG:4326")

# Load Host Maps & Sanitize
anem_res <- readRDS("Rdata/anemENMs.RDS")
host_maps <- anem_res$maps
if (inherits(host_maps, "PackedSpatRaster")) host_maps <- terra::unwrap(host_maps)
names(host_maps) <- gsub(" ", "_", names(host_maps))
names(host_maps) <- gsub("\\.", "_", names(host_maps))

amph_occ <- read_csv("data/amph_occ_env_final_dataset.csv", show_col_types = FALSE)
int_mat <- read.csv("data/interaction_matrix.csv", row.names = 1)
colnames(int_mat) <- gsub("\\.", "_", colnames(int_mat))
rownames(int_mat) <- gsub("\\.", "_", rownames(int_mat))

species_list <- unique(amph_occ$species)

# --- PREPARE FOR PARALLEL (CRITICAL STEP) ---
# Wrap rasters so they survive the trip to the worker nodes
env_rast_packed <- terra::wrap(env_rast)
host_maps_packed <- terra::wrap(host_maps)

# --- Parallel Cluster Setup ---
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

# Export objects explicitly to workers
# Note: functions are exported, but wrapped rasters are passed via arguments or export
clusterExport(cl, c("run_species_model", "get_biotic_layer", "species_list", 
                    "amph_occ", "int_mat", "OUT_DIR", 
                    "env_rast_packed", "host_maps_packed"))

clusterEvalQ(cl, {
  library(terra)
  library(ENMeval)
  library(dplyr)
})

cat("Starting Parallel Loop over", N_ITERATIONS, "iterations...\n")

# --- The Loop ---
results <- foreach(iter = 1:N_ITERATIONS, 
                   .packages = c('terra', 'dplyr', 'ENMeval'), 
                   .errorhandling = "pass") %dopar% { # "pass" lets us see errors without crashing everything
                     
                     # UNWRAP INSIDE THE WORKER
                     env_rast_worker <- terra::unwrap(env_rast_packed)
                     host_maps_worker <- terra::unwrap(host_maps_packed)
                     
                     # Create a log file for this worker/iteration to see progress
                     log_file <- file.path(OUT_DIR, paste0("log_iter_", iter, ".txt"))
                     cat("Starting Iteration", iter, "\n", file = log_file, append = TRUE)
                     
                     for (sp in species_list) {
                       sp_clean <- gsub(" ", "_", sp)
                       
                       # Define Output Filenames
                       f_env <- file.path(OUT_DIR, paste0(sp_clean, "_env_iter_", iter, ".rds"))
                       f_bio <- file.path(OUT_DIR, paste0(sp_clean, "_biotic_iter_", iter, ".rds"))
                       f_cmb <- file.path(OUT_DIR, paste0(sp_clean, "_combined_iter_", iter, ".rds"))
                       
                       # --- MODEL A: ENV ONLY ---
                       if (!file.exists(f_env)) {
                         cat("  Running Env:", sp, "\n", file = log_file, append = TRUE)
                         res <- run_species_model(sp, amph_occ, env_rast_worker, "env", seed = iter)
                         if(!is.null(res)) saveRDS(res$stats, f_env)
                       }
                       
                       # Prepare Biotic Layer
                       bio_layer <- get_biotic_layer(sp, int_mat, host_maps_worker)
                       
                       if (!is.null(bio_layer)) {
                         # --- MODEL B: BIOTIC ONLY ---
                         if (!file.exists(f_bio)) {
                           cat("  Running Biotic:", sp, "\n", file = log_file, append = TRUE)
                           res <- run_species_model(sp, amph_occ, env_rast_worker, "biotic_only", 
                                                    seed = iter, host_map_stack = bio_layer)
                           if(!is.null(res)) saveRDS(res$stats, f_bio)
                         }
                         
                         # --- MODEL C: COMBINED ---
                         if (!file.exists(f_cmb)) {
                           cat("  Running Combined:", sp, "\n", file = log_file, append = TRUE)
                           res <- run_species_model(sp, amph_occ, env_rast_worker, "combined", 
                                                    seed = iter, host_map_stack = bio_layer)
                           if(!is.null(res)) saveRDS(res$stats, f_cmb)
                         }
                       }
                     }
                     return(iter)
                   }

stopCluster(cl)
cat("Iterations complete. Check 'outputs/models/iterations' for RDS files and logs.\n")
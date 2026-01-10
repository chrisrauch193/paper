# scripts/00_config.R
# ------------------------------------------------------------------------------
# GLOBAL CONFIGURATION
# ------------------------------------------------------------------------------

# --- RUN SETTINGS ---
RUN_ID           <- "test_run"   # "test_run" or "final_run"
N_CORES          <- 4            # 30 for final
WIPE_PREDICTIONS <- TRUE         # Force re-run?

# --- METHODOLOGY TOGGLES ---
USE_SPATIAL_THINNING  <- TRUE    # Recommended: TRUE
USE_SPATIAL_TUNING    <- TRUE    # Recommended: TRUE (BlockCV)
BG_SAMPLING_METHOD    <- "paper_exact" # Options: "paper_exact", "nearest_neighbor", "random"

# --- MODEL PARAMETERS ---
if (RUN_ID == "final_run") {
  N_HOST_BOOT <- 10
  N_FISH_BOOT <- 40
  TUNE_ARGS   <- list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5))
} else {
  N_HOST_BOOT <- 2
  N_FISH_BOOT <- 2
  TUNE_ARGS   <- list(fc = c("L", "LQ", "H"), rm = c(1, 2, 5))
}

# --- SPECIES & VARS ---
TARGET_HOSTS <- c("Entacmaea_quadricolor", "Heteractis_magnifica")
TARGET_FISH  <- c("Amphiprion_clarkii", "Amphiprion_frenatus")
# STATIC_VARS  <- c("rugosity", "bathymetry", "slope", "aspect") 
STATIC_VARS  <- c("rugosity") 


# --- PATHS ---
OUTPUT_ROOT <- file.path("outputs", RUN_ID)
ENV_PATH    <- "data/final_env_stack.tif"
FUT_DIR     <- "data/env/future_pca"
# scripts/00_config.R
# ------------------------------------------------------------------------------
# GLOBAL CONFIGURATION
# ------------------------------------------------------------------------------

# --- RUN SETTINGS ---
RUN_ID           <- "final_run"   # "test_run" or "final_run"
N_CORES          <- 24            # 30 for final
WIPE_PREDICTIONS <- FALSE         # Force re-run?

# --- METHODOLOGY TOGGLES ---
USE_SPATIAL_THINNING  <- TRUE      # Recommended: TRUE
USE_SPATIAL_TUNING    <- TRUE      # Recommended: TRUE (ENMeval block partitions for tuning)
BG_SAMPLING_METHOD    <- "paper_exact" # Options: "paper_exact", "nearest_neighbor", "random"

# --- MODEL PARAMETERS ---
if (RUN_ID == "final_run") {
  N_HOST_BOOT <- 10
  N_FISH_BOOT <- 40
  TUNE_ARGS   <- list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5))
} else {
  N_HOST_BOOT <- 3
  N_FISH_BOOT <- 2
  TUNE_ARGS   <- list(fc = c("L", "LQ", "H"), rm = c(1, 2, 5))
}

# --- SPECIES & VARS ---
TARGET_HOSTS <- c("Entacmaea_quadricolor", "Heteractis_magnifica")
TARGET_FISH  <- c("Amphiprion_clarkii", "Amphiprion_frenatus")

# comment out to run all valid spp
TARGET_HOSTS <- NULL
TARGET_FISH  <- NULL

# Static predictors that should be included wherever appropriate
STATIC_VARS  <- c("rugosity")

# IMPORTANT: explicit PCs to use for environmental predictors
ENV_PCS <- c("PC1", "PC2", "PC3", "PC4")

# Explicit predictor sets (PC5 is intentionally excluded)
ENV_PREDICTORS <- unique(c(ENV_PCS, STATIC_VARS))  # PC1-4 + rugosity

# Background bias "env distance" space (keep small & interpretable)
BG_ENV_VARS <- c("PC1", "PC2")  # only these are used for env-distance in bias background

# --- EVALUATION / CV SETTINGS (FINAL RUN DEFENSIBILITY) ---
# These settings govern your *reported* metrics in bootstrap evaluation.
CV_METHOD      <- "block_quadrant"  # "block_quadrant" (4 spatial blocks) or "random"
CV_FOLDS       <- 4L
CV_ASSIGNMENT  <- "rotate"          # "rotate" or "random"
STRICT_BLOCKCV <- TRUE              # if TRUE and blockCV cannot be formed, error (no silent optimism)

MIN_TEST_PRES  <- 5L
MIN_TEST_BG    <- 200L
EVAL_BG_MAX    <- 5000L

# Bootstrap behavior (training resample)
BOOTSTRAP_PRES_WITH_REPLACEMENT <- TRUE
BOOTSTRAP_BG_WITH_REPLACEMENT   <- TRUE

# Extra outputs
DO_PERM_IMPORTANCE <- FALSE
BOYCE_N_BINS        <- 20L

# --- BACKGROUND SAMPLING (defensible + explicit) ---
BG_N_BG        <- 10000L
BG_CAND_MULT   <- 3L
BG_ALPHA       <- 0.5
BG_GEO_METRIC  <- "haversine"  # "haversine" meters (recommended); "degrees" (not recommended)
BG_BUFFER_M    <- 1000000      # used only for random background

# --- PATHS ---
OUTPUT_ROOT <- file.path("outputs", RUN_ID)
ENV_PATH    <- "data/final_env_stack.tif"
FUT_DIR     <- "data/env/future_pca"

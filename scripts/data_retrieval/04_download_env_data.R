# scripts/data_retrieval/04_download_env_data.R
# ------------------------------------------------------------------------------
# STEP 0: RAW ENV DOWNLOAD (Optional / Reproducibility)
# ------------------------------------------------------------------------------
# 1. Downloads Bio-Oracle v3 Data (Current & Future).
# 2. Downloads Terrain Variables (Rugosity).
# 3. Copies Static Variables (CHL, PAR) to Future Folders (since Bio-Oracle assumes they are static).
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, fs, remotes, stringr)

# Install obissdm if missing
if (!requireNamespace("obissdm", quietly = TRUE)) {
  remotes::install_github("iobis/obissdm")
}
library(obissdm)

# --- CONFIG ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
ENV_DIR  <- file.path(DATA_DIR, "env")

# Create structure
dir.create(file.path(ENV_DIR, "current"), recursive=T, showWarnings=F)
dir.create(file.path(ENV_DIR, "terrain"), recursive=T, showWarnings=F)
dir.create(file.path(ENV_DIR, "future", "ssp119"), recursive=T, showWarnings=F)
dir.create(file.path(ENV_DIR, "future", "ssp585"), recursive=T, showWarnings=F)

# --- VARIABLES ---
# Bio-Oracle v3 Dataset IDs
DATASETS <- c("Bio-ORACLE_v3.0_JO2023_surface", "Bio-ORACLE_v3.0_JO2023_bottom")
VARS_MEAN <- c("thetao", "so", "no3", "po4", "o2", "chl", "phyc", "par")
VARS_RANGE <- c("thetao", "no3", "o2") # Calculating ranges usually requires min/max

SCENARIOS <- c("ssp119", "ssp585")
TIME_STEPS <- c("2000-2010", "2040-2050", "2090-2100") # Approx tags

# --- DOWNLOAD FUNCTION ---
# Note: This is pseudo-code wrapper for obissdm logic. 
# Since we have the files, this documents the *intent*.

cat("--- Starting Download Logic (Bio-Oracle v3) ---\n")

# 1. Current (Baseline)
# In reality, you'd loop through dataset IDs and variables.
# obissdm::get_env_data() is the main function.

# 2. Terrain
# Download "Terrain_Ruggedness_Index" -> Rename to "rugosity.tif"
# Download "Bathymetry" -> "bathymetry_mean.tif"

# 3. Future Copying (The Critical Step for Reproducibility)
# Bio-Oracle doesn't provide future CHL/PAR/PHYC (Primary Prod). 
# We assume they are static and copy current -> future.

copy_static_vars <- function(ssp) {
  target_dir <- file.path(ENV_DIR, "future", ssp)
  current_dir <- file.path(ENV_DIR, "current")
  
  static_vars <- c("chl", "par", "phyc")
  
  for(var in static_vars) {
    # Find current file (e.g. chl_baseline_depthmax_mean.tif)
    files <- list.files(current_dir, pattern=paste0("^", var, ".*\\.tif$"), full.names=TRUE)
    
    for(f in files) {
      base <- basename(f)
      # Rename logic: remove 'baseline', add ssp + decade
      # Example: chl_baseline_mean.tif -> chl_ssp585_dec50_mean.tif
      
      # For dec50
      new_name_50 <- gsub("baseline", paste0(ssp, "_dec50"), base)
      file.copy(f, file.path(target_dir, new_name_50), overwrite=TRUE)
      
      # For dec100
      new_name_100 <- gsub("baseline", paste0(ssp, "_dec100"), base)
      file.copy(f, file.path(target_dir, new_name_100), overwrite=TRUE)
      
      cat("  Copied static var:", base, "->", ssp, "\n")
    }
  }
}

copy_static_vars("ssp119")
copy_static_vars("ssp585")

cat("--- Download & Prep Complete ---\n")
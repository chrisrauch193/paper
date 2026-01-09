# scripts/00_prepare_final_env_dataset.R
# ------------------------------------------------------------------------------
# MASTER ENV PREP: CLUSTER-SAFE VERSION
# ------------------------------------------------------------------------------
# 1. Manual 0-360 Rotation (Fixes 'unknown option: left' error).
# 2. Selects Paper's Tropical Provinces + Temperate Expansion Zones.
# 3. Trains PCA on Climate Only.
# 4. Exports Training CSV (Current) and TIFs (Future).
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, stringr, sf, tools)

# --- CONFIGURATION ------------------------------------------------------------
BASE_DIR    <- getwd() 
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
TERRAIN_DIR <- file.path(DATA_DIR, "env", "terrain")
FUTURE_RAW  <- file.path(DATA_DIR, "env", "future")

# Outputs
OUTPUT_CSV     <- file.path(DATA_DIR, "selected_environmental_variables.csv")
FUTURE_PCA_DIR <- file.path(DATA_DIR, "env", "future_pca")
dir.create(FUTURE_PCA_DIR, recursive = TRUE, showWarnings = FALSE)

# Input Files
MEOW_SHP      <- file.path(DATA_DIR, "marine_regions.shp") 
RUGOSITY_FILE <- "rugosity.tif"

# Variables
CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

N_PCA_AXES <- 5 
SEED       <- 42

# --- HELPER: ROBUST ROTATION --------------------------------------------------
# Function to convert -180..180 to 0..360 manually
# This avoids version conflicts with terra::rotate(left=FALSE)
rotate_to_360 <- function(r) {
  # If already 0-360, return as is
  if (xmin(r) >= 0) return(r)
  
  # Split into West (-180 to 0) and East (0 to 180)
  west <- terra::crop(r, terra::ext(-180, 0, -90, 90))
  east <- terra::crop(r, terra::ext(0, 180, -90, 90))
  
  # Shift West to the right side (180 to 360)
  # using shift() which is safer than setting extent manually
  west_shifted <- terra::shift(west, dx=360)
  
  # Merge East then West
  r_rotated <- terra::merge(east, west_shifted)
  
  # Fix extent precision issues
  terra::ext(r_rotated) <- c(0, 360, ymin(r), ymax(r))
  
  return(r_rotated)
}

# --- 1. PREPARE STUDY AREA (MEOW MASK) ----------------------------------------
cat("--- Preparing Study Area (MEOW) ---\n")

if(!file.exists(MEOW_SHP)) stop("MEOW Shapefile missing at: ", MEOW_SHP)

# Load MEOW
meow <- st_read(MEOW_SHP, quiet=TRUE)

# A. SHIFT LONGITUDE (0-360)
meow <- st_shift_longitude(meow)

# B. FILTER REGIONS
# Paper Guy's Tropical Set + Expansion Zones
target_provinces <- c(
  18:41, 58, 55, 9,           # Tropical Core
  49, 50, 51,                 # Temp Australasia
  10, 11, 12,                 # Warm Temp NW Pacific
  52, 53                      # Warm Temp Southern Africa
)

study_area <- meow %>% filter(PROV_CODE %in% target_provinces)

# --- 2. LOAD & PROCESS ENVIRONMENT --------------------------------------------
cat("--- Loading & Shifting Environmental Layers ---\n")

# A. Climate Stack
clim_list <- list()
for(v in CLIMATE_VARS) {
  f <- file.path(RAW_ENV_DIR, paste0(v, ".tif"))
  if(!file.exists(f)) stop("Missing: ", f)
  r <- terra::rast(f)
  
  # MANUAL ROTATE
  r <- rotate_to_360(r)
  
  names(r) <- v
  clim_list[[v]] <- r
}
clim_stack <- terra::rast(clim_list)

# B. Rugosity
rug_path <- file.path(TERRAIN_DIR, RUGOSITY_FILE)
if(!file.exists(rug_path)) stop("Rugosity missing")
rug <- terra::rast(rug_path)
rug <- rotate_to_360(rug) # Manual Rotate
names(rug) <- "rugosity"

# Align Rugosity to Climate Grid
if(!terra::compareGeom(rug, clim_stack, stopOnError=FALSE)) {
  cat("Resampling Rugosity...\n")
  rug <- terra::resample(rug, clim_stack, method="bilinear")
}

# --- 3. CREATE FINAL MASK -----------------------------------------------------
cat("--- Rasterizing Mask ---\n")

# A. MEOW Mask (Polygon to Raster)
meow_mask <- terra::rasterize(terra::vect(study_area), clim_stack, field=1)

# B. Bathymetry Mask (-200m)
bathy_path <- file.path(TERRAIN_DIR, "bathymetry_mean.tif")
if(!file.exists(bathy_path)) stop("Bathymetry missing")
bathy <- terra::rast(bathy_path)
bathy <- rotate_to_360(bathy) # Manual Rotate

if(!terra::compareGeom(bathy, clim_stack, stopOnError=FALSE)) {
  bathy <- terra::resample(bathy, clim_stack, method="near")
}
depth_mask <- (bathy > -200)
depth_mask[depth_mask == 0] <- NA

# Combine
final_mask <- meow_mask * depth_mask

# Apply Mask
clim_stack <- terra::mask(clim_stack, final_mask)
rug        <- terra::mask(rug, final_mask)

# Crop to data extent (removes empty ocean)
clim_stack <- terra::trim(clim_stack)
rug        <- terra::crop(rug, clim_stack)
final_mask <- terra::crop(final_mask, clim_stack) 

cat("Final Dimensions:", paste(dim(clim_stack), collapse="x"), "\n")

# --- 4. PCA & CURRENT EXPORT --------------------------------------------------
cat("--- Training PCA ---\n")
set.seed(SEED)
samp <- terra::spatSample(clim_stack, size=50000, method="random", na.rm=TRUE, xy=FALSE)
samp <- samp[complete.cases(samp), ]
pca_model <- prcomp(samp, center=TRUE, scale.=TRUE)
saveRDS(pca_model, file.path(DATA_DIR, "env", "pca_model.rds"))

cat("--- Saving Current CSV ---\n")
pca_curr <- terra::predict(clim_stack, pca_model, index=1:N_PCA_AXES)
names(pca_curr) <- paste0("PC", 1:N_PCA_AXES)

final_curr <- c(pca_curr, rug)
df_curr <- terra::as.data.frame(final_curr, xy=TRUE, na.rm=TRUE)
df_curr <- as.data.frame(lapply(df_curr, as.numeric)) # Force Numeric

write_csv(df_curr, OUTPUT_CSV)
cat("Saved:", OUTPUT_CSV, "\n")

# --- 5. FUTURE PROJECTIONS ----------------------------------------------------
cat("--- Processing Future Scenarios ---\n")
scenarios <- list.dirs(FUTURE_RAW, full.names=FALSE, recursive=FALSE)

for(scen in scenarios) {
  if(scen=="") next
  files <- list.files(file.path(FUTURE_RAW, scen), pattern="\\.tif$", full.names=TRUE)
  time_steps <- unique(str_extract(basename(files), "2050|2100|dec50|dec100")) 
  time_steps <- time_steps[!is.na(time_steps)]
  
  for(ts in time_steps) {
    fut_list <- list()
    missing <- FALSE
    for(v in CLIMATE_VARS) {
      v_simple <- str_remove(v, "_baseline")
      matches <- grep(v_simple, files, value=TRUE)
      matches <- grep(ts, matches, value=TRUE)
      if(length(matches) >= 1) {
        r <- terra::rast(matches[1])
        
        # MANUAL ROTATE FUTURE
        r <- rotate_to_360(r)
        
        names(r) <- v 
        fut_list[[v]] <- r
      } else { missing <- TRUE }
    }
    
    if(!missing) {
      fut_stack <- terra::rast(fut_list)
      
      # Apply the Master Mask
      if(!terra::compareGeom(fut_stack, final_mask, stopOnError=FALSE)) {
        fut_stack <- terra::resample(fut_stack, final_mask, method="bilinear")
      }
      fut_stack <- terra::mask(fut_stack, final_mask)
      
      # PCA Projection
      pca_fut <- terra::predict(fut_stack, pca_model, index=1:N_PCA_AXES)
      names(pca_fut) <- paste0("PC", 1:N_PCA_AXES)
      
      # Add Static Rugosity
      rug_fut <- terra::resample(rug, pca_fut, method="near")
      names(rug_fut) <- "rugosity"
      
      final_fut <- c(pca_fut, rug_fut)
      
      terra::writeRaster(final_fut, file.path(FUTURE_PCA_DIR, paste0(scen, "_", ts, ".tif")), overwrite=TRUE)
      cat("    Saved:", paste0(scen, "_", ts, ".tif"), "\n")
    }
  }
}
cat("DONE.\n")
# scripts/00_prepare_final_env_dataset.R
# ------------------------------------------------------------------------------
# MASTER ENV PREP: WCMC DYNAMIC SELECTION + EXPANSION ZONES
# ------------------------------------------------------------------------------
# 1. Loads MEOW & WCMC Shapefiles (0-360 Shift).
# 2. Dynamically selects MEOW regions containing WCMC Reefs.
# 3. Adds Warm-Temperate Expansion Zones for 2100 projections.
# 4. Masks Bio-Oracle data by Depth (-200m) and these Regions.
# 5. Trains PCA (Climate Only) & Exports Final Datasets.
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, stringr, sf, tools)

# --- CONFIGURATION ------------------------------------------------------------
BASE_DIR    <- getwd() 
DATA_DIR    <- file.path(BASE_DIR, "data")
SHP_DIR     <- file.path(DATA_DIR, "shapefiles") # New Shapefile path
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
TERRAIN_DIR <- file.path(DATA_DIR, "env", "terrain")
FUTURE_RAW  <- file.path(DATA_DIR, "env", "future")

# Outputs
OUTPUT_CSV     <- file.path(DATA_DIR, "selected_environmental_variables.csv")
FUTURE_PCA_DIR <- file.path(DATA_DIR, "env", "future_pca")
dir.create(FUTURE_PCA_DIR, recursive = TRUE, showWarnings = FALSE)

# Input Files
MEOW_FILE <- file.path(SHP_DIR, "meow_ecos.shp")
WCMC_FILE <- file.path(SHP_DIR, "WCMC008_CoralReef2018_Py_v4_1.shp")
RUGOSITY_FILE <- "rugosity.tif" # Assumed in TERRAIN_DIR

# Variables
CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

N_PCA_AXES <- 5 
SEED       <- 42

# --- HELPER: ROBUST RASTER ROTATION (0-360) -----------------------------------
rotate_to_360 <- function(r) {
  if (xmin(r) >= 0) return(r)
  west <- terra::crop(r, terra::ext(-180, 0, -90, 90))
  east <- terra::crop(r, terra::ext(0, 180, -90, 90))
  west_shifted <- terra::shift(west, dx=360)
  r_rotated <- terra::merge(east, west_shifted)
  terra::ext(r_rotated) <- c(0, 360, ymin(r), ymax(r))
  return(r_rotated)
}

# --- 1. DEFINE STUDY AREA (WCMC + MEOW + INDO-PACIFIC FILTER) -----------------
cat("--- Defining Study Area ---\n")

if(!file.exists(MEOW_FILE) || !file.exists(WCMC_FILE)) stop("Missing Shapefiles!")

# A. Load Vectors
meow <- st_read(MEOW_FILE, quiet=TRUE)
wcmc <- st_read(WCMC_FILE, quiet=TRUE)

# B. Shift Longitude to 0-360 & Fix Topology
cat("   Shifting Longitudes & Fixing Topologies...\n")
meow <- st_make_valid(st_shift_longitude(meow))
wcmc <- st_make_valid(st_shift_longitude(wcmc))

# C. DEFINE TARGET REALMS (Crucial Fix)
# We strictly select only Indo-Pacific Realms to exclude the Atlantic/Caribbean.
target_realms <- c(
  "Western Indo-Pacific", 
  "Central Indo-Pacific", 
  "Eastern Indo-Pacific",
  "Temperate Australasia",      # Expansion Zone (South)
  "Temperate Northern Pacific", # Expansion Zone (North)
  "Temperate Southern Africa"   # Expansion Zone (West)
)

# Filter MEOW to these Realms FIRST
meow_ip <- meow %>% filter(REALM %in% target_realms)

# D. Dynamic Selection: Intersect Filtered MEOW with Coral Reefs
cat("   Intersecting Indo-Pacific MEOW with Coral Reefs...\n")
# Finds regions within our target realms that contain reefs
meow_tropical <- st_filter(meow_ip, wcmc)
tropical_prov_codes <- unique(meow_tropical$PROV_CODE)

# E. Add Neighbors (Expansion Zones)
# We explicitly ensure these neighbors are included even if they don't have reefs yet
expansion_codes <- c(49, 50, 51, 10, 11, 12, 52, 53)

final_prov_codes <- unique(c(tropical_prov_codes, expansion_codes))

# Final Study Polygon
study_area <- meow %>% filter(PROV_CODE %in% final_prov_codes)

cat("   Selected", nrow(study_area), "MEOW Ecoregions in the Indo-Pacific.\n")

# --- 2. LOAD & PROCESS ENV LAYERS ---------------------------------------------
cat("--- Loading Environmental Layers ---\n")

# A. Climate Stack
clim_list <- list()
for(v in CLIMATE_VARS) {
  f <- file.path(RAW_ENV_DIR, paste0(v, ".tif"))
  if(!file.exists(f)) stop("Missing: ", f)
  r <- terra::rast(f)
  r <- rotate_to_360(r) # Align to 0-360
  names(r) <- v
  clim_list[[v]] <- r
}
clim_stack <- terra::rast(clim_list)

# B. Rugosity
rug_path <- file.path(TERRAIN_DIR, RUGOSITY_FILE)
if(!file.exists(rug_path)) stop("Rugosity missing")
rug <- terra::rast(rug_path)
rug <- rotate_to_360(rug)
names(rug) <- "rugosity"

# Align Rugosity
if(!terra::compareGeom(rug, clim_stack, stopOnError=FALSE)) {
  rug <- terra::resample(rug, clim_stack, method="bilinear")
}

# --- 3. CREATE FINAL MASK -----------------------------------------------------
cat("--- Rasterizing Study Mask ---\n")

# A. MEOW Mask
meow_mask <- terra::rasterize(terra::vect(study_area), clim_stack, field=1)

# B. Bathymetry Mask (-200m)
bathy_path <- file.path(TERRAIN_DIR, "bathymetry_mean.tif")
if(!file.exists(bathy_path)) stop("Bathymetry missing")
bathy <- terra::rast(bathy_path)
bathy <- rotate_to_360(bathy)

if(!terra::compareGeom(bathy, clim_stack, stopOnError=FALSE)) {
  bathy <- terra::resample(bathy, clim_stack, method="near")
}
depth_mask <- (bathy > -200)
depth_mask[depth_mask == 0] <- NA

# Combine & Crop
final_mask <- meow_mask * depth_mask
clim_stack <- terra::mask(clim_stack, final_mask)
rug        <- terra::mask(rug, final_mask)

# Trim empty space
clim_stack <- terra::trim(clim_stack)
rug        <- terra::crop(rug, clim_stack)
final_mask <- terra::crop(final_mask, clim_stack)

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
        r <- rotate_to_360(r) # Align Future too
        names(r) <- v 
        fut_list[[v]] <- r
      } else { missing <- TRUE }
    }
    
    if(!missing) {
      fut_stack <- terra::rast(fut_list)
      
      # Apply Master Mask
      if(!terra::compareGeom(fut_stack, final_mask, stopOnError=FALSE)) {
        fut_stack <- terra::resample(fut_stack, final_mask, method="bilinear")
      }
      fut_stack <- terra::mask(fut_stack, final_mask)
      
      pca_fut <- terra::predict(fut_stack, pca_model, index=1:N_PCA_AXES)
      names(pca_fut) <- paste0("PC", 1:N_PCA_AXES)
      
      rug_fut <- terra::resample(rug, pca_fut, method="near")
      names(rug_fut) <- "rugosity"
      
      final_fut <- c(pca_fut, rug_fut)
      
      terra::writeRaster(final_fut, file.path(FUTURE_PCA_DIR, paste0(scen, "_", ts, ".tif")), overwrite=TRUE)
      cat("    Saved:", paste0(scen, "_", ts, ".tif"), "\n")
    }
  }
}
cat("DONE.\n")
# scripts/0b_env_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 2: RAW ENV PREPARATION (Dual Mode)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf)

BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
TERRAIN_DIR <- file.path(DATA_DIR, "env", "terrain")
SHP_DIR     <- file.path(DATA_DIR, "shapefiles")
CORAL_SHP   <- file.path(SHP_DIR, "WCMC008_CoralReef2018_Py_v4_1.shp")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# 1. SETUP PATHS BASED ON MODE
if (PIPELINE_MODE == "REPLICATION") {
  REGION_SHP <- file.path(DATA_DIR, "marine_regions_strict.shp")
  OUT_RDS    <- file.path(DATA_DIR, "env_stack_intermediate_strict.rds")
} else {
  REGION_SHP <- file.path(DATA_DIR, "marine_regions.shp")
  OUT_RDS    <- file.path(DATA_DIR, "env_stack_intermediate.rds")
}

# 2. LOAD RASTERS (Standard logic)
CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

rotate_to_360 <- function(r) {
  if (xmin(r) >= 0) return(r)
  west <- crop(r, ext(-180, 0, -90, 90)); east <- crop(r, ext(0, 180, -90, 90))
  west_s <- shift(west, dx=360); m <- merge(east, west_s)
  ext(m) <- c(0, 360, ymin(r), ymax(r))
  return(m)
}

cat("Loading Rasters...\n")
clim_list <- list()
for(v in CLIMATE_VARS) clim_list[[v]] <- rotate_to_360(rast(file.path(RAW_ENV_DIR, paste0(v, ".tif"))))
clim_stack <- rast(clim_list)

rug <- rotate_to_360(rast(file.path(TERRAIN_DIR, "rugosity.tif")))
names(rug) <- "rugosity"
bathy <- rotate_to_360(rast(file.path(TERRAIN_DIR, "bathymetry_mean.tif")))

# 3. MASKING LOGIC
regions <- st_read(REGION_SHP, quiet=TRUE)
region_mask <- rasterize(vect(regions), clim_stack, field="PROV_CODE")

if(!compareGeom(rug, clim_stack, stopOnError=FALSE)) rug <- resample(rug, clim_stack)
if(!compareGeom(bathy, clim_stack, stopOnError=FALSE)) bathy <- resample(bathy, clim_stack)

if (PIPELINE_MODE == "REPLICATION") {
  cat("  REPLICATION: Hard Crop (30-240) + Union Mask\n")
  
  # A. Physio Mask
  physio_mask <- (bathy > -200) * (clim_stack[["thetao_baseline_depthmax_mean"]] > 20)
  physio_mask[physio_mask == 0] <- NA
  
  # B. Coral Mask
  coral <- st_read(CORAL_SHP, quiet=TRUE)
  coral_mask <- rasterize(vect(coral), clim_stack, field=1, background=NA, touches=TRUE) %>% rotate_to_360()
  
  # C. Union
  final_mask <- (classify(physio_mask, cbind(NA,0)) + classify(coral_mask, cbind(NA,0))) > 0
  final_mask[final_mask == 0] <- NA
  final_mask <- final_mask * region_mask
  
  # D. HARD CROP (Matches Paper Guy Exact Frame)
  crop_ext <- ext(30, 240, -40, 40)
  clim_stack <- crop(clim_stack, crop_ext)
  final_mask <- crop(final_mask, crop_ext)
  rug        <- crop(rug, crop_ext)
  region_mask <- crop(region_mask, crop_ext)
  
} else {
  cat("  EXPANSION: Hybrid Mask + Full Extent\n")
  final_mask <- region_mask * (bathy > -1000)
  
  # Trim only empty space, no hard geographic crop
  clim_stack <- mask(clim_stack, final_mask) %>% trim()
  final_mask <- crop(final_mask, clim_stack)
  rug        <- crop(rug, clim_stack) %>% mask(final_mask)
  region_mask <- crop(region_mask, clim_stack) %>% mask(final_mask)
}

final_mask[final_mask==0] <- NA
clim_stack <- mask(clim_stack, final_mask)

# 4. SAVE
saveRDS(list(clim = terra::wrap(clim_stack), rug = terra::wrap(rug), mask = terra::wrap(final_mask), region_mask = terra::wrap(region_mask)), OUT_RDS)
cat("Saved RDS to:", OUT_RDS, "\n")
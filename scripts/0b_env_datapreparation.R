# scripts/0b_env_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 2: RAW ENV PREPARATION (Exact Framing Fix)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
# "REPLICATION" = Paper Guy's Logic (Hard Crop 30-240, Union Mask)
# "EXPANSION"   = Thesis Logic (Hybrid Mask, Full Range)
PIPELINE_MODE <- "REPLICATION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf)

BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
TERRAIN_DIR <- file.path(DATA_DIR, "env", "terrain")
SHP_DIR     <- file.path(DATA_DIR, "shapefiles")
REGION_SHP  <- file.path(DATA_DIR, "marine_regions.shp")
CORAL_SHP   <- file.path(SHP_DIR, "WCMC008_CoralReef2018_Py_v4_1.shp")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

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

# 1. Load Climate
cat("Processing Climate Layers...\n")
clim_list <- list()
for(v in CLIMATE_VARS) {
  r <- rast(file.path(RAW_ENV_DIR, paste0(v, ".tif")))
  clim_list[[v]] <- rotate_to_360(r)
}
clim_stack <- rast(clim_list)

# 2. Load Rugosity & Bathy
rug <- rotate_to_360(rast(file.path(TERRAIN_DIR, "rugosity.tif")))
names(rug) <- "rugosity"
bathy <- rotate_to_360(rast(file.path(TERRAIN_DIR, "bathymetry_mean.tif")))

# 3. CREATE MASKS & CROP
regions <- st_read(REGION_SHP, quiet=TRUE)
region_mask <- rasterize(vect(regions), clim_stack, field="PROV_CODE")

if(!compareGeom(rug, clim_stack, stopOnError=FALSE)) rug <- resample(rug, clim_stack)
if(!compareGeom(bathy, clim_stack, stopOnError=FALSE)) bathy <- resample(bathy, clim_stack)

if (PIPELINE_MODE == "REPLICATION") {
  cat("  REPLICATION MODE: Applying Hard Crop (30-240) & Union Mask...\n")
  
  # A. Physiological Mask (Depth > -200 & Temp > 20)
  depth_condition <- (bathy > -200)
  temp_condition  <- (clim_stack[["thetao_baseline_depthmax_mean"]] > 20)
  physio_mask <- depth_condition * temp_condition
  physio_mask[physio_mask == 0] <- NA
  
  # B. Coral Mask
  if(file.exists(CORAL_SHP)) {
    cat("    Rasterizing Coral...\n")
    coral <- st_read(CORAL_SHP, quiet=TRUE)
    raw_template <- rast(file.path(RAW_ENV_DIR, paste0(CLIMATE_VARS[1], ".tif")))
    
    # touches=TRUE captures all pixels touching reef
    coral_mask_raw <- rasterize(vect(coral), raw_template, field=1, background=NA, touches=TRUE)
    coral_mask <- rotate_to_360(coral_mask_raw)
    if(!compareGeom(coral_mask, clim_stack, stopOnError=FALSE)) {
      coral_mask <- resample(coral_mask, clim_stack, method="near")
    }
    
    # C. Union Logic
    cat("    Combining Masks...\n")
    p_filled <- classify(physio_mask, cbind(NA, 0))
    c_filled <- classify(coral_mask, cbind(NA, 0))
    combined <- p_filled + c_filled
    final_mask <- combined > 0
    final_mask[final_mask == 0] <- NA
    
    final_mask <- final_mask * region_mask
    
  } else { stop("WCMC Coral Shapefile missing!") }
  
  # !!! THE FRAMING FIX !!! 
  # Apply Hard Crop exactly like Paper Guy (30 to 240 Longitude)
  cat("    Applying Exact Framing Crop (30, 240, -40, 40)...\n")
  crop_ext <- ext(30, 240, -40, 40)
  
  # Crop everything to this exact box
  clim_stack <- crop(clim_stack, crop_ext)
  final_mask <- crop(final_mask, crop_ext)
  rug        <- crop(rug, crop_ext)
  region_mask <- crop(region_mask, crop_ext)
  
} else {
  # --- EXPANSION MODE ---
  cat("  EXPANSION MODE: Hybrid Mask...\n")
  depth_mask <- (bathy > -1000)
  final_mask <- region_mask * depth_mask
  
  # Just trim empty space
  clim_stack <- mask(clim_stack, final_mask)
  clim_stack <- trim(clim_stack)
  
  # Align others
  final_mask <- crop(final_mask, clim_stack)
  rug        <- crop(rug, clim_stack); rug <- mask(rug, final_mask)
  region_mask <- crop(region_mask, clim_stack); region_mask <- mask(region_mask, final_mask)
}

final_mask[final_mask==0] <- NA

# 4. FINAL APPLY
cat("Applying Final Masks...\n")
clim_stack <- mask(clim_stack, final_mask)

cat("Final Dimensions:", paste(dim(clim_stack), collapse="x"), "\n")

saveRDS(list(
  clim = terra::wrap(clim_stack), 
  rug  = terra::wrap(rug), 
  mask = terra::wrap(final_mask), 
  region_mask = terra::wrap(region_mask)
), file.path(DATA_DIR, "env_stack_intermediate.rds"))

cat("Saved Intermediate Stack (Wrapped).\n")
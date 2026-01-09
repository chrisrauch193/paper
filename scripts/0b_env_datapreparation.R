# scripts/0b_env_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 2: RAW ENV PREPARATION (With Master Switch)
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

# 3. CREATE MASKS BASED ON MODE
regions <- st_read(REGION_SHP, quiet=TRUE)
region_mask <- rasterize(vect(regions), clim_stack, field="PROV_CODE")

if(!compareGeom(rug, clim_stack, stopOnError=FALSE)) rug <- resample(rug, clim_stack)
if(!compareGeom(bathy, clim_stack, stopOnError=FALSE)) bathy <- resample(bathy, clim_stack)

# --- SWITCH LOGIC FOR MASKING ---
if (PIPELINE_MODE == "REPLICATION") {
  # A. DEPTH: Strict -200m (Paper Guy Standard)
  depth_mask <- (bathy > -200)
  
  # B. BIOTIC: Strict Coral Reef Intersection
  if(file.exists(CORAL_SHP)) {
    cat("  REPLICATION MODE: Masking to WCMC Coral Reefs...\n")
    
    # 1. Load Coral (Keep in original -180/180 to avoid geometry errors)
    coral <- st_read(CORAL_SHP, quiet=TRUE)
    
    # 2. Load a temporary raw template (original -180/180 env layer)
    # We use this to rasterize the coral safely before rotating
    raw_template <- rast(file.path(RAW_ENV_DIR, paste0(CLIMATE_VARS[1], ".tif")))
    
    # 3. Rasterize Coral on -180/180 grid
    cat("    Rasterizing Coral Polygons (this may take a moment)...\n")
    coral_mask_raw <- rasterize(vect(coral), raw_template, field=1, background=NA)
    
    # 4. Rotate the Mask to 0-360 (Matches clim_stack)
    cat("    Rotating Coral Mask to 0-360...\n")
    coral_mask <- rotate_to_360(coral_mask_raw)
    
    # 5. Align with Climate Stack (Fix slight grid shifts if any)
    if(!compareGeom(coral_mask, clim_stack, stopOnError=FALSE)) {
      coral_mask <- resample(coral_mask, clim_stack, method="near")
    }
    
    # Combined Mask: Region (Tropical) + Depth + Coral
    final_mask <- region_mask * depth_mask * coral_mask
    
  } else {
    stop("WCMC Coral Shapefile missing! Cannot run Replication Mode.")
  }
  
} else {
  # --- EXPANSION MODE ---
  # A. DEPTH: Relaxed -1000m (To catch Hawaii/Marquesas islands)
  depth_mask <- (bathy > -1000)
  
  # B. BIOTIC: Hybrid Mask (Potential Habitat)
  cat("  EXPANSION MODE: Masking to Potential Habitat (MEOW + Depth)...\n")
  final_mask <- region_mask * depth_mask
}

# Set 0s to NA
final_mask[final_mask==0] <- NA

# 4. APPLY MASKS & SAVE
cat("Applying Final Masks and Trimming...\n")
clim_stack <- mask(clim_stack, final_mask)
clim_stack <- trim(clim_stack)

final_mask  <- crop(final_mask, clim_stack)
rug         <- crop(rug, clim_stack); rug <- mask(rug, final_mask)
region_mask <- crop(region_mask, clim_stack); region_mask <- mask(region_mask, final_mask)

cat("Final Dimensions:", paste(dim(clim_stack), collapse="x"), "\n")

saveRDS(list(
  clim = terra::wrap(clim_stack), 
  rug  = terra::wrap(rug), 
  mask = terra::wrap(final_mask), 
  region_mask = terra::wrap(region_mask)
), file.path(DATA_DIR, "env_stack_intermediate.rds"))

cat("Saved Intermediate Stack (Wrapped).\n")
# scripts/0d.env_var_selection.R
# ------------------------------------------------------------------------------
# STEP 3: PCA TRANSFORMATION & FUTURE PROJECTION
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf, stringr)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
FUTURE_RAW  <- file.path(DATA_DIR, "env", "future")
REGION_SHP  <- file.path(DATA_DIR, "marine_regions.shp")
N_PCA_AXES  <- 5

rotate_to_360 <- function(r) {
  if (xmin(r) >= 0) return(r)
  west <- crop(r, ext(-180, 0, -90, 90)); east <- crop(r, ext(0, 180, -90, 90))
  west_s <- shift(west, dx=360)
  m <- merge(east, west_s)
  ext(m) <- c(0, 360, ymin(r), ymax(r))
  return(m)
}

# --- 1. LOAD & UNWRAP ---
cat("Loading Intermediate Data...\n")
dat <- readRDS(file.path(DATA_DIR, "env_stack_intermediate.rds"))
clim_stack  <- terra::unwrap(dat$clim)
rug         <- terra::unwrap(dat$rug)
final_mask  <- terra::unwrap(dat$mask)
region_mask <- terra::unwrap(dat$region_mask)

# --- 2. TRAIN PCA ---
cat("Training PCA...\n")
set.seed(42)
samp <- spatSample(clim_stack, 50000, method="random", na.rm=TRUE, xy=FALSE)
samp <- samp[complete.cases(samp),]
pca_model <- prcomp(samp, center=TRUE, scale.=TRUE)
saveRDS(pca_model, file.path(DATA_DIR, "env", "pca_model.rds"))

# --- 3. EXPORT CURRENT ENV ---
cat("Predicting Current PCA...\n")
pca_map <- predict(clim_stack, pca_model, index=1:N_PCA_AXES)
names(pca_map) <- paste0("PC", 1:N_PCA_AXES)
final_stack <- c(pca_map, rug)

df_out <- as.data.frame(final_stack, xy=TRUE, na.rm=TRUE)
df_out <- as.data.frame(lapply(df_out, as.numeric))
write.csv(df_out, file.path(DATA_DIR, "selected_environmental_variables.csv"), row.names=FALSE)
cat("Saved: selected_environmental_variables.csv\n")

# --- 4. EXPORT MARINE REGIONS (EXACT NAMING FIX) ---
cat("Generating marine_regions.csv...\n")
regions_df <- as.data.frame(region_mask, xy=TRUE, na.rm=TRUE)
names(regions_df)[3] <- "PROV_CODE"

regions_sf <- st_read(REGION_SHP, quiet=TRUE)
meta <- st_drop_geometry(regions_sf) %>% dplyr::select(PROV_CODE, PROVINCE, REALM)

final_regions <- left_join(regions_df, meta, by="PROV_CODE") %>%
  dplyr::select(x, y, province=PROVINCE, realm=REALM) %>%
  # !!! CRITICAL FIX !!!
  # Use [^[:alnum:]] to replace slashes AND hyphens with underscores
  mutate(
    province = gsub("[^[:alnum:]]", "_", province), 
    realm = gsub("[^[:alnum:]]", "_", realm)
  )

write.csv(final_regions, file.path(DATA_DIR, "marine_regions.csv"), row.names=FALSE)
cat("Saved: marine_regions.csv (Naming matched to Paper Guy)\n")

# --- 5. FUTURE PROJECTIONS ---
cat("Processing Futures...\n")
FUTURE_PCA_DIR <- file.path(DATA_DIR, "env", "future_pca")
dir.create(FUTURE_PCA_DIR, recursive=T, showWarnings=F)
CLIMATE_VARS <- names(clim_stack)

scenarios <- list.dirs(FUTURE_RAW, full.names=F, recursive=F)
for(scen in scenarios) {
  if(scen=="") next
  files <- list.files(file.path(FUTURE_RAW, scen), pattern="\\.tif$", full.names=T)
  time_steps <- unique(str_extract(basename(files), "2050|2100|dec50|dec100"))
  time_steps <- time_steps[!is.na(time_steps)]
  
  for(ts in time_steps) {
    fut_list <- list(); missing <- FALSE
    for(v in CLIMATE_VARS) {
      v_simp <- str_remove(v, "_baseline")
      match <- grep(v_simp, files, value=T)
      match <- grep(ts, match, value=T)
      if(length(match)>=1) {
        r <- rotate_to_360(rast(match[1]))
        names(r) <- v; fut_list[[v]] <- r
      } else { missing <- TRUE }
    }
    if(!missing) {
      stk <- rast(fut_list)
      stk <- resample(stk, clim_stack)
      stk <- mask(stk, final_mask)
      pca_f <- predict(stk, pca_model, index=1:N_PCA_AXES)
      names(pca_f) <- paste0("PC", 1:N_PCA_AXES)
      final_f <- c(pca_f, rug)
      terra::writeRaster(final_f, file.path(FUTURE_PCA_DIR, paste0(scen, "_", ts, ".tif")), overwrite=TRUE)
      cat("  Saved:", paste0(scen, "_", ts, ".tif"), "\n")
    }
  }
}
cat("DONE.\n")
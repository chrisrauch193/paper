# scripts/0d.env_var_selection.R
# ------------------------------------------------------------------------------
# STEP 3: VARIABLE SELECTION & RASTER SAVING (Dual Output)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
# "REPLICATION" = Raw Variables (Paper Guy Style)
# "EXPANSION"   = PCA Axes (PC1-PC5)
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf, FactoMineR, factoextra)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")

# --- OUTPUT PATHS ---
OUT_ENV_CSV <- file.path(DATA_DIR, "selected_environmental_variables.csv") # For Plots
OUT_ENV_TIF <- file.path(DATA_DIR, "final_env_stack.tif")                # For SDM Pipeline
OUT_REGIONS <- file.path(DATA_DIR, "marine_regions.csv")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# 1. SETUP INPUT PATHS BASED ON MODE
if (PIPELINE_MODE == "REPLICATION") {
  IN_RDS     <- file.path(DATA_DIR, "env_stack_intermediate_strict.rds")
  REGION_SHP <- file.path(DATA_DIR, "marine_regions_strict.shp")
} else {
  IN_RDS     <- file.path(DATA_DIR, "env_stack_intermediate.rds")
  REGION_SHP <- file.path(DATA_DIR, "marine_regions.shp")
}

# 2. LOAD INTERMEDIATE DATA
dat <- readRDS(IN_RDS)
clim_stack <- terra::unwrap(dat$clim)
rug        <- terra::unwrap(dat$rug)

# 3. PROCESSING
if (PIPELINE_MODE == "REPLICATION") {
  cat("  Logic: Variable Selection (Paper Guy)\n")
  
  # A. Convert to DF
  env.vars_df <- as.data.frame(clim_stack, xy=TRUE, na.rm=TRUE)
  
  # B. PCA for Selection Only
  pca_data <- env.vars_df %>% dplyr::select(-x, -y)
  pca_res <- PCA(scale(pca_data), graph = FALSE)
  
  # C. Select Top Contributors
  eig <- as.data.frame(get_eigenvalue(pca_res))
  stick <- rev(cumsum(rev(1/(1:nrow(eig)))))
  ncp <- max(2, sum(eig$eigenvalue >= stick[1:nrow(eig)]))
  
  loadings <- abs(pca_res$var$coord[, 1:ncp])
  sel_vars <- unique(as.vector(apply(loadings, 2, function(x) names(sort(x, decreasing=T)[1:3]))))
  cat("  Selected Vars:", paste(sel_vars, collapse=", "), "\n")
  
  # D. Create Final Stack (Subset)
  # We select the specific layers from the original stack + rugosity
  sel_stack <- clim_stack[[sel_vars]]
  # Resample rugosity to match exactly if needed, though they should match from 0b
  if(!compareGeom(sel_stack, rug, stopOnError=FALSE)) rug <- resample(rug, sel_stack)
  
  final_stack <- c(sel_stack, rug)
  names(final_stack) <- c(sel_vars, "rugosity")
  
} else {
  cat("  Logic: PCA Projection (Thesis Standard)\n")
  
  # A. Train PCA on Sample
  set.seed(42)
  samp <- terra::spatSample(clim_stack, 50000, method="random", na.rm=TRUE, xy=FALSE)
  samp <- samp[complete.cases(samp),]
  pca_model <- prcomp(samp, center=TRUE, scale.=TRUE)
  
  # Save Model
  saveRDS(pca_model, file.path(DATA_DIR, "env", "pca_model.rds"))
  
  # B. Project Current Climate to PCs
  pca_map <- terra::predict(clim_stack, pca_model, index=1:5)
  names(pca_map) <- paste0("PC", 1:5)
  
  # C. Create Final Stack
  if(!compareGeom(pca_map, rug, stopOnError=FALSE)) rug <- resample(rug, pca_map)
  final_stack <- c(pca_map, rug)
  names(final_stack) <- c(paste0("PC", 1:5), "rugosity")
}

# 4. SAVE OUTPUTS (Dual Save)

# A. Save TIFF (Robust Input for SDM)
terra::writeRaster(final_stack, OUT_ENV_TIF, overwrite=TRUE)
cat("Saved Final Raster Stack to:", OUT_ENV_TIF, "\n")

# B. Save CSV (For Plots/Legacy)
df_out <- as.data.frame(final_stack, xy=TRUE, na.rm=TRUE)
write.csv(df_out, OUT_ENV_CSV, row.names=FALSE)
cat("Saved Env CSV to:         ", OUT_ENV_CSV, "\n")

# 5. GENERATE REGION MAP (RESCUE LOGIC)
cat("Mapping Points to Regions...\n")
env_sf <- st_as_sf(df_out[,c("x","y")], coords=c("x","y"), crs=4326)
regions_sf <- st_read(REGION_SHP, quiet=TRUE)
env_joined <- st_join(env_sf, regions_sf["PROVINCE"], left=TRUE)

missing <- which(is.na(env_joined$PROVINCE))
if(length(missing) > 0) {
  cat("  Rescuing", length(missing), "points...\n")
  nearest <- st_nearest_feature(env_sf[missing,], regions_sf)
  env_joined$PROVINCE[missing] <- regions_sf$PROVINCE[nearest]
}

# Fix Metadata
meta <- st_drop_geometry(regions_sf) %>% dplyr::select(PROVINCE, REALM) %>% distinct()

final_regions <- env_joined %>% 
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% st_drop_geometry() %>%
  left_join(meta, by="PROVINCE") %>%
  mutate(province = gsub("[^[:alnum:]]", "_", PROVINCE), realm = gsub("[^[:alnum:]]", "_", REALM)) %>%
  dplyr::select(x, y, province, realm)

write.csv(final_regions, OUT_REGIONS, row.names=FALSE)
cat("Saved Region Map to:", OUT_REGIONS, "\n")
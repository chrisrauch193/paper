# scripts/0d.env_var_selection.R
# ------------------------------------------------------------------------------
# STEP 3: VARIABLE SELECTION (Namespace Conflict Fixed)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf, FactoMineR, factoextra)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# 1. SETUP PATHS
if (PIPELINE_MODE == "REPLICATION") {
  IN_RDS     <- file.path(DATA_DIR, "env_stack_intermediate_strict.rds")
  OUT_ENV    <- file.path(DATA_DIR, "selected_environmental_variables_strict.csv")
  OUT_REGIONS <- file.path(DATA_DIR, "marine_regions_strict.csv")
  REGION_SHP <- file.path(DATA_DIR, "marine_regions_strict.shp")
} else {
  IN_RDS     <- file.path(DATA_DIR, "env_stack_intermediate.rds")
  OUT_ENV    <- file.path(DATA_DIR, "selected_environmental_variables.csv")
  OUT_REGIONS <- file.path(DATA_DIR, "marine_regions.csv")
  REGION_SHP <- file.path(DATA_DIR, "marine_regions.shp")
}

# 2. LOAD
dat <- readRDS(IN_RDS)
clim_stack <- terra::unwrap(dat$clim)
rug        <- terra::unwrap(dat$rug)

# 3. PROCESSING
if (PIPELINE_MODE == "REPLICATION") {
  # Paper Guy's Logic (Raw Variables)
  env.vars_df <- as.data.frame(clim_stack, xy=TRUE, na.rm=TRUE)
  
  # !!! FIX: Use dplyr::select to avoid terra conflict !!!
  pca_data <- env.vars_df %>% dplyr::select(-x, -y)
  
  pca_res <- PCA(scale(pca_data), graph = FALSE)
  
  eig <- as.data.frame(get_eigenvalue(pca_res))
  stick <- rev(cumsum(rev(1/(1:nrow(eig)))))
  ncp <- max(2, sum(eig$eigenvalue >= stick[1:nrow(eig)]))
  
  loadings <- abs(pca_res$var$coord[, 1:ncp])
  sel_vars <- unique(as.vector(apply(loadings, 2, function(x) names(sort(x, decreasing=T)[1:3]))))
  
  df_out <- cbind(env.vars_df[,c("x","y")], env.vars_df[, sel_vars])
  df_out$rugosity <- terra::extract(rug, df_out[,1:2])[,2]
  
} else {
  # Thesis Logic (PCA Scores)
  # For simplicity of pipeline structure, using raw vars or placeholder logic here
  # If you have specific PCA projection logic for expansion, ensure dplyr::select is used there too
  df_out <- as.data.frame(c(clim_stack, rug), xy=TRUE, na.rm=TRUE)
}

# 4. SAVE ENV
write.csv(df_out, OUT_ENV, row.names=FALSE)
cat("Saved Env Vars to:", OUT_ENV, "\n")

# 5. GENERATE REGION MAP (RESCUE)
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

# !!! FIX: Use dplyr::select !!!
meta <- st_drop_geometry(regions_sf) %>% dplyr::select(PROVINCE, REALM) %>% distinct()

final_regions <- env_joined %>% 
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>% st_drop_geometry() %>%
  left_join(meta, by="PROVINCE") %>%
  mutate(province = gsub("[^[:alnum:]]", "_", PROVINCE), realm = gsub("[^[:alnum:]]", "_", REALM)) %>%
  dplyr::select(x, y, province, realm)

write.csv(final_regions, OUT_REGIONS, row.names=FALSE)
cat("Saved Region Map to:", OUT_REGIONS, "\n")
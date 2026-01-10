# scripts/0d.env_var_selection.R
# ------------------------------------------------------------------------------
# STEP 3: VARIABLE SELECTION & REGION MAPPING (Crash Fixed)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
PIPELINE_MODE <- "REPLICATION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, sf, stringr, FactoMineR, factoextra, missMDA)

BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
FUTURE_RAW  <- file.path(DATA_DIR, "env", "future")
REGION_SHP  <- file.path(DATA_DIR, "marine_regions.shp")
N_PCA_AXES  <- 5

rotate_to_360 <- function(r) {
  if (xmin(r) >= 0) return(r)
  west <- crop(r, ext(-180, 0, -90, 90)); east <- crop(r, ext(0, 180, -90, 90))
  west_s <- shift(west, dx=360); m <- merge(east, west_s)
  ext(m) <- c(0, 360, ymin(r), ymax(r))
  return(m)
}

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# --- 1. LOAD DATA ---
cat("Loading Data...\n")
dat <- readRDS(file.path(DATA_DIR, "env_stack_intermediate.rds"))
clim_stack  <- terra::unwrap(dat$clim)
rug         <- terra::unwrap(dat$rug)

# --- 2. EXECUTE LOGIC BASED ON MODE ---

if (PIPELINE_MODE == "REPLICATION") {
  cat(" executing Paper Guy's Variable Selection...\n")
  
  # 1. Convert to DF (na.rm=TRUE drops land)
  env.vars_df <- as.data.frame(clim_stack, xy=TRUE, na.rm=TRUE)
  if(nrow(env.vars_df) == 0) stop("Error: Dataframe empty. Check masks.")
  
  # 2. PCA
  cat("  Running PCA...\n")
  pca_data <- env.vars_df %>% dplyr::select(-x, -y)
  pca_res <- PCA(scale(pca_data), graph = FALSE)
  
  # !!! FIX: Convert Matrix to DF !!!
  eig_values <- as.data.frame(get_eigenvalue(pca_res))
  
  stick_values <- rev(cumsum(rev(1/(1:ncol(pca_data)))))
  selected_ncp <- sum(eig_values$eigenvalue >= stick_values[1:length(eig_values$eigenvalue)])
  if(selected_ncp < 2) selected_ncp <- 2
  
  # 3. Select Variables
  loading_matrix <- abs(pca_res$var$coord[, 1:selected_ncp])
  selected_vars <- apply(loading_matrix, 2, function(x) names(sort(x, decreasing = TRUE)[1:3]))
  selected_vars <- unique(as.vector(selected_vars))
  cat("  Selected Vars:", paste(selected_vars, collapse=", "), "\n")
  
  # 4. Export Raw
  final_var.sel <- cbind(env.vars_df[,c("x","y")], env.vars_df[, selected_vars])
  rug_df <- as.data.frame(rug, xy=TRUE, na.rm=TRUE)
  df_out <- left_join(final_var.sel, rug_df %>% dplyr::select(x, y, rugosity), by=c("x","y"))
  
  write.csv(df_out, file.path(DATA_DIR, "selected_environmental_variables.csv"), row.names=FALSE)
  cat("Saved: selected_environmental_variables.csv (Raw Vars)\n")
  
} else {
  # --- EXPANSION MODE ---
  cat(" executing Thesis PCA Projection...\n")
  set.seed(42)
  samp <- spatSample(clim_stack, 50000, method="random", na.rm=TRUE, xy=FALSE)
  samp <- samp[complete.cases(samp),]
  pca_model <- prcomp(samp, center=TRUE, scale.=TRUE)
  saveRDS(pca_model, file.path(DATA_DIR, "env", "pca_model.rds"))
  
  pca_map <- predict(clim_stack, pca_model, index=1:N_PCA_AXES)
  names(pca_map) <- paste0("PC", 1:N_PCA_AXES)
  final_stack <- c(pca_map, rug)
  
  df_out <- as.data.frame(final_stack, xy=TRUE, na.rm=TRUE)
  df_out <- as.data.frame(lapply(df_out, as.numeric))
  write.csv(df_out, file.path(DATA_DIR, "selected_environmental_variables.csv"), row.names=FALSE)
  cat("Saved: selected_environmental_variables.csv (PC Scores)\n")
  
  # Futures
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
        stk <- resample(stk, clim_stack); stk <- mask(stk, terra::unwrap(dat$mask))
        pca_f <- predict(stk, pca_model, index=1:N_PCA_AXES)
        names(pca_f) <- paste0("PC", 1:N_PCA_AXES)
        final_f <- c(pca_f, rug)
        terra::writeRaster(final_f, file.path(FUTURE_PCA_DIR, paste0(scen, "_", ts, ".tif")), overwrite=TRUE)
        cat("  Saved:", paste0(scen, "_", ts, ".tif"), "\n")
      }
    }
  }
}

# --- 3. GENERATE MARINE REGIONS ---
cat("Generating marine_regions.csv with Nearest Feature Rescue...\n")

env_pts <- df_out %>% dplyr::select(x, y) 
env_sf <- st_as_sf(env_pts, coords = c("x", "y"), crs = 4326)
regions_sf <- st_read(REGION_SHP, quiet=TRUE) %>% st_make_valid()

env_joined <- st_join(env_sf, regions_sf["PROVINCE"], left = TRUE)

missing_idx <- which(is.na(env_joined$PROVINCE))
if (length(missing_idx) > 0) {
  cat("  Rescuing", length(missing_idx), "points outside polygons using Nearest Feature...\n")
  nearest_idx <- st_nearest_feature(env_joined[missing_idx, ], regions_sf)
  env_joined$PROVINCE[missing_idx] <- regions_sf$PROVINCE[nearest_idx]
}

meta <- st_drop_geometry(regions_sf) %>% dplyr::select(PROVINCE, REALM) %>% distinct()
final_regions <- env_joined %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  left_join(meta, by="PROVINCE") %>%
  dplyr::select(x, y, province=PROVINCE, realm=REALM) %>%
  mutate(
    province = gsub("[^[:alnum:]]", "_", province), 
    realm = gsub("[^[:alnum:]]", "_", realm)
  )

write.csv(final_regions, file.path(DATA_DIR, "marine_regions.csv"), row.names=FALSE)
cat("Saved: marine_regions.csv\n")
cat("DONE.\n")
# scripts/0a_occ_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 4: OCCURRENCE CLEANING & MERGING (Replicates 0a logic)
# ------------------------------------------------------------------------------
# 1. Loads Pre-downloaded Occurrences (from your retrieval script).
# 2. Filters points by the Marine Regions defined in Step 1/2.
# 3. Attaches Environmental Data (PC1-5 + Rugosity) to points.
# 4. Exports final datasets for modeling.
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, terra, sf)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")

# --- 1. LOAD ENV & REGIONS ---
cat("Loading Environmental Data & Regions...\n")

# Load Output from 0d (PCA Env Data)
env_df <- read_csv(file.path(DATA_DIR, "selected_environmental_variables.csv"), show_col_types=F)
# Convert to Raster for extraction
env_rast <- rast(env_df, type="xyz", crs="EPSG:4326")

# Load Marine Regions (Map)
mr_df <- read_csv(file.path(DATA_DIR, "marine_regions.csv"), show_col_types=F)
# Create a raster of valid regions (1=Valid, NA=Invalid)
# FIX: Explicitly use dplyr::select to avoid conflict with terra::select
mr_rast <- rast(mr_df %>% dplyr::select(x, y, province) %>% mutate(val=1), type="xyz", crs="EPSG:4326")

# --- 2. PROCESSING FUNCTION ---
process_files <- function(input_dir, output_file) {
  files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)
  if(length(files) == 0) {
    cat("  No files found in:", input_dir, "\n")
    return(NULL)
  }
  
  master_df <- data.frame()
  
  for(f in files) {
    d <- read_csv(f, show_col_types=F)
    
    # Standardize Column Names (Handle decimals vs x/y)
    if("decimalLongitude" %in% names(d)) d <- d %>% rename(x=decimalLongitude, y=decimalLatitude)
    if(!"species" %in% names(d)) d$species <- tools::file_path_sans_ext(basename(f))
    
    # Keep only relevant columns
    d <- d %>% dplyr::select(x, y, species)
    
    # 2. SHIFT LONGITUDE (0-360)
    # The map is 0-360. If points are -180 (e.g. Hawaii), shift them.
    d$x <- ifelse(d$x < 0, d$x + 360, d$x)
    
    # 3. SPATIAL FILTERING (Must fall in valid Marine Region)
    # Extract values from the Region Raster
    vals <- terra::extract(mr_rast, d[, c("x", "y")])
    d <- d[!is.na(vals$val), ] # Keep only points in valid regions
    
    # 4. ATTACH ENV DATA & THIN
    # Extract Env Data (PC1..5 + Rugosity)
    env_vals <- terra::extract(env_rast, d[, c("x", "y")])
    
    # Combine (Drop the 'ID' column from extract)
    d_full <- cbind(d, env_vals[,-1]) 
    
    # Remove points with missing Env data (e.g. on land)
    d_full <- d_full[complete.cases(d_full), ] 
    
    # GRID THINNING: Remove duplicates falling in the same grid cell
    # Since we attached env data, points in the same pixel will have identical env values.
    # We remove duplicates based on x, y, and species.
    d_full <- distinct(d_full, x, y, species, .keep_all=TRUE)
    
    master_df <- bind_rows(master_df, d_full)
  }
  
  # Format Final Output (Matches Paper Guy's Format: species, x, y, env...)
  # Actually, his input to 'ENmodel' usually just needs species, x, y. 
  # But attaching env data here is good for checking. 
  # Let's save the version with JUST x,y,species for strict replication if needed,
  # or keep env if his script expects it. 
  # Looking at his 1_reg_models.R: "amph_occ <- read.csv(...)" and "env <- read.csv(...)"
  # He loads them separately. So this file should probably just be x,y,species.
  
  final_out <- master_df %>% dplyr::select(x, y, species)
  
  write.csv(final_out, output_file, row.names=FALSE)
  cat("Saved:", output_file, "with", nrow(final_out), "records.\n")
}

# --- 3. RUN PROCESSING ---

cat("Processing Anemones...\n")
# Input: Your downloaded CSVs (e.g. 289169.csv)
process_files(file.path(DATA_DIR, "occurrence", "anemone"), 
              file.path(DATA_DIR, "anem_occ_env_final_dataset.csv"))

cat("Processing Clownfish...\n")
process_files(file.path(DATA_DIR, "occurrence", "anemonefish"), 
              file.path(DATA_DIR, "amph_occ_env_final_dataset.csv"))

cat("DONE.\n")
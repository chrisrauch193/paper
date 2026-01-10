# scripts/0a_occ_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 4: OCCURRENCE CLEANING (Fixed Namespace Conflict)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, sf, tools)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")

# --- UNIFIED OUTPUT NAMES ---
ANEM_OUT  <- file.path(DATA_DIR, "anem_occ_env_final_dataset.csv")
CLOWN_OUT <- file.path(DATA_DIR, "amph_occ_env_final_dataset.csv")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# 1. SELECT SHAPEFILE
if (PIPELINE_MODE == "REPLICATION") {
  REGION_SHP <- file.path(DATA_DIR, "marine_regions_strict.shp")
} else {
  REGION_SHP <- file.path(DATA_DIR, "marine_regions.shp")
}

# 2. LOAD REGIONS
regions_sf <- st_read(REGION_SHP, quiet=TRUE) %>% st_make_valid() %>% st_union()

# 3. PROCESSING FUNCTION
process_files <- function(input_dir, output_file) {
  files <- list.files(input_dir, pattern="\\.csv$", full.names=TRUE)
  if(length(files)==0) return(NULL)
  
  master <- data.frame()
  for(f in files) {
    d <- read_csv(f, show_col_types=F)
    
    # Normalize names
    if("decimalLongitude" %in% names(d)) d <- rename(d, x=decimalLongitude, y=decimalLatitude)
    if(!"species" %in% names(d)) d$species <- file_path_sans_ext(basename(f))
    
    # !!! FIX: Use dplyr::select to avoid terra conflict !!!
    d <- d %>% dplyr::select(x, y, species) %>% filter(complete.cases(.))
    
    # Shift 0-360
    d$x <- ifelse(d$x < 0, d$x + 360, d$x)
    
    # Spatial Filter
    pts <- st_as_sf(d, coords=c("x","y"), crs=4326)
    inside <- st_intersects(pts, regions_sf, sparse=FALSE)
    
    # Save valid points
    d_clean <- d[as.vector(inside),] %>% distinct(x, y, species, .keep_all=T)
    master <- bind_rows(master, d_clean)
  }
  write.csv(master, output_file, row.names=FALSE)
  cat("Saved:", output_file, "with", nrow(master), "records.\n")
}

cat("Processing Anemones...\n")
process_files(file.path(DATA_DIR, "occurrence", "anemone"), ANEM_OUT)

cat("Processing Clownfish...\n")
process_files(file.path(DATA_DIR, "occurrence", "anemonefish"), CLOWN_OUT)
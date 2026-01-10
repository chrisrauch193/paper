# scripts/0a_occ_datapreparation.R
# ------------------------------------------------------------------------------
# STEP 4: OCCURRENCE CLEANING (Polygons instead of Grid)
# ------------------------------------------------------------------------------
# This script matches the Paper Guy's logic: 
# It filters points that fall inside the target MEOW Provinces (Polygons).
# It does NOT filter by the strict Coral/Depth Grid yet (that happens in the Model).
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, sf, tools)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
REGION_SHP <- file.path(DATA_DIR, "marine_regions.shp")

# --- 1. LOAD REGIONS (SHAPEFILE) ---
cat("Loading Marine Regions Shapefile...\n")
if(!file.exists(REGION_SHP)) stop("Marine Regions Shapefile missing. Run 0c first.")

# Load and ensure valid geometry
regions_sf <- st_read(REGION_SHP, quiet=TRUE) %>% 
  st_make_valid() %>%
  st_union() # Union all provinces into one big "Allowed Area" polygon

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
    
    # Standardize Column Names
    if("decimalLongitude" %in% names(d)) d <- d %>% rename(x=decimalLongitude, y=decimalLatitude)
    if(!"species" %in% names(d)) d$species <- tools::file_path_sans_ext(basename(f))
    
    # Keep strictly X, Y, Species (Matches Paper Guy's Format)
    d <- d %>% dplyr::select(x, y, species) %>% filter(complete.cases(.))
    
    # 1. SHIFT LONGITUDE (0-360)
    d$x <- ifelse(d$x < 0, d$x + 360, d$x)
    
    # 2. SPATIAL FILTERING (Polygon Intersection)
    # Convert points to SF object
    pts_sf <- st_as_sf(d, coords = c("x", "y"), crs = 4326, remove = FALSE)
    
    # Check intersection with the Allowed Region Polygon
    # sparse=FALSE returns a matrix of TRUE/FALSE
    inside <- st_intersects(pts_sf, regions_sf, sparse = FALSE)
    
    # Keep only points inside
    d_clean <- d[as.vector(inside), ]
    
    # 3. THINNING (Exact Duplicates)
    d_clean <- distinct(d_clean, x, y, species, .keep_all=TRUE)
    
    master_df <- bind_rows(master_df, d_clean)
  }
  
  # Final Save (x, y, species ONLY)
  write.csv(master_df, output_file, row.names=FALSE)
  cat("Saved:", output_file, "with", nrow(master_df), "records.\n")
}

# --- 3. RUN PROCESSING ---

cat("Processing Anemones...\n")
process_files(file.path(DATA_DIR, "occurrences", "anemone"), 
              file.path(DATA_DIR, "anem_occ_env_final_dataset.csv"))

cat("Processing Clownfish...\n")
process_files(file.path(DATA_DIR, "occurrences", "anemonefish"), 
              file.path(DATA_DIR, "amph_occ_env_final_dataset.csv"))

cat("DONE.\n")
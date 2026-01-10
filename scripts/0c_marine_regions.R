# scripts/0c_marine_regions.R
# ------------------------------------------------------------------------------
# STEP 1: PREPARE MARINE REGIONS (Final Scientific Selection)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
# "REPLICATION" = Strict Match (27 Provs, No Easter Island) -> Saves *_strict files
# "EXPANSION"   = Validated Thesis Scope (Replication + SE Australia ONLY) -> Saves standard files
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(sf, dplyr, readr, terra)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
SHP_DIR  <- file.path(DATA_DIR, "shapefiles")
MEOW_SHP <- file.path(SHP_DIR, "meow_ecos.shp")
PAPER_FILE <- file.path(DATA_DIR, "paper_guy_meow_ecos_df.csv")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# 1. DEFINE FILE NAMES BASED ON MODE
if (PIPELINE_MODE == "REPLICATION") {
  OUT_SHP  <- file.path(DATA_DIR, "marine_regions_strict.shp")
  OUT_CSV  <- file.path(DATA_DIR, "meow_ecos_df_strict.csv")
} else {
  OUT_SHP  <- file.path(DATA_DIR, "marine_regions.shp")
  OUT_CSV  <- file.path(DATA_DIR, "meow_ecos_df.csv")
}

# 2. LOAD & FILTER
regions <- st_read(MEOW_SHP, quiet = TRUE)

# --- DEFINITIONS ---
# Paper Guy: 9 (Japan), 18-41 (Tropical), 55 (E. Aus), 58 (W. Aus)
# Note: 42 (Easter Island) is EXCLUDED.
strict_codes <- c(9, 18:41, 55, 58)

# --- EXPANSION LOGIC (VALIDATED BY DEEP RESEARCH) ---
# Valid Expansion:
# 54 = Southeast Australian Shelf (Realized niche extension for A. latezonatus)

# Rejected False Positives (Commented Out):
# 51 = Agulhas (Too cold, no anemones)
# 53 = Northern New Zealand (No resident populations, only vagrants)
# 57 = Southwest Australian Shelf (Too cold, kelp dominated)
# 56 = Southern New Zealand (Subantarctic, way too cold)

expansion_codes <- c(strict_codes, 54) 

# Apply Filter
if (PIPELINE_MODE == "REPLICATION") {
  target_codes <- strict_codes
} else {
  target_codes <- expansion_codes
}

regions_filtered <- regions %>% filter(PROV_CODE %in% target_codes)

cat("Selected", nrow(regions_filtered), "Provinces.\n")
if(PIPELINE_MODE == "EXPANSION") cat("Includes SE Australian Shelf (Code 54).\n")

# 3. PROCESS GEOMETRY
regions_shifted <- st_shift_longitude(regions_filtered) %>% st_make_valid()

prov_data <- regions_shifted %>%
  group_by(PROVINCE) %>%
  summarise(
    PROV_CODE = first(PROV_CODE),
    REALM = first(REALM),
    geometry = st_union(geometry)
  ) %>%
  ungroup() %>% st_make_valid()

# 4. SAVE OUTPUTS
st_write(prov_data, OUT_SHP, delete_layer = TRUE, quiet=TRUE)

# Generate CSV
provs_geom <- prov_data %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
provs_geom$id <- 1:nrow(provs_geom)
coords <- st_coordinates(provs_geom) %>% as.data.frame()
colnames(coords) <- c("long", "lat", "feature", "id")
meow_df <- left_join(coords, st_drop_geometry(provs_geom), by="id")

write.csv(meow_df, OUT_CSV, row.names=FALSE)
cat("Saved Shapefile to:", OUT_SHP, "\n")
cat("Saved CSV to:      ", OUT_CSV, "\n")

# --- 5. VALIDATION SUMMARY ---
if(file.exists(PAPER_FILE)) {
  cat("\n--- COMPARISON VS PAPER GUY ---\n")
  paper_df <- read_csv(PAPER_FILE, show_col_types=F)
  paper_provs <- unique(paper_df$PROVINCE)
  my_provs <- unique(prov_data$PROVINCE)
  
  diffs <- setdiff(my_provs, paper_provs)
  if(length(diffs)==0) cat("Result: PERFECT MATCH with Paper Guy.\n")
  else {
    cat("Result: Expansion Provinces added:\n")
    print(diffs)
  }
}
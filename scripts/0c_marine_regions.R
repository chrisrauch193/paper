# scripts/0c_marine_regions.R
# ------------------------------------------------------------------------------
# STEP 1: PREPARE MARINE REGIONS (Robust Name-Based Selection)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
# "REPLICATION" = Tropical Indo-Pacific Only (Paper Guy's Scope)
# "EXPANSION"   = Tropical + Warm Temperate Expansion Zones (Thesis Scope)
PIPELINE_MODE <- "EXPANSION" 

if(!require("pacman")) install.packages("pacman")
pacman::p_load(sf, dplyr, readr, terra)

# --- CONFIG ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
SHP_DIR  <- file.path(DATA_DIR, "shapefiles")
MEOW_SHP <- file.path(SHP_DIR, "meow_ecos.shp")

cat("--- RUNNING IN", PIPELINE_MODE, "MODE ---\n")

# --- 1. LOAD SHAPEFILE ---
cat("Loading MEOW Shapefile...\n")
regions <- st_read(MEOW_SHP, quiet = TRUE)

# --- 2. DEFINE REGIONS BASED ON MODE ---

if (PIPELINE_MODE == "REPLICATION") {
  # --- PAPER GUY'S EXACT SCOPE ---
  # We select the Tropical Realms by NAME to avoid ID mismatches.
  # This covers his "18:41, 58, 55, 9" range accurately.
  cat("Selecting Tropical Core (Paper Guy's Scope)...\n")
  
  target_realms <- c(
    "Central Indo-Pacific",
    "Eastern Indo-Pacific",
    "Western Indo-Pacific"
  )
  
  # Filter by Realm
  regions_filtered <- regions %>% filter(REALM %in% target_realms)
  
  # Note: He might exclude Hawaii (Prov 40?). 
  # If you need to strictly exclude Hawaii like some clownfish studies:
  # regions_filtered <- regions_filtered %>% filter(PROVINCE != "Hawaii")
  
} else {
  # --- THESIS EXPANSION SCOPE ---
  cat("Selecting Tropical Core + Expansion Zones...\n")
  
  tropical_realms <- c("Central Indo-Pacific", "Eastern Indo-Pacific", "Western Indo-Pacific")
  
  expansion_provinces <- c(
    "Warm Temperate Northwest Pacific", # Japan (Kuroshio)
    "Southeast Australian Shelf",       # E. Australia (EAC)
    "Southwest Australian Shelf",       # W. Australia (Leeuwin)
    "West Central Australian Shelf",    # W. Australia (Ningaloo)
    "East Central Australian Shelf",    # Transition
    "Agulhas",                          # S. Africa
    "Northern New Zealand"              # NZ
  )
  
  regions_filtered <- regions %>% 
    filter(REALM %in% tropical_realms | PROVINCE %in% expansion_provinces)
}

cat("Selected", nrow(regions_filtered), "Provinces.\n")

# --- 3. SHIFT LONGITUDE (0-360) ---
cat("Shifting Longitude to 0-360...\n")
regions_shifted <- st_shift_longitude(regions_filtered)
regions_shifted <- st_make_valid(regions_shifted)

# --- 4. DISSOLVE BY PROVINCE ---
prov_data <- regions_shifted %>%
  group_by(PROVINCE) %>%
  summarise(
    PROV_CODE = first(PROV_CODE),
    REALM = first(REALM),
    RLM_CODE = first(RLM_CODE),
    Lat_Zone = first(Lat_Zone),
    geometry = st_union(geometry)
  ) %>%
  ungroup() %>%
  st_make_valid()

# --- 5. EXPORT ---
st_write(prov_data, file.path(DATA_DIR, "marine_regions.shp"), delete_layer = TRUE, quiet=TRUE)

provs_geom <- prov_data %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
provs_geom$id <- 1:nrow(provs_geom)
coords <- st_coordinates(provs_geom) %>% as.data.frame()
colnames(coords) <- c("long", "lat", "feature", "id")
meow_df <- left_join(coords, st_drop_geometry(provs_geom), by="id")

write.csv(meow_df, file.path(DATA_DIR, "meow_ecos_df.csv"), row.names=FALSE)
cat("Saved: meow_ecos_df.csv\n")
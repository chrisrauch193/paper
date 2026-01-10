# scripts/0c_marine_regions.R
# ------------------------------------------------------------------------------
# STEP 1: PREPARE MARINE REGIONS (Hard-Coded Lists)
# ------------------------------------------------------------------------------

# !!! MASTER SWITCH !!!
# "REPLICATION" = Matches Paper Guy (Tropical + Specific Temperate)
# "EXPANSION"   = Thesis Scope (Tropical + Wider Expansion Zones)
PIPELINE_MODE <- "REPLICATION" 

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

# --- 2. DEFINE REGIONS (Hard-Coded Raw Names) ---
# Note: These names match the raw shapefile BEFORE any underscores are added.

if (PIPELINE_MODE == "REPLICATION") {
  cat("Selecting Paper Guy Scope...\n")
  
  # A. Tropical Core
  tropical_realms <- c(
    "Central Indo-Pacific",
    "Eastern Indo-Pacific",
    "Western Indo-Pacific"
  )
  
  # B. Specific "Temperate" Additions (Paper Guy includes these)
  force_include <- c(
    "West Central Australian Shelf",    # Ningaloo
    "East Central Australian Shelf",    # Southern GBR
    "Warm Temperate Northwest Pacific", # Southern Japan
    "Hawaii", 
    "Easter Island",
    "Southeast Polynesia"
  )
  
  regions_filtered <- regions %>% 
    filter(REALM %in% tropical_realms | PROVINCE %in% force_include)
  
} else {
  cat("Selecting Thesis Expansion Scope...\n")
  
  # A. Tropical Core
  tropical_realms <- c("Central Indo-Pacific", "Eastern Indo-Pacific", "Western Indo-Pacific")
  
  # B. Thesis Expansion Zones
  expansion_provinces <- c(
    "Warm Temperate Northwest Pacific", # Japan
    "Southeast Australian Shelf",       # E. Australia
    "Southwest Australian Shelf",       # W. Australia
    "West Central Australian Shelf",    # W. Australia
    "East Central Australian Shelf",    # Transition
    "Agulhas",                          # S. Africa
    "Northern New Zealand"              # NZ
  )
  
  regions_filtered <- regions %>% 
    filter(REALM %in% tropical_realms | PROVINCE %in% expansion_provinces)
}

cat("Selected", nrow(regions_filtered), "Provinces.\n")

# --- 3. CLEAN UP & EXPORT ---
# Shift to 0-360
regions_shifted <- st_shift_longitude(regions_filtered)
regions_shifted <- st_make_valid(regions_shifted)

# Dissolve by Province
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

# Save Shapefile (Raw Names)
st_write(prov_data, file.path(DATA_DIR, "marine_regions.shp"), delete_layer = TRUE, quiet=TRUE)

# Export CSV for Reference (Raw Names)
provs_geom <- prov_data %>% st_cast("MULTIPOLYGON") %>% st_cast("POLYGON")
provs_geom$id <- 1:nrow(provs_geom)
coords <- st_coordinates(provs_geom) %>% as.data.frame()
colnames(coords) <- c("long", "lat", "feature", "id")
meow_df <- left_join(coords, st_drop_geometry(provs_geom), by="id")

write.csv(meow_df, file.path(DATA_DIR, "meow_ecos_df.csv"), row.names=FALSE)
cat("Saved: meow_ecos_df.csv\n")
# scripts/data_retrieval/00_download_and_clean_occurrences.R
#-------------------------------------------------------------------------------
# MASTER OCCURRENCE SCRIPT (Platinum Version v3 - Fixed Columns)
# 1. Downloads from OBIS & GBIF using synchronized species lists.
# 2. Cleans coordinates (duplicates, zeros, depth).
# 3. FILTERS by MEOW PROVINCES (Exact replication of Jiménez et al. logic).
# 4. Saves individual files as "[AphiaID].csv" with 'decimalLongitude/Latitude'.
# 5. Aggregates all into Master CSVs with 'x/y' for paper comparison.
#-------------------------------------------------------------------------------

# --- 1. Setup ---
if (!exists("config")) source("scripts/config.R")
pacman::p_load(dplyr, readr, robis, rgbif, sf, terra, CoordinateCleaner, stringr, fs)

# Define Output Directories (Using Config Paths)
dir.create(config$anemone_occurrence_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$anemonefish_occurrence_dir, recursive = TRUE, showWarnings = FALSE)

# Load Species Lists
anemone_list <- read_csv(config$anemone_species_list_file, show_col_types = FALSE)
fish_list <- read_csv(config$anemonefish_species_list_file, show_col_types = FALSE)

# Combine for processing loop
all_species <- bind_rows(
  anemone_list %>% mutate(Group = "Anemone"),
  fish_list %>% mutate(Group = "Anemonefish")
)

# --- 2. Load MEOW Shapefile ---
meow_path <- config$meow_provinces_shapefile # Use config path if available
if(is.null(meow_path)) meow_path <- file.path(config$base_dir, "data", "shapefiles", "meow_ecos.shp")

if(file.exists(meow_path)) {
  cat("Loading MEOW shapefile for regional filtering...\n")
  meow_sf <- sf::st_read(meow_path, quiet = TRUE)
  meow_sf <- sf::st_make_valid(meow_sf)
} else {
  stop("MEOW shapefile not found at: ", meow_path)
}

# --- 3. DEFINE REGIONAL LOGIC (From Jiménez et al.) ---

# A. General Indo-Pacific Filter
valid_realms_general <- c(
  "Western Indo-Pacific", "Central Indo-Pacific", "Eastern Indo-Pacific",
  "Temperate Australasia", "Temperate Northern Pacific"
)
general_poly <- meow_sf %>% filter(REALM %in% valid_realms_general) %>% sf::st_union()

# B. Specific Clownfish Province Lists (Hardcoded from Paper)
clownfish_regions <- list(
  "Amphiprion_akallopisos" = c("Western Indian Ocean", "West and South Indian Shelf", "Central Indian Ocean Islands", "Andaman", "Java Transitional"),
  "Amphiprion_clarkii" = c("Somali/Arabian", "West and South Indian Shelf", "Central Indian Ocean Islands", "Andaman", "Java Transitional", "Sunda Shelf", "South China Sea", "West Central Australian Shelf", "Northwest Australian Shelf", "Western Coral Triangle", "Warm Temperate Northwest Pacific", "South Kuroshio", "Sahul Shelf", "Tropical Northwestern Pacific", "Northeast Australian Shelf", "Eastern Coral Triangle"),
  "Amphiprion_ephippium" = c("Andaman"),
  "Amphiprion_frenatus" = c("Andaman", "Java Transitional", "Western Coral Triangle", "Sunda Shelf", "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific", "Tropical Northwestern Pacific"),
  "Amphiprion_ocellaris" = c("Andaman", "Java Transitional", "Northwest Australian Shelf", "Western Coral Triangle", "Sunda Shelf", "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific", "Tropical Northwestern Pacific", "Sahul Shelf"),
  "Amphiprion_sebae" = c("Somali/Arabian", "West and South Indian Shelf", "South China Sea", "Andaman", "Central Indian Ocean Islands", "Sunda Shelf", "Java Transitional", "Western Coral Triangle"),
  "Amphiprion_bicinctus" = c("Somali/Arabian", "Red Sea and Gulf of Aden"),
  "Amphiprion_chagosensis" = c("Central Indian Ocean Islands"),
  "Amphiprion_nigripes" = c("Central Indian Ocean Islands", "West and South Indian Shelf"),
  "Amphiprion_chrysopterus" = c("Western Coral Triangle", "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific", "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", "Tropical Southwestern Pacific", "Central Polynesia", "Southeast Polynesia"),
  "Amphiprion_akindynos" = c("East Central Australian Shelf", "Northeast Australian Shelf", "Sahul Shelf", "Tropical Southwestern Pacific"),
  "Amphiprion_latezonatus" = c("East Central Australian Shelf", "Tropical Southwestern Pacific", "Lord Howe and Norfolk Islands"),
  "Amphiprion_melanopus" = c("Java Transitional", "Western Coral Triangle", "Sunda Shelf", "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific", "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", "East Central Australian Shelf", "Tropical Southwestern Pacific", "Central Polynesia"),
  "Amphiprion_leucokranos" = c("Eastern Coral Triangle"),
  "Amphiprion_percula" = c("Western Coral Triangle", "Marshall, Gilbert and Ellis Islands", "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", "East Central Australian Shelf", "Southeast Polynesia", "Tropical Southwestern Pacific"),
  "Amphiprion_perideraion" = c("Andaman", "Java Transitional", "West Central Australian Shelf", "Northwest Australian Shelf", "Western Coral Triangle", "Sunda Shelf", "South China Sea", "South Kuroshio", "Warm Temperate Northwest Pacific", "Marshall, Gilbert and Ellis Islands", "Tropical Northwestern Pacific", "Eastern Coral Triangle", "Sahul Shelf", "Northeast Australian Shelf", "East Central Australian Shelf", "Tropical Southwestern Pacific", "Central Polynesia"),
  "Amphiprion_polymnus" = c("Andaman", "Java Transitional", "Northwest Australian Shelf", "Northeast Australian Shelf", "Sahul Shelf", "Eastern Coral Triangle", "South China Sea", "Sunda Shelf", "Western Coral Triangle", "Tropical Northwestern Pacific", "South Kuroshio", "Warm Temperate Northwest Pacific"),
  "Amphiprion_sandaracinos" = c("Eastern Coral Triangle", "South China Sea", "Sunda Shelf", "Western Coral Triangle", "Java Transitional", "South Kuroshio", "Northwest Australian Shelf", "Sahul Shelf"),
  "Premnas_biaculeatus" = c("Eastern Coral Triangle", "Java Transitional", "Northeast Australian Shelf", "Sahul Shelf", "Sunda Shelf", "Western Coral Triangle", "Tropical Southwestern Pacific", "Tropical Northwestern Pacific", "South Kuroshio", "Andaman", "Northwest Australian Shelf"),
  "Amphiprion_mccullochi" = c("Lord Howe and Norfolk Islands"),
  "Amphiprion_rubrocinctus" = c("Northwest Australian Shelf", "Sahul Shelf", "West Central Australian Shelf"),
  "Amphiprion_omanensis" = c("Somali/Arabian"),
  "Amphiprion_barberi" = c("Tropical Southwestern Pacific", "Central Polynesia", "Southeast Polynesia"),
  "Amphiprion_tricinctus" = c("Tropical Northwestern Pacific", "Marshall, Gilbert and Ellis Islands"),
  "Amphiprion_allardi" = c("Western Indian Ocean"),
  "Amphiprion_chrysogaster" = c("Western Indian Ocean"),
  "Amphiprion_fuscocaudatus" = c("Western Indian Ocean"),
  "Amphiprion_pacificus" = c("Eastern Coral Triangle", "Tropical Southwestern Pacific"),
  "Amphiprion_latifasciatus" = c("Western Indian Ocean")
)

# --- 4. The Processing Function ---

target_cols <- c("decimalLongitude", "decimalLatitude", "eventDate", "basisOfRecord", "datasetName", "occurrenceStatus", "depth")

process_species <- function(sp_row) {
  sp_name <- sp_row$scientificName
  aphia_id <- sp_row$AphiaID
  group <- sp_row$Group
  paper_key <- gsub(" ", "_", sp_name)
  
  # Determine Output Directory based on Group
  target_dir <- if(group == "Anemone") config$anemone_occurrence_dir else config$anemonefish_occurrence_dir
  
  # Construct Filename using AphiaID
  clean_filename <- file.path(target_dir, paste0(aphia_id, ".csv"))
  
  if(file.exists(clean_filename) && !config$force_rerun$occurrence_download) {
    cat("  [SKIP] ", sp_name, " (ID:", aphia_id, ")\n"); return(NULL)
  }
  
  cat("\nProcessing:", sp_name, "(", group, ")\n")
  
  # --- A. DOWNLOAD ---
  
  # 1. OBIS
  obis_data <- tryCatch({
    d <- robis::occurrence(taxonid = aphia_id)
    if (!is.null(d) && nrow(d) > 0) {
      d %>% dplyr::select(any_of(target_cols)) %>% mutate(source = "OBIS")
    } else { NULL }
  }, error = function(e) { cat("    OBIS Error:", e$message, "\n"); NULL })
  
  # 2. GBIF
  gbif_data <- tryCatch({
    backbone <- rgbif::name_backbone(name = sp_name)
    if(!is.null(backbone$usageKey)) {
      d <- rgbif::occ_data(taxonKey = backbone$usageKey, hasCoordinate = TRUE, limit = 10000)$data
      if (!is.null(d) && nrow(d) > 0) {
        d %>% dplyr::select(any_of(target_cols)) %>% mutate(source = "GBIF")
      } else { NULL }
    } else { NULL }
  }, error = function(e) { cat("    GBIF Error:", e$message, "\n"); NULL })
  
  raw_df <- bind_rows(obis_data, gbif_data)
  
  if(nrow(raw_df) == 0) { cat("    ! No records found.\n"); return(NULL) }
  raw_df$species <- paper_key 
  
  # --- B. BASIC CLEANING ---
  clean_df <- raw_df %>%
    filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
    cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", species = "species", verbose = FALSE) %>%
    cc_zero(lon = "decimalLongitude", lat = "decimalLatitude", verbose = FALSE)
  
  # --- C. SPECIAL ANEMONE LOGIC ---
  if (sp_name %in% c("Macrodactyla doreensis", "Heteractis malu")) {
    clean_df <- clean_df %>% filter(decimalLongitude >= 100)
    cat("    Applied Special Anemone Filter (<100 Lon removed).\n")
  }
  
  # --- D. DEPTH FILTERING ---
  if("depth" %in% names(clean_df)) {
    clean_df <- clean_df %>% filter(is.na(depth) | (depth >= 0 & depth <= 200))
  }
  
  # --- E. REGIONAL FILTERING ---
  filter_poly <- NULL
  if (group == "Anemonefish" && paper_key %in% names(clownfish_regions)) {
    target_provinces <- clownfish_regions[[paper_key]]
    cat("    Applying Specific Filter:", length(target_provinces), "Provinces\n")
    filter_poly <- meow_sf %>% filter(PROVINCE %in% target_provinces) %>% sf::st_union()
  } else {
    cat("    Applying General Indo-Pacific Realm Filter\n")
    filter_poly <- general_poly
  }
  
  if(nrow(clean_df) > 0) {
    pts_sf <- sf::st_as_sf(clean_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)
    valid_indices <- sf::st_intersects(pts_sf, filter_poly, sparse = FALSE)
    final_df <- clean_df[as.vector(valid_indices), ]
    
    cat("    Final Count:", nrow(final_df), "(dropped", nrow(clean_df) - nrow(final_df), "points)\n")
    
    # --- F. SAVE (Keeping decimalLongitude/Latitude for SDM script) ---
    if(nrow(final_df) > 0) {
      write_csv(final_df, clean_filename)
    } else {
      cat("    ! All points filtered out.\n")
    }
  }
}

# Run Loop
for(i in 1:nrow(all_species)) {
  process_species(all_species[i, ])
}

# --- 5. CONSOLIDATE FINAL DATASETS (Comparison/Validation) ---
cat("\n--- Aggregating Final Datasets (Saving as x/y for paper comparison) ---\n")

# A. Consolidate Anemonefish
amph_files <- fs::dir_ls(config$anemonefish_occurrence_dir, glob = "*.csv")
if(length(amph_files) > 0) {
  # Read files and rename cols to x/y ONLY for this master file
  amph_master <- purrr::map_dfr(amph_files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    if(all(c("decimalLongitude","decimalLatitude","species") %in% names(d))) {
      d %>% 
        rename(x = decimalLongitude, y = decimalLatitude) %>%
        select(x, y, species)
    } else { return(NULL) }
  })
  write_csv(amph_master, file.path(config$base_dir, "data", "amph_occ_env_final_dataset.csv"))
  cat("Saved Master Anemonefish File: data/amph_occ_env_final_dataset.csv\n")
}

# B. Consolidate Anemones
anem_files <- fs::dir_ls(config$anemone_occurrence_dir, glob = "*.csv")
if(length(anem_files) > 0) {
  anem_master <- purrr::map_dfr(anem_files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    if(all(c("decimalLongitude","decimalLatitude","species") %in% names(d))) {
      d %>% 
        rename(x = decimalLongitude, y = decimalLatitude) %>%
        select(x, y, species)
    } else { return(NULL) }
  })
  write_csv(anem_master, file.path(config$base_dir, "data", "anem_occ_env_final_dataset.csv"))
  cat("Saved Master Anemone File: data/anem_occ_env_final_dataset.csv\n")
}

cat("\n--- PIPELINE COMPLETE ---\n")
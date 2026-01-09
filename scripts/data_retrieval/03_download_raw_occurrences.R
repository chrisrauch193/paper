# scripts/data_retrieval/03_download_raw_occurrences.R
# ------------------------------------------------------------------------------
# STEP 0: RAW OCCURRENCE DOWNLOAD (Optional / Reproducibility)
# ------------------------------------------------------------------------------
# 1. Downloads occurrence data from OBIS & GBIF.
# 2. Performs basic cleaning (coordinates, duplicates).
# 3. Saves individual CSVs to 'data/occurrences/anemone' & 'anemonefish'.
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, robis, rgbif, CoordinateCleaner, fs, stringr)

# --- CONFIG ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_ANEM <- file.path(DATA_DIR, "occurrences", "anemone")
OUT_FISH <- file.path(DATA_DIR, "occurrences", "anemonefish")

dir.create(OUT_ANEM, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FISH, recursive = TRUE, showWarnings = FALSE)

# --- SPECIES LISTS ---
# Define species here directly to be self-contained
anemone_species <- c(
  "Heteractis magnifica", "Heteractis malu", "Heteractis crispa",  
  "Cryptodendrum adhaesivum", "Macrodactyla doreensis", "Stichodactyla mertensii",
  "Stichodactyla gigantea", "Heteractis aurora", "Stichodactyla haddoni", "Entacmaea quadricolor"
)

fish_species <- c(
  "Amphiprion akallopisos", "Amphiprion clarkii", "Amphiprion ephippium", 
  "Amphiprion frenatus", "Amphiprion ocellaris", "Amphiprion sebae", 
  "Amphiprion bicinctus", "Amphiprion chagosensis", "Amphiprion nigripes", 
  "Amphiprion chrysopterus", "Amphiprion akindynos", "Amphiprion latezonatus", 
  "Amphiprion melanopus", "Amphiprion leucokranos", "Amphiprion percula", 
  "Amphiprion perideraion", "Amphiprion polymnus", "Amphiprion sandaracinos", 
  "Premnas biaculeatus", "Amphiprion mccullochi", "Amphiprion rubrocinctus", 
  "Amphiprion omanensis", "Amphiprion barberi", "Amphiprion tricinctus", 
  "Amphiprion allardi", "Amphiprion chrysogaster", "Amphiprion fuscocaudatus", 
  "Amphiprion pacificus", "Amphiprion latifasciatus"
)

target_cols <- c("decimalLongitude", "decimalLatitude", "eventDate", "basisOfRecord", "datasetName", "occurrenceStatus", "depth")

# --- PROCESSING FUNCTION ---
download_species <- function(sp_name, group) {
  
  clean_name <- gsub(" ", "_", sp_name)
  target_dir <- if(group == "Anemone") OUT_ANEM else OUT_FISH
  outfile <- file.path(target_dir, paste0(clean_name, ".csv")) # Save by Name, not AphiaID for simplicity here
  
  if(file.exists(outfile)) {
    cat("  [SKIP] ", sp_name, " (Exists)\n")
    return(NULL)
  }
  
  cat("\nProcessing:", sp_name, "...\n")
  
  # 1. OBIS
  obis_data <- tryCatch({
    robis::occurrence(scientificname = sp_name) %>% 
      dplyr::select(any_of(target_cols)) %>% 
      mutate(source = "OBIS")
  }, error = function(e) { NULL })
  
  # 2. GBIF
  gbif_data <- tryCatch({
    key <- rgbif::name_backbone(name = sp_name)$usageKey
    if(!is.null(key)) {
      rgbif::occ_data(taxonKey = key, hasCoordinate = TRUE, limit = 5000)$data %>%
        dplyr::select(any_of(target_cols)) %>% 
        mutate(source = "GBIF")
    } else { NULL }
  }, error = function(e) { NULL })
  
  # Combine
  df <- bind_rows(obis_data, gbif_data)
  
  if(nrow(df) == 0) {
    cat("  ! No records found.\n")
    return(NULL)
  }
  
  # Basic Cleaning
  df_clean <- df %>%
    filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
    cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", species = "species", verbose = FALSE) %>%
    cc_zero(lon = "decimalLongitude", lat = "decimalLatitude", verbose = FALSE)
  
  # Special Filters
  if("depth" %in% names(df_clean)) {
    df_clean <- df_clean %>% filter(is.na(depth) | (depth >= 0 & depth <= 200))
  }
  if (sp_name %in% c("Macrodactyla doreensis", "Heteractis malu")) {
    df_clean <- df_clean %>% filter(decimalLongitude >= 100)
  }
  
  # Save
  df_clean$species <- clean_name
  write_csv(df_clean, outfile)
  cat("  Saved:", nrow(df_clean), "records.\n")
}

# --- EXECUTE ---
for(sp in anemone_species) download_species(sp, "Anemone")
for(sp in fish_species) download_species(sp, "Fish")

cat("\n--- RAW DOWNLOAD COMPLETE ---\n")
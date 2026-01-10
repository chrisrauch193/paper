# scripts/99_nuclear_debug.R
# ------------------------------------------------------------------------------
# NUCLEAR DEBUGGER: Single Thread, Inline Functions, Verbose Prints
# ------------------------------------------------------------------------------

# 1. CLEAN SLATE
rm(list=ls())
gc()

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, maxnet, ENMeval, dismo)

# --- CONFIGURATION ---
DEBUG_SPECIES <- "Heteractis_magnifica"
ENV_PATH      <- "data/final_env_stack.tif" 
OCC_PATH      <- "data/anem_occ_env_final_dataset.csv"

cat("=== NUCLEAR DEBUGGER STARTED ===\n")

# --- 1. DEFINE WORKER LOCALLY (To guarantee code version) ---
run_debug_worker <- function(occ_df, r_stack) {
  
  cat("\n[WORKER] Starting...\n")
  
  # A. ANONYMIZE (The Fix)
  cat("[WORKER] Original Raster Names:", paste(names(r_stack), collapse=", "), "\n")
  safe_names <- paste0("v", sprintf("%02d", 1:nlyr(r_stack)))
  names(r_stack) <- safe_names
  cat("[WORKER] Renamed Raster To:    ", paste(names(r_stack), collapse=", "), "\n")
  
  # B. DATA PREP
  # Vectorize points
  coords <- as.matrix(occ_df[, c("x","y")])
  v_occ <- vect(coords, crs="EPSG:4326")
  
  # Buffer & Mask
  cat("[WORKER] Creating Buffer...\n")
  v_buff <- aggregate(buffer(v_occ, width=1000000))
  r_mask <- mask(crop(r_stack, v_buff), v_buff)
  
  # Sample Background
  cat("[WORKER] Sampling Background...\n")
  bg_coords <- spatSample(r_mask, size=1000, method="random", na.rm=TRUE, xy=TRUE, values=FALSE)
  
  # Extract Data
  cat("[WORKER] Extracting Data...\n")
  bg_data <- extract(r_stack, bg_coords)
  pres_data <- extract(r_stack, coords)
  
  # CLEANING
  if("ID" %in% names(bg_data)) bg_data$ID <- NULL
  if("ID" %in% names(pres_data)) pres_data$ID <- NULL
  
  # Force Numeric DataFrame
  bg_data <- as.data.frame(lapply(bg_data, as.numeric))
  pres_data <- as.data.frame(lapply(pres_data, as.numeric))
  
  # Remove NAs
  bg_data <- na.omit(bg_data)
  pres_data <- na.omit(pres_data)
  
  cat("[WORKER] Training Data Columns:", paste(names(bg_data), collapse=", "), "\n")
  cat("[WORKER] N Presences:", nrow(pres_data), "| N Background:", nrow(bg_data), "\n")
  
  # C. MODELING
  # Combine
  p_vec <- c(rep(1, nrow(pres_data)), rep(0, nrow(bg_data)))
  data_df <- rbind(pres_data, bg_data)
  
  cat("[WORKER] Fitting Maxnet...\n")
  tryCatch({
    # We explicitly print the formula maxnet will use
    print(head(data_df))
    
    mod <- maxnet::maxnet(p_vec, data_df, 
                          maxnet.formula(p_vec, data_df, classes="lq"), 
                          regmult=1)
    
    cat("[WORKER] Model Fit SUCCESS!\n")
    
    # Predict
    cat("[WORKER] Projecting...\n")
    pred_map <- predict(r_stack, mod, type="logistic", na.rm=TRUE)
    plot(pred_map, main="Debug Prediction")
    cat("[WORKER] Prediction SUCCESS!\n")
    
  }, error = function(e) {
    cat("\n[WORKER] CRASH IN MAXNET/PREDICT:\n")
    print(e)
  })
}

# --- 2. EXECUTION ---

# Load Data
cat("\n[MAIN] Loading Env...\n")
if(!file.exists(ENV_PATH)) stop("TIF missing")
r <- terra::rast(ENV_PATH)
print(r)

cat("\n[MAIN] Loading Occ...\n")
occ <- read_csv(OCC_PATH, show_col_types=F) %>% filter(species == DEBUG_SPECIES)
print(nrow(occ))

# Run
run_debug_worker(occ, r)

cat("\n=== DEBUGGER FINISHED ===\n")
# scripts/sdm_runs/06d_run_sdm_anemonefish_combined.R
#-------------------------------------------------------------------------------
# Run Combined SDMs for Anemonefish Species using Environmental Predictors (PCA)
# and a Biotic Predictor (Max Host Suitability) (SDMtune Workflow).
# Includes Parallel Processing, Logging, Progress, and Species-Specific Logs.
# Includes SAC thinning + OBIS background.
# Saves outputs to target structure defined in config.R (v7 - Tweaking Integrated).
# Assumes 06a produced host predictions (_pca suffix) for ALL scenarios.
#-------------------------------------------------------------------------------
cat("--- Running Script 06d: Run Anemonefish Combined (Env PCA + Biotic) SDMs (v7 - Tweaking Integrated) ---\n")

# --- 1. Setup: Load Config FIRST ---
# if (file.exists("scripts/config.R")) {
#   source("scripts/config.R")
#   if (!exists("config") || !is.list(config)) {
#     stop("FATAL: 'config' list object not found or invalid after sourcing scripts/config.R")
#   }
# } else {
#   stop("FATAL: Configuration file 'scripts/config.R' not found.")
# }

# --- 2. Load Required Packages ---
if (!require("pacman")) install.packages("pacman", dependencies = TRUE)
pacman::p_load(terra, sf, dplyr, readr, SDMtune, tools, stringr, log4r, future, furrr, progressr)

# --- 3. Source Helper Functions ---
source(file.path(config$helpers_dir, "logging_setup.R"))
source(file.path(config$helpers_dir, "env_processing_helpers.R"))
source(file.path(config$helpers_dir, "sdm_modeling_helpers.R")) # Needs updated version

# --- 4. Setup Logging ---
logger <- setup_logger(log_file = config$log_file_path,
                       log_level = config$log_level,
                       append = config$log_append,
                       log_to_console = config$log_to_console,
                       console_level = config$log_console_level)

log4r::info(logger, "--- Starting Script 06d: Run Anemonefish Combined (Env PCA + Biotic) SDMs (v7) ---")

# --- 5. Define Group Specifics & Predictor Type ---
group_name <- "anemonefish"
species_list_file <- config$anemonefish_species_list_file
occurrence_dir <- config$anemonefish_occurrence_dir
host_group_name <- "anemone"

if(!config$use_pca_predictors) { log4r::fatal(logger, "Script 06d requires config$use_pca_predictors=TRUE."); stop("This script requires PCA env predictors.") }
predictor_type_suffix <- "_combined_pca"
host_predictor_type_suffix <- "_pca" # Suffix of host predictions to load (from 06a)

log4r::info(logger, paste("--- Processing Group:", group_name, "---"))
log4r::info(logger, paste("--- Using Predictors: Env PCA Components + Max Host Suitability", host_predictor_type_suffix ,"---"))

# --- 6. Load Environmental Predictor Paths (PCA) ---
pca_paths_rds <- config$pca_raster_paths_rds_path
if (is.null(pca_paths_rds) || !is.character(pca_paths_rds) || nchar(pca_paths_rds) == 0 || !file.exists(pca_paths_rds)) { log4r::fatal(logger, paste("PCA paths RDS missing/invalid:", pca_paths_rds %||% "NULL")); stop("PCA paths RDS invalid.") }
env_predictor_paths_list <- tryCatch(readRDS(pca_paths_rds), error = function(e) { log4r::fatal(logger, paste("Failed load PCA paths RDS:", e$message)); NULL })
if (!is.list(env_predictor_paths_list) || length(env_predictor_paths_list) == 0) { log4r::fatal(logger, "PCA paths list empty/invalid."); stop("PCA paths list invalid.") }
scenarios_to_process <- names(env_predictor_paths_list)
log4r::info(logger, paste("Loaded Env PCA raster paths for scenarios:", paste(scenarios_to_process, collapse=", ")))

# --- 7. Load Processed Host Association Data ---
association_file_path <- config$anemone_fish_association_file
if(!file.exists(association_file_path)) { log4r::fatal(logger, "Association file missing."); stop("Association file missing.") }
associations_df <- tryCatch(readr::read_csv(association_file_path, show_col_types = FALSE), error = function(e){ log4r::fatal(logger, paste("Failed load associations:", e$message)); NULL })
if(is.null(associations_df)) stop("Failed to load associations file.")
log4r::info(logger, paste("Loaded", nrow(associations_df), "fish-host associations."))

# --- 8. Get File Paths for ALL Host Anemone Predictions ---
# (Keep existing logic - already checks target location)
log4r::info(logger, "Gathering file paths for ALL potential host anemone predictions...")
host_species_list_file <- config$anemone_species_list_file
if (!file.exists(host_species_list_file)) stop("Anemone list missing.")
host_species_df <- readr::read_csv(host_species_list_file, show_col_types = FALSE)

all_host_prediction_paths <- list()
missing_host_preds_overall <- FALSE
for (scenario in scenarios_to_process) {
  all_host_prediction_paths[[scenario]] <- list(); missing_files_scenario = 0
  log4r::debug(logger, paste("Checking host prediction paths for scenario:", scenario))
  for (i in 1:nrow(host_species_df)) {
    host_sci_name <- host_species_df$scientificName[i]
    host_sci_name_sanitized <- gsub(" ", "_", host_sci_name)
    host_pred_file_path <- construct_prediction_filename(host_sci_name_sanitized, scenario, host_predictor_type_suffix, config) # Use helper with _pca suffix
    if (!is.null(host_pred_file_path) && file.exists(host_pred_file_path)) {
      all_host_prediction_paths[[scenario]][[host_sci_name_sanitized]] <- host_pred_file_path
    } else {
      log4r::warn(logger, paste("Host prediction file missing in target location:", basename(host_pred_file_path %||% paste(host_sci_name_sanitized, scenario, host_predictor_type_suffix)), "(Scenario:", scenario, ")"))
      missing_files_scenario = missing_files_scenario + 1
    }
  }
  if(length(all_host_prediction_paths[[scenario]]) == 0 && missing_files_scenario > 0){ log4r::error(logger, paste("FATAL: No host predictions found for ANY species in target location for scenario:", scenario, ". Ensure 06a ran successfully.")); missing_host_preds_overall <- TRUE }
  else if (missing_files_scenario > 0) { log4r::warn(logger, paste("Found", length(all_host_prediction_paths[[scenario]]), "host files, but", missing_files_scenario, "missing for scenario:", scenario)) }
  else { log4r::info(logger, paste("Found all", length(all_host_prediction_paths[[scenario]]), "host files for scenario:", scenario)) }
}
if(missing_host_preds_overall) stop("Halting script because required host prediction files are missing for one or more scenarios. Check logs and ensure 06a completed successfully.")


# --- 9. Create Intermediate Output Dirs ---
# (Keep existing logic)
base_intermediate_model_path <- config$models_dir_intermediate
base_intermediate_results_path <- config$results_dir_intermediate
base_species_log_path <- config$species_log_dir
if (any(sapply(list(base_intermediate_model_path, base_intermediate_results_path, base_species_log_path, group_name, predictor_type_suffix), function(x) is.null(x) || !is.character(x) || nchar(x) == 0))) { stop("FATAL: Base path/group/suffix invalid.") }
intermediate_models_dir <- file.path(base_intermediate_model_path, paste0(group_name, predictor_type_suffix))
intermediate_results_dir <- file.path(base_intermediate_results_path, paste0(group_name, predictor_type_suffix))
tryCatch({dir.create(intermediate_models_dir, recursive=T, showW=F);dir.create(intermediate_results_dir, recursive=T, showW=F);dir.create(base_species_log_path, recursive=T, showW=F);log4r::debug(logger, paste("Intermediate dirs checked/created:", group_name, predictor_type_suffix))}, error=function(e){log4r::fatal(logger,paste("Failed create intermediate dirs:",e$message));stop("Dir creation failed.")})

# --- 10. Load Species List ---
tryCatch({ species_df <- readr::read_csv(species_list_file, show_col_types = FALSE) }, error = function(e) { log4r::fatal(logger, paste("Failed load species list:", e$message)); stop("Species list failed.") })
log4r::info(logger, paste("Loaded", nrow(species_df), "species from", basename(species_list_file)))

# --- 11. Define Function to Process Single Species (Combined Model v7 - Tweaking Integrated) ---
process_species_sdm_combined <- function(species_row, config, env_predictor_paths_list, all_host_prediction_paths, associations_df,
                                         group_name, predictor_type_suffix, occurrence_dir,
                                         tuning_scenario = "current") {
  
  species_name <- species_row$scientificName
  species_name_sanitized <- gsub(" ", "_", species_name)
  species_aphia_id <- species_row$AphiaID
  
  species_log_file <- file.path(config$species_log_dir, paste0(species_name_sanitized, predictor_type_suffix, "_detail.log"))
  slog <- function(level, ...) { msg <- paste(Sys.time(), paste0("[",level,"]"), paste0("[", species_name, "]"), paste0(..., collapse = " ")); cat(msg, "\n", file = species_log_file, append = TRUE) }
  slog("INFO", "--- Starting Combined processing (Script 06d v7 - Tweaking Integrated) ---")
  
  # --- Define Intermediate Paths ---
  intermediate_models_dir <- file.path(config$models_dir_intermediate, paste0(group_name, predictor_type_suffix))
  intermediate_results_dir <- file.path(config$results_dir_intermediate, paste0(group_name, predictor_type_suffix))
  final_model_file <- file.path(intermediate_models_dir, paste0("sdm_model_", species_name_sanitized, predictor_type_suffix, ".rds"))
  tuning_rds_file <- file.path(intermediate_results_dir, paste0("sdm_tuning_", species_name_sanitized, predictor_type_suffix, "_object.rds"))
  
  # --- Get Associated Hosts ---
  associated_hosts <- associations_df %>% filter(AnemonefishScientificName == species_name) %>% pull(AssociatedAnemoneScientificName)
  associated_hosts_sanitized <- gsub(" ", "_", associated_hosts)
  if (length(associated_hosts_sanitized) == 0) { msg <- paste0("Skipping: No associated hosts found."); slog("WARN", msg); return(list(status = "skipped_no_hosts", species = species_name, occurrence_count = NA, message = msg)) }
  slog("DEBUG", paste("Associated hosts:", paste(associated_hosts_sanitized, collapse=", ")))
  
  # --- Prepare GLOBAL Combined Predictor Stack for Tuning Scenario ---
  slog("DEBUG", "Preparing GLOBAL COMBINED predictor stack for tuning/base scenario:", tuning_scenario)
  tuning_predictor_stack_global <- NULL
  tryCatch({
    tuning_env_path <- env_predictor_paths_list[[tuning_scenario]]
    if (is.null(tuning_env_path) || !file.exists(tuning_env_path)) stop(paste("Env (PCA) stack path missing for tuning:", tuning_env_path %||% "NULL"))
    tuning_env_stack <- terra::rast(tuning_env_path)
    reference_geom_tuning <- tuning_env_stack[[1]] # Use first PCA layer as reference
    
    host_paths_for_tuning <- all_host_prediction_paths[[tuning_scenario]][names(all_host_prediction_paths[[tuning_scenario]]) %in% associated_hosts_sanitized]
    if (length(host_paths_for_tuning) == 0) stop("No host prediction paths found for associated hosts in tuning scenario.")
    if (length(host_paths_for_tuning) < length(associated_hosts_sanitized)) { slog("WARN", "Missing paths for some associated hosts in tuning scenario.") }
    
    host_rasters_tuning_list <- list()
    for(host_name in names(host_paths_for_tuning)){
      host_path <- host_paths_for_tuning[[host_name]]
      host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN", paste("Failed load host", basename(host_path),":", e$message)); NULL})
      if(!is.null(host_rast)){
        if (!terra::compareGeom(reference_geom_tuning, host_rast, stopOnError=FALSE, res=TRUE)) { slog("DEBUG", paste("Resampling host", basename(host_path), "for tuning.")); host_rast <- tryCatch(terra::resample(host_rast, reference_geom_tuning, method="bilinear"), error = function(e) { slog("ERROR", paste("Failed resample host", basename(host_path),":", e$message)); NULL }) }
        if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_tuning_list[[host_name]] <- host_rast }
      }
    }
    if (length(host_rasters_tuning_list) == 0) stop("Failed load/resample ANY required host rasters for tuning.")
    
    host_stack_tuning <- terra::rast(host_rasters_tuning_list)
    max_host_suitability_tuning <- terra::app(host_stack_tuning, fun = "max", na.rm = TRUE); names(max_host_suitability_tuning) <- "host_suitability_max"
    slog("DEBUG", "Max host suitability layer created for tuning scenario.")
    
    tuning_predictor_stack_global <- c(tuning_env_stack, max_host_suitability_tuning) # Combine Env PCA + Max Host
    slog("INFO", paste("GLOBAL Tuning predictor stack created:", paste(names(tuning_predictor_stack_global), collapse=", ")))
    rm(tuning_env_stack, host_stack_tuning, max_host_suitability_tuning, host_rasters_tuning_list); gc()
    
  }, error = function(e) {
    msg <- paste0("Skipping: Failed to prepare GLOBAL combined tuning stack: ", e$message); slog("ERROR", msg)
    if (exists("tuning_env_stack")) rm(tuning_env_stack); if (exists("host_stack_tuning")) rm(host_stack_tuning); if (exists("max_host_suitability_tuning")) rm(max_host_suitability_tuning); if (exists("host_rasters_tuning_list")) rm(host_rasters_tuning_list); gc()
    return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = msg))
  })
  if (is.null(tuning_predictor_stack_global)) { return(list(status = "error_tuning_predictors", species = species_name, occurrence_count = NA, message = paste0("Failed to create GLOBAL tuning_predictor_stack."))) }
  raster_crs_terra <- terra::crs(tuning_predictor_stack_global)
  
  # --- Load and Clean Initial Occurrences (NO CELL THINNING HERE) ---
  config_for_occ_load <- config
  config_for_occ_load$thinning_method <- NULL
  slog("DEBUG", "Loading/cleaning initial occurrences (before SAC thinning).")
  occ_data_result_raw <- load_clean_individual_occ_coords(species_aphia_id, occurrence_dir, config_for_occ_load, logger=NULL, species_log_file=species_log_file)
  if (is.null(occ_data_result_raw) || is.null(occ_data_result_raw$coords) || occ_data_result_raw$count < config$min_occurrences_sdm) {
    msg <- paste0("Skipping: Insufficient occurrences before SAC thinning (", occ_data_result_raw$count %||% 0, ")."); slog("WARN", msg); return(list(status = "skipped_occurrences_raw", species = species_name, occurrence_count = occ_data_result_raw$count %||% 0, message = msg))
  }
  occs_coords_raw <- occ_data_result_raw$coords
  occurrence_count_raw <- occ_data_result_raw$count
  slog("INFO", paste("Occurrence count after basic cleaning:", occurrence_count_raw))
  
  # --- <<< Spatial Autocorrelation (SAC) Thinning >>> ---
  occs_coords_thinned <- occs_coords_raw
  occurrence_count_after_thinning <- occurrence_count_raw
  thinning_dist_applied <- NA
  
  sac_thin_result <- thin_occurrences_by_sac(
    occs_coords = occs_coords_raw,
    predictor_stack = tuning_predictor_stack_global, # Use GLOBAL combined stack
    config = config,
    logger = NULL,
    species_log_file = species_log_file
  )
  
  if (!is.null(sac_thin_result)) {
    occs_coords_thinned <- sac_thin_result$coords_thinned
    occurrence_count_after_thinning <- sac_thin_result$n_thinned
    thinning_dist_applied <- sac_thin_result$thinning_distance_km
    slog("INFO", paste("Occurrence count after SAC thinning:", occurrence_count_after_thinning, "(Distance:", thinning_dist_applied %||% "NA", "km)"))
    if(occurrence_count_after_thinning < config$min_occurrences_sdm) {
      msg <- paste0("Skipping: Insufficient occurrences after SAC thinning (", occurrence_count_after_thinning, ")."); slog("WARN", msg); return(list(status = "skipped_occurrences_sac", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
  } else { slog("WARN", "SAC thinning function returned NULL or failed. Proceeding with unthinned data.") }
  
  # --- Use thinned coordinates ---
  occs_coords <- occs_coords_thinned
  occs_sf_clean <- sf::st_as_sf(as.data.frame(occs_coords), coords=c("longitude","latitude"), crs=config$occurrence_crs)
  occs_sf_clean$AphiaID <- species_aphia_id
  if(sf::st_crs(occs_sf_clean) != sf::st_crs(raster_crs_terra)){
    occs_sf_clean <- sf::st_transform(occs_sf_clean, crs=sf::st_crs(raster_crs_terra))
  }
  
  # --- Check Existing Model ---
  final_model <- NULL; tuning_output <- NULL; load_existing_model <- FALSE
  train_swd <- NULL; species_specific_stack <- NULL
  
  if (!config$force_rerun$run_biotic_sdms && file.exists(final_model_file)) { # Check biotic flag
    slog("INFO", "Final model file exists. Attempting to load...")
    final_model <- tryCatch(readRDS(final_model_file), error=function(e) { slog("ERROR", "Failed to load existing model:", e$message); NULL})
    if(!is.null(final_model) && inherits(final_model, "SDMmodel")){
      slog("INFO", "Successfully loaded existing final model. Skipping tuning/training.")
      load_existing_model <- TRUE
      slog("DEBUG", "Regenerating background points and species stack for loaded model VI/Eval...")
      background_return_eval <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file, seed_offset = config$global_seed)
      if(!is.null(background_return_eval) && !is.null(background_return_eval$background_points) && !is.null(background_return_eval$species_specific_stack)) {
        background_points_eval <- background_return_eval$background_points
        species_specific_stack <- background_return_eval$species_specific_stack
        train_swd <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points_eval, env = species_specific_stack, verbose = FALSE) }, error = function(e) { slog("WARN", "Failed prepareSWD for VI/Eval on loaded model:", e$message); NULL })
        rm(background_points_eval); gc()
      } else { slog("WARN", "Failed background/stack generation for VI/Eval on loaded model.") }
      if(file.exists(tuning_rds_file)) {
        tuning_output <- tryCatch(readRDS(tuning_rds_file), error=function(e){slog("WARN","Could not load tuning RDS."); NULL})
      } else { slog("WARN", "Tuning RDS file not found, cannot log test metrics for loaded model.") }
      if(is.null(tuning_output)) { tuning_output <- NULL }
    } else {
      msg <- paste0("Existing final model file invalid/failed load. Will re-run."); slog("WARN", msg);
      final_model <- NULL; load_existing_model <- FALSE
    }
  }
  
  # --- Tuning & Training (if not loading existing) ---
  if (!load_existing_model) {
    slog("INFO", "Final model not found/invalid or rerun forced. Proceeding with tuning/training.")
    slog("DEBUG", "Generating background points and species stack for training...")
    background_return <- generate_sdm_background_obis(occs_sf_clean, tuning_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file, seed_offset = config$global_seed)
    if (is.null(background_return) || is.null(background_return$background_points) || is.null(background_return$species_specific_stack)) {
      msg <- paste0("Skipping: Failed background point/stack generation for training."); slog("ERROR", msg); return(list(status = "error_background", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    background_points <- background_return$background_points
    species_specific_stack <- background_return$species_specific_stack # Use this masked Combined stack
    slog("DEBUG", "Background points generated and stack masked for training.")
    
    train_swd <- tryCatch({ SDMtune::prepareSWD(species = species_name, p = occs_coords, a = background_points, env = species_specific_stack, verbose = FALSE) }, error = function(e) { slog("ERROR", paste("Failed prepareSWD for tuning/training:", e$message)); NULL })
    if (is.null(train_swd)) { rm(background_points, species_specific_stack); gc(); return(list(status = "error_swd_preparation", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("SWD prep failed."))) }
    slog("DEBUG", "SWD object prepared for tuning/training.")
    
    # --- <<< SPLIT SWD INTO TRAINING AND TESTING FOR INDEPENDENT EVALUATION >>> ---
    slog("INFO", "Splitting SWD data into training and testing sets (e.g., 80/20)...")
    test_proportion <- config$sdm_test_proportion %||% 0.2 # Add sdm_test_proportion to config, default 0.2
    
    # Ensure reproducibility of the split
    split_seed <- config$global_seed + species_aphia_id # Species-specific seed for the split
    
    swd_parts <- tryCatch({
      SDMtune::trainValTest(train_swd, test = test_proportion, only_presence = FALSE, seed = split_seed)
    }, error = function(e) {
      slog("ERROR", paste("Failed to split SWD data:", e$message)); NULL
    })
    
    if (is.null(swd_parts) || length(swd_parts) != 2) {
      msg <- "Data splitting into train/test failed. Halting species processing."
      slog("ERROR", msg)
      rm(background_points_full, species_specific_stack, initial_swd_data); gc()
      return(list(status = "error_swd_split", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    } else {
      train_swd <- swd_parts[[1]]
      test_swd  <- swd_parts[[2]]
      slog("INFO", paste("SWD data split: Train (", nrow(train_swd@data), " total points), Test (", nrow(test_swd@data), " total points)"))
      slog("DEBUG", paste("  Train Presences:", sum(train_swd@pa == 1), "Absences:", sum(train_swd@pa == 0)))
      slog("DEBUG", paste("  Test Presences:", sum(test_swd@pa == 1), "Absences:", sum(test_swd@pa == 0)))
      if (nrow(train_swd@data) < config$min_occurrences_sdm * 2 || sum(train_swd@pa == 1) < config$min_occurrences_sdm ) { 
        slog("WARN", "Training set size or presence count after split might be too small (Total:", nrow(train_swd@data), "Pres:", sum(train_swd@pa==1), "). Model performance may be affected.")
        # Decide if this should be a fatal error
        # msg <- "Training set too small after split."; slog("ERROR", msg)
        # return(list(status = "error_train_set_too_small", species = species_name, message = msg))
      }
      if (nrow(test_swd@data) == 0 || sum(test_swd@pa == 1) == 0 || sum(test_swd@pa == 0) == 0) {
        slog("WARN", "Test set is empty or lacks presences/absences. Independent test metrics will be NA.")
        test_swd <- NULL # Set to NULL if unusable
      }
    }
    # --- <<< END SWD SPLIT >>> ---
    
    slog("DEBUG", "Creating spatial folds...")
    spatial_folds <- create_spatial_cv_folds_simplified(train_swd, species_specific_stack, config, logger=NULL, species_log_file)
    if(is.null(spatial_folds)) {
      msg <- "Failed to create spatial folds."; slog("ERROR", msg); rm(background_points, species_specific_stack, train_swd); gc(); return(list(status = "error_cv_folds", species=species_name, occurrence_count = occurrence_count_after_thinning, message=msg))
    }
    
    slog("INFO", "Starting hyperparameter tuning.")
    tuning_output <- run_sdm_tuning_scv(occs_coords, species_specific_stack, background_points, config, logger=NULL, species_name, species_log_file=species_log_file)
    if (is.null(tuning_output) || is.null(attr(tuning_output, "best_hypers"))) {
      msg <- paste0("Skipping: Tuning failed."); slog("ERROR", msg); rm(background_points, species_specific_stack, train_swd, spatial_folds); gc(); return(list(status = "error_tuning", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    best_hypers <- attr(tuning_output, "best_hypers")
    
    # if(!save_tuning_results(tuning_output, species_name_sanitized, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file, group_name=group_name)) {
    #   rm(background_points, species_specific_stack, train_swd, spatial_folds, tuning_output); gc(); return(list(status = "error_saving_tuning_results", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save tuning results.")))
    # }
    
    slog("INFO", "Starting final model training.")
    final_model <- train_final_sdm(occs_coords, species_specific_stack, background_points, best_hypers, config, logger=NULL, species_name, species_log_file=species_log_file)
    if(is.null(final_model)){
      msg <- paste0("Skipping: Training failed."); slog("ERROR", msg); rm(background_points, species_specific_stack, train_swd, spatial_folds, tuning_output); gc(); return(list(status = "error_training", species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
    slog("INFO", "Final model training complete.")
    
    if(!save_final_model(final_model, species_name_sanitized, predictor_type_suffix, group_name, config, logger=NULL, species_log_file=species_log_file)) {
      rm(background_points, species_specific_stack, train_swd, spatial_folds, tuning_output, final_model); gc(); return(list(status = "error_saving_model", species = species_name, occurrence_count = occurrence_count_after_thinning, message = paste0("Failed save final model.")))
    }
    rm(background_points, spatial_folds); gc()
  }
  
  # --- Log Metrics ---
  if (!is.null(final_model) && !is.null(species_specific_stack)) {
    log_final_model_metrics(
      final_model = final_model, train_swd = train_swd, test_swd = test_swd,
      tuning_predictor_stack = species_specific_stack, tuning_output = tuning_output,
      species_name_sanitized = species_name_sanitized, group_name = group_name,
      predictor_type_suffix = predictor_type_suffix, config = config, logger = NULL,
      species_log_file = species_log_file )
  } else { slog("WARN", "Final model or species-specific stack object is NULL, cannot log metrics.") }
  
  # --- Variable Importance ---
  if (!is.null(final_model)) {
    if(!is.null(train_swd)){
      slog("INFO", "Calculating and saving variable importance...")
      vi_success <- calculate_and_save_vi(
        final_model = final_model, training_swd = train_swd,
        species_name_sanitized = species_name_sanitized, group_name = group_name,
        predictor_type_suffix = predictor_type_suffix, config = config, logger = NULL,
        species_log_file = species_log_file)
      if(!vi_success) slog("WARN", "Variable importance calculation/saving failed.")
    } else { slog("WARN", "Skipping variable importance: training SWD data is missing.") }
  } else { slog("WARN", "Skipping variable importance: final model object is unavailable.") }
  
  # --- Prediction Loop ---
  # Requires GLOBAL predictor stacks
  predictions_made = 0; prediction_errors = 0
  scenarios_to_predict <- names(env_predictor_paths_list)
  if (config$do_final_prediction) {
    slog("INFO", paste("Starting predictions for", length(scenarios_to_predict), "scenarios."))
    
    if (!is.null(final_model)) {
      for (pred_scenario in scenarios_to_predict) {
        slog("DEBUG", paste("  Predicting scenario:", pred_scenario))
        target_pred_file_path <- construct_prediction_filename(species_name_sanitized, pred_scenario, predictor_type_suffix, config)
        if (is.null(target_pred_file_path)) { slog("WARN", "  Could not construct prediction filename. Skipping."); prediction_errors <- prediction_errors + 1; next }
        rerun_flag <- config$force_rerun$run_biotic_sdms # Check biotic flag
        if (!rerun_flag && file.exists(target_pred_file_path)) { slog("DEBUG", "  Prediction exists. Skipping."); next }
        
        # --- Prepare GLOBAL Stack for Prediction Scenario ---
        pred_predictor_stack_global <- NULL
        tryCatch({
          pred_env_path <- env_predictor_paths_list[[pred_scenario]]
          if (is.null(pred_env_path) || !file.exists(pred_env_path)) stop(paste("Env PCA stack missing:", pred_scenario))
          pred_env_stack <- terra::rast(pred_env_path) # Load Env PCA
          reference_geom_pred <- pred_env_stack[[1]]
          
          host_paths_for_pred <- all_host_prediction_paths[[pred_scenario]][names(all_host_prediction_paths[[pred_scenario]]) %in% associated_hosts_sanitized]
          if (length(host_paths_for_pred) == 0) stop(paste("No host paths found:", pred_scenario))
          
          host_rasters_pred_list <- list()
          for(host_name in names(host_paths_for_pred)){
            host_path <- host_paths_for_pred[[host_name]]
            host_rast <- tryCatch(terra::rast(host_path), error=function(e){slog("WARN",paste("Failed load host",basename(host_path),":",e$message)); NULL})
            if(!is.null(host_rast)){
              if (!terra::compareGeom(reference_geom_pred, host_rast, stopOnError=FALSE, res=TRUE)) { slog("DEBUG", paste("Resampling host", basename(host_path))); host_rast <- tryCatch(terra::resample(host_rast, reference_geom_pred, method="bilinear"), error=function(e){slog("ERROR",paste("Failed resample",basename(host_path),":",e$message)); NULL}) }
              if(!is.null(host_rast)){ names(host_rast) <- host_name; host_rasters_pred_list[[host_name]] <- host_rast }
            }
          }
          if (length(host_rasters_pred_list) == 0) stop(paste("Failed load/resample hosts:", pred_scenario))
          
          host_stack_scenario <- terra::rast(host_rasters_pred_list)
          max_host_suitability_scenario <- terra::app(host_stack_scenario, fun = "max", na.rm = TRUE); names(max_host_suitability_scenario) <- "host_suitability_max"
          pred_predictor_stack_global <- c(pred_env_stack, max_host_suitability_scenario) # Assign to outer scope
          slog("DEBUG", paste("  GLOBAL Predictor stack for prediction:", paste(names(pred_predictor_stack_global), collapse=", ")))
          rm(pred_env_stack, host_stack_scenario, max_host_suitability_scenario, host_rasters_pred_list); gc()
          
        }, error = function(e) {
          slog("ERROR", paste("  Failed preparing GLOBAL combined stack for prediction:", pred_scenario, "Error:", e$message))
          if (exists("pred_env_stack")) rm(pred_env_stack); if (exists("host_stack_scenario")) rm(host_stack_scenario); if (exists("max_host_suitability_scenario")) rm(max_host_suitability_scenario); if (exists("host_rasters_pred_list")) rm(host_rasters_pred_list)
          gc(); pred_predictor_stack_global <- NULL
        })
        # --- End GLOBAL Stack Prep ---
        
        if(is.null(pred_predictor_stack_global)) { prediction_errors <- prediction_errors + 1; next }
        
        # Predict & Save
        prediction_output <- predict_sdm_suitability(final_model, pred_predictor_stack_global, config, logger=NULL, species_log_file=species_log_file)
        if (inherits(prediction_output, "SpatRaster")) {
          save_success <- save_sdm_prediction(prediction_output, species_name_sanitized, pred_scenario, predictor_type_suffix, config, logger=NULL, species_log_file=species_log_file)
          if(save_success) predictions_made <- predictions_made + 1 else prediction_errors <- prediction_errors + 1
          rm(prediction_output); gc()
        } else { slog("WARN", paste("  Prediction failed:", prediction_output)); prediction_errors <- prediction_errors + 1 }
        rm(pred_predictor_stack_global); gc()
      } # End prediction scenario loop
    } else {
      prediction_errors <- length(scenarios_to_predict); msg <- paste0("Skipping predictions: Final model unavailable."); slog("ERROR", msg)
      status <- if(load_existing_model) "error_loading_model" else "error_training"; return(list(status = status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = msg))
    }
  }
  
  # --- Prepare return status ---
  final_status <- "success"; status_message <- paste0("Finished Combined SDM. Occs:", occurrence_count_after_thinning, ". Preds attempted:", length(scenarios_to_predict), ". Made:", predictions_made, ". Errors/Skipped:", prediction_errors, ".")
  if (prediction_errors > 0 && predictions_made == 0) { final_status <- "error_prediction_all"; slog("ERROR", status_message) }
  else if (prediction_errors > 0) { final_status <- "success_with_pred_errors"; slog("WARN", status_message) }
  else { slog("INFO", status_message) }
  
  # --- Clean up ---
  if (!is.null(species_specific_stack)) rm(species_specific_stack)
  if (!is.null(tuning_predictor_stack_global)) rm(tuning_predictor_stack_global)
  if (!is.null(occs_coords)) rm(occs_coords)
  if (!is.null(train_swd)) rm(train_swd)
  if (!load_existing_model && !is.null(tuning_output) && inherits(tuning_output, "SDMtune")) rm(tuning_output)
  if (exists("final_model", inherits=FALSE)) rm(final_model)
  gc()
  
  return(list(status = final_status, species = species_name, occurrence_count = occurrence_count_after_thinning, message = status_message))
} # End process_species_sdm_combined function


# --- 12. Setup Parallel Backend & Run ---
# (Keep existing logic)
if (config$use_parallel && config$num_cores > 1) { log4r::info(logger, paste("Setting up parallel backend with", config$num_cores, "cores.")); gc(full=TRUE); future::plan(future::multisession, workers = config$num_cores, gc = TRUE) } else { log4r::info(logger, "Running sequentially."); future::plan(future::sequential) }
progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")

log4r::info(logger, paste("Starting Combined SDM processing for", nrow(species_df), group_name, "species..."))

results_list <- progressr::with_progress({
  furrr::future_map(1:nrow(species_df), ~{
    process_species_sdm_combined(species_row = species_df[.x, ], config = config, env_predictor_paths_list = env_predictor_paths_list, all_host_prediction_paths = all_host_prediction_paths, associations_df = associations_df,
                                 group_name = group_name, predictor_type_suffix = predictor_type_suffix, occurrence_dir = occurrence_dir, tuning_scenario = "current")
  }, .options = furrr_options(seed = TRUE, scheduling = 1.0))
})

log4r::info(logger, "Parallel/sequential processing complete.")

# --- 13. Process Results ---
# (Keep existing summary logic)
occurrence_counts <- list(); success_count <- 0; error_count <- 0; skipped_count <- 0; partial_success_count <- 0
log4r::info(logger, paste("--- Processing Results Summary (", group_name, predictor_type_suffix, "Models) ---"))
for (res in results_list) {
  if (is.null(res)) { error_count <- error_count + 1; log4r::error(logger, "Received NULL result."); next }
  log_level_func <- switch(res$status, success=log4r::info, success_with_pred_errors=log4r::warn, skipped_occurrences_raw=log4r::warn, skipped_occurrences_sac=log4r::warn, skipped_no_hosts=log4r::warn, log4r::error)
  log_level_func(logger, res$message)
  occurrence_counts[[res$species]] <- if(!is.null(res$occurrence_count)) res$occurrence_count else NA
  if (res$status == "success") success_count <- success_count + 1
  else if (res$status == "success_with_pred_errors") partial_success_count <- partial_success_count + 1
  else if (grepl("skipped", res$status)) skipped_count <- skipped_count + 1
  else error_count <- error_count + 1
}
log4r::info(logger, paste("--- Overall Summary (", group_name, predictor_type_suffix, "Models) ---"))
log4r::info(logger, paste("Total Species Targets:", length(results_list)))
log4r::info(logger, paste("Fully Successful Runs:", success_count))
log4r::info(logger, paste("Partially Successful Runs:", partial_success_count))
log4r::info(logger, paste("Skipped:", skipped_count))
log4r::info(logger, paste("Errors:", error_count))
if (length(occurrence_counts) > 0) {
  occ_count_df <- data.frame( Species = names(occurrence_counts), OccurrenceCountAfterThinning = unlist(occurrence_counts)) %>% dplyr::arrange(Species)
  occ_count_df_print <- occ_count_df; occ_count_df_print$OccurrenceCountAfterThinning[is.na(occ_count_df_print$OccurrenceCountAfterThinning)] <- "N/A"
  occ_count_file <- file.path(config$log_dir_base, paste0("occurrence_counts_", group_name, predictor_type_suffix, ".csv"))
  tryCatch({ readr::write_csv(occ_count_df, occ_count_file); log4r::info(logger, paste("Occurrence counts saved:", occ_count_file))}, error = function(e) { log4r::error(logger, paste("Failed save counts:", e$message)) })
  cat("\n--- Final Occurrence Counts (", group_name, predictor_type_suffix, ", after final thinning) ---\n"); print(occ_count_df_print)
} else { log4r::warn(logger, "No occurrence counts were recorded.") }

future::plan(future::sequential); gc(full=TRUE)
log4r::info(logger, "--- Script 06d (v7 - Tweaking Integrated) finished. ---")
#-------------------------------------------------------------------------------
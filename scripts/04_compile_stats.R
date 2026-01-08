# scripts/04_compile_stats.R
# Aggregates the 50 iterations into a single CSV for analysis

library(dplyr)
library(readr)
library(stringr)

res_dir <- "outputs/models/iterations"

# --- FIX: Only list .rds files (ignore .txt logs) ---
files <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)

if(length(files) == 0) stop("No .rds files found in ", res_dir)

results_df <- data.frame()

cat("Compiling", length(files), "model result files...\n")

for (f in files) {
  # Parse filename for metadata
  fname <- basename(f)
  
  # Filename structure: Species_Name_type_iter_X.rds
  # We split by "_iter_" to separate the species/model part from the iteration part
  parts <- strsplit(fname, "_iter_")[[1]]
  
  if(length(parts) < 2) {
    warning("Skipping malformed filename: ", fname)
    next
  }
  
  sp_model <- parts[1] # e.g. Amphiprion_clarkii_env
  iter_part <- parts[2] # e.g. 1.rds
  iter <- as.numeric(gsub(".rds", "", iter_part))
  
  # Determine type logic
  if (grepl("_env$", sp_model)) {
    type <- "env_only"
    sp <- gsub("_env$", "", sp_model)
  } else if (grepl("_biotic$", sp_model)) {
    type <- "biotic_only"
    sp <- gsub("_biotic$", "", sp_model)
  } else if (grepl("_combined$", sp_model)) {
    type <- "combined"
    sp <- gsub("_combined$", "", sp_model)
  } else {
    warning("Unknown model type in filename: ", fname)
    next
  }
  
  # Read Data safely
  tryCatch({
    dat <- readRDS(f)
    
    # Check if dat is valid (it should be a dataframe with metrics)
    if(is.null(dat) || nrow(dat) == 0) {
      warning("Empty data in file: ", fname)
    } else {
      # Extract metrics (adapt column names if your ENMeval version differs)
      # Usually: auc.val.avg, cbi.val.avg, AICc
      row_data <- data.frame(
        species = sp,
        model_type = type,
        iteration = iter,
        AUC_test = if("auc.val.avg" %in% names(dat)) dat$auc.val.avg else NA,
        CBI_test = if("cbi.val.avg" %in% names(dat)) dat$cbi.val.avg else NA,
        AICc = if("AICc" %in% names(dat)) dat$AICc else NA
      )
      
      results_df <- bind_rows(results_df, row_data)
    }
  }, error = function(e) {
    warning("Failed to read ", fname, ": ", e$message)
  })
}

# Save final table
out_path <- "outputs/tables/model_comparison_50iter.csv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_csv(results_df, out_path)

cat("Success! Compiled stats saved to:", out_path, "\n")
cat("Total Rows:", nrow(results_df), "\n")
print(head(results_df))
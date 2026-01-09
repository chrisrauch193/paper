# scripts/99_env_justification.R
# ------------------------------------------------------------------------------
# ENVIRONMENTAL JUSTIFICATION & VIF ANALYSIS
# ------------------------------------------------------------------------------
# 1. Runs PCA on Current Climate Variables.
# 2. Generates Scree Plot (Variance Explained) to justify choosing 5 Axes.
# 3. Generates Variable Contribution Plots (Interpretation of PC1/PC2).
# 4. Checks VIF of the Final Stack (PC1-5 + Rugosity).
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(terra, dplyr, readr, ggplot2, factoextra, usdm, corrplot)

# --- CONFIG ---
BASE_DIR    <- getwd()
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_ENV_DIR <- file.path(DATA_DIR, "env", "current")
OUTPUT_PLOT <- file.path(BASE_DIR, "outputs", "justification_plots")
dir.create(OUTPUT_PLOT, recursive = TRUE, showWarnings = FALSE)

CLIMATE_VARS <- c("sws_baseline_depthsurf_mean", "so_baseline_depthmax_mean", 
                  "thetao_baseline_depthmax_mean", "no3_baseline_depthmax_mean", 
                  "no3_baseline_depthmax_range", "chl_baseline_depthmax_mean", 
                  "o2_baseline_depthmax_range", "phyc_baseline_depthmax_mean")

# --- 1. LOAD DATA ---
cat("Loading Data...\n")
clim_list <- list()
for(v in CLIMATE_VARS) {
  r <- terra::rast(file.path(RAW_ENV_DIR, paste0(v, ".tif")))
  names(r) <- v
  clim_list[[v]] <- r
}
clim_stack <- terra::rast(clim_list)

# Sample for PCA (50k points)
set.seed(42)
samp <- terra::spatSample(clim_stack, size=50000, method="random", na.rm=TRUE, xy=FALSE)
samp <- samp[complete.cases(samp), ]

# --- 2. RUN PCA ---
cat("Running PCA...\n")
pca_res <- prcomp(samp, center=TRUE, scale.=TRUE)

# --- 3. JUSTIFICATION PLOTS ---

# A. Scree Plot (Variance Explained)
# This proves why you chose 5 axes (usually >90% variance)
p1 <- fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "Variance Explained by PC Axes", 
       subtitle = "Justification for selecting top 5 PCs")
ggsave(file.path(OUTPUT_PLOT, "01_scree_plot.png"), p1, width=8, height=6)

# B. Variable Contributions (What is PC1?)
# This helps you explain "PC1 represents Temperature/Nutrients"
p2 <- fviz_pca_var(pca_res, col.var = "contrib", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   repel = TRUE) +
  labs(title = "Variable Contributions to PC1 & PC2")
ggsave(file.path(OUTPUT_PLOT, "02_variable_contributions.png"), p2, width=8, height=8)

# --- 4. VIF ANALYSIS (The "Defensible" Check) ---
cat("Calculating VIF for Final Predictors...\n")

# Predict PC axes for the sample
pca_scores <- predict(pca_res, samp)
df_pca <- as.data.frame(pca_scores[, 1:5]) # Take top 5

# Add Dummy Rugosity (Random noise for check, or load real if needed)
# Since we know Rugosity is spatially independent, we can simulate or load real.
# Loading real is better:
rug_r <- terra::rast(file.path(DATA_DIR, "env", "terrain", "rugosity.tif"))
# We need to extract rugosity at the SAME points as the sample
# This requires coordinates, but our 'samp' lost them.
# Let's re-sample specifically for VIF.

vif_samp <- terra::spatSample(c(clim_stack, rug_r), size=10000, method="random", na.rm=TRUE, xy=FALSE)
vif_samp <- vif_samp[complete.cases(vif_samp), ]

# Transform Climate cols to PCs
clim_only <- vif_samp[, CLIMATE_VARS]
pca_vif   <- predict(pca_res, clim_only)[, 1:5] # PC1-5

# Combine PC1-5 + Rugosity
final_predictors <- data.frame(pca_vif, rugosity = vif_samp[, "rugosity"])

# Calculate VIF
vif_res <- usdm::vif(final_predictors)
print(vif_res)

# Save VIF Table
write.csv(vif_res, file.path(OUTPUT_PLOT, "03_final_vif_scores.csv"), row.names=FALSE)

cat("--- FINISHED ---\n")
cat("Check 'outputs/justification_plots/' for your thesis figures.\n")
# scripts/01_functions_core.R
# ------------------------------------------------------------------------------
# CORE UTILITIES (Thinning, Background, Logging, Status helpers)
# ------------------------------------------------------------------------------

write_log <- function(path, msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", timestamp, "] ", msg, "\n"), file = path, append = TRUE)
}

thin_occurrences <- function(occ_df, env_rast) {
  cells <- terra::cellFromXY(env_rast, as.matrix(occ_df[, c("x", "y")]))
  dups <- duplicated(cells)
  occ_df[!dups, ]
}

# ---- Geo distance helpers (meters) ----
.haversine_m <- function(lon1, lat1, lon2, lat2) {
  R <- 6371000.0
  to_rad <- pi / 180.0
  lon1 <- lon1 * to_rad; lat1 <- lat1 * to_rad
  lon2 <- lon2 * to_rad; lat2 <- lat2 * to_rad
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  2 * R * atan2(sqrt(a), sqrt(pmax(1e-15, 1 - a)))
}

# ---- Spatial fold assignment (simple, stable, defensible) ----
make_block_quadrant_folds <- function(coords_xy, ref_xy = coords_xy) {
  xy <- as.matrix(coords_xy[, 1:2])
  ref <- as.matrix(ref_xy[, 1:2])
  mx <- stats::median(ref[, 1], na.rm = TRUE)
  my <- stats::median(ref[, 2], na.rm = TRUE)
  
  fx <- ifelse(xy[, 1] <= mx, 1L, 2L)
  fy <- ifelse(xy[, 2] <= my, 0L, 2L)
  as.integer(fx + fy)  # 1..4
}

assign_spatial_folds <- function(pres_xy, bg_xy, method = "block_quadrant", k = 4L) {
  method <- tolower(method)
  if (method == "random") {
    set.seed(1)
    pres_fold <- sample(seq_len(k), nrow(pres_xy), replace = TRUE)
    bg_fold   <- sample(seq_len(k), nrow(bg_xy),   replace = TRUE)
    return(list(pres_fold = pres_fold, bg_fold = bg_fold, method = "random", k = k))
  }
  
  if (k != 4L) stop("block_quadrant requires k=4")
  pres_fold <- make_block_quadrant_folds(pres_xy, ref_xy = pres_xy)
  bg_fold   <- make_block_quadrant_folds(bg_xy,   ref_xy = pres_xy)  # use presence medians
  list(pres_fold = pres_fold, bg_fold = bg_fold, method = "block_quadrant", k = 4L)
}

# ---- Boyce-like index from predicted values ----
boyce_index_values <- function(pred_pres, pred_bg, n_bins = 20L) {
  p <- as.numeric(pred_pres); a <- as.numeric(pred_bg)
  p <- p[is.finite(p)]
  a <- a[is.finite(a)]
  if (length(p) < 5 || length(a) < 50) return(NA_real_)
  
  rng <- range(c(p, a), na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) <= 0) return(NA_real_)
  
  bins <- seq(rng[1], rng[2], length.out = n_bins + 1L)
  hp  <- graphics::hist(p, breaks = bins, plot = FALSE)$counts
  ha  <- graphics::hist(a, breaks = bins, plot = FALSE)$counts
  
  exp <- ha / sum(ha)
  obs <- hp / sum(hp)
  valid <- is.finite(exp) & is.finite(obs) & exp > 0
  if (sum(valid) < 3) return(NA_real_)
  
  ratio <- obs[valid] / exp[valid]
  mids  <- (bins[-1] + bins[-length(bins)]) / 2
  mids  <- mids[valid]
  
  suppressWarnings(stats::cor(mids, ratio, method = "spearman"))
}

# ------------------------------------------------------------------------------
# Background sampling
# ------------------------------------------------------------------------------

get_bias_corrected_background <- function(
    occ_coords, env_stack,
    n_bg = 10000, alpha = 0.5, method = "paper_exact",
    env_cols = NULL, cand_mult = 3L, geo_metric = "haversine", seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  if (inherits(occ_coords, "data.frame")) occ <- as.matrix(occ_coords[, c("x", "y")])
  else occ <- as.matrix(occ_coords)
  
  if (is.null(env_cols)) env_cols <- names(env_stack)[1:2]
  env_cols <- intersect(env_cols, names(env_stack))
  if (length(env_cols) < 1) stop("No valid env_cols found for bias background.")
  
  occ_env <- terra::extract(env_stack[[env_cols]], occ)
  occ_env <- na.omit(cbind(x = occ[, 1], y = occ[, 2], occ_env))
  if (nrow(occ_env) == 0) return(NULL)
  
  candidates <- terra::spatSample(
    env_stack[[env_cols]],
    size = as.integer(n_bg * cand_mult),
    method = "random", na.rm = TRUE, xy = TRUE, values = TRUE
  )
  candidates <- na.omit(candidates)
  if (nrow(candidates) < n_bg) return(NULL)
  
  cand_xy <- as.matrix(candidates[, c("x", "y")])
  occ_xy  <- as.matrix(occ_env[, c("x", "y")])
  
  cand_env <- as.matrix(candidates[, env_cols, drop = FALSE])
  occ_envv <- as.matrix(occ_env[, env_cols, drop = FALSE])
  
  geo_metric <- tolower(geo_metric)
  if (geo_metric == "haversine") {
    dist_geo <- apply(cand_xy, 1, function(pt) {
      d <- .haversine_m(pt[1], pt[2], occ_xy[, 1], occ_xy[, 2])
      if (method == "paper_exact") mean(d) else min(d)
    })
  } else {
    dist_geo <- apply(cand_xy, 1, function(pt) {
      d <- sqrt((occ_xy[, 1] - pt[1])^2 + (occ_xy[, 2] - pt[2])^2)
      if (method == "paper_exact") mean(d) else min(d)
    })
  }
  
  dist_env <- apply(cand_env, 1, function(pt) {
    d <- sqrt(rowSums((occ_envv - matrix(pt, nrow = nrow(occ_envv), ncol = ncol(occ_envv), byrow = TRUE))^2))
    if (method == "paper_exact") mean(d) else min(d)
  })
  
  d_geo_norm <- dist_geo / max(dist_geo, na.rm = TRUE)
  d_env_norm <- dist_env / max(dist_env, na.rm = TRUE)
  
  sampling_prob <- 1 - ((1 - d_geo_norm)^alpha) * ((1 - d_env_norm)^(1 - alpha))
  sampling_prob[!is.finite(sampling_prob)] <- 0
  
  selected_idx <- sample(seq_len(nrow(candidates)), size = n_bg, prob = sampling_prob, replace = FALSE)
  candidates[selected_idx, c("x", "y")]
}

get_random_background <- function(occ_coords, env_stack, n_bg = 10000, buffer_m = 1000000) {
  if (inherits(occ_coords, "data.frame")) coords_mat <- as.matrix(occ_coords[, c("x", "y")])
  else coords_mat <- as.matrix(occ_coords)
  
  vect_occ <- terra::vect(coords_mat, crs = "EPSG:4326", type = "points")
  
  # Buffer in meters by projecting to EPSG:3857
  vect_occ_m  <- terra::project(vect_occ, "EPSG:3857")
  vect_buff_m <- terra::aggregate(terra::buffer(vect_occ_m, width = buffer_m))
  vect_buff   <- terra::project(vect_buff_m, "EPSG:4326")
  
  env_crop <- terra::crop(env_stack, vect_buff)
  env_mask <- terra::mask(env_crop, vect_buff)
  
  terra::spatSample(env_mask, size = n_bg, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
}

# ------------------------------------------------------------------------------
# Biotic layer (host-weighted)
# ------------------------------------------------------------------------------

get_biotic_layer <- function(fish_sp, host_stack, int_mat, debug_path = NULL) {
  fish_clean <- gsub(" ", "_", fish_sp)
  if (!fish_clean %in% rownames(int_mat)) {
    if (!is.null(debug_path)) write_log(debug_path, paste("DEBUG BIOTIC:", fish_clean, "not in matrix."))
    return(NULL)
  }
  
  row_idx <- which(rownames(int_mat) == fish_clean)
  w_vals <- as.numeric(int_mat[row_idx, ])
  names(w_vals) <- colnames(int_mat)
  weights <- w_vals[w_vals > 0]
  
  available_hosts <- intersect(names(weights), names(host_stack))
  if (length(available_hosts) == 0) {
    if (!is.null(debug_path)) {
      write_log(debug_path, paste("DEBUG BIOTIC FAIL:", fish_clean))
      write_log(debug_path, paste("  > Need:", paste(names(weights), collapse = ", ")))
      write_log(debug_path, paste("  > Have:", paste(names(host_stack), collapse = ", ")))
    }
    return(NULL)
  }
  
  sub_stack <- host_stack[[available_hosts]]
  sub_weights <- as.numeric(weights[available_hosts])
  biotic_layer <- terra::weighted.mean(sub_stack, w = sub_weights, na.rm = TRUE)
  names(biotic_layer) <- "biotic_suitability"
  biotic_layer
}

# ------------------------------------------------------------------------------
# STATUS HELPERS (for robust skipping)
# ------------------------------------------------------------------------------

read_status_file <- function(path) {
  if (!file.exists(path)) return(list(status = "NONE", n = 0L))
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(txt) == 0) return(list(status = "NONE", n = 0L))
  line <- txt[1]
  
  if (grepl("^SKIP_NO_BIOTIC", line)) return(list(status = "SKIP_NO_BIOTIC", n = 0L))
  
  n <- NA_integer_
  m <- regmatches(line, regexpr("n=\\d+", line))
  if (length(m) > 0 && nzchar(m)) n <- as.integer(sub("n=", "", m))
  
  if (grepl("^OK", line)) return(list(status = "OK", n = n))
  list(status = "UNKNOWN", n = n)
}

write_status_ok <- function(path, n) {
  writeLines(paste0("OK n=", as.integer(n)), path)
}

write_status_progress <- function(path, completed, target) {
  writeLines(paste0("PROGRESS n=", as.integer(completed), "/", as.integer(target)), path)
}

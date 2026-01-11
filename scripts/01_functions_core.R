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
  return(occ_df[!dups, ])
}

get_bias_corrected_background <- function(occ_coords, env_stack, n_bg=10000, alpha=0.5, method="paper_exact") {
  if(inherits(occ_coords, "data.frame")) occ <- as.matrix(occ_coords[, c("x", "y")])
  else occ <- as.matrix(occ_coords)
  
  env_cols <- names(env_stack)[1:2]
  occ_env <- terra::extract(env_stack[[env_cols]], occ)
  occ_env <- na.omit(cbind(occ, occ_env))
  if(nrow(occ_env) == 0) return(NULL)
  
  candidates <- terra::spatSample(env_stack, size=n_bg*3, method="random", na.rm=TRUE, xy=TRUE, values=TRUE)
  candidates <- na.omit(candidates)
  
  cand_xy <- as.matrix(candidates[, c("x", "y")])
  occ_xy  <- as.matrix(occ_env[, c("x", "y")])
  cand_env <- as.matrix(candidates[, env_cols])
  occ_env_vals <- as.matrix(occ_env[, env_cols])
  
  if(method == "paper_exact") {
    dist_geo <- apply(cand_xy, 1, function(pt) mean(sqrt((occ_xy[,1]-pt[1])^2 + (occ_xy[,2]-pt[2])^2)))
    dist_env <- apply(cand_env, 1, function(pt) mean(sqrt((occ_env_vals[,1]-pt[1])^2 + (occ_env_vals[,2]-pt[2])^2)))
  } else {
    dist_geo <- apply(cand_xy, 1, function(pt) min(sqrt((occ_xy[,1]-pt[1])^2 + (occ_xy[,2]-pt[2])^2)))
    dist_env <- apply(cand_env, 1, function(pt) min(sqrt((occ_env_vals[,1]-pt[1])^2 + (occ_env_vals[,2]-pt[2])^2)))
  }
  
  d_geo_norm <- dist_geo / max(dist_geo, na.rm=TRUE)
  d_env_norm <- dist_env / max(dist_env, na.rm=TRUE)
  sampling_prob <- 1 - ((1 - d_geo_norm)^alpha) * ((1 - d_env_norm)^(1 - alpha))
  
  if(nrow(candidates) > n_bg) {
    selected_idx <- sample(1:nrow(candidates), size=n_bg, prob=sampling_prob, replace=FALSE)
    bg_final <- candidates[selected_idx, c("x", "y")]
  } else {
    bg_final <- candidates[, c("x", "y")]
  }
  return(bg_final)
}

get_random_background <- function(occ_coords, env_stack, n_bg=10000) {
  if(inherits(occ_coords, "data.frame")) coords_mat <- as.matrix(occ_coords[, c("x", "y")])
  else coords_mat <- as.matrix(occ_coords)
  vect_occ <- terra::vect(coords_mat, crs="EPSG:4326", type="points")
  vect_buff <- terra::aggregate(terra::buffer(vect_occ, width=1000000))
  env_crop <- terra::crop(env_stack, vect_buff)
  env_mask <- terra::mask(env_crop, vect_buff)
  bg_coords <- terra::spatSample(env_mask, size=n_bg, method="random", na.rm=TRUE, xy=TRUE, values=FALSE)
  return(bg_coords)
}

get_biotic_layer <- function(fish_sp, host_stack, int_mat, debug_path=NULL) {
  fish_clean <- gsub(" ", "_", fish_sp)
  if(!fish_clean %in% rownames(int_mat)) {
    if(!is.null(debug_path)) write_log(debug_path, paste("DEBUG BIOTIC:", fish_clean, "not in matrix."))
    return(NULL)
  }
  
  row_idx <- which(rownames(int_mat) == fish_clean)
  w_vals <- as.numeric(int_mat[row_idx, ])
  names(w_vals) <- colnames(int_mat)
  weights <- w_vals[w_vals > 0]
  
  available_hosts <- intersect(names(weights), names(host_stack))
  if(length(available_hosts) == 0) {
    if(!is.null(debug_path)) {
      write_log(debug_path, paste("DEBUG BIOTIC FAIL:", fish_clean))
      write_log(debug_path, paste("  > Need:", paste(names(weights), collapse=", ")))
      write_log(debug_path, paste("  > Have:", paste(names(host_stack), collapse=", ")))
    }
    return(NULL)
  }
  
  sub_stack <- host_stack[[available_hosts]]
  sub_weights <- as.numeric(weights[available_hosts])
  biotic_layer <- terra::weighted.mean(sub_stack, w=sub_weights, na.rm=TRUE)
  names(biotic_layer) <- "biotic_suitability"
  return(biotic_layer)
}

# ------------------------------------------------------------------------------
# STATUS HELPERS (for robust skipping)
# ------------------------------------------------------------------------------

read_status_file <- function(path) {
  if(!file.exists(path)) return(list(status="NONE", n=0L))
  txt <- tryCatch(readLines(path, warn=FALSE), error=function(e) character(0))
  if(length(txt) == 0) return(list(status="NONE", n=0L))
  line <- txt[1]
  
  if(grepl("^SKIP_NO_BIOTIC", line)) return(list(status="SKIP_NO_BIOTIC", n=0L))
  
  n <- NA_integer_
  m <- regmatches(line, regexpr("n=\\d+", line))
  if(length(m) > 0 && nzchar(m)) n <- as.integer(sub("n=", "", m))
  
  if(grepl("^OK", line)) return(list(status="OK", n=n))
  return(list(status="UNKNOWN", n=n))
}

write_status_ok <- function(path, n) {
  writeLines(paste0("OK n=", as.integer(n)), path)
}

write_status_progress <- function(path, completed, target) {
  writeLines(paste0("PROGRESS n=", as.integer(completed), "/", as.integer(target)), path)
}

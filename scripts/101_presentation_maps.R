# scripts/101_presentation_maps.R
# ------------------------------------------------------------------------------
# PRESENTATION MAPS: BASELINE vs. EXPANSION (Zoomed Out & Greyed Logic)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, readr, ggrepel, viridis, patchwork, stringr)

# --- CONFIG ---
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs", "figures")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)

# FILES
FILE_STRICT <- file.path(DATA_DIR, "meow_ecos_df_strict.csv") # Baseline
FILE_EXPAND <- file.path(DATA_DIR, "meow_ecos_df.csv")        # Expansion

# --- 1. DATA PREP FUNCTION ---
prepare_presentation_data <- function() {
  
  if(!file.exists(FILE_EXPAND) || !file.exists(FILE_STRICT)) {
    stop("Missing CSV files. Please run pipeline 0c first.")
  }
  
  # Load Data
  df_expand <- read_csv(FILE_EXPAND, show_col_types = FALSE)
  df_strict <- read_csv(FILE_STRICT, show_col_types = FALSE)
  
  # Identify Strict Provinces
  strict_provs <- unique(df_strict$PROVINCE)
  
  # Standardize Expansion Data
  colnames(df_expand) <- tolower(colnames(df_expand))
  df <- df_expand %>% 
    rename(long = long, lat = lat, PROVINCE = province, REALM = realm, group = id)
  
  # Shift 0-360
  if(min(df$long, na.rm=TRUE) < 0) {
    df <- df %>% mutate(long = ifelse(long < 0, long + 360, long))
  }
  
  # Calculate Centroids & Labels
  centroids <- df %>% 
    group_by(PROVINCE, REALM) %>% 
    summarise(X_cent = mean(long), Y_cent = mean(lat), .groups="drop") %>% 
    arrange(X_cent) %>% 
    mutate(
      Map_Index = row_number(),
      Clean_Name = str_trunc(PROVINCE, 25, "right"),
      # Tag if Strict or New
      Is_Strict = PROVINCE %in% strict_provs
    )
  
  # Join back
  df_final <- df %>% 
    left_join(centroids %>% dplyr::select(PROVINCE, Map_Index, Is_Strict), by="PROVINCE")
  
  return(list(polys = df_final, labs = centroids))
}

# --- 2. PLOTTING FUNCTION ---
create_slide <- function(data_list, title, mode = "FULL") {
  # mode = "GREY_EXPANSION" or "FULL"
  
  polys <- data_list$polys
  labs  <- data_list$labs
  
  # --- COLOR LOGIC ---
  # Define base colors
  realm_cols <- c(
    "Central Indo-Pacific" = "#E6AB02", 
    "Eastern Indo-Pacific" = "#A6761D",
    "Western Indo-Pacific" = "#66A61E", 
    "Temperate Northern Pacific" = "#7570B3",
    "Temperate Australasia" = "#E7298A", 
    "Temperate Southern Africa" = "#D95F02",
    "EXCLUDED" = "grey85" # The Ghost Color
  )
  
  # Modify categories based on Mode
  if(mode == "GREY_EXPANSION") {
    # If it's NOT strict, call it "EXCLUDED" to make it grey
    polys <- polys %>% mutate(fill_cat = ifelse(Is_Strict, REALM, "EXCLUDED"))
    labs  <- labs  %>% mutate(fill_cat = ifelse(Is_Strict, REALM, "EXCLUDED"))
  } else {
    # Full Color
    polys$fill_cat <- polys$REALM
    labs$fill_cat  <- labs$REALM
  }
  
  # --- MAP ---
  p_map <- ggplot() +
    borders("world2", colour="gray90", fill="gray95", size=0.1) +
    
    geom_polygon(data = polys, aes(x=long, y=lat, group=group, fill=fill_cat),
                 color = "white", size = 0.1, alpha = 0.9) +
    
    geom_text_repel(data = labs, aes(x=X_cent, y=Y_cent, label=Map_Index),
                    size=3, fontface="bold", min.segment.length=0, 
                    box.padding=0.4, segment.size=0.3, segment.color="gray50",
                    bg.color="white", bg.r=0.15) +
    
    scale_fill_manual(values = realm_cols) +
    
    # *** ZOOMED OUT FRAMING ***
    # 15 to 260 covers South Africa to Easter Island buffer
    # -50 to 50 covers NZ to Japan
    coord_fixed(xlim = c(15, 260), ylim = c(-50, 50), expand = FALSE) +
    
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      panel.background = element_rect(fill = "#E0F3F8", color = NA),
      axis.text = element_blank(), 
      panel.grid = element_blank()
    )
  
  # --- KEY ---
  n_rows <- ceiling(nrow(labs) / 4)
  labs$col_group <- rep(1:4, each = n_rows, length.out = nrow(labs))
  
  p_key <- ggplot(labs, aes(x = 0, y = -Map_Index, label = paste0(Map_Index, ". ", Clean_Name))) +
    geom_text(aes(color = fill_cat, alpha = ifelse(fill_cat=="EXCLUDED", 0.5, 1)), 
              hjust = 0, size = 3, fontface = "bold") +
    
    scale_color_manual(values = realm_cols) +
    scale_alpha_identity() +
    
    facet_wrap(~col_group, scales = "free") +
    theme_void() +
    theme(legend.position = "none", strip.text = element_blank()) +
    scale_x_continuous(expand = c(0, 0.1))
  
  return(p_map / p_key + plot_layout(heights = c(4, 1)))
}

# --- 3. GENERATE ---
cat("Preparing Data...\n")
data <- prepare_presentation_data()

cat("Generating Slide 1: Baseline (New Regions Greyed)...\n")
s1 <- create_slide(data, "Study Extent: Strict Replication (Baseline)", mode = "GREY_EXPANSION")
ggsave(file.path(OUTPUT_DIR, "Slide_1_Baseline_Greyed.png"), s1, width = 14, height = 10, dpi = 300)

cat("Generating Slide 2: Full Expansion (All Colored)...\n")
s2 <- create_slide(data, "Study Extent: Proposed Expansion", mode = "FULL")
ggsave(file.path(OUTPUT_DIR, "Slide_2_Expansion_Full.png"), s2, width = 14, height = 10, dpi = 300)

print(s1)
print(s2)
cat("DONE.\n")
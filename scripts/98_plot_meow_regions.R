# scripts/98_plot_meow_regions.R
# ------------------------------------------------------------------------------
# FIGURE S1: GLOBAL DEBUG (No Cropping, No Zooming)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(sf, ggplot2, dplyr, readr, ggrepel, viridis, patchwork, gridExtra, stringr)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
SHP_DIR  <- file.path(DATA_DIR, "shapefiles")
MEOW_SHP <- file.path(SHP_DIR, "meow_ecos.shp")
PAPER_GUY_CSV <- file.path(DATA_DIR, "paper_guy_meow_ecos_df.csv")

dir.create(file.path(BASE_DIR, "outputs", "figures"), recursive=TRUE, showWarnings=FALSE)

cat("Loading Data...\n")
meow_raw <- st_read(MEOW_SHP, quiet = TRUE)
meow_sf <- st_shift_longitude(meow_raw) # 0-360 Shift

# --- 1. DEFINITIONS ---
paper_codes <- c(9, 18:41, 55, 58) 
thesis_codes <- c(paper_codes, 51, 53, 57, 56, 54)

# --- 2. PREPARE SF DATA ---
prepare_data <- function(sf_data, codes) {
  subset <- sf_data %>% 
    filter(PROV_CODE %in% codes) %>% 
    st_make_valid()
  
  centroids <- subset %>% 
    st_centroid() %>% st_coordinates() %>% as.data.frame() %>% 
    rename(X_cent = X, Y_cent = Y)
  
  subset %>% 
    bind_cols(centroids) %>% 
    arrange(X_cent) %>% 
    mutate(
      Map_Index = row_number(),
      Clean_Name = str_trunc(PROVINCE, 25, "right"), 
      Label_Text = paste0(Map_Index, ". ", Clean_Name)
    ) %>%
    mutate(REALM = factor(REALM))
}

data_rep <- prepare_data(meow_sf, paper_codes)
data_exp <- prepare_data(meow_sf, thesis_codes)

# --- 3. GLOBAL MAP FUNCTION ---
plot_global <- function(data, title) {
  
  cols <- c("Central Indo-Pacific"="#E6AB02", "Eastern Indo-Pacific"="#A6761D",
            "Western Indo-Pacific"="#66A61E", "Temperate Northern Pacific"="#7570B3",
            "Temperate Australasia"="#E7298A", "Temperate Southern Africa"="#D95F02")
  
  labs <- data %>% st_drop_geometry() %>% dplyr::select(X_cent, Y_cent, Map_Index) %>% distinct()
  
  ggplot() +
    # "world2" is explicitly designed for 0-360 Pacific Centered maps
    borders("world2", colour="gray90", fill="gray95", size=0.1) +
    
    geom_sf(data=data, aes(fill=REALM), color="white", size=0.1, alpha=0.9) +
    
    geom_text_repel(data=labs, aes(x=X_cent, y=Y_cent, label=Map_Index),
                    size=3, fontface="bold", min.segment.length=0, 
                    box.padding=0.5, max.overlaps = Inf) +
    
    scale_fill_manual(values=cols) +
    
    # !!! GLOBAL VIEW: No cropping, just let it show the whole 0-360 world !!!
    coord_sf(xlim = c(0, 360), ylim = c(-90, 90), expand = FALSE) +
    
    labs(title=title, x=NULL, y=NULL) +
    theme_minimal() +
    theme(legend.position="none", plot.title=element_text(face="bold", hjust=0.5))
}

plot_key <- function(data) {
  key <- data %>% st_drop_geometry() %>% 
    dplyr::select(Map_Index, Clean_Name, REALM) %>% 
    distinct(Map_Index, .keep_all=TRUE) %>% arrange(Map_Index)
  
  n_rows <- ceiling(nrow(key)/4)
  key$col <- rep(1:4, each=n_rows, length.out=nrow(key))
  
  ggplot(key, aes(x=0, y=-Map_Index, label=paste0(Map_Index, ". ", Clean_Name))) +
    geom_text(aes(color=REALM), hjust=0, size=2.8, fontface="bold") +
    scale_color_manual(values=c("Central Indo-Pacific"="#E6AB02", "Eastern Indo-Pacific"="#A6761D",
                                "Western Indo-Pacific"="#66A61E", "Temperate Northern Pacific"="#7570B3",
                                "Temperate Australasia"="#E7298A", "Temperate Southern Africa"="#D95F02")) +
    facet_wrap(~col, scales="free") +
    theme_void() + theme(legend.position="none", strip.text=element_blank())
}

# --- 4. EXECUTE ---

cat("Generating Global Replication Plot...\n")
p1 <- plot_global(data_rep, "Global View: Replication (Strict 27)")
k1 <- plot_key(data_rep)
f1 <- p1 / k1 + plot_layout(heights=c(3, 1))
ggsave(file.path(BASE_DIR, "outputs", "figures", "FigS1_Global_Replication.png"), f1, width=14, height=10)

cat("Generating Global Expansion Plot...\n")
p2 <- plot_global(data_exp, "Global View: Thesis Expansion")
k2 <- plot_key(data_exp)
f2 <- p2 / k2 + plot_layout(heights=c(3, 1))
ggsave(file.path(BASE_DIR, "outputs", "figures", "FigS1_Global_Expansion.png"), f2, width=14, height=10)

print(f1)
print(f2)
cat("DONE.\n")
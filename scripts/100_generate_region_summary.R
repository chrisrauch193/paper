# scripts/100_generate_region_summary.R
# ------------------------------------------------------------------------------
# FINAL STEP: GENERATE CLEAN LLM PROMPT (Robust Version)
# ------------------------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, stringr)

BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# --- 1. LOAD DATA ---
# We compare YOUR current configuration vs. PAPER GUY
FILE_USER <- file.path(DATA_DIR, "meow_ecos_df.csv")
FILE_PAPER <- file.path(DATA_DIR, "paper_guy_meow_ecos_df.csv")

get_hierarchy <- function(filepath) {
  if(!file.exists(filepath)) return(NULL)
  read_csv(filepath, show_col_types = FALSE) %>% 
    dplyr::select(REALM, PROVINCE) %>% 
    distinct() %>% 
    arrange(REALM, PROVINCE)
}

user_df  <- get_hierarchy(FILE_USER)
paper_df <- get_hierarchy(FILE_PAPER)

if(is.null(user_df) || is.null(paper_df)) stop("Missing CSV files. Run pipelines 0c first.")

# --- 2. ANALYZE DIFFERENCES ---
# What extra regions do you have that he doesn't?
added_provs <- setdiff(user_df$PROVINCE, paper_df$PROVINCE)
missing_provs <- setdiff(paper_df$PROVINCE, user_df$PROVINCE)

# --- 3. BUILD TEXT CONTENT ---
# We build a list of strings, then write them all at once.
txt <- c()
add_line <- function(line) { txt <<- c(txt, line) }

add_line("I am conducting a biogeographical study on Clownfish and Sea Anemones.")
add_line("I need to validate my study extent ('Expanded') against a baseline 'Strict Replication'.")
add_line("")
add_line("### DATASET 1: STRICT BASELINE (Paper Guy)")
add_line(paste("Total Provinces:", nrow(paper_df)))
add_line("Contains strictly Tropical Indo-Pacific regions.")
add_line("Regions:")

current_realm <- ""
for(i in 1:nrow(paper_df)) {
  r <- paper_df$REALM[i]
  p <- paper_df$PROVINCE[i]
  if(r != current_realm) {
    add_line(paste0("\n**Realm: ", r, "**"))
    current_realm <- r
  }
  add_line(paste0("- ", p))
}

add_line("\n--------------------------------------------------------------------\n")

add_line("### DATASET 2: MY PROPOSED EXPANSION")
add_line(paste("Total Provinces:", nrow(user_df)))
add_line("This dataset includes the Baseline plus the following NEW additions:")
add_line("")

if(length(added_provs) > 0) {
  added_info <- user_df %>% filter(PROVINCE %in% added_provs)
  for(i in 1:nrow(added_info)) {
    add_line(paste0("- ", added_info$PROVINCE[i], " (Realm: ", added_info$REALM[i], ")"))
  }
} else {
  add_line("[No new provinces added - This dataset is identical to Baseline]")
}

if(length(missing_provs) > 0) {
  add_line("")
  add_line("**NOTE: I have EXCLUDED these Baseline provinces:**")
  add_line(paste("-", missing_provs, collapse="\n"))
}

add_line("\n--------------------------------------------------------------------\n")

add_line("### RESEARCH QUESTION FOR LLM:")
add_line("Given the goal of modelling the realized niche and future distribution of Amphiprioninae (Clownfish)")
add_line("and their host Sea Anemones under climate change:")
add_line("1. Scientifically, do the 'NEW additions' contain viable populations of clownfish or host anemones?")
add_line("2. Specifically evaluate 'Agulhas', 'Northern New Zealand', and 'Southeast Australian Shelf'.")
add_line("   - Do these represent valid range shifts?")
add_line("   - Or are they false positives that should be removed?")
add_line("3. Are there any other critical transition zones I am missing?")

# --- 4. WRITE TO FILE ---
out_file <- file.path(OUTPUT_DIR, "prompt_for_deep_research.txt")
writeLines(txt, out_file)

cat("----------------------------------------------------------------\n")
cat("SUCCESS! Clean prompt saved to:\n")
cat(out_file, "\n")
cat("----------------------------------------------------------------\n")
cat("You can now open that file and copy-paste it into the AI.\n")
# Locate shared paths for Rscript or source().
.egfr_script <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.egfr_script)) {
  .egfr_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(.egfr_arg)) stop("Run with Rscript or source this file.")
  .egfr_script <- sub("^--file=", "", .egfr_arg[[1]])
}
.egfr_src <- dirname(normalizePath(.egfr_script, mustWork = TRUE))
while (!file.exists(file.path(.egfr_src, "common", "02_project_paths.R"))) {
  .egfr_parent <- dirname(.egfr_src)
  if (.egfr_parent == .egfr_src) stop("Cannot locate src/common/02_project_paths.R")
  .egfr_src <- .egfr_parent
}
source(file.path(.egfr_src, "common", "02_project_paths.R"))
use_working_directory(".")

library(UpSetR)
library(dplyr)
library(tidyr)
library(stringr)

# --- 1. DATA PREP ---
data <- read.csv(project_path("outputs/standard_master_table/EGFR_Arm2_Master_Standardized.csv"))

tools_cols <- c("AlphaMissense_Standard", "VARITY_Standard", "SuSPect_Standard", "MetaLR_Standard", "MetaSVM_Standard", "MetaSNP_Standard", "PHDSNPg_Standard")
tools_clean <- c("AlphaMissense", "VARITY", "SuSPect", "MetaLR", "MetaSVM", "MetaSNP", "PHDSNPg")

df_bin <- data %>% select(all_of(tools_cols))
colnames(df_bin) <- tools_clean

df_bin <- df_bin %>%
  mutate(across(everything(), ~ ifelse(.x == "Deleterious", 1, 0))) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  mutate(across(everything(), as.integer)) %>%
  as.data.frame()

# --- 2. PLOTTING ---
output_file <- project_path("outputs/figures/EGFR_Arm2_Disease_Tendency_UpSet.png")

# THE FIX: Increased height from 1600 to 2400 to stop vertical squishing
png(filename = output_file, width = 3600, height = 2400, res = 300)

suppressWarnings({
  upset(df_bin, 
        sets = rev(tools_clean),
        nintersects = 20, 
        order.by = "freq",
        
        main.bar.color = "#E41A1C", 
        sets.bar.color = "#4DAF4A", 
        
        matrix.color   = "#000000", 
        shade.color    = "#E5E5E5", 
        shade.alpha    = 0.6,       
        
        mb.ratio = c(0.55, 0.45), 
        
        point.size = 3.5, 
        line.size = 1.0, 
        
        show.numbers = "yes",
        number.angles = 0,
        text.scale = 2.5,
  )
})

dev.off()
cat("Success! Taller canvas applied to separate chart and matrix. Saved to", output_file, "\n")
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
use_working_directory("outputs/figures")

# --- 0. DEPENDENCIES ---
# Set default CRAN mirror so the script doesn't pause to ask you for one
options(repos = c(CRAN = "https://cloud.r-project.org"))

# List of required packages
required_packages <- c("ggplot2", "cowplot", "magick")

# Check which ones are missing and install them
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages)
}

cat("All dependencies loaded successfully!\n")

# --- 1. SET DIRECTORIES & AUTO-CROP FUNCTION ---

library(ggplot2)
library(cowplot)
library(magick)

# Corrected for Docker internal volume path
base_dir <- project_path("outputs/figures/docking")

# Load your 3 Erlotinib images
cat("Loading Erlotinib images...\n")
img_A <- image_trim(image_read(file.path(base_dir, "WT_TKD_Erlotinib.png")))
img_B <- image_trim(image_read(file.path(base_dir, "V774M_Erlotinib.png")))
img_C <- image_trim(image_read(file.path(base_dir, "L861Q_Erlotinib.png")))

# The Sniper Grid Generator
generate_grid_overlay <- function(img) {
  ggdraw() + 
    draw_image(img) + 
    # Draw horizontal and vertical grid lines every 0.1 units
    geom_hline(yintercept = seq(0, 1, by = 0.1), color = "red", alpha = 0.3) +
    geom_vline(xintercept = seq(0, 1, by = 0.1), color = "red", alpha = 0.3) +
    # Add the text coordinates so you can read them
    geom_text(data = expand.grid(x = seq(0.1, 0.9, 0.1), y = seq(0.1, 0.9, 0.1)),
              aes(x = x, y = y, label = paste0(x, ",", y)), 
              size = 4, color = "blue", fontface = "bold") +
    theme(plot.margin = margin(20, 20, 20, 20))
}

# Generate and save a calibrated version of each image
cat("Generating sniper grids...\n")
ggsave(file.path(base_dir, "CALIBRATE_WT_Erlotinib.png"), plot = generate_grid_overlay(img_A), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_V774M_Erlotinib.png"), plot = generate_grid_overlay(img_B), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_L861Q_Erlotinib.png"), plot = generate_grid_overlay(img_C), width = 8, height = 6)

cat("Calibration images saved successfully! Open them on your Windows host to map your exact arrow coordinates.\n")
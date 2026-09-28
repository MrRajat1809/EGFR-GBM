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

library(ggplot2)
library(cowplot)
library(magick)

base_dir <- project_path("outputs/figures/docking")

# Load your images
img_A <- image_trim(image_read(file.path(base_dir, "WT_Cetuximab.png")))
img_B <- image_trim(image_read(file.path(base_dir, "C240Y_Cetuximab.png")))
img_C <- image_trim(image_read(file.path(base_dir, "G598V_Cetuximab.png")))
img_D <- image_trim(image_read(file.path(base_dir, "WT_Nimotuzumab.png")))
img_E <- image_trim(image_read(file.path(base_dir, "C240Y_Nimotuzumab.png")))
img_F <- image_trim(image_read(file.path(base_dir, "G598V_Nimotuzumab.png")))

# The Sniper Grid Generator (Fixed for grob conversion)
generate_grid_overlay <- function(img) {
  ggdraw() + 
    draw_image(img) + # This is the fix!
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
ggsave(file.path(base_dir, "CALIBRATE_A.png"), plot = generate_grid_overlay(img_A), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_B.png"), plot = generate_grid_overlay(img_B), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_C.png"), plot = generate_grid_overlay(img_C), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_D.png"), plot = generate_grid_overlay(img_D), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_E.png"), plot = generate_grid_overlay(img_E), width = 8, height = 6)
ggsave(file.path(base_dir, "CALIBRATE_F.png"), plot = generate_grid_overlay(img_F), width = 8, height = 6)

cat("Calibration images saved! Open them to read your exact arrow coordinates.\n")
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
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!require("png", quietly = TRUE)) install.packages("png") 

library(ggplot2)
library(cowplot)
library(png)
library(grid) # Base R package, requires no installation

# --- 1. SET DIRECTORIES AND LOAD THE IMAGES ---
cat("Loading high-res figures natively...\n")

base_dir <- project_path("outputs/figures")

# Helper function to read a PNG and convert it to a Grid Object (Grob)
# This completely bypasses the need for the 'magick' package!
get_png_grob <- function(file_name) {
  img <- readPNG(file.path(base_dir, file_name))
  rasterGrob(img, interpolate = TRUE)
}

imgA <- get_png_grob("Fig6A.png")
imgB <- get_png_grob("Fig6B.png")
imgC <- get_png_grob("Fig6C.png")
imgD <- get_png_grob("Fig6D.png")
imgE <- get_png_grob("Fig6E.png")
imgF <- get_png_grob("Fig6F.png")

# --- 2. BUILD THE NESTED GRID ---
cat("Combining panels and applying labels...\n")

# Step 2a: Stack C and D vertically to form the left column of the bottom section
col_CD <- plot_grid(
  imgC, imgD, 
  ncol = 1, 
  labels = c("C", "D"), 
  label_size = 24, 
  label_fontface = "bold"
)

# Step 2b: Create the bottom row combining the C/D column, E, and F
bottom_row <- plot_grid(
  col_CD, imgE, imgF, 
  ncol = 3, 
  labels = c("", "E", "F"), 
  label_size = 24, 
  label_fontface = "bold",
  rel_widths = c(1, 1, 1) 
)

# Step 2c: Combine A, B, and the bottom row into the final master panel
master_panel <- plot_grid(
  imgA, imgB, bottom_row, 
  ncol = 1, 
  labels = c("A", "B", ""), 
  label_size = 24, 
  label_fontface = "bold",
  rel_heights = c(1, 1, 1.2) 
)

# --- 3. EXPORT THE MASTER FIGURE ---
output_file <- file.path(base_dir, "Fig6_Master_Panel_Final.png")

# Saving as a massive, publication-ready 300 DPI composite
ggsave(output_file, plot = master_panel, width = 20, height = 24, dpi = 300, bg = "white")

cat("Success! Master panel saved to", output_file, "\n")

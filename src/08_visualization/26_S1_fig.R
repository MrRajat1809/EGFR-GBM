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

Sys.setenv(MAGICK_MEMORY_LIMIT = "4GB")
Sys.setenv(MAGICK_MAP_LIMIT = "4GB")
Sys.setenv(MAGICK_DISK_LIMIT = "8GB")
Sys.setenv(MAGICK_THREAD_LIMIT = "1")

library(magick)
library(rsvg)
library(grid)

base_dir <- project_path("outputs/figures/supplementary")

# --- 1. SET PARAMETERS & LOAD IMAGES ---
cat("Loading native images...\n")

# A4 Width at 300 DPI is ~2480 pixels. 
# With 150px margins on both left and right (300px total), 
# the internal image width should be exactly 2180 to fit perfectly.
target_width <- 2180 

# Load images and force a pure white background
svg_img <- image_read_svg(file.path(base_dir, "NetSurfPro_wild_type.svg"), width = target_width)
svg_img <- image_background(svg_img, "white")

png1 <- image_read(file.path(base_dir, "SOPMA_wild_type.png"))
png1 <- image_background(png1, "white")

png2 <- image_read(file.path(base_dir, "PSIPRED_wild_type.png"))
png2 <- image_background(png2, "white")

# --- 2. SCALE & TRIM ---
cat("Scaling and cleaning edges...\n")

png1_scaled <- image_scale(png1, as.character(target_width))
png2_scaled <- image_scale(png2, as.character(target_width))

# Custom function to shave vertical pixels using image_crop
shave_vertical <- function(img, pixels) {
  info <- image_info(img)
  new_height <- info$height - (pixels * 2)
  geom <- paste0(info$width, "x", new_height, "+0+", pixels)
  return(image_crop(img, geom))
}

# Shave 2 pixels off top and bottom to destroy bounding boxes
png1_scaled <- shave_vertical(png1_scaled, 2)
png2_scaled <- shave_vertical(png2_scaled, 2)
svg_img     <- shave_vertical(svg_img, 2)

# --- 3. COMBINE & ADD MARGINS ---
cat("Combining figures and adding page margins...\n")

# Stack the images flush against each other
final_img <- image_append(c(png1_scaled, png2_scaled, svg_img), stack = TRUE)

# Record the exact heights of the sub-images
h1 <- image_info(png1_scaled)$height
h2 <- image_info(png2_scaled)$height

# Add the 150px white canvas margin
margin_px <- 150
final_img <- image_border(final_img, "white", paste0(margin_px, "x", margin_px))

# --- 4. ANNOTATE IN THE MARGINS ---
cat("Placing labels cleanly in the white space...\n")

text_color <- "black"
# Slightly increased font size since the overall image is wider now
label_size <- 55 

final_img <- image_annotate(final_img, "A", size = label_size, color = text_color, gravity = "northwest", location = paste0("+50+", margin_px + 20))
final_img <- image_annotate(final_img, "B", size = label_size, color = text_color, gravity = "northwest", location = paste0("+50+", margin_px + h1 + 20))
final_img <- image_annotate(final_img, "C", size = label_size, color = text_color, gravity = "northwest", location = paste0("+50+", margin_px + h1 + h2 + 20))

# --- 5. EXPORT WITH EXACT WIDTH & DYNAMIC LENGTH ---
cat("Writing PDF with A4 width and dynamic length to prevent squishing...\n")

output_path <- file.path(base_dir, "S1_fig.pdf")

info <- image_info(final_img)

# A4 paper width is exactly 8.27 inches
a4_width_inches <- 8.27 

# Calculate the exact proportional height based on the final stacked image
aspect_ratio <- info$height / info$width
dynamic_height_inches <- a4_width_inches * aspect_ratio

# Open PDF with locked A4 width and free-flowing height
pdf(file = output_path, title = "S1_fig", width = a4_width_inches, height = dynamic_height_inches)

# Draw the image filling 100% of the newly calculated canvas
grid.raster(as.raster(final_img), width = unit(1, "npc"), height = unit(1, "npc"))

dev.off()

cat("Done! A4 width achieved, dynamic length preserved. Saved at:", output_path, "\n")
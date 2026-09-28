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
library(magick)
library(ggplot2)
library(grid)

cat("Generating strict coordinate legend...\n")

# --- 1. GENERATE BULLETPROOF LEGEND ---
legend_plot <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  
  # Left Column: Confidence Lines
  annotate("segment", x = 0.05, xend = 0.15, y = 0.70, yend = 0.70, color = "grey40", linewidth = 1.5) +
  annotate("text", x = 0.17, y = 0.70, label = "High Confidence (score \u2265 0.900)", hjust = 0, size = 5, fontface = "bold") +
  
  annotate("segment", x = 0.05, xend = 0.15, y = 0.30, yend = 0.30, color = "grey60", linewidth = 1, linetype = "dashed") +
  annotate("text", x = 0.17, y = 0.30, label = "Medium Confidence", hjust = 0, size = 5) +
  
  # Right Column: Biological Clusters
  annotate("point", x = 0.65, y = 0.85, color = "firebrick", size = 6) + 
  annotate("text", x = 0.68, y = 0.85, label = "Core Receptor Complex", hjust = 0, size = 5) +
  
  annotate("point", x = 0.65, y = 0.50, color = "royalblue", size = 6) + 
  annotate("text", x = 0.68, y = 0.50, label = "Downstream Modulators", hjust = 0, size = 5) +
  
  annotate("point", x = 0.65, y = 0.15, color = "forestgreen", size = 6) + 
  annotate("text", x = 0.68, y = 0.15, label = "Regulatory Elements", hjust = 0, size = 5)

ggsave("temp_legend.png", plot = legend_plot, width = 8, height = 2, dpi = 300, bg = "white")

# --- 2. LOAD & CLEAN IMAGES VIA MAGICK ---
cat("Loading and cleaning images via Magick...\n")

img_A <- image_read("GeneMania.png") |> image_background("white")
img_B_net <- image_read("STRING.png") |> image_background("white")
img_B_leg <- image_read("temp_legend.png") |> image_background("white")

shave_edges <- function(img, pixels) {
  info <- image_info(img)
  geom <- paste0(info$width - (pixels*2), "x", info$height - (pixels*2), "+", pixels, "+", pixels)
  return(image_crop(img, geom))
}

img_A <- shave_edges(img_A, 2)
img_B_net <- shave_edges(img_B_net, 2)

# --- 3. SCALE & COMBINE ---
cat("Stacking panels...\n")

target_col_width <- 1800 

img_A_scaled <- image_scale(img_A, as.character(target_col_width))
img_B_net_scaled <- image_scale(img_B_net, as.character(target_col_width))
img_B_leg_scaled <- image_scale(img_B_leg, as.character(target_col_width))

panel_B <- image_append(c(img_B_net_scaled, img_B_leg_scaled), stack = TRUE)

max_height <- max(image_info(img_A_scaled)$height, image_info(panel_B)$height)
img_A_final <- image_extent(img_A_scaled, paste0(target_col_width, "x", max_height), gravity = "center", color = "white")
panel_B_final <- image_extent(panel_B, paste0(target_col_width, "x", max_height), gravity = "center", color = "white")

final_img <- image_append(c(img_A_final, panel_B_final), stack = FALSE)

# --- 4. TIGHT MARGINS & ANNOTATIONS ---
# Dropped from the massive 150px to a tight 75px
margin_px <- 75
final_img <- image_border(final_img, "white", paste0(margin_px, "x", margin_px))

label_size <- 65
final_img <- image_annotate(final_img, "A", size = label_size, color = "black", weight = 700, 
                            gravity = "northwest", location = paste0("+", margin_px - 20, "+", margin_px - 20))
final_img <- image_annotate(final_img, "B", size = label_size, color = "black", weight = 700, 
                            gravity = "northwest", location = paste0("+", margin_px + target_col_width, "+", margin_px - 20))

# --- 5. EXPORT DYNAMICALLY WRAPPED PDF ---
cat("Calculating tight PDF dimensions...\n")

info <- image_info(final_img)

# We fix the width to a standard high-res size (12 inches) 
# and let the exact pixel aspect ratio dictate the height.
target_width_inches <- 12
aspect_ratio <- info$height / info$width
dynamic_height_inches <- target_width_inches * aspect_ratio

cat("Exporting...\n")

# paper = "special" ensures the canvas matches our dynamic math perfectly
pdf(file = "S2_fig.pdf", title = "S2_fig", width = target_width_inches, height = dynamic_height_inches, paper = "special", bg = "white")
grid.raster(as.raster(final_img), width = unit(1, "npc"), height = unit(1, "npc"))
dev.off()

# Clean up temp file
unlink("temp_legend.png")

cat("Done! Image exported with a tight, dynamic canvas and small margins.\n")
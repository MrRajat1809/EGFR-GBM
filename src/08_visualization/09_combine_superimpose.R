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
if (!require("magick", quietly = TRUE)) install.packages("magick")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(ggplot2)
library(cowplot)
library(magick)

# --- 1. SET DIRECTORIES, LOAD & CROP ---
cat("Loading and cropping 8 high-res ChimeraX panels...\n")

# Docking panels in the shared output directory.
base_dir <- project_path("outputs/figures/docking")

load_and_crop <- function(filename) {
  img_path <- file.path(base_dir, filename)
  if (!file.exists(img_path)) stop(paste("File not found:", img_path))
  img <- image_read(img_path)
  img_cropped <- image_trim(img) 
  return(ggdraw() + draw_image(img_cropped))
}

# Load Macro images
m1 <- load_and_crop("Panel1_ECD_vs_C240Y.png")
m2 <- load_and_crop("Panel2_ECD_vs_G598V.png")
m3 <- load_and_crop("Panel3_TKD_vs_V774M.png")
m4 <- load_and_crop("Panel4_TKD_vs_L861Q.png")

# Load Micro (Zoom) images
z1 <- load_and_crop("Panel1_ECD_vs_C240Y_zoom.png")
z2 <- load_and_crop("Panel2_ECD_vs_G598V_zoom.png")
z3 <- load_and_crop("Panel3_TKD_vs_V774M_zoom.png")
z4 <- load_and_crop("Panel4_TKD_vs_L861Q_zoom.png")

# --- 2. ADD CURVED ARROWS & TEXT TO ZOOM IMAGES ---
cat("Annotating zoom images with text and curved arrows...\n")

# Helper function to dynamically draw the text and curved arrow over the image
annotate_zoom <- function(zoom_img, top_text, bottom_text) {
  ggdraw(zoom_img) +
    # Top text (Amino acid + position)
    draw_text(top_text, x = 0.5, y = 0.95, size = 25, vjust = 1, fontface = "plain") +
    # Bottom text (Mutated amino acid)
    draw_text(bottom_text, x = 0.75, y = 0.05, size = 25, vjust = 0, fontface = "plain") +
    # Curved arrow swooping from top-left to bottom-right
    geom_curve(aes(x = 0.35, y = 0.85, xend = 0.65, yend = 0.12),
               curvature = 0.4, linewidth = 1.2, color = "black",
               arrow = arrow(length = unit(0.04, "npc"), type = "closed"))
}

# Apply annotations (Coordinates inside the function might need slight tweaking depending on your exact protein crop)
z1_anno <- annotate_zoom(z1, "Cysteine (240th position)", "Tyrosine")
z2_anno <- annotate_zoom(z2, "Glycine (598th position)", "Valine")
z3_anno <- annotate_zoom(z3, "Valine (774th position)", "Methionine")
z4_anno <- annotate_zoom(z4, "Leucine (861st position)", "Glutamine")

# --- 3. FUSE MACRO & ANNOTATED MICRO INTO SINGLE PANELS ---
cat("Fusing Macro and Micro pairs...\n")

padding <- theme(plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt"))

panel_A <- plot_grid(m1, z1_anno, ncol = 2, rel_widths = c(5, 1)) + padding
panel_B <- plot_grid(m2, z2_anno, ncol = 2, rel_widths = c(5, 1)) + padding
panel_C <- plot_grid(m3, z3_anno, ncol = 2, rel_widths = c(5, 1)) + padding
panel_D <- plot_grid(m4, z4_anno, ncol = 2, rel_widths = c(5, 1)) + padding

# --- 4. BUILD THE FINAL MASTER GRID & ADD LEGEND ---
cat("Building final A, B, C, D grid with Legend...\n")

master_panel <- plot_grid(
  panel_A, panel_B, 
  panel_C, panel_D,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"), 
  label_size = 32, label_fontface = "bold", align = "vh"
)

# Dummy ggplot for the legend
legend_data <- data.frame(
  Type = factor(
    c("Wild type", "C240Y mutant", "G598V mutant", "V774M mutant", "L861Q mutant"),
    levels = c("Wild type", "C240Y mutant", "G598V mutant", "V774M mutant", "L861Q mutant")
  ), x = 1, y = 1
)

legend_colors <- c("Wild type" = "#4DAF4A", 
                   "C240Y mutant" = "salmon", 
                   "G598V mutant" = "orange", 
                   "V774M mutant" = "red", 
                   "L861Q mutant" = "blue")

dummy_plot <- ggplot(legend_data, aes(x = x, y = y, fill = Type)) +
  geom_point(shape = 22, size = 12, color = "black") +
  scale_fill_manual(values = legend_colors, name = NULL) +
  theme_minimal() +
  theme(
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), panel.grid = element_blank(),
    legend.position = "left",                     
    legend.text = element_text(size = 25),             # Font size set to 25
    legend.key.size = unit(1.0, "cm"),
    legend.spacing.y = unit(1.0, "cm")               # Added spacing between entries
  ) +
  guides(fill = guide_legend(byrow = TRUE))           # Required to make spacing.y work

extracted_legend <- get_legend(dummy_plot)

# Append the legend FIRST (left side), then the image grid
final_figure <- plot_grid(
  extracted_legend, master_panel, 
  ncol = 2, 
  rel_widths = c(1, 4.5) # Legend gets 1 part, images get 4.5 parts
)

# --- 5. EXPORT FOR HALF-PAGE JOURNAL SIZE ---
output_file <- file.path(base_dir, "Figure_1_Final_HalfPage.png")

ggsave(output_file, plot = final_figure, width = 36, height = 18, dpi = 300, bg = "white")

cat("Success! Polished master panel with arrows and left legend saved to", output_file, "\n")

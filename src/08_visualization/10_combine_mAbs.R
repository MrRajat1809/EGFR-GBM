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

# --- 1. SET DIRECTORIES & AUTO-CROP FUNCTION ---
cat("Loading and cropping 6 high-res ChimeraX docking panels...\n")
base_dir <- project_path("outputs/figures/docking")

load_and_crop <- function(filename) {
  img_path <- file.path(base_dir, filename)
  if (!file.exists(img_path)) stop(paste("File not found:", img_path))
  
  img <- image_read(img_path)
  img_cropped <- image_trim(img) 
  return(ggdraw() + draw_image(img_cropped))
}

img_A <- load_and_crop("WT_Cetuximab.png")
img_B <- load_and_crop("C240Y_Cetuximab.png")
img_C <- load_and_crop("G598V_Cetuximab.png")
img_D <- load_and_crop("WT_Nimotuzumab.png")
img_E <- load_and_crop("C240Y_Nimotuzumab.png")
img_F <- load_and_crop("G598V_Nimotuzumab.png")

# --- 2. ADVANCED ANNOTATION (PER-PANEL CONTROL) ---
cat("Annotating panels with custom sweeping arrows...\n")

# Updated function: includes independent curvature control, shorter arrowheads, and panel margins
annotate_custom <- function(img, left_text, right_text, l_xend, l_yend, r_xend, r_yend, l_curve = -0.2, r_curve = 0.2) {
  ggdraw(img) +
    # Left text and arrow (Start points shifted slightly left to avoid overlap)
    draw_text(left_text, x = 0.15, y = 0.95, size = 25, fontface = "plain", color = "black") +
    geom_curve(aes(x = 0.15, y = 0.90, xend = l_xend, yend = l_yend),
               curvature = l_curve, linewidth = 1.0, color = "black",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    # Right text and arrow (Start points shifted slightly right to avoid overlap)
    draw_text(right_text, x = 0.85, y = 0.05, size = 25, fontface = "plain", color = "black") +
    geom_curve(aes(x = 0.85, y = 0.10, xend = r_xend, yend = r_yend),
               curvature = r_curve, linewidth = 1.0, color = "black",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    # Consistent whitespace between panels
    theme(plot.margin = margin(20, 20, 20, 20)) 
}

# Fine-tuned end points and curvatures specifically for each image's geometry
panel_E <- annotate_custom(img_A, "Cetuximab", "WT",    l_xend=0.20, l_yend=0.60, r_xend=0.80, r_yend=0.20)
panel_F <- annotate_custom(img_B, "Cetuximab", "C240Y", l_xend=0.20, l_yend=0.60, r_xend=0.80, r_yend=0.20)
panel_G <- annotate_custom(img_C, "Cetuximab", "G598V", l_xend=0.20, l_yend=0.60, r_xend=0.80, r_yend=0.20)

panel_H <- annotate_custom(img_D, "Nimotuzumab", "WT",    l_xend=0.20, l_yend=0.60, r_xend=0.80, r_yend=0.20)
panel_I <- annotate_custom(img_E, "Nimotuzumab", "C240Y", l_xend=0.20, l_yend=0.60, r_xend=0.80, r_yend=0.20)
panel_J <- annotate_custom(img_F, "Nimotuzumab", "G598V", l_xend=0.20, l_yend=0.40, r_xend=0.80, r_yend=0.20)

# --- 3. BUILD THE FINAL 2x3 MASTER GRID ---
cat("Building final A-F grid...\n")

master_panel <- plot_grid(
  panel_E, panel_F, panel_G,
  panel_H, panel_I, panel_J,
  ncol = 3, nrow = 2,
  labels = c("E", "F", "G", "H", "I", "J"), 
  label_size = 25, label_fontface = "bold", align = "vh"
)

# --- 4. BUILD THE COLOR LEGEND ---
cat("Generating legend...\n")

legend_data <- data.frame(
  Type = factor(
    c("Cetuximab", "Nimotuzumab", "Wild-type (WT)", "C240Y Mutant", "G598V Mutant"),
    levels = c("Cetuximab", "Nimotuzumab", "Wild-type (WT)", "C240Y Mutant", "G598V Mutant")
  ), x = 1, y = 1
)

legend_colors <- c("Cetuximab" = "tan", 
                   "Nimotuzumab" = "sienna", 
                   "Wild-type (WT)" = "darkseagreen", 
                   "C240Y Mutant" = "cadetblue", 
                   "G598V Mutant" = "slategray")

dummy_plot <- ggplot(legend_data, aes(x = x, y = y, fill = Type)) +
  geom_point(shape = 22, size = 16, color = "black") + 
  scale_fill_manual(values = legend_colors, name = "Protein Chains") +
  theme_minimal() +
  theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 25, face = "bold"), 
    legend.text = element_text(size = 25),                 
    legend.key.size = unit(1.8, "cm"),                     
    legend.spacing.y = unit(0.8, "cm"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 100)     
  ) +
  guides(fill = guide_legend(byrow = TRUE))

extracted_legend <- get_legend(dummy_plot)

# --- 5. BUILD THE HORIZONTAL GRADIENT LEGEND ---
dummy_gradient <- ggplot(data.frame(x = 1, y = 1, val = c(-1, 1)), aes(x = x, y = y, fill = val)) +
  geom_point() +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    breaks = c(-1, 1),
    labels = c("Donor", "Acceptor"),
    name = "H-bonds",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(6, "cm"),  # Controls the length of the horizontal line
      barheight = unit(0.8, "cm"), # Controls the thickness of the line
      ticks = FALSE
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25),
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0), # Adds space between the top legend and this one
    plot.margin = margin(t = 0, r = 0, b = 0, l = 100)
  )

extracted_gradient <- get_legend(dummy_gradient)

# Stack the two legends vertically
combined_legends <- plot_grid(extracted_legend, extracted_gradient, ncol = 1, rel_heights = c(3, 1))

# Combine the image grid with the stacked legends
final_figure <- plot_grid(
  master_panel, combined_legends, 
  ncol = 2, 
  rel_widths = c(7, 1.5) 
)

# --- 6. EXPORT FINAL FIGURE ---
output_file <- file.path(base_dir, "Docking_Master_Grid_Perfected.png")

# Maintained the large canvas size
ggsave(output_file, plot = final_figure, width = 34, height = 18, dpi = 300, bg = "white")

cat("Success! Perfected master panel with dual legends saved to", output_file, "\n")
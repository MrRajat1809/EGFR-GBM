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

required_packages <- c("ggplot2", "cowplot", "magick")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  cat("Installing missing packages:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages)
}

library(ggplot2)
library(cowplot)
library(magick)

# --- 1. SET DIRECTORIES & AUTO-CROP FUNCTION ---
cat("Loading and cropping 3 high-res Erlotinib docking panels...\n")
base_dir <- project_path("outputs/figures/docking")

load_and_crop <- function(filename) {
  img_path <- file.path(base_dir, filename)
  if (!file.exists(img_path)) stop(paste("File not found:", img_path))
  
  img <- image_read(img_path)
  img_cropped <- image_trim(img) 
  return(ggdraw() + draw_image(img_cropped))
}

img_A <- load_and_crop("WT_TKD_Erlotinib.png")
img_B <- load_and_crop("V774M_Erlotinib.png")
img_C <- load_and_crop("L861Q_Erlotinib.png")

# --- 2. ADVANCED ANNOTATION (PER-PANEL CONTROL) ---
cat("Annotating panels with custom sweeping arrows...\n")

annotate_custom <- function(img, left_text, right_text, l_xend, l_yend, r_xend, r_yend, l_curve = -0.2, r_curve = 0.2) {
  ggdraw(img) +
    # Left text and arrow (Drug Name)
    draw_text(left_text, x = 0.15, y = 0.95, size = 25, fontface = "plain", color = "black") +
    geom_curve(aes(x = 0.15, y = 0.90, xend = l_xend, yend = l_yend),
               curvature = l_curve, linewidth = 1.0, color = "black",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    # Right text and arrow (Protein State)
    draw_text(right_text, x = 0.85, y = 0.05, size = 25, fontface = "plain", color = "black") +
    geom_curve(aes(x = 0.85, y = 0.10, xend = r_xend, yend = r_yend),
               curvature = r_curve, linewidth = 1.0, color = "black",
               arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    # Consistent whitespace between panels
    theme(plot.margin = margin(20, 20, 20, 20)) 
}

# ==============================================================================
# INSERT YOUR CALIBRATED COORDINATES HERE
# Replace the 0.5 placeholders with the (x,y) numbers from your sniper grids!
# ==============================================================================
panel_A <- annotate_custom(img_A, "Erlotinib", "WT",    
                           l_xend=0.3, l_yend=0.6, r_xend=0.6, r_yend=0.3)

panel_B <- annotate_custom(img_B, "Erlotinib", "V774M", 
                           l_xend=0.3, l_yend=0.6, r_xend=0.6, r_yend=0.3)

panel_C <- annotate_custom(img_C, "Erlotinib", "L861Q", 
                           l_xend=0.3, l_yend=0.6, r_xend=0.6, r_yend=0.3)
# ==============================================================================

# --- 3. BUILD THE FINAL 1x3 MASTER GRID ---
cat("Building final A-C grid...\n")

master_panel <- plot_grid(
  panel_A, panel_B, panel_C,
  ncol = 3, nrow = 1,
  labels = c("K", "L", "M"), 
  label_size = 25, label_fontface = "bold", align = "vh"
)

# --- 4. BUILD THE COLOR LEGEND ---
cat("Generating legends...\n")

legend_data <- data.frame(
  Type = factor(
    c("Erlotinib", "Wild-type (WT)", "V774M Mutant", "L861Q Mutant"),
    levels = c("Erlotinib", "Wild-type (WT)", "V774M Mutant", "L861Q Mutant")
  ), x = 1, y = 1
)

legend_colors <- c("Erlotinib" = "skyblue", 
                   "Wild-type (WT)" = "lightgreen", 
                   "V774M Mutant" = "thistle", 
                   "L861Q Mutant" = "wheat")

dummy_plot <- ggplot(legend_data, aes(x = x, y = y, fill = Type)) +
  geom_point(shape = 22, size = 16, color = "black") + 
  scale_fill_manual(values = legend_colors, name = "Protein & Ligand") +
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
      barwidth = unit(6, "cm"),  
      barheight = unit(0.8, "cm"), 
      ticks = FALSE
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25),
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0), 
    plot.margin = margin(t = 0, r = 0, b = 0, l = 100)
  )

extracted_gradient <- get_legend(dummy_gradient)

# Stack the two legends vertically
combined_legends <- plot_grid(extracted_legend, extracted_gradient, ncol = 1, rel_heights = c(3, 1))

# Combine the image grid with the stacked legends
# Uses 6:2 ratio to account for the wider 1x3 aspect ratio
final_figure <- plot_grid(
  master_panel, combined_legends, 
  ncol = 2, 
  rel_widths = c(6.5, 2) 
)

# --- 6. EXPORT FINAL FIGURE ---
output_file <- file.path(base_dir, "TKD_Erlotinib_Master_Grid.png")

# Width 30 for 3 panels side-by-side, height 12 for a balanced single row
ggsave(output_file, plot = final_figure, width = 30, height = 12, dpi = 300, bg = "white")

cat("Success! Perfected TKD master panel with dual legends saved to", output_file, "\n")
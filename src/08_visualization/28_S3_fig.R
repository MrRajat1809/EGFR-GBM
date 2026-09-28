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
# install.packages(c("magick", "ggplot2", "cowplot"))

library(magick)
library(ggplot2)
library(cowplot)

# --- 1. PREPARE THE 3D IMAGE PANEL ---
cat("Processing 3D structure image...\n")

# Load image, crop empty transparent space, and give it a MUTED BLACK background
img_3d <- image_read("EGFR_All_Clusters.png") |>
  image_trim() |>
  image_background("grey15")

# Convert to a ggplot object
p_image <- ggdraw() + 
  draw_image(img_3d, scale = 0.95) + 
  theme(
    plot.background = element_rect(fill = "grey15", color = NA),
    panel.border = element_rect(color = "black", linewidth = 2, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# --- 2. PREPARE DATA FOR 1D SEQUENCE PLOT ---
cat("Generating 1D sequence plot...\n")

mutations <- data.frame(
  res = c(
    # Cluster 1
    62, 63, 108, 222, 223, 229, 240, 252, 254, 256, 263, 289, 304,
    # Cluster 2
    270, 596, 598, 608, 620, 624, 628, 636,
    # Cluster 3
    774, 861, 1068
  ),
  cluster = factor(c(
    rep("Cluster 1", 13),
    rep("Cluster 2", 8),
    rep("Cluster 3", 3)
  ), levels = c("Cluster 1", "Cluster 2", "Cluster 3"))
)

cluster_colors <- c(
  "Cluster 1" = "firebrick",
  "Cluster 2" = "cornflowerblue", 
  "Cluster 3" = "sienna"
)

# --- 3. BUILD 1D SEQUENCE PLOT ---
p_seq <- ggplot() +
  # Draw the main EGFR backbone as a LIGHT GREEN rectangle
  annotate("rect", xmin = 1, xmax = 1210, ymin = -0.5, ymax = 0.5, 
           fill = "lightgreen", color = "black", linewidth = 0.5) +
  
  # Add the mutations as brightly colored vertical bars across the backbone
  geom_segment(data = mutations, 
               aes(x = res, xend = res, y = -0.55, yend = 0.55, color = cluster), 
               linewidth = 1.5) +
  
  # Lock the X-axis to the exact length of the protein
  scale_x_continuous(breaks = seq(0, 1210, by = 200), limits = c(-10, 1220), expand = c(0,0)) +
  scale_y_continuous(limits = c(-2, 2)) + 
  
  # Apply colors and style the legend
  scale_color_manual(name = "Mutation Clusters", values = cluster_colors) +
  
  # Clean up the theme
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm")
  ) +
  labs(x = "EGFR Amino Acid Position")

# --- 4. ASSEMBLE WITH 4-SIDED MARGINS ---
cat("Assembling final figure...\n")

# Stack the plots first
stacked_plots <- plot_grid(
  p_image, 
  p_seq, 
  ncol = 1, 
  rel_heights = c(1, 0.25)
)

# Wrap the entire stacked grid in a new drawing layer to force an outer margin
final_page <- ggdraw(stacked_plots) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(40, 40, 40, 40)
  )

# --- 5. EXPORT WITH METADATA ---
cat("Exporting to PDF with internal metadata...\n")

output_file <- "S3_fig.pdf"

# Using the base pdf() function allows us to inject the "title" directly into the file's metadata
pdf(file = output_file, title = "S3_fig", width = 10, height = 10, bg = "white")
print(final_page) # final_page is a ggplot object, so we print it to the PDF device
dev.off()

cat("Done! Image successfully generated as", output_file, "\n")
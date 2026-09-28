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

library(magick)
library(ggplot2)
library(cowplot)

cat("Loading and trimming images...\n")
# 1. Read and instantly shave off invisible white borders
fig1 <- image_trim(image_read("docking/Fig1.png"))
fig2 <- image_trim(image_read("docking/Fig2.png"))
fig3 <- image_trim(image_read("docking/Fig3.png"))

cat("Calculating exact dynamic aspect ratios...\n")
# 2. Get pixel dimensions of the trimmed images
info1 <- image_info(fig1)
info2 <- image_info(fig2)
info3 <- image_info(fig3)

# Calculate Height-to-Width ratio for each panel
ar1 <- info1$height / info1$width
ar2 <- info2$height / info2$width
ar3 <- info3$height / info3$width

cat("Converting to plot objects...\n")
p1 <- ggdraw() + draw_image(fig1)
p2 <- ggdraw() + draw_image(fig2)
p3 <- ggdraw() + draw_image(fig3)

cat("Stacking with normalized relative heights...\n")
# 3. rel_heights ensures the 1-row image isn't given the same vertical space as a 2-row image.
# This completely eliminates the gaps between the rows.
combined_plot <- plot_grid(
  p1, p2, p3, 
  ncol = 1, 
  rel_heights = c(ar1, ar2, ar3)
)

cat("Exporting to disk with exact canvas sizing...\n")
# 4. We set a fixed width (8 inches for a standard journal page)
base_width <- 8.0 

# We mathematically calculate the exact height needed so absolutely no extra white space is generated
total_height <- base_width * (ar1 + ar2 + ar3)

# Save it. 
ggsave("docking/Combined_Figures.png", plot = combined_plot, 
       width = base_width, height = total_height, dpi = 300, bg = "white", limitsize = FALSE)

cat("Success! Images stacked without distortion or gaps.\n")
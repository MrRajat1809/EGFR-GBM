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
use_working_directory(".")

# --- 0. DEPENDENCIES ---
# Force secure HTTPS mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!require("magick", quietly = TRUE)) install.packages("magick")

library(ggplot2)
library(cowplot)
library(magick)

# --- 1. LOAD THE IMAGES ---
cat("Loading high-res figures...\n")
img1 <- ggdraw() + draw_image(project_path("outputs/figures/EGFR_Arm1_Pathogenicity_UpSet.png"))
img2 <- ggdraw() + draw_image(project_path("outputs/figures/EGFR_Arm1_Pathogenicity_Rangeplot.png"))
img3 <- ggdraw() + draw_image(project_path("outputs/figures/EGFR_Arm2_Disease_Tendency_UpSet.png"))
img4 <- ggdraw() + draw_image(project_path("outputs/figures/EGFR_Arm2_Disease_Tendency_Rangeplot.png"))

# --- 2. BUILD THE GRID ---
cat("Combining panels and applying labels...\n")
master_panel <- plot_grid(
  img1, img2, img3, img4, 
  labels = c("A", "B", "C", "D"), 
  label_size = 32,                    
  label_fontface = "bold",
  ncol = 2,                           
  align = 'hv'                        
)

# --- 3. EXPORT THE MASTER FIGURE ---
output_file <- project_path("outputs/figures/EGFR_Figure1_Master_Panel.png")

ggsave(output_file, plot = master_panel, width = 20, height = 16, dpi = 300, bg = "white")

cat("Success! Master panel saved with equal areas to", output_file, "\n")
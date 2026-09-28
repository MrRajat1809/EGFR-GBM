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

# =============================================================================
# Figure 5 Master Assembly Script
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!require("png", quietly = TRUE)) install.packages("png")

library(ggplot2)
library(cowplot)
library(png)
library(grid)

cat("Starting final Figure 5 assembly...\n")

base_dir <- getwd()

# =============================================================================
# 1. IMAGE TRIMMING FUNCTIONS
# =============================================================================

# Reduced pad from 20 to 5 to pull margins closer without totally removing them
trim_gentle <- function(img, tol = 0.995, pad = 5) {
  non_white_rows <- apply(img, 1, function(r) any(r < tol))
  non_white_cols <- apply(img, 2, function(c) any(c < tol))

  r <- which(non_white_rows)
  c <- which(non_white_cols)

  r1 <- max(min(r) - pad, 1)
  r2 <- min(max(r) + pad, dim(img)[1])
  c1 <- max(min(c) - pad, 1)
  c2 <- min(max(c) + pad, dim(img)[2])

  img[r1:r2, c1:c2, , drop = FALSE]
}

# Tight trimming for D panel remains the same
trim_tight <- function(img, tol = 0.995, pad = 5) {
  non_white_rows <- apply(img, 1, function(r) any(r < tol))
  non_white_cols <- apply(img, 2, function(c) any(c < tol))

  r <- which(non_white_rows)
  c <- which(non_white_cols)

  r1 <- max(min(r) - pad, 1)
  r2 <- min(max(r) + pad, dim(img)[1])
  c1 <- max(min(c) - pad, 1)
  c2 <- min(max(c) + pad, dim(img)[2])

  img[r1:r2, c1:c2, , drop = FALSE]
}

# =============================================================================
# 2. IMAGE LOADING FUNCTION
# =============================================================================

get_png_grob <- function(file_name, trim_fun) {
  path <- file.path(base_dir, file_name)
  if (!file.exists(path)) stop(paste("Missing file:", path))
  
  cat("Loading:", file_name, "\n")
  img <- readPNG(path)
  img_trim <- trim_fun(img)
  rasterGrob(img_trim, interpolate = TRUE)
}

# =============================================================================
# 3. LOAD PANELS
# =============================================================================

imgA <- get_png_grob("Fig5A.png", trim_gentle)
imgB <- get_png_grob("Fig5B.png", trim_gentle)
imgC <- get_png_grob("Fig5C.png", trim_gentle)
imgD <- get_png_grob("Fig5D.png", trim_tight)

cat("All panels loaded successfully.\n")

# =============================================================================
# 4. BUILD PANEL LAYOUT
# =============================================================================

# (panel_BC remains exactly the same as before)
panel_BC <- plot_grid(
  imgB,
  imgC,
  ncol = 2,
  align = "h",
  labels = c("B", "C"),
  label_size = 26,
  label_fontface = "bold",
  rel_widths = c(1.05, 1)
)

master_panel <- plot_grid(
  imgA,
  panel_BC,
  imgD,
  ncol = 1,
  labels = c("A", "", "D"),
  label_size = 26,
  label_fontface = "bold",
  
  # RESTORED PANEL A HEIGHT
  # Panel A gets '1.15' again so the heatmap doesn't squish.
  # Panel D stays at a lower value ('0.5') to prevent bottom whitespace.
  rel_heights = c(1.15, 1.8, 0.5) 
)

cat("Layout finalized.\n")


# =============================================================================
# 5. EXPORT FINAL FIGURE
# =============================================================================

output_file <- file.path(base_dir, "Fig5_Master_Panel_FINAL_300dpi.png")

ggsave(
  filename = output_file,
  plot = master_panel,
  width = 12,
  
  # Bumped height slightly to 11.5 to accommodate Panel A's restored size
  # (Still much better than the original 14, which caused the huge gaps)
  height = 11.5, 
  
  dpi = 300,
  bg = "white"
)

cat("\nSUCCESS\n")
cat("Figure saved to:\n")
cat(output_file, "\n")

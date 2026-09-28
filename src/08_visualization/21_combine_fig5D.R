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

# --- DEPENDENCIES -------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("png", quietly = TRUE)) install.packages("png")
if (!require("abind", quietly = TRUE)) install.packages("abind")

library(png)
library(abind)
library(grid)

cat("Merging D1 and D2 with controlled micro-spacing...\n")

base_dir <- getwd()

# --- READ IMAGES --------------------------------------------------------------

img1 <- readPNG("Fig5D1.png")   # chromosome
img2 <- readPNG("Fig5D2.png")   # enhancer

# --- TRIM TRUE EMPTY ROWS ONLY ------------------------------------------------

trim_rows <- function(img, tol = 0.999) {

  keep <- apply(img, 1, function(r)
    any(r < tol)
  )

  img[which(keep), , , drop = FALSE]
}

img1 <- trim_rows(img1)
img2 <- trim_rows(img2)

# --- MATCH WIDTHS -------------------------------------------------------------

w1 <- dim(img1)[2]
w2 <- dim(img2)[2]

max_w <- max(w1, w2)

pad_width_center <- function(img, target_w) {

  h <- dim(img)[1]
  w <- dim(img)[2]
  ch <- dim(img)[3]

  if (w == target_w) return(img)

  total_pad <- target_w - w

  left_pad <- floor(total_pad / 2)
  right_pad <- ceiling(total_pad / 2)

  pad_left <- array(
    1,
    dim = c(h, left_pad, ch)
  )

  pad_right <- array(
    1,
    dim = c(h, right_pad, ch)
  )

  abind(
    pad_left,
    img,
    pad_right,
    along = 2
  )
}

img1 <- pad_width_center(img1, max_w)
img2 <- pad_width_center(img2, max_w)

# --- ADD CONTROLLED MICRO GAP -------------------------------------------------

# 1 point ≈ 1/72 inch
# At 300 DPI → ~4 pixels
# So 3 points ≈ 12 pixels

gap_pixels <- 10   # tiny visual breathing space

gap <- array(
  1,
  dim = c(gap_pixels, max_w, dim(img1)[3])
)

# --- STACK (STRICT ORDER) -----------------------------------------------------

merged_img <- abind(
  img1,
  gap,
  img2,
  along = 1
)

# --- EXPORT -------------------------------------------------------------------

output_file <- file.path(base_dir, "Fig5D.png")

png(
  filename = output_file,
  width = 12,
  height = 3,
  units = "in",
  res = 300
)

grid.newpage()
grid.raster(merged_img)

dev.off()

cat("SUCCESS\n")
cat("Created:", output_file, "\n")
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

# Set the working directory to the exact subfolder with the underscore
setwd(project_path("outputs/figures/discovery_studio"))

cat("Loading and processing 9 docking images...\n")

# --- 1. IMAGE PROCESSING FUNCTION ---
prep_panel <- function(file_name, panel_title) {
  # Check if file exists to prevent errors
  if(!file.exists(file_name)) {
     stop(paste("Cannot find file:", file_name, "- Check the name and directory!"))
  }
    
  img <- image_read(file_name) |>
    image_trim() |>
    image_background("white")
  
  p <- ggdraw() + 
    draw_image(img, scale = 0.90) + 
    # Changed fontface from "bold" to "plain"
    draw_label(panel_title, x = 0.5, y = 0.98, vjust = 1, fontface = "plain", size = 16) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      # Removed the border around the images
      panel.border = element_blank(), 
      plot.margin = margin(5, 5, 5, 5)
    )
  return(p)
}

# --- 2. LOAD ALL 9 PANELS (Using exact _Best.png names) ---
# ROW 1: Cetuximab
p_wt_cetux    <- prep_panel("WT_Cetuximab_Best.png", "WT - Cetuximab Complex")
p_c240y_cetux <- prep_panel("C240Y_Cetuximab_Best.png", "C240Y - Cetuximab Complex")
p_g598v_cetux <- prep_panel("G598V_Cetuximab_Best.png", "G598V - Cetuximab Complex")

# ROW 2: Nimotuzumab
p_wt_nimo    <- prep_panel("WT_Nimotuzumab_Best.png", "WT - Nimotuzumab Complex")
p_c240y_nimo <- prep_panel("C240Y_Nimotuzumab_Best.png", "C240Y - Nimotuzumab Complex")
p_g598v_nimo <- prep_panel("G598V_Nimotuzumab_Best.png", "G598V - Nimotuzumab Complex")

# ROW 3: Erlotinib (Note the WT_TKD name)
p_wt_erlo    <- prep_panel("WT_TKD_Erlotinib_Best.png", "WT - Erlotinib Complex")
p_v774m_erlo <- prep_panel("V774M_Erlotinib_Best.png", "V774M - Erlotinib Complex")
p_l861q_erlo <- prep_panel("L861Q_Erlotinib_Best.png", "L861Q - Erlotinib Complex")


# --- 3. ASSEMBLE THE 3x3 GRID ---
cat("Weaving panels into a 3x3 matrix...\n")

master_grid <- plot_grid(
  p_wt_cetux, p_c240y_cetux, p_g598v_cetux,
  p_wt_nimo, p_c240y_nimo, p_g598v_nimo,
  p_wt_erlo, p_v774m_erlo, p_l861q_erlo,
  ncol = 3,
  nrow = 3,
  labels = "AUTO", 
  label_size = 24,
  align = "hv"
)

# Add outer margin
final_figure <- ggdraw(master_grid) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(30, 30, 30, 30)
  )

# --- 4. EXPORT AS PUBLICATION-GRADE PDF ---
cat("Exporting final figure...\n")

# Save it back up one level in the main figures folder
output_file <- project_path("outputs/figures/S4_fig.pdf")

# 15x15 inches keeps the 3x3 squares perfectly proportioned
pdf(file = output_file, title = "S4_fig", width = 15, height = 15, bg = "white")
print(final_figure)
dev.off()

cat("Success! Master figure saved as", output_file, "\n")
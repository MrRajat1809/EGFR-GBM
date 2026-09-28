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

library(ggplot2)
library(dplyr)

read_xvg <- function(filepath, system_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  colnames(df) <- c("Time", "RMSD")
  df$System <- system_name
  return(df)
}

# 1. Load Data
wt_rmsd <- read_xvg("../gromacs/WT_Erlotinib/rmsd_protein_200.xvg", "Wild Type")
v774m_rmsd <- read_xvg("../gromacs/V774M_Erlotinib/rmsd_protein_200.xvg", "V774M Mutant")
l861q_rmsd <- read_xvg("../gromacs/L861Q_Erlotinib/rmsd_protein_200.xvg", "L861Q Mutant")

master_rmsd <- bind_rows(wt_rmsd, v774m_rmsd, l861q_rmsd)

# 2. Set Factor Levels (V774M=Red, L861Q=Orange, WT=Green)
master_rmsd$System <- factor(master_rmsd$System, 
                             levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

# 3. Generate Plot
rmsd_plot <- ggplot(master_rmsd, aes(x = Time, y = RMSD, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("V774M Mutant" = "#E41A1C", 
                                "L861Q Mutant" = "#FF7F00", 
                                "Wild Type" = "#4DAF4A")) +
  # Expand = c(0,0) removes the "floating" gap at the origin
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45)) +
  theme_light(base_size = 14) +
  theme(
    legend.position = "none", # <-- Legend completely removed
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.3),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    panel.border = element_rect(color = "#d1d1d1", fill = NA, linewidth = 0.5)
  ) +
  labs(title = "Root Mean Square Deviation (RMSD)", 
       x = "Time (ns)", 
       y = "RMSD (nm)")

# 4. Save as high-resolution PNG
ggsave("RMSD_HighRes.png", plot = rmsd_plot, width = 8, height = 6, units = "in", dpi = 300)
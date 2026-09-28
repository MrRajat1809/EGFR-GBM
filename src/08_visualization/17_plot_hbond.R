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

# Custom function to safely read GROMACS .xvg files
read_xvg_hbond <- function(filepath, system_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  # H-bond outputs have Time and Number of H-bonds in the first two columns
  df <- df[, 1:2]
  colnames(df) <- c("Time", "Hbonds")
  df$System <- system_name
  return(df)
}

# 1. Load Data
wt_hbond <- read_xvg_hbond("../gromacs/WT_Erlotinib/hbnum_200.xvg", "Wild Type")
v774m_hbond <- read_xvg_hbond("../gromacs/V774M_Erlotinib/hbnum_200.xvg", "V774M Mutant")
l861q_hbond <- read_xvg_hbond("../gromacs/L861Q_Erlotinib/hbnum_200.xvg", "L861Q Mutant")

master_hbond <- bind_rows(wt_hbond, v774m_hbond, l861q_hbond)

# 2. Set Factor Levels (WT first, then Mutants, to match standard publication formatting)
master_hbond$System <- factor(master_hbond$System, 
                              levels = c("Wild Type", "V774M Mutant", "L861Q Mutant"))

# 3. Generate Plot
hbond_plot <- ggplot(master_hbond, aes(x = Time, y = Hbonds, color = System)) +
  # Using geom_step because H-bonds are discrete integers (0, 1, 2, etc.)
  geom_step(alpha = 0.8, linewidth = 0.5) +
  
  scale_color_manual(values = c("Wild Type" = "#4DAF4A", 
                                "V774M Mutant" = "#E41A1C", 
                                "L861Q Mutant" = "#FF7F00")) +
  
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201)) + 
  # Setting Y-axis to dynamically fit the max number of bonds, usually 3 or 4
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(master_hbond$Hbonds) + 1)) +
  
  theme_light(base_size = 14) +
  theme(
    legend.position = "none", # Legend removed since facets are labeled
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.3),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    panel.border = element_rect(color = "#d1d1d1", fill = NA, linewidth = 0.5),
    # Styling the facet labels
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(face = "bold", color = "black", size = 12)
  ) +
  
  # Stack the plots vertically (3 rows, 1 column) so they don't overlap
  facet_wrap(~System, ncol = 1) +
  
  labs(title = "Protein-Ligand Hydrogen Bond Dynamics", 
       x = "Time (ns)", 
       y = "Number of H-Bonds")

# 4. Save as high-resolution PNG (Height increased to 8 to accommodate the 3 stacked panels)
ggsave("HBOND_HighRes.png", plot = hbond_plot, width = 8, height = 8, units = "in", dpi = 300)
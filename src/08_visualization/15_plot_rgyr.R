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

read_xvg_rgyr <- function(filepath, system_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  df <- df[, 1:2]
  colnames(df) <- c("Time", "Rg")
  
  # Convert ps to ns so it matches the 0-200 X-axis perfectly
  df$Time <- df$Time / 1000 
  
  df$System <- system_name
  return(df)
}

wt_rgyr <- read_xvg_rgyr("../gromacs/WT_Erlotinib/rgyr_200.xvg", "Wild Type")
v774m_rgyr <- read_xvg_rgyr("../gromacs/V774M_Erlotinib/rgyr_200.xvg", "V774M Mutant")
l861q_rgyr <- read_xvg_rgyr("../gromacs/L861Q_Erlotinib/rgyr_200.xvg", "L861Q Mutant")

master_rgyr <- bind_rows(wt_rgyr, v774m_rgyr, l861q_rgyr)
master_rgyr$System <- factor(master_rgyr$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

rgyr_plot <- ggplot(master_rgyr, aes(x = Time, y = Rg, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("V774M Mutant" = "#E41A1C", "L861Q Mutant" = "#FF7F00", "Wild Type" = "#4DAF4A")) +
  # Removed limits to prevent dropping data, kept expand=c(0,0) for the clean corners
  scale_x_continuous(expand = c(0, 0)) +
  theme_light(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.3),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    panel.border = element_rect(color = "#d1d1d1", fill = NA, linewidth = 0.5)
  ) +
  labs(title = "Radius of Gyration (RGYR)", x = "Time (ns)", y = "RGYR (nm)")

ggsave("RGYR_HighRes.png", plot = rgyr_plot, width = 8, height = 6, units = "in", dpi = 300)
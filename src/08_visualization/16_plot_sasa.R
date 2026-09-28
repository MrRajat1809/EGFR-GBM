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

read_xvg_sasa <- function(filepath, system_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  df <- df[, 1:2]
  colnames(df) <- c("Time", "SASA")
  
  # Convert ps to ns
  df$Time <- df$Time / 1000 
  
  df$System <- system_name
  return(df)
}

wt_sasa <- read_xvg_sasa("../gromacs/WT_Erlotinib/sasa_200.xvg", "Wild Type")
v774m_sasa <- read_xvg_sasa("../gromacs/V774M_Erlotinib/sasa_200.xvg", "V774M Mutant")
l861q_sasa <- read_xvg_sasa("../gromacs/L861Q_Erlotinib/sasa_200.xvg", "L861Q Mutant")

master_sasa <- bind_rows(wt_sasa, v774m_sasa, l861q_sasa)
master_sasa$System <- factor(master_sasa$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

sasa_plot <- ggplot(master_sasa, aes(x = Time, y = SASA, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("V774M Mutant" = "#E41A1C", "L861Q Mutant" = "#FF7F00", "Wild Type" = "#4DAF4A")) +
  # Removed limits to prevent dropping data
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
  labs(title = "Solvent Accessible Surface Area (SASA)", x = "Time (ns)", y = expression(bold(paste("Area (", nm^2, ")"))))

ggsave("SASA_HighRes.png", plot = sasa_plot, width = 8, height = 6, units = "in", dpi = 300)
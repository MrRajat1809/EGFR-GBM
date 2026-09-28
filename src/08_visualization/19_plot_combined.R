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

# Load essential libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# 1. DATA READING FUNCTIONS
# ==========================================

read_xvg_standard <- function(filepath, system_name, col_names) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)[, 1:2]
  colnames(df) <- col_names
  df$System <- system_name
  return(df)
}

read_xvg_timeconv <- function(filepath, system_name, col_names) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)[, 1:2]
  colnames(df) <- col_names
  df[,1] <- df[,1] / 1000 # Convert ps to ns
  df$System <- system_name
  return(df)
}

# ==========================================
# 2. LOAD ALL DATA
# ==========================================

# A. RMSD
wt_rmsd <- read_xvg_standard("../gromacs/WT_Erlotinib/rmsd_protein_200.xvg", "Wild Type", c("Time", "RMSD"))
v774m_rmsd <- read_xvg_standard("../gromacs/V774M_Erlotinib/rmsd_protein_200.xvg", "V774M Mutant", c("Time", "RMSD"))
l861q_rmsd <- read_xvg_standard("../gromacs/L861Q_Erlotinib/rmsd_protein_200.xvg", "L861Q Mutant", c("Time", "RMSD"))
master_rmsd <- bind_rows(wt_rmsd, v774m_rmsd, l861q_rmsd)
master_rmsd$System <- factor(master_rmsd$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

# B. RMSF
wt_rmsf <- read_xvg_standard("../gromacs/WT_Erlotinib/rmsf_200.xvg", "Wild Type", c("Residue", "RMSF"))
v774m_rmsf <- read_xvg_standard("../gromacs/V774M_Erlotinib/rmsf_200.xvg", "V774M Mutant", c("Residue", "RMSF"))
l861q_rmsf <- read_xvg_standard("../gromacs/L861Q_Erlotinib/rmsf_200.xvg", "L861Q Mutant", c("Residue", "RMSF"))
master_rmsf <- bind_rows(wt_rmsf, v774m_rmsf, l861q_rmsf)
master_rmsf$System <- factor(master_rmsf$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

# C. RGYR
wt_rgyr <- read_xvg_timeconv("../gromacs/WT_Erlotinib/rgyr_200.xvg", "Wild Type", c("Time", "Rg"))
v774m_rgyr <- read_xvg_timeconv("../gromacs/V774M_Erlotinib/rgyr_200.xvg", "V774M Mutant", c("Time", "Rg"))
l861q_rgyr <- read_xvg_timeconv("../gromacs/L861Q_Erlotinib/rgyr_200.xvg", "L861Q Mutant", c("Time", "Rg"))
master_rgyr <- bind_rows(wt_rgyr, v774m_rgyr, l861q_rgyr)
master_rgyr$System <- factor(master_rgyr$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

# D. SASA
wt_sasa <- read_xvg_timeconv("../gromacs/WT_Erlotinib/sasa_200.xvg", "Wild Type", c("Time", "SASA"))
v774m_sasa <- read_xvg_timeconv("../gromacs/V774M_Erlotinib/sasa_200.xvg", "V774M Mutant", c("Time", "SASA"))
l861q_sasa <- read_xvg_timeconv("../gromacs/L861Q_Erlotinib/sasa_200.xvg", "L861Q Mutant", c("Time", "SASA"))
master_sasa <- bind_rows(wt_sasa, v774m_sasa, l861q_sasa)
master_sasa$System <- factor(master_sasa$System, levels = c("V774M Mutant", "L861Q Mutant", "Wild Type"))

# ==========================================
# 3. DEFINE COMMON AESTHETICS
# ==========================================
custom_colors <- c("V774M Mutant" = "#E41A1C", "L861Q Mutant" = "#FF7F00", "Wild Type" = "#4DAF4A")

# Base size set to 15, plain text for axis labels
shared_theme <- theme_light(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.3),
    axis.title = element_text(face = "plain", color = "black"), # Only Plain Text
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0), # Keep Bold
    panel.border = element_rect(color = "#d1d1d1", fill = NA, linewidth = 0.5)
  )

# ==========================================
# 4. BUILD INDIVIDUAL PLOTS
# ==========================================

plot_a <- ggplot(master_rmsd, aes(x = Time, y = RMSD, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.45)) +
  shared_theme +
  labs(title = "A)   Root Mean Square Deviation (RMSD)", x = "Time (ns)", y = "RMSD (nm)")

plot_b <- ggplot(master_rmsf, aes(x = Residue, y = RMSF, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  shared_theme +
  labs(title = "B)   Root Mean Square Fluctuation (RMSF)", x = "Residue", y = "RMSF (nm)")

plot_c <- ggplot(master_rgyr, aes(x = Time, y = Rg, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201)) +
  shared_theme +
  labs(title = "C)   Radius of Gyration (RGYR)", x = "Time (ns)", y = "RGYR (nm)")

# Removed the bold() styling from the expression
plot_d <- ggplot(master_sasa, aes(x = Time, y = SASA, color = System)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 201)) +
  shared_theme +
  labs(title = "D)   Solvent Accessible Surface Area (SASA)", x = "Time (ns)", y = expression(paste("Area (", nm^2, ")")))

# ==========================================
# 5. ASSEMBLE AND SAVE WITH PATCHWORK
# ==========================================

combined_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.text = element_text(face = "plain", size = 15), # Plain text
    legend.key.width = unit(1, "cm") # Reduced the breadth of the color lines
  )

ggsave("Combined_MD_Analysis_HighRes.png", plot = combined_plot, width = 16, height = 10, units = "in", dpi = 300)
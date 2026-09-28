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

# Load the libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# 1. Universal function to read GROMACS .xvg files
read_metric <- function(filepath, system_name, metric_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  df <- df[, 1:2]
  colnames(df) <- c("Index", "Value")
  df$System <- system_name
  df$Metric <- metric_name
  return(df)
}

# 2. Re-order to match your reference image (Mutants first, WT last)
systems <- c("V774M", "L861Q", "Wild Type") 
my_colors <- c("V774M" = "#e41a1c", "L861Q" = "#ff7f00", "Wild Type" = "#4daf4a")

# 3. Load ALL the data and apply the new ordering
wt_rmsd <- read_metric("../gromacs/WT_Erlotinib/rmsd_protein_200.xvg", "Wild Type", "RMSD")
v7_rmsd <- read_metric("../gromacs/V774M_Erlotinib/rmsd_protein_200.xvg", "V774M", "RMSD")
l8_rmsd <- read_metric("../gromacs/L861Q_Erlotinib/rmsd_protein_200.xvg", "L861Q", "RMSD")
df_rmsd <- bind_rows(wt_rmsd, v7_rmsd, l8_rmsd)
df_rmsd$System <- factor(df_rmsd$System, levels = systems)

wt_rmsf <- read_metric("../gromacs/WT_Erlotinib/rmsf_200.xvg", "Wild Type", "RMSF")
v7_rmsf <- read_metric("../gromacs/V774M_Erlotinib/rmsf_200.xvg", "V774M", "RMSF")
l8_rmsf <- read_metric("../gromacs/L861Q_Erlotinib/rmsf_200.xvg", "L861Q", "RMSF")
df_rmsf <- bind_rows(wt_rmsf, v7_rmsf, l8_rmsf)
df_rmsf$System <- factor(df_rmsf$System, levels = systems)

wt_rgyr <- read_metric("../gromacs/WT_Erlotinib/rgyr_200.xvg", "Wild Type", "RGYR")
v7_rgyr <- read_metric("../gromacs/V774M_Erlotinib/rgyr_200.xvg", "V774M", "RGYR")
l8_rgyr <- read_metric("../gromacs/L861Q_Erlotinib/rgyr_200.xvg", "L861Q", "RGYR")
df_rgyr <- bind_rows(wt_rgyr, v7_rgyr, l8_rgyr)
df_rgyr$System <- factor(df_rgyr$System, levels = systems)

wt_sasa <- read_metric("../gromacs/WT_Erlotinib/sasa_200.xvg", "Wild Type", "SASA")
v7_sasa <- read_metric("../gromacs/V774M_Erlotinib/sasa_200.xvg", "V774M", "SASA")
l8_sasa <- read_metric("../gromacs/L861Q_Erlotinib/sasa_200.xvg", "L861Q", "SASA")
df_sasa <- bind_rows(wt_sasa, v7_sasa, l8_sasa)
df_sasa$System <- factor(df_sasa$System, levels = systems)

# 4. Define the statistical comparisons
my_comparisons <- list(
  c("V774M", "L861Q"),
  c("V774M", "Wild Type"),
  c("L861Q", "Wild Type")
)

# 5. Build the simpler, cleaner plotting function
create_bar_plot <- function(data, y_label, title) {
  ggbarplot(data, x = "System", y = "Value",
            fill = "System", color = "black",
            palette = my_colors,
            add = "mean_sd",
            error.plot = "errorbar",
            add.params = list(width = 0.25)) + 
    
    # Text inside the bar (Pushed further down via vjust = 4.5 to avoid error bar overlap)
    stat_summary(fun = mean, geom = "text",
                 aes(label = round(after_stat(y), 3)),
                 vjust = 4.5, size = 5.5, color = "black", fontface = "plain") +
    
    # Statistical brackets configured to reach down and touch the error bars
    stat_compare_means(comparisons = my_comparisons, 
                       label = "p.signif", 
                       method = "t.test",
                       step.increase = 0.08, # Controls spacing between stacked brackets
                       tip.length = 0.04) +  # Lengthens the vertical legs to touch the error bars
    
    # Theme configuration
    theme_bw(base_size = 15) + # theme_bw provides the faint background grid naturally
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "plain", color = "black"),
      axis.text.x = element_text(face = "plain", color = "black", margin = margin(t = 5)),
      axis.text.y = element_text(face = "plain", color = "black"),
      plot.title = element_text(face = "plain", size = 16, hjust = 0, margin = margin(b = 10)),
      
      # Removing the harsh solid axis lines
      axis.line = element_blank(),
      
      # Softening the grid and borders to match the reference image aesthetics
      panel.grid.major = element_line(color = "#eaeaea", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#eaeaea", fill = NA, linewidth = 0.8)
    ) +
    labs(y = y_label, title = title)
}

# 6. Generate the four individual plots
p_rmsd <- create_bar_plot(df_rmsd, "Mean RMSD (nm)", "A)   Root Mean Square Deviation")
p_rmsf <- create_bar_plot(df_rmsf, "Mean RMSF (nm)", "B)   Root Mean Square Fluctuation")
p_rgyr <- create_bar_plot(df_rgyr, "Mean RGYR (nm)", "C)   Radius of Gyration")
p_sasa <- create_bar_plot(df_sasa, expression(paste("Mean SASA (", nm^2, ")")), "D)   Solvent Accessible Surface Area")

# 7. Stitch them together into the 2x2 grid
master_figure <- ggarrange(p_rmsd, p_rmsf, p_rgyr, p_sasa, 
                           ncol = 2, nrow = 2)

# 8. Export directly to PDF
ggsave("Combined_MD_BarPlots_Master.pdf", plot = master_figure, width = 12, height = 10)
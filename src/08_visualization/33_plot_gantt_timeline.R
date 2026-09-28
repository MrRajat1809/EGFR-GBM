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

# Set the working directory to your Docker mount point
setwd(project_path("outputs/figures"))
cat("Maximizing plot area and resetting font size to 12...\n")

# 1. Define the Timeline Data
project_data <- data.frame(
  Task = factor(c(
    "Data Retrieval & Pre-filtering",
    "Coding Mutations Analysis",
    "Molecular Docking",
    "Molecular Dynamics (200ns)",
    "Non-Coding Mutations Analysis",
    "Final Review & Viva-Voce Prep"
  ), levels = rev(c(
    "Data Retrieval & Pre-filtering",
    "Coding Mutations Analysis",
    "Molecular Docking",
    "Molecular Dynamics (200ns)",
    "Non-Coding Mutations Analysis",
    "Final Review & Viva-Voce Prep"
  ))),
  Start = as.Date(c("2026-01-26", "2026-02-05", "2026-02-28", "2026-03-15", "2026-03-30", "2026-04-18")),
  End = as.Date(c("2026-02-04", "2026-02-27", "2026-03-14", "2026-03-30", "2026-04-17", "2026-04-26")),
  Phase = c("Data Preparation", "Structural & Functional Analysis", "Structural & Functional Analysis", "Simulation", "Regulatory Analysis", "Conclusion")
)

# 2. Define Specific Milestone
viva_date <- as.Date("2026-04-25")

cat("Building the aesthetics...\n")

# 3. Build the Gantt Chart
p_gantt <- ggplot(project_data, aes(x = Start, xend = End, y = Task, yend = Task, color = Phase)) +
  geom_segment(linewidth = 8) +
  
  # Highlight the Final Viva-Voce date
  geom_vline(xintercept = as.numeric(viva_date), linetype = "dashed", color = "#C2185B", linewidth = 1) +
  # Annotation size reduced to match base_size 12
  annotate("text", x = viva_date - 1, y = 6.4, label = "Final Viva-Voce\n(April 25)", color = "#C2185B", hjust = 1, size = 4, lineheight = 0.8) +
  
  scale_color_manual(values = c(
    "Data Preparation" = "#4A90E2", 
    "Structural & Functional Analysis" = "#F5A623", 
    "Simulation" = "#7ED321", 
    "Regulatory Analysis" = "#9013FE",
    "Conclusion" = "#607D8B"
  )) +
  
  # Set font size exactly to 12
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(t = 10),
    
    plot.title = element_text(hjust = 0.5, face = "plain", margin = margin(b = 15)),
    
    # Dramatically reduced dead space so the plot expands
    plot.margin = margin(t = 20, r = 50, b = 20, l = 20),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  ) +
  
  coord_cartesian(clip = "off") +
  
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5)) +
  
  labs(
    title = "Project Timeline: Structural & Functional Effects of EGFR nsSNPs",
    x = "Timeline (Jan 26 - Apr 26)",
    y = NULL,
    color = "Project Phase"
  ) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d", limits = c(as.Date("2026-01-24"), as.Date("2026-04-28")))

cat("Saving Maximized High-Res Gantt Chart (300 DPI)...\n")

# 4. Save the Output (Wide, horizontal canvas to let the bars stretch)
png("EGFR_Project_GanttChart_MaxPlot_300DPI.png", width = 2800, height = 1200, res = 300)
print(p_gantt)
dev.off()

cat("\n--- MAXIMIZED GANTT CHART COMPLETED AND SAVED! ---\n")
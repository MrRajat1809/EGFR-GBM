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

# Install required packages (if not already present)
req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs, repos = "http://cran.us.r-project.org")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Define the Graphviz flowchart for the Massive MD Snaking Layout
flowchart <- grViz("
digraph molecular_dynamics_snake {
  
  # Graph-level aesthetics
  graph [layout = dot, rankdir = TB, nodesep = 0.5, ranksep = 0.9]
  
  # Base Node styling - increased font and margins for the massive canvas
  node [shape = box, style = 'filled, rounded', fontname = 'Helvetica-Bold', 
        fontsize = 16, color = '#333333', penwidth = 1.5, margin = '0.3,0.2']

  # --- Define Colored Nodes by Phase ---
  
  # Phase 1: System Preparation (Light Blue)
  node [fillcolor = '#E3F2FD']
  1 [label = '1. Best Docking Pose\\nSelection']
  2 [label = '2. Protein Topology\\n(AMBER99SB-ILDN)']
  3 [label = '3. Ligand Parameterization\\n(ANTECHAMBER)']
  4 [label = '4. Simulation Box Creation\\n(Cubic, 1.0 nm Cutoff)']
  5 [label = '5. Water Solvation\\n(TIP3P Model)']
  6 [label = '6. Ion Addition & Neut.\\n(0.15M Concentration)']
  
  # Phase 2: Minimization & Equilibration (Light Yellow)
  node [fillcolor = '#FFF9C4']
  7 [label = '7. Energy Minimization\\n(Steepest Descent)']
  8 [label = '8. NVT Equilibration\\n(310 K, V-rescale)']
  9 [label = '9. NPT Equilibration\\n(1 bar, Parrinello-Rahman)']
  10 [label = '10. Physics Parameters\\n(PME, LINCS, Verlet)']
  
  # Phase 3: Production MD (Light Green)
  node [fillcolor = '#E8F5E9']
  11 [label = '11. Production MD Run\\n(200 ns, 2 fs step)']
  12 [label = '12. Trajectory Mapping\\n& Processing']
  
  # Phase 4: Analysis & Validation (Light Purple)
  node [fillcolor = '#F3E5F5']
  13 [label = '13. Structural Analysis\\n(RMSD, RMSF, Rg, SASA, H-bonds)']
  14 [label = '14. Statistical Validation\\n(Welch T-Test, ANOVA)']
  15 [label = '15. Data Visualization\\n(R Studio & ggplot2)']

  # --- ROW 1 (Left to Right) ---
  { 
    rank = same;
    edge [color = '#444444', penwidth = 3, arrowsize = 1.2, dir=forward];
    1 -> 2 -> 3 -> 4 -> 5;
  }
  
  # --- THE BEND 1 (Drop down from 5 to 6) ---
  edge [color = '#444444', penwidth = 3, arrowsize = 1.2, dir=forward];
  5 -> 6 [weight=10];

  # --- ROW 2 (Right to Left Visually) ---
  # dir=back makes the arrows point leftwards to continue the S-shape
  { 
    rank = same;
    edge [color = '#444444', penwidth = 3, arrowsize = 1.2, dir=back];
    10 -> 9 -> 8 -> 7 -> 6;
  }
  
  # --- THE BEND 2 (Drop down from 10 to 11) ---
  edge [color = '#444444', penwidth = 3, arrowsize = 1.2, dir=forward];
  10 -> 11 [weight=10];

  # --- ROW 3 (Left to Right) ---
  { 
    rank = same;
    edge [color = '#444444', penwidth = 3, arrowsize = 1.2, dir=forward];
    11 -> 12 -> 13 -> 14 -> 15;
  }

  # --- INVISIBLE ALIGNMENT EDGES ---
  # This tricks Graphviz into keeping the 3 rows perfectly aligned in 5 columns
  edge [style=invis, weight=100];
  1 -> 10 -> 11;
  2 -> 9 -> 12;
  3 -> 8 -> 13;
  4 -> 7 -> 14;
  6 -> 15; # 5 to 6 is already a visible structural edge, so we link 6 to 15
}
")

# 1. Export the DiagrammeR object to raw SVG code
svg_code <- export_svg(flowchart)

# 2. Render the SVG to a Massive High-Resolution PNG 
# 4800 pixels guarantees extreme sharpness for big canvases
rsvg_png(charToRaw(svg_code), file = "MD_Simulation_FinalBoss.png", width = 4800)

cat("Success: Massive MD flowchart generated as 'MD_Simulation_FinalBoss.png'\n")
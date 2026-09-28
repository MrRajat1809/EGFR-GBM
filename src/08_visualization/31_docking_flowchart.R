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

# Define the Graphviz flowchart for a Colored Snaking Layout
flowchart <- grViz("
digraph molecular_docking_snake_colored {
  
  # Graph-level aesthetics
  graph [layout = dot, rankdir = TB, nodesep = 0.4, ranksep = 0.8]
  
  # Base Node styling
  node [shape = box, style = 'filled, rounded', fontname = 'Helvetica-Bold', 
        fontsize = 14, color = '#333333', penwidth = 1.2, margin = '0.2,0.15']

  # --- Define Colored Nodes by Phase ---
  
  # Phase 1: Preparation (Light Blue)
  node [fillcolor = '#E3F2FD']
  1 [label = '1. Best Model\\nSelection']
  2 [label = '2. Drug\\nSelection']
  3 [label = '3. Structure Cleaning\\n& Chain Joining']
  4 [label = '4. Creating AIR File\\n(Active/Passive)']
  
  # Phase 2: it0 Rigid Body (Light Yellow)
  node [fillcolor = '#FFF9C4']
  5 [label = '5. Topology Gen.\\n(topoaa)']
  6 [label = '6. Rigid-body\\nDocking (rigidbody)']
  7 [label = '7. Score Models\\n(caprieval)']
  8 [label = '8. Select Top 200\\n(seletop)']
  
  # Phase 3: it1 Flexible Refinement (Light Green)
  node [fillcolor = '#E8F5E9']
  9 [label = '9. Semi-flexible\\nRefinement (flexref)']
  10 [label = '10. Score Refined\\nModels (caprieval)']
  11 [label = '11. Select Top 100\\n(seletop)']
  
  # Phase 4: Clustering & Analysis (Light Purple)
  node [fillcolor = '#F3E5F5']
  12 [label = '12. FCC Clustering\\n(clustfcc)']
  13 [label = '13. Pick Top Clusters\\n(seletopclusts)']
  14 [label = '14. Final Cluster\\nEvaluation (caprieval)']
  15 [label = '15. Interface Contact\\nMap (contactmap)']

  # --- TOP ROW (Left to Right) ---
  { 
    rank = same;
    edge [color = '#666666', penwidth = 2.5, arrowsize = 1.2, dir=forward];
    1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8;
  }
  
  # --- THE BEND (Drop down from 8 to 9) ---
  edge [color = '#666666', penwidth = 2.5, arrowsize = 1.2, dir=forward];
  8 -> 9 [weight=10];

  # --- BOTTOM ROW (Right to Left Visually) ---
  { 
    rank = same;
    edge [color = '#666666', penwidth = 2.5, arrowsize = 1.2, dir=back];
    15 -> 14 -> 13 -> 12 -> 11 -> 10 -> 9;
  }

  # --- INVISIBLE ALIGNMENT EDGES ---
  edge [style=invis, weight=100];
  2 -> 15;
  3 -> 14;
  4 -> 13;
  5 -> 12;
  6 -> 11;
  7 -> 10;
}
")

# 1. Export the DiagrammeR object to raw SVG code
svg_code <- export_svg(flowchart)

# 2. Render the SVG to a High-Resolution PNG 
rsvg_png(charToRaw(svg_code), file = "Molecular_Docking_Snake_Colored.png", width = 3600)

cat("Success: Colored flowchart generated as high-resolution 'Molecular_Docking_Snake_Colored.png'\n")
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

# Install required packages
req_pkgs <- c("DiagrammeR", "DiagrammeRsvg", "rsvg")
new_pkgs <- req_pkgs[!(req_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs, repos = "http://cran.us.r-project.org")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Define the Graphviz flowchart
flowchart <- grViz("
digraph homology_modelling {
  
  # Graph-level aesthetics
  graph [layout = dot, rankdir = TB, nodesep = 0.8, ranksep = 0.6, compound = true]
  
  # Default Node styling 
  node [shape = box, style = 'filled, rounded', fillcolor = '#FFFFFF', 
        fontname = 'Helvetica-Bold', fontsize = 12, 
        color = '#333333', penwidth = 1.2, margin = '0.2,0.1']
  
  # Default Edge styling 
  edge [color = '#666666', penwidth = 2, arrowsize = 1.2]

  # --- Define Background Phases (Clusters) ---
  subgraph cluster_phase1 {
    style = 'filled, rounded'
    fillcolor = '#E0F7FA' 
    color = 'transparent'
    1; 2;
  }
  
  subgraph cluster_phase2 {
    style = 'filled, rounded'
    fillcolor = '#F3E5F5' 
    color = 'transparent'
    3; 4;
  }
  
  subgraph cluster_phase3 {
    style = 'filled, rounded'
    fillcolor = '#FFF3E0' 
    color = 'transparent'
    5; 5_sub;
  }
  
  subgraph cluster_phase4 {
    style = 'filled, rounded'
    fillcolor = '#E8F5E9' 
    color = 'transparent'
    6;
  }

  # --- Define Nodes ---
  1 [label = '1. Sequence Retrieval']
  2 [label = '2. Template Search']
  3 [label = '3. Model Building']
  4 [label = '4. Loop Modelling']
  5 [label = '5. Structural Assessment']
  6 [label = '6. Superimpose']

  # Callout node styling 
  node [shape = note, fontname = 'Helvetica', fontsize = 10, fillcolor = '#FFFFFF']
  5_sub [label = 'Structural Assessment includes:\\l  • SWISS Model Structural Assessment\\l  • SAVESv6.1\\l']

  # --- Define Edges ---
  1 -> 2 -> 3 -> 4 -> 5 -> 6

  # Force the callout to be on the same horizontal level as node 5
  {rank = same; 5; 5_sub;}
  
  # Invisible edge or styled edge to connect the callout
  5 -> 5_sub [style = dashed, arrowhead = none, penwidth = 1.5, color = '#888888']
}
")

# 1. Export the DiagrammeR object to raw SVG code
svg_code <- export_svg(flowchart)

# 2. Render the SVG to a High-Resolution PNG
# Setting width to 2400 pixels yields roughly 300 DPI for a standard publication column/page width.
rsvg_png(charToRaw(svg_code), file = "Homology_Modelling_Flowchart.png", width = 2400)

cat("Success: Flowchart generated as high-resolution 'Homology_Modelling_Flowchart.png'\n")
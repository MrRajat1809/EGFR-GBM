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

# sequence_plot.R
library(grid)

# 1. Open a PNG graphics device (High Quality)
png("fasta_annotated_super_enhancer.png", width = 14, height = 12, units = "in", res = 300)
grid.newpage()

# 2. Custom function to draw sequence with dynamic arrow text
draw_dna <- function(header, seq, y_pos, hl_start = NA, hl_end = NA, hl_color = "black", underline_idx = 16, do_arrow = FALSE, arrow_text = "") {
  
  # Draw the sequence header
  grid.text(header, x = 0.5, y = y_pos + 0.035, just = "center", 
            gp = gpar(fontsize = 18, fontfamily = "mono", fontface = "italic", col = "gray20"))
  
  # Split the sequence
  chars <- strsplit(seq, "")[[1]]
  n <- length(chars)
  
  # Calculate X spacing
  spacing <- 0.022 
  start_x <- 0.5 - (n / 2) * spacing + (spacing / 2)
  
  for (i in 1:n) {
    current_col <- "black"
    if (!is.na(hl_start) && i >= hl_start && i <= hl_end) {
      current_col <- hl_color
    }
    
    cx <- start_x + (i - 1) * spacing
    
    # Draw character
    grid.text(chars[i], x = cx, y = y_pos, just = "center", 
              gp = gpar(fontsize = 28, fontfamily = "mono", fontface = "bold", lwd = 2, col = current_col))
    
    # Draw underline
    if (!is.na(underline_idx) && i == underline_idx) {
      grid.lines(x = unit(c(cx - spacing*0.4, cx + spacing*0.4), "npc"),
                 y = unit(c(y_pos - 0.015, y_pos - 0.015), "npc"),
                 gp = gpar(col = current_col, lwd = 4))
    }
  }
  
  # Draw the annotation arrow and custom text
  if (do_arrow) {
    mid_i <- (hl_start + hl_end) / 2
    mid_x <- start_x + (mid_i - 1) * spacing
    
    # Arrow
    grid.lines(x = unit(c(mid_x, mid_x), "npc"),
               y = unit(c(y_pos - 0.10, y_pos - 0.03), "npc"),
               arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
               gp = gpar(col = "red", fill = "red", lwd = 3))
    
    # Dynamic Label
    grid.text(arrow_text,
              x = mid_x, y = y_pos - 0.12, just = "top",
              gp = gpar(fontsize = 20, fontfamily = "mono", fontface = "bold", col = "red", lineheight = 1.2))
  }
}

# 3. Draw the sequences (Spacing adjusted for 5 rows)

# --- LOCUS 1: The OLIG2 Docking Bay ---
draw_dna(header = ">chr7_55259524_T_A_WildType",
         seq = "TTGGGCTGGCCAAACTGCTGGGTGCGGAAGA",
         y_pos = 0.90, hl_start = 15, hl_end = 20, hl_color = "green4", underline_idx = 16)

draw_dna(header = ">chr7_55259524_T_A_Mutated (Forward Strand)",
         seq = "TTGGGCTGGCCAAACAGCTGGGTGCGGAAGA",
         y_pos = 0.74, hl_start = 15, hl_end = 20, hl_color = "red", underline_idx = 16, 
         do_arrow = TRUE, arrow_text = "Transcription Factor\n(OLIG2 motif)")

# --- LOCUS 2: The SOX2 Docking Bay ---
draw_dna(header = ">chr7_55259509_T_G_WildType",
         seq = "TCAAGATCACAGATTTTGGGCTGGCCAAACT",
         y_pos = 0.50, underline_idx = 16)

# Mutated Forward (Highlighting the structural break G)
draw_dna(header = ">chr7_55259509_T_G_Mutated (Forward Strand)",
         seq = "TCAAGATCACAGATTGTGGGCTGGCCAAACT",
         y_pos = 0.34, hl_start = 16, hl_end = 16, hl_color = "red", underline_idx = 16)

# Mutated Reverse Complement (Highlighting the ACAAT motif)
draw_dna(header = ">chr7_55259509_T_G_Mutated (Reverse Complement)",
         seq = "AGTTTGGCCAGCCCACAATCTGTGATCTTGA",
         y_pos = 0.18, hl_start = 15, hl_end = 19, hl_color = "red", underline_idx = 16, 
         do_arrow = TRUE, arrow_text = "Transcription Factor\n(SOX2 motif)")

# 4. Close and save
dev.off()

cat("Success! The objective, publication-ready plot 'fasta_annotated_super_enhancer.png' has been saved.\n")

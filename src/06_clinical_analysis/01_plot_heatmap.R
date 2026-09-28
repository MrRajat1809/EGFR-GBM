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

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(ggplot2)
library(dplyr)
library(pheatmap)
library(grid)      
library(gtable)    

setwd(project_path("outputs/figures"))
cat("Starting Pan-Cancer Heatmap Generation (Ultra-Wide, Anti-Cropping Edition)...\n")
gene_target <- "EGFR"

# =========================================================================
# 1. CALCULATE EXPRESSION SCORES DIRECTLY FROM LOCAL MASTER FILES
# =========================================================================
cat("1. Extracting Pan-Cancer Expression from local EB++ file...\n")
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
clin_data <- tryCatch(read.delim(project_path("data/raw/clinical/Survival_SupplementalTable_S1_20171025_xena_sp"), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data) || is.null(clin_data)) stop("CRITICAL ERROR: Pan-Cancer Xena files missing!")

gene_col <- colnames(expr_data)[1]
egfr_row <- expr_data[grepl(paste0("^", gene_target, "$|^", gene_target, "\\|"), expr_data[[gene_col]], ignore.case=TRUE), ]
expr_vals <- as.numeric(egfr_row[1, -1])
df_expr <- data.frame(sample = colnames(expr_data)[-1], Expression = expr_vals, stringsAsFactors=FALSE)

# Isolate Primary Tumors ("01")
df_expr$Sample_Type <- substr(df_expr$sample, 14, 15)
df_expr <- df_expr[df_expr$Sample_Type == "01", ]
df_expr$Patient_ID <- substr(df_expr$sample, 1, 12)
df_expr <- df_expr[!duplicated(df_expr$Patient_ID), ]

clin_sub <- clin_data[, c("sample", "cancer type abbreviation")]
colnames(clin_sub) <- c("Patient_ID", "Cancer")
clin_sub$Patient_ID <- substr(clin_sub$Patient_ID, 1, 12)
clin_sub <- clin_sub[!duplicated(clin_sub$Patient_ID), ]

merged <- merge(df_expr, clin_sub, by="Patient_ID")

# Calculate Median Expression per Cancer and scale it to a Z-Score
expr_summary <- merged %>%
  group_by(Cancer) %>%
  summarize(Median_Expr = median(Expression, na.rm = TRUE)) %>%
  mutate(Expression_Score = as.numeric(scale(Median_Expr)))

# =========================================================================
# 2. LOAD & PROCESS COMBINED MOLECULAR SCORE FROM CBIOPORTAL
# =========================================================================
cat("2. Loading and processing cBioPortal Molecular Data...\n")
cbio_file <- project_path("data/raw/clinical/cancer_types_summary.txt") 

if (!file.exists(cbio_file)) stop("CRITICAL ERROR: cancer_types_summary.txt not found in folder!")

cbio <- read.delim(cbio_file, stringsAsFactors=FALSE, check.names=FALSE)

study_col <- grep("Study", colnames(cbio), ignore.case = TRUE, value = TRUE)[1]
freq_col <- grep("Freq|Percent", colnames(cbio), ignore.case = TRUE, value = TRUE)[1]

cbio$Mol_Freq <- as.numeric(gsub("%", "", cbio[[freq_col]]))

cbio_sum <- cbio %>%
  group_by(!!sym(study_col)) %>%
  summarize(Total_Freq = sum(Mol_Freq, na.rm = TRUE))

cbio_sum$Clean_Study <- trimws(gsub("\\(TCGA.*\\)", "", cbio_sum[[study_col]]))

tcga_map <- c(
  "Acute Myeloid Leukemia" = "LAML", "Adrenocortical Carcinoma" = "ACC",
  "Bladder Urothelial Carcinoma" = "BLCA", "Brain Lower Grade Glioma" = "LGG",
  "Breast Invasive Carcinoma" = "BRCA", "Cervical Squamous Cell Carcinoma" = "CESC",
  "Cholangiocarcinoma" = "CHOL", "Colorectal Adenocarcinoma" = "COAD",
  "Esophageal Adenocarcinoma" = "ESCA", "Glioblastoma Multiforme" = "GBM",
  "Head and Neck Squamous Cell Carcinoma" = "HNSC", "Kidney Chromophobe" = "KICH",
  "Kidney Renal Clear Cell Carcinoma" = "KIRC", "Kidney Renal Papillary Cell Carcinoma" = "KIRP",
  "Liver Hepatocellular Carcinoma" = "LIHC", "Lung Adenocarcinoma" = "LUAD",
  "Lung Squamous Cell Carcinoma" = "LUSC", "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "DLBC",
  "Mesothelioma" = "MESO", "Ovarian Serous Cystadenocarcinoma" = "OV",
  "Pancreatic Adenocarcinoma" = "PAAD", "Pheochromocytoma and Paraganglioma" = "PCPG",
  "Prostate Adenocarcinoma" = "PRAD", "Sarcoma" = "SARC",
  "Skin Cutaneous Melanoma" = "SKCM", "Stomach Adenocarcinoma" = "STAD",
  "Testicular Germ Cell Tumors" = "TGCT", "Thymoma" = "THYM",
  "Thyroid Carcinoma" = "THCA", "Uterine Carcinosarcoma" = "UCS",
  "Uterine Corpus Endometrial Carcinoma" = "UCEC", "Uveal Melanoma" = "UVM"
)

cbio_sum$Cancer <- tcga_map[cbio_sum$Clean_Study]
expr_summary <- merge(expr_summary, cbio_sum[, c("Cancer", "Total_Freq")], by = "Cancer", all.x = TRUE)

expr_summary$Total_Freq[is.na(expr_summary$Total_Freq)] <- 0

if ("COAD" %in% expr_summary$Cancer && "READ" %in% expr_summary$Cancer) {
  if (expr_summary$Total_Freq[expr_summary$Cancer == "READ"] == 0) {
    expr_summary$Total_Freq[expr_summary$Cancer == "READ"] <- expr_summary$Total_Freq[expr_summary$Cancer == "COAD"]
  }
}

expr_summary$Molecular_Score <- as.numeric(scale(expr_summary$Total_Freq))

# =========================================================================
# 3. BUILD THE PLOT MATRIX
# =========================================================================
cat("3. Building the Heatmap Matrix...\n")
cancers <- expr_summary$Cancer

plot_matrix <- rbind(
  expr_summary$Expression_Score,
  expr_summary$Molecular_Score
)
colnames(plot_matrix) <- cancers
rownames(plot_matrix) <- c("Expression score", "Combined molecular score")

plot_matrix[plot_matrix > 4] <- 4
plot_matrix[plot_matrix < -4] <- -4

# =========================================================================
# 4. RENDER THE HEATMAP & MODIFY GTABLE
# =========================================================================
cat("4. Rendering the Hierarchical Heatmap...\n")

palette_colors <- colorRampPalette(c("blue", "white", "red"))(100)

ph <- pheatmap(
  plot_matrix,
  color = palette_colors,
  breaks = seq(-4, 4, length.out = 101),
  cluster_rows = FALSE,      
  cluster_cols = TRUE,       
  treeheight_col = 80,       
  show_rownames = TRUE,      
  show_colnames = TRUE,
  fontsize_col = 14,
  fontsize_row = 16,
  angle_col = 90,            
  cellwidth = 30,            
  cellheight = 30,           
  border_color = "black",    # Set inner lines base color to black
  legend = FALSE,            
  main = NA,                 
  silent = TRUE 
)

# --- 4.1 MOVE ROW LABELS TO FAR LEFT ---
rn_idx <- which(ph$gtable$layout$name == "row_names")
mat_idx <- which(ph$gtable$layout$name == "matrix")

# Enforce exactly 1pt thickness on inner heatmap lines
ph$gtable$grobs[[mat_idx]] <- editGrob(ph$gtable$grobs[[mat_idx]], gp = gpar(col = "black", lwd = 1))

rn_col_orig <- ph$gtable$layout$l[rn_idx]
mat_col_orig <- ph$gtable$layout$l[mat_idx]
rn_width <- ph$gtable$widths[rn_col_orig]

ph$gtable <- gtable_add_cols(ph$gtable, rn_width, pos = mat_col_orig - 1)

new_rn_col <- mat_col_orig
ph$gtable$layout$l[rn_idx] <- new_rn_col
ph$gtable$layout$r[rn_idx] <- new_rn_col

shifted_old_rn_col <- rn_col_orig + 1
ph$gtable$widths[shifted_old_rn_col] <- unit(0, "cm")

rn_grob <- ph$gtable$grobs[[rn_idx]]
rn_grob$x <- unit(1, "npc") - unit(5, "points") 
rn_grob$hjust <- 1 
ph$gtable$grobs[[rn_idx]] <- rn_grob

mat_pos <- ph$gtable$layout[ph$gtable$layout$name == "matrix", ]
nc <- ncol(plot_matrix)
nr <- nrow(plot_matrix)


# --- 4.2 ADD HIGHLIGHT BOX AROUND HIGHEST COMBINED MOLECULAR SCORE ---
mol_scores <- plot_matrix["Combined molecular score", ]
max_col_name <- names(mol_scores)[which.max(mol_scores)]

orig_col_idx <- which(colnames(plot_matrix) == max_col_name)
plot_col_idx <- which(ph$tree_col$order == orig_col_idx)

hl_rect <- rectGrob(
  x = (plot_col_idx - 0.5) / nc,
  y = (nr - 2 + 0.5) / nr,   
  width = 1 / nc,
  height = 1 / nr,
  gp = gpar(col = "blue", fill = NA, lwd = 4, lty = "dashed") 
)

ph$gtable <- gtable_add_grob(
  ph$gtable, hl_rect, 
  t = mat_pos$t, l = mat_pos$l, b = mat_pos$b, r = mat_pos$r, 
  name = "highlight_box"
)


# --- 4.3 DENDROGRAM LINE THICKNESS (lwd = 2) ---
tree_idx <- which(ph$gtable$layout$name == "col_tree")
if(length(tree_idx) > 0) {
  ph$gtable$grobs[[tree_idx]] <- editGrob(ph$gtable$grobs[[tree_idx]], gp = gpar(lwd = 2))
}


# --- 4.4 BOLD 2-POINT PERIPHERAL BORDER AROUND MATRIX ---
peri_border <- rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2))
ph$gtable <- gtable_add_grob(
  ph$gtable, peri_border, 
  t = mat_pos$t, l = mat_pos$l, b = mat_pos$b, r = mat_pos$r, 
  name = "peripheral_border"
)


# --- 4.5 ADD CUSTOM HORIZONTAL LEGEND TO BOTTOM-LEFT ---
# Add a new row at the bottom of the plot
ph$gtable <- gtable_add_rows(ph$gtable, heights = unit(1.5, "inches"), pos = -1)
bottom_row <- nrow(ph$gtable)

# Build legend aligned entirely to the right edge of its bounding box (which will be the heatmap edge)
leg_grob <- gTree(children = gList(
  textGrob("Normalized Score", x = unit(1, "npc"), y = 0.65, just = "right", gp = gpar(fontsize = 16, fontface = "plain")),
  rasterGrob(matrix(palette_colors, nrow=1), x = unit(1, "npc"), y = 0.4, just = "right", width = unit(2.5, "inches"), height = unit(0.2, "inches"), interpolate=FALSE),
  rectGrob(x = unit(1, "npc"), y = 0.4, just = "right", width = unit(2.5, "inches"), height = unit(0.2, "inches"), gp = gpar(col = "black", fill = NA, lwd = 1)),
  textGrob("4", x = unit(1, "npc"), y = 0.15, just = "center", gp = gpar(fontsize = 14, fontface = "plain")),
  textGrob("0", x = unit(1, "npc") - unit(1.25, "inches"), y = 0.15, just = "center", gp = gpar(fontsize = 14, fontface = "plain")),
  textGrob("-4", x = unit(1, "npc") - unit(2.5, "inches"), y = 0.15, just = "center", gp = gpar(fontsize = 14, fontface = "plain"))
))

# Insert the legend in the bottom row, directly inside the row labels column
ph$gtable <- gtable_add_grob(
  ph$gtable, leg_grob, 
  t = bottom_row, l = new_rn_col, b = bottom_row, r = new_rn_col, 
  name = "horiz_legend_bottom_left", clip = "off"
)


# --- 4.6 EXPORT IMAGE ---
png("Fig5C_Custom_Heatmap_FINAL_300DPI.png", width = 6600, height = 3000, res = 300)
grid.newpage()
grid.draw(ph$gtable)
dev.off()

cat("\n--- PANEL 5C (HEATMAP) UNCLIPPED & PERFECTED! ---\n")

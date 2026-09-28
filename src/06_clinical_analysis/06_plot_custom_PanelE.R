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
library(ggpubr)

setwd(project_path("outputs/figures"))
cat("Starting the Master Panel E (Compact Portrait Edition)...\n")
gene_target <- "EGFR"

cat("1. Loading Gold-Standard Pan-Cancer Data...\n")
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
clin_data <- tryCatch(read.delim(project_path("data/raw/clinical/Survival_SupplementalTable_S1_20171025_xena_sp"), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data) || is.null(clin_data)) {
  stop("CRITICAL ERROR: Could not find the Pan-Cancer Xena files!")
}

cat("2. Extracting EGFR & Strictly Filtering for Primary Tumors...\n")
gene_col <- colnames(expr_data)[1]
egfr_row <- expr_data[grepl(paste0("^", gene_target, "$|^", gene_target, "\\|"), expr_data[[gene_col]], ignore.case=TRUE), ]
if(nrow(egfr_row) == 0) stop("CRITICAL ERROR: EGFR not found.")

expr_vals <- as.numeric(egfr_row[1, -1])
expr_samples <- colnames(expr_data)[-1]
df_expr <- data.frame(sample = expr_samples, Expression = expr_vals, stringsAsFactors=FALSE)

# STRICT PRIMARY TUMOR FILTER ("01")
df_expr$Sample_Type <- substr(df_expr$sample, 14, 15)
df_expr <- df_expr[df_expr$Sample_Type == "01", ]
df_expr$Patient_ID <- substr(df_expr$sample, 1, 12)
df_expr <- df_expr[!duplicated(df_expr$Patient_ID), ]

cat("3. Merging Curated TCGA-CDR Clinical Data...\n")
clin_sub <- clin_data[, c("sample", "cancer type abbreviation", "histological_grade", "histological_type")]
colnames(clin_sub) <- c("Patient_ID", "Cancer", "Grade", "Histology")
clin_sub$Patient_ID <- substr(clin_sub$Patient_ID, 1, 12)
clin_sub <- clin_sub[!duplicated(clin_sub$Patient_ID), ]

merged <- merge(df_expr, clin_sub, by="Patient_ID")
merged <- merged[merged$Cancer %in% c("GBM", "LGG"), ]
merged$Grade <- ifelse(merged$Cancer == "GBM", "G4", merged$Grade)

# Calculate global Y-axis limits so both plots match perfectly
min_expr <- min(merged$Expression, na.rm = TRUE)
max_expr <- max(merged$Expression, na.rm = TRUE)

custom_colors <- c("#FF8C00", "#9ACD32", "#4169E1", "#43A047", "#E53935", "#8E24AA")
plot_list <- list()

cat("4. Applying Compact Portrait Visuals...\n")

# --- PLOT 1: GRADE ---
plot_data_g <- merged[!is.na(merged$Grade) & merged$Grade != "" & !grepl("not reported|unknown", merged$Grade, ignore.case=TRUE), ]

if(nrow(plot_data_g) > 0) {
  plot_data_g$Category <- case_when(
    grepl("G2|Grade II\\b", plot_data_g$Grade, ignore.case=TRUE) ~ "G2",
    grepl("G3|Grade III\\b", plot_data_g$Grade, ignore.case=TRUE) ~ "G3",
    grepl("G4|Grade IV\\b", plot_data_g$Grade, ignore.case=TRUE) ~ "G4",
    TRUE ~ "Other"
  )
  plot_data_g <- plot_data_g[plot_data_g$Category != "Other", ]
  plot_data_g$Category <- factor(plot_data_g$Category, levels = c("G2", "G3", "G4"))
  
  p_grade <- ggplot(plot_data_g, aes(x = Category, y = Expression)) +
    stat_boxplot(geom = "errorbar", width = 0.2, linetype = "dashed", color = "black") +
    geom_boxplot(linetype = "dashed", color = "black", fill = NA, outlier.shape = 16, outlier.size = 1) +
    geom_boxplot(aes(fill = Category), color = "transparent", outlier.shape = NA) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 1) +
    scale_fill_manual(values = custom_colors) +
    
    # NEW: Title at the top, lock Y-axis scale
    ggtitle("Grade") +
    scale_y_continuous(limits = c(min_expr, max_expr)) +
    
    theme_classic(base_size = 12) +
    theme(
      # NEW: Angle 90 (perpendicular), vjust 0.5 centers it to the tick mark
      axis.text.x = element_text(color = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.x = element_blank(), # Removed bottom title
      axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 10)),
      plot.title = element_text(hjust = 0.5, size = 14, face = "plain", margin = margin(b = 15)),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 5, b = 10, l = 10)
    ) +
    labs(y = paste("TPM of", gene_target))
    
  plot_list[["Grade"]] <- p_grade
}

# --- PLOT 2: HISTOLOGY ---
plot_data_h <- merged[!is.na(merged$Histology) & merged$Histology != "" & !grepl("not reported|unknown", merged$Histology, ignore.case=TRUE), ]

if(nrow(plot_data_h) > 0) {
  plot_data_h$Category <- gsub(",.*", "", plot_data_h$Histology) 
  plot_data_h$Category <- as.factor(plot_data_h$Category)
  
  p_hist <- ggplot(plot_data_h, aes(x = Category, y = Expression)) +
    stat_boxplot(geom = "errorbar", width = 0.2, linetype = "dashed", color = "black") +
    geom_boxplot(linetype = "dashed", color = "black", fill = NA, outlier.shape = 16, outlier.size = 1) +
    geom_boxplot(aes(fill = Category), color = "transparent", outlier.shape = NA) +
    stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", linewidth = 1) +
    scale_fill_manual(values = custom_colors) +
    
    # NEW: Title at the top, lock Y-axis scale
    ggtitle("Histology") +
    scale_y_continuous(limits = c(min_expr, max_expr)) +
    
    theme_classic(base_size = 12) +
    theme(
      # NEW: Angle 90 (perpendicular)
      axis.text.x = element_text(color = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.x = element_blank(), # Removed bottom title
      axis.title.y = element_blank(), # Removed Y title so it stitches closely to Plot 1
      plot.title = element_text(hjust = 0.5, size = 14, face = "plain", margin = margin(b = 15)),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 10, b = 10, l = 5)
    ) 
    
  plot_list[["Histology"]] <- p_hist
}

cat("5. Saving the High-Res Portrait Image (300 DPI)...\n")
if(length(plot_list) > 0) {
  # Arrange them side-by-side, perfectly aligned
  final_plot <- ggarrange(plotlist = plot_list, ncol = length(plot_list), nrow = 1, align = "h")
  
  # NEW: TALL Portrait Dimensions (Width is smaller than Height)
  png("Fig6E_Custom_ClinicalBoxplots_Portrait_300DPI.png", width = 1400, height = 1800, res = 300)
  print(final_plot)
  dev.off()
  cat("\n--- PANEL E COMPACT PORTRAIT PERFECTED! ---\n")
} else {
  cat("\nERROR: No valid Grade/Histology data could be extracted.\n")
}

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

# Install the dual-scale package if it is missing
if (!require("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale", repos = "http://cran.us.r-project.org")
}

library(ggplot2)
library(dplyr)
library(ggnewscale)

setwd(project_path("outputs/figures"))
cat("Starting Final Polish for Panel D (Local MD Anderson Data)...\n")
gene_target <- "EGFR"

cat("1. Loading Local Xena Pan-Cancer Expression Data...\n")
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data)) stop("CRITICAL ERROR: Could not find the Pan-Cancer Xena file!")

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

cat("2. Loading Official MD Anderson ESTIMATE Scores...\n")
# Read the provided txt files directly
gbm_est <- tryCatch(read.delim(project_path("data/raw/clinical/GBM_estimate.txt"), stringsAsFactors=FALSE), error=function(e) NULL)
lgg_est <- tryCatch(read.delim(project_path("data/raw/clinical/LGG_estimate.txt"), stringsAsFactors=FALSE), error=function(e) NULL)

if(is.null(gbm_est) || is.null(lgg_est)) {
  stop("CRITICAL ERROR: Could not find GBM_estimate.txt or LGG_estimate.txt!")
}

clean_est <- function(df, cancer_name) {
  # Grab the exact columns from the MD Anderson files safely
  data.frame(
    Patient_ID = substr(df$ID, 1, 12),
    ESTIMATE = as.numeric(df$ESTIMATE_score),
    Immune = as.numeric(df$Immune_score),
    Stromal = as.numeric(df$Stromal_score),
    Cancer = cancer_name,
    stringsAsFactors = FALSE
  )
}

gbm_clean <- clean_est(gbm_est, "GBM")
lgg_clean <- clean_est(lgg_est, "LGG")

meta_combined <- rbind(gbm_clean, lgg_clean)
# Deduplicate just in case
meta_combined <- meta_combined[!duplicated(meta_combined$Patient_ID), ]

cat("3. Merging and Calculating EXACT Statistical Math...\n")
merged <- merge(df_expr, meta_combined, by="Patient_ID")

# Drop any patients missing ANY score to ensure strict 1:1 math
merged <- merged[!is.na(merged$Expression) & !is.na(merged$ESTIMATE) & !is.na(merged$Immune) & !is.na(merged$Stromal), ]

res_list <- list()
for (canc in c("GBM", "LGG")) {
  sub_df <- merged[merged$Cancer == canc, ]
  
  if (nrow(sub_df) > 10) {
    # THE MATH FIX: Calculate exact Pearson correlation using the TRUE patient N
    actual_n <- nrow(sub_df)
    
    cor_est <- cor.test(sub_df$Expression, sub_df$ESTIMATE, method = "pearson")
    cor_imm <- cor.test(sub_df$Expression, sub_df$Immune, method = "pearson")
    cor_str <- cor.test(sub_df$Expression, sub_df$Stromal, method = "pearson")
    
    res_list[[length(res_list) + 1]] <- data.frame(Cancer = canc, Score = "ESTIMATE", R = cor_est$estimate, Pval = cor_est$p.value, N = actual_n)
    res_list[[length(res_list) + 1]] <- data.frame(Cancer = canc, Score = "Immune", R = cor_imm$estimate, Pval = cor_imm$p.value, N = actual_n)
    res_list[[length(res_list) + 1]] <- data.frame(Cancer = canc, Score = "Stromal", R = cor_str$estimate, Pval = cor_str$p.value, N = actual_n)
  }
}

df_plot <- do.call(rbind, res_list)
rownames(df_plot) <- NULL

cat("4. Formatting P-Values for the Split Heatmap...\n")
df_plot$Pval[df_plot$Pval < 1e-16] <- 1e-16
df_plot$NegLogP <- -log10(df_plot$Pval)

df_plot$Score <- factor(df_plot$Score, levels = c("ESTIMATE", "Immune", "Stromal"))
df_plot$Cancer <- factor(df_plot$Cancer, levels = c("LGG", "GBM")) 

df_plot$x <- as.numeric(df_plot$Score)
df_plot$y <- as.numeric(df_plot$Cancer)

cat("5. Generating Split-Triangle Coordinates...\n")
poly_cor <- data.frame()
poly_pval <- data.frame()

for(i in 1:nrow(df_plot)) {
  cx <- df_plot$x[i]
  cy <- df_plot$y[i]
  
  poly_cor <- rbind(poly_cor, data.frame(
    id = i, R = df_plot$R[i],
    x = c(cx - 0.5, cx + 0.5, cx + 0.5),
    y = c(cy - 0.5, cy - 0.5, cy + 0.5)
  ))
  
  poly_pval <- rbind(poly_pval, data.frame(
    id = i, NegLogP = df_plot$NegLogP[i],
    x = c(cx - 0.5, cx - 0.5, cx + 0.5),
    y = c(cy - 0.5, cy + 0.5, cy + 0.5)
  ))
}

cat("6. Building the Scientifically Accurate Split-Cell Heatmap...\n")
pD_custom <- ggplot() +
  
  geom_polygon(data = poly_pval, aes(x = x, y = y, group = id, fill = NegLogP), color = "black", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "#E65100", name = "P-value\n(-log10)") + 
  new_scale_fill() +
  
  geom_polygon(data = poly_cor, aes(x = x, y = y, group = id, fill = R), color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "#005b96", mid = "white", high = "#C2185B", midpoint = 0, limit = c(-1, 1), name = "Correlation (R)") +
  
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 2.55, ymax = 2.75, fill = "#1E88E5", color="black", linewidth=0.5) + 
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = 2.55, ymax = 2.75, fill = "#FF8C00", color="black", linewidth=0.5) + 
  annotate("rect", xmin = 2.55, xmax = 3.45, ymin = 2.55, ymax = 2.75, fill = "#43A047", color="black", linewidth=0.5) + 
  
  annotate("text", x = 1, y = 2.92, label = "ESTIMATE", size = 3.17, color = "black") +
  annotate("text", x = 2, y = 2.92, label = "Immune", size = 3.17, color = "black") +
  annotate("text", x = 3, y = 2.92, label = "Stromal", size = 3.17, color = "black") +
  
  scale_x_continuous(limits = c(0.5, 3.5)) +
  scale_y_continuous(breaks = 1:2, labels = levels(df_plot$Cancer), limits = c(0.5, 3.1)) +
  
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(color = "black", size = 9, face = "plain"), 
    legend.position = "right",
    legend.title = element_text(size = 9, face = "plain"), 
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  coord_fixed(ratio = 1, clip = "off") 

cat("7. Saving the High-Res Image (300 DPI)...\n")
png("Fig6D_Custom_SplitHeatmap_Scientific_300DPI.png", width = 1600, height = 1000, res = 300)
print(pD_custom)
dev.off()

cat("\n--- PANEL D SCIENTIFICALLY PERFECTED! ---\n")
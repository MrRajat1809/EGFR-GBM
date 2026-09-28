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
library(survival)

setwd(project_path("outputs/figures"))
cat("Starting Final Polish for Panel C (Scientific Forest Plot)...\n")
gene_target <- "EGFR"

cat("1. Loading Local Xena Pan-Cancer Data...\n")
# Load Batch-Corrected Expression and TCGA-CDR Clinical/Survival Data
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
clin_data <- tryCatch(read.delim(project_path("data/raw/clinical/Survival_SupplementalTable_S1_20171025_xena_sp"), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data) || is.null(clin_data)) {
  stop("CRITICAL ERROR: Could not find the Pan-Cancer Xena files!")
}

cat("2. Extracting EGFR and Formatting...\n")
gene_col <- colnames(expr_data)[1]
egfr_row <- expr_data[grepl(paste0("^", gene_target, "$|^", gene_target, "\\|"), expr_data[[gene_col]], ignore.case=TRUE), ]

if(nrow(egfr_row) == 0) stop("CRITICAL ERROR: EGFR not found in the expression file.")

expr_vals <- as.numeric(egfr_row[1, -1])
expr_samples <- colnames(expr_data)[-1]
df_expr <- data.frame(sample = expr_samples, Expression = expr_vals, stringsAsFactors=FALSE)

# THE SCIENCE FIX 1: Strictly isolate Primary Solid Tumors ("01")
df_expr$Sample_Type <- substr(df_expr$sample, 14, 15)
df_expr <- df_expr[df_expr$Sample_Type == "01", ]
df_expr$Patient_ID <- substr(df_expr$sample, 1, 12)
df_expr <- df_expr[!duplicated(df_expr$Patient_ID), ]

cat("3. Merging Clinical Data and Calculating Hazard Ratios...\n")
# Extract Cancer Type and Overall Survival (OS, OS.time) from TCGA-CDR
cancer_col <- "cancer type abbreviation"
if(!(cancer_col %in% colnames(clin_data))) cancer_col <- "type" # fallback

clin_clean <- clin_data[, c("sample", cancer_col, "OS", "OS.time")]
colnames(clin_clean) <- c("Patient_ID", "Cancer", "OS", "OS.time")
clin_clean$Patient_ID <- substr(clin_clean$Patient_ID, 1, 12)
clin_clean <- clin_clean[!duplicated(clin_clean$Patient_ID), ]

merged <- merge(df_expr, clin_clean, by="Patient_ID", all.x=TRUE)
merged <- merged[!is.na(merged$OS) & !is.na(merged$OS.time) & !is.na(merged$Expression) & merged$OS.time > 0, ]

# Convert days to months
merged$time_months <- merged$OS.time / 30.4167

target_cancers <- c("GBM", "LGG")
res_list <- list()

# Run the strict Cox Proportional-Hazards Regression
for (cancer in target_cancers) {
  cancer_df <- merged[merged$Cancer == cancer, ]
  if(nrow(cancer_df) > 20) { 
    cox_model <- tryCatch(coxph(Surv(time_months, OS) ~ Expression, data = cancer_df), error=function(e) NULL)
    
    if(!is.null(cox_model)) {
      s <- summary(cox_model)
      res_list[[cancer]] <- data.frame(
        Cancer = cancer,
        HR = s$conf.int[1, 1],
        Lower = s$conf.int[1, 3],
        Upper = s$conf.int[1, 4],
        Pval = s$coefficients[1, 5]
      )
    }
  }
}

df_sub <- do.call(rbind, res_list)

cat("4. Formatting the Statistical Text...\n")
df_sub$HR_text <- sprintf("%.2f (%.2f-%.2f)", df_sub$HR, df_sub$Lower, df_sub$Upper)
df_sub$Pval_text <- ifelse(df_sub$Pval < 0.001, "<0.001", sprintf("%.3f", df_sub$Pval))
df_sub$Cancer <- factor(df_sub$Cancer, levels = c("LGG", "GBM")) 

cat("5. Building the Polished 300 DPI Visual...\n")

# Calculate the perfect explicit X-coordinates so nothing gets cut off
x_min_val <- min(df_sub$Lower, 1) - 0.3 
x_max_upper <- max(df_sub$Upper)

# Force the text to live INSIDE the coordinate grid, not in the margins
x_hr_pos <- x_max_upper + (x_max_upper * 0.25)
x_pval_pos <- x_max_upper + (x_max_upper * 0.70)
x_plot_end <- x_pval_pos + (x_max_upper * 0.15) 

pC_custom <- ggplot(df_sub, aes(y = Cancer, x = HR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.15, color = "#005b96", linewidth = 1) +
  geom_point(size = 4, color = "#C2185B") +
  
  # Cancer labels
  geom_text(aes(x = x_min_val, label = Cancer), size = 3.17, hjust = 0, color="black") +
  
  # Statistical text plotted precisely on our calculated X-coordinates
  geom_text(aes(x = x_hr_pos, label = HR_text), size = 3.17, hjust = 0, color="black") +
  geom_text(aes(x = x_pval_pos, label = Pval_text), size = 3.17, hjust = 0, color="black") +
  
  # Plain Text Headers aligned perfectly above the text
  annotate("text", x = x_min_val, y = 2.4, label = "Cancer", hjust = 0, size = 3.17) +
  annotate("text", x = x_hr_pos, y = 2.4, label = "Hazard Ratio (95% CI)", hjust = 0, size = 3.17) +
  annotate("text", x = x_pval_pos, y = 2.4, label = "p-value", hjust = 0, size = 3.17) +
  
  theme_classic(base_size = 9) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10) 
  ) +
  # x_plot_end guarantees the grid reaches past the p-value column
  coord_cartesian(clip = "off", ylim = c(1, 2.2), xlim = c(x_min_val, x_plot_end)) + 
  labs(x = "Hazard Ratio of EGFR", y = "")

cat("6. Saving High-Res Image (300 DPI)...\n")
# Dimensions strictly doubled (1400x700) to keep identical physical size but push it to 300 DPI
png("Fig6C_Custom_Forest_Polished_300DPI.png", width = 1400, height = 700, res = 300)
print(pC_custom)
dev.off()

rm(df_sub, pC_custom, res_list)
gc()
cat("\n--- PANEL C SCIENTIFICALLY PERFECTED! ---\n")

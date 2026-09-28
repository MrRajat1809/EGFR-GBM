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

setwd(project_path("outputs/figures"))
cat("Starting Final Polish for Panel B (Scientific Xena Edition)...\n")
gene_target <- "EGFR"

cat("1. Loading Local Xena Pan-Cancer Data...\n")
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
clin_data <- tryCatch(read.delim(project_path("data/raw/clinical/Survival_SupplementalTable_S1_20171025_xena_sp"), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
immune_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/Subtype_Immune_Model_Based.txt.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data) || is.null(clin_data) || is.null(immune_data)) {
  stop("CRITICAL ERROR: Could not find one of the three Xena files! Make sure Subtype_Immune_Model_Based.txt.gz is downloaded.")
}

cat("2. Extracting EGFR and Formatting...\n")
gene_col <- colnames(expr_data)[1]
egfr_row <- expr_data[grepl(paste0("^", gene_target, "$|^", gene_target, "\\|"), expr_data[[gene_col]], ignore.case=TRUE), ]

if(nrow(egfr_row) == 0) stop("CRITICAL ERROR: EGFR not found in the expression file.")

expr_vals <- as.numeric(egfr_row[1, -1])
expr_samples <- colnames(expr_data)[-1]
df_expr <- data.frame(sample = expr_samples, Expression = expr_vals, stringsAsFactors=FALSE)

# THE SCIENCE FIX 1: Thorsson Immune Subtypes are ONLY valid for Primary Tumors. 
# We explicitly strip out Normal tissue and Metastases.
df_expr$Sample_Type <- substr(df_expr$sample, 14, 15)
df_expr <- df_expr[df_expr$Sample_Type == "01", ]
df_expr$Patient_ID <- substr(df_expr$sample, 1, 12)
df_expr <- df_expr[!duplicated(df_expr$Patient_ID), ]

cat("3. Merging Clinical Cancer Abbreviations...\n")
clin_sub <- clin_data[, c("sample", "cancer type abbreviation")]
colnames(clin_sub) <- c("Patient_ID", "Cancer")
clin_sub$Patient_ID <- substr(clin_sub$Patient_ID, 1, 12)
clin_sub <- clin_sub[!duplicated(clin_sub$Patient_ID), ]

merged <- merge(df_expr, clin_sub, by="Patient_ID", all.x=TRUE)
merged <- merged[!is.na(merged$Cancer), ]

cat("4. Formatting Thorsson Immune Subtypes...\n")
colnames(immune_data)[1] <- "Raw_Barcode"
colnames(immune_data)[2] <- "Raw_Subtype"
immune_data$Patient_ID <- substr(immune_data$Raw_Barcode, 1, 12)
immune_data <- immune_data[!duplicated(immune_data$Patient_ID), ]

immune_data <- immune_data %>%
  mutate(
    Subtype_Code = case_when(
      grepl("C1", Raw_Subtype) ~ "Wound healing (C1)",
      grepl("C2", Raw_Subtype) ~ "IFN-y dominant (C2)",
      grepl("C3", Raw_Subtype) ~ "Inflammatory (C3)",
      grepl("C4", Raw_Subtype) ~ "Lymphocyte depleted (C4)",
      grepl("C5", Raw_Subtype) ~ "Immunologically quiet (C5)",
      grepl("C6", Raw_Subtype) ~ "TGF-b dominant (C6)",
      TRUE ~ "N/A"
    )
  )

full_data <- merge(merged, immune_data[, c("Patient_ID", "Subtype_Code")], by="Patient_ID", all.x=TRUE)
full_data$Cluster <- ifelse(is.na(full_data$Subtype_Code), "N/A", full_data$Subtype_Code)

cat("5. SCIENTIFIC FIX 2: Calculating Cancer-Specific Expression Ratios...\n")
plot_data <- full_data %>%
  group_by(Cancer) %>%
  # Calculate the median expression for EACH cancer independently
  mutate(Cancer_Median = median(Expression, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Cancer, Cluster) %>%
  summarise(
    Count = n(),
    # A patient is only "Overexpressed" if they beat their OWN tissue's median baseline
    Overexpressed = sum(Expression > Cancer_Median, na.rm = TRUE),
    Ratio = Overexpressed / Count,
    .groups = "drop"
  ) %>%
  filter(Count > 0)

plot_data$Cluster <- factor(plot_data$Cluster, levels = c(
  "N/A", "TGF-b dominant (C6)", "Lymphocyte depleted (C4)",
  "Immunologically quiet (C5)", "Wound healing (C1)",
  "IFN-y dominant (C2)", "Inflammatory (C3)"
))

# Ensure X-axis is perfectly alphabetical left-to-right
plot_data$Cancer <- factor(plot_data$Cancer, levels = sort(unique(plot_data$Cancer)))

cat("6. Building the Polished 300 DPI Bubble Plot...\n")
pB <- ggplot(plot_data, aes(x = Cancer, y = Cluster)) +
  geom_point(aes(size = Count, color = Ratio)) +
  scale_color_gradient(low = "#B0C4DE", high = "#0066CC", 
                       name = "Ratio of overexpressed samples in each cluster") +
  scale_size_continuous(range = c(1.5, 6), 
                        name = "Number of samples in each cluster for EGFR") +
  guides(size = guide_legend(override.aes = list(color = "#0066CC"))) +
  theme_bw(base_size = 9) + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major = element_line(color = "grey90", linetype = "solid"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(margin = margin(r = 10)),
    # Add bottom margin so the X-axis label breathes
    axis.title.x = element_text(margin = margin(t = 15))
  ) +
  labs(x = "Cancer Types", y = "Clusters")

cancers_ordered <- levels(plot_data$Cancer)
gbm_pos <- which(cancers_ordered == "GBM")
lgg_pos <- which(cancers_ordered == "LGG")

if(length(gbm_pos) > 0) pB <- pB + annotate("rect", xmin = gbm_pos - 0.5, xmax = gbm_pos + 0.5, ymin = 0.5, ymax = 7.5, color = "darkblue", fill = NA, linetype = "dashed")
if(length(lgg_pos) > 0) pB <- pB + annotate("rect", xmin = lgg_pos - 0.5, xmax = lgg_pos + 0.5, ymin = 0.5, ymax = 7.5, color = "darkblue", fill = NA, linetype = "dashed")

cat("7. Saving High-Res Image (300 DPI)...\n")
# Dimensions strictly doubled to preserve exact visual proportions but yield 300 DPI
png("Fig6B_Custom_ImmuneBubble_Polished_300DPI.png", width = 3000, height = 1000, res = 300)
print(pB)
dev.off()

rm(full_data, plot_data, pB, immune_data)
gc()
cat("\n--- PANEL B SCIENTIFICALLY PERFECTED! ---\n")
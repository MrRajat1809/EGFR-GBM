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
cat("Starting Final Aesthetic Polish for Panel A (Scientific Xena Edition)...\n")

gene_target <- "EGFR"

cat("1. Loading Local Xena Pan-Cancer Data...\n")
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

cat("3. Defining True Normals and Tumors...\n")
df_expr$Sample_Type <- substr(df_expr$sample, 14, 15)
df_expr <- df_expr[df_expr$Sample_Type %in% c("01", "11"), ]
df_expr$Group <- ifelse(df_expr$Sample_Type == "11", "Normal", "Tumor")

df_expr$Patient_ID <- substr(df_expr$sample, 1, 12)
df_expr <- df_expr[!duplicated(df_expr$sample), ]

cat("4. Merging Clinical Metadata...\n")
clin_sub <- clin_data[, c("sample", "cancer type abbreviation", "ajcc_pathologic_tumor_stage", "clinical_stage", "histological_grade")]
colnames(clin_sub) <- c("Patient_ID", "Cancer", "path_stage", "clin_stage", "grade")
clin_sub$Patient_ID <- substr(clin_sub$Patient_ID, 1, 12)
clin_sub <- clin_sub[!duplicated(clin_sub$Patient_ID), ]

merged <- merge(df_expr, clin_sub, by="Patient_ID", all.x=TRUE)
merged <- merged[!is.na(merged$Cancer) & merged$Cancer != "LAML", ]

cat("5. Scientifically Rigorous Staging (WHO Grading + AJCC)...\n")
merged$stage_consol <- ifelse(!is.na(merged$path_stage) & merged$path_stage != "" & merged$path_stage != "[Not Applicable]", merged$path_stage, merged$clin_stage)

full_data <- merged %>%
  mutate(
    Stage_Cat = case_when(
      Group == "Normal" ~ "Normal",
      Cancer == "GBM" & Group == "Tumor" ~ "Advanced",
      Cancer == "LGG" & Group == "Tumor" & grepl("G2|Grade II\\b", grade, ignore.case=TRUE) ~ "Early",
      Cancer == "LGG" & Group == "Tumor" & grepl("G3|Grade III\\b|G4|Grade IV\\b", grade, ignore.case=TRUE) ~ "Advanced",
      grepl("Stage I$|Stage I[A-C]$|Stage II$|Stage II[A-C]$", stage_consol, ignore.case = TRUE) ~ "Early",
      grepl("Stage III|Stage IV", stage_consol, ignore.case = TRUE) ~ "Advanced",
      TRUE ~ "No staging information"
    )
  )

full_data$Stage_Cat <- factor(full_data$Stage_Cat, levels = c("Early", "Advanced", "Normal", "No staging information"))

cat("6. Downsampling Dots...\n")
set.seed(42)
dots_data <- full_data %>%
  group_by(Cancer, Stage_Cat) %>%
  sample_n(size = min(n(), 50)) %>% 
  ungroup()

cat("7. Building the aesthetically perfect ggplot with Span Arrow...\n")
custom_colors <- c("Early" = "#E6B800",       
                   "Advanced" = "#C2185B",    
                   "Normal" = "#005b96",      
                   "No staging information" = "#9E9E9E") 

# Dynamically calculate the number of cancers and a safe Y-position for the arrow
num_cancers <- length(unique(full_data$Cancer))
arrow_y <- min(full_data$Expression, na.rm = TRUE) - 2.5

pA <- ggplot() +
  geom_boxplot(data = full_data, aes(x = Cancer, y = Expression), 
               outlier.shape = NA, fill = NA, color = "grey70", width = 0.6) +
               
  geom_jitter(data = dots_data, aes(x = Cancer, y = Expression, color = Stage_Cat), 
              position = position_jitter(width = 0.15), size = 0.6, alpha = 0.8) +
              
  scale_color_manual(values = custom_colors, name = "Stages") +
  
           
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "left",
    legend.title = element_text(margin = margin(b = 5)),
    # NEW: Increased bottom margin to prevent clipping the arrow
    plot.margin = margin(t = 20, r = 20, b = 45, l = 20),
    # NEW: Pushed the title down slightly more to clear the arrow
    axis.title.x = element_text(margin = margin(t = 25))
  ) +
  # NEW: Allow drawing outside the standard plot boundaries
  coord_cartesian(clip = "off") +
  labs(x = "Cancer Types", y = paste0(gene_target, " expression")) +
  guides(color = guide_legend(override.aes = list(size = 4)))

cancers_ordered <- sort(unique(full_data$Cancer))
gbm_pos <- which(cancers_ordered == "GBM")
lgg_pos <- which(cancers_ordered == "LGG")

if(length(gbm_pos) > 0) pA <- pA + annotate("rect", xmin = gbm_pos - 0.5, xmax = gbm_pos + 0.5, ymin = -Inf, ymax = Inf, color = "darkblue", fill = NA, linetype = "dashed")
if(length(lgg_pos) > 0) pA <- pA + annotate("rect", xmin = lgg_pos - 0.5, xmax = lgg_pos + 0.5, ymin = -Inf, ymax = Inf, color = "darkblue", fill = NA, linetype = "dashed")

cat("8. Saving High-Res Image (300 DPI)...\n")
png("Fig6A_Custom_PanStage_Polished_300DPI.png", width = 2800, height = 900, res = 300)
print(pA)
dev.off()

rm(full_data, dots_data, pA)
gc()

cat("\n--- PANEL A SCIENTIFICALLY PERFECTED! ---\n")

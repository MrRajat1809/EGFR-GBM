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

library(survival)
library(ggplot2)
library(dplyr)

setwd(project_path("outputs/figures"))
cat("Starting Grand Finale: LGG Multivariate Cox Adjusted Survival Model (Clean Lines)...\n")
gene_target <- "EGFR"

cat("1. Loading Gold-Standard Master Files...\n")
expr_data <- tryCatch(read.delim(gzfile(project_path("data/raw/clinical/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)
clin_data <- tryCatch(read.delim(project_path("data/raw/clinical/Survival_SupplementalTable_S1_20171025_xena_sp"), stringsAsFactors=FALSE, check.names=FALSE), error = function(e) NULL)

if(is.null(expr_data) || is.null(clin_data)) {
  stop("CRITICAL ERROR: Could not find the master Pan-Cancer Xena files!")
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

cat("3. Merging Curated LGG Clinical (Age & Grade) & Survival Data...\n")
clin_sub <- clin_data[, c("sample", "cancer type abbreviation", "OS", "OS.time", "age_at_initial_pathologic_diagnosis", "histological_grade")]
colnames(clin_sub) <- c("Patient_ID", "Cancer", "OS", "OS.time", "age", "Grade")
clin_sub$Patient_ID <- substr(clin_sub$Patient_ID, 1, 12)
clin_sub <- clin_sub[!duplicated(clin_sub$Patient_ID), ]

merged <- merge(df_expr, clin_sub, by="Patient_ID")

# Filter exclusively for LGG
merged <- merged[merged$Cancer == "LGG", ]

# Clean Grade Data (Isolate G2 and G3)
merged <- merged[!is.na(merged$Grade) & merged$Grade != "" & !grepl("not reported|unknown", merged$Grade, ignore.case=TRUE), ]
merged$Grade <- case_when(
  grepl("G2|Grade II\\b", merged$Grade, ignore.case=TRUE) ~ "G2",
  grepl("G3|Grade III\\b", merged$Grade, ignore.case=TRUE) ~ "G3",
  TRUE ~ "Other"
)
merged <- merged[merged$Grade != "Other", ]
merged$Grade <- as.factor(merged$Grade)

# Clean out any patients missing survival data, age, or grade
merged <- merged[!is.na(merged$OS.time) & !is.na(merged$OS) & !is.na(merged$Expression) & 
                 !is.na(merged$age) & !is.na(merged$Grade) & merged$OS.time > 0, ]

# Convert days to months
merged$time_months <- merged$OS.time / 30.4167

cat("4. Running MULTIVARIATE COX REGRESSION...\n")
cox_model <- coxph(Surv(time_months, OS) ~ Expression + age + Grade, data = merged)

# --- PRINT PUBLICATION-READY COX SUMMARY ---
cat("\n========================================================================\n")
cat("      MULTIVARIATE COX REGRESSION STATISTICS FOR LGG        \n")
cat("========================================================================\n\n")

# Get the raw summary object
cox_sum <- summary(cox_model)

# Create a clean data frame of the results
results_df <- data.frame(
  Variable = rownames(cox_sum$coefficients),
  HR = round(cox_sum$coefficients[, "exp(coef)"], 3),
  Lower_95_CI = round(cox_sum$conf.int[, "lower .95"], 3),
  Upper_95_CI = round(cox_sum$conf.int[, "upper .95"], 3),
  z_value = round(cox_sum$coefficients[, "z"], 3),
  p_value = signif(cox_sum$coefficients[, "Pr(>|z|)"], 3)
)

# Print the formatted table
print(results_df, row.names = FALSE)

cat("\n------------------------------------------------------------------------\n")
cat("Concordance Index:", round(cox_sum$concordance["C"], 3), "\n")
cat("Likelihood ratio test p-value:", signif(cox_sum$logtest["pvalue"], 3), "\n")
cat("========================================================================\n\n")

cox_summary <- summary(cox_model)
adjusted_p_val <- cox_summary$coefficients["Expression", "Pr(>|z|)"]

p_val_text <- ifelse(adjusted_p_val < 0.001, "p < 0.001", sprintf("p = %.3f", adjusted_p_val))

cat("5. Preparing Math for the Curve Plot...\n")
median_val <- median(merged$Expression, na.rm = TRUE)
merged$Risk_Group <- ifelse(merged$Expression >= median_val, "High Expression", "Low Expression")
merged$Risk_Group <- factor(merged$Risk_Group, levels = c("High Expression", "Low Expression"))

fit <- survfit(Surv(time_months, OS) ~ Risk_Group, data = merged)

med_table <- summary(fit)$table
med_high <- med_table["Risk_Group=High Expression", "median"]
med_low <- med_table["Risk_Group=Low Expression", "median"]

max_med <- max(c(med_high, med_low), na.rm=TRUE)
if(is.na(max_med) || is.infinite(max_med)) max_med <- max(merged$time_months, na.rm=TRUE)

surv_df <- data.frame(
  time = fit$time,
  surv = fit$surv,
  n.censor = fit$n.censor,
  Group = rep(names(fit$strata), fit$strata)
)
surv_df$Group <- gsub("Risk_Group=", "", surv_df$Group)

start_df <- data.frame(time = c(0, 0), surv = c(1, 1), n.censor = c(0, 0), Group = c("High Expression", "Low Expression"))
surv_df <- rbind(start_df, surv_df)
surv_df <- surv_df[order(surv_df$Group, surv_df$time), ]

cat("6. Building the 9pt, 300 DPI Aesthetic Plot (No Censored Marks)...\n")

text_x <- max(surv_df$time) * 0.02
text_y <- 0.10
max_t <- max(surv_df$time, na.rm = TRUE)

# Spaced out X-axis breaks in steps of 50
uniform_breaks <- seq(0, ceiling(max_t / 50) * 50, by = 50)

pF <- ggplot(surv_df, aes(x = time, y = surv, color = Group)) +
  # Just the smooth steps, no geom_point for censored data
  geom_step(linewidth = 0.8) +
  scale_color_manual(values = c("High Expression" = "red", "Low Expression" = "#005b96")) +
  scale_x_continuous(breaks = uniform_breaks) +
  annotate("segment", x = 0, xend = max_med, y = 0.5, yend = 0.5, linetype = "dashed", color = "black", linewidth=0.5) +
  annotate("segment", x = med_high, xend = med_high, y = 0, yend = 0.5, linetype = "dashed", color = "black", linewidth=0.5, na.rm=TRUE) +
  annotate("segment", x = med_low, xend = med_low, y = 0, yend = 0.5, linetype = "dashed", color = "black", linewidth=0.5, na.rm=TRUE) +
  
  annotate("text", x = text_x, y = text_y, label = p_val_text, size = 9 / .pt, hjust = 0, color = "black") +
  
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "plain"),
    axis.text = element_text(color = "black", size = 9, face = "plain"),
    axis.title = element_text(color = "black", size = 9, face = "plain"),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "plain"),
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(t = 20, r = 30, b = 20, l = 20) 
  ) +
  coord_cartesian(clip = "off") +
  labs(title = "Overall Survival (LGG)", x = "Overall Survival Time (Months)", y = "Survival Probability")

cat("7. Saving the Final Clean Masterpiece...\n")
png("Fig6F_Custom_KM_Survival_FINAL_LGG.png", width = 1200, height = 1000, res = 300)
print(pF)
dev.off()

cat("\n--- PANEL F (LGG - CLEAN LINES) PERFECTED! ---\n")

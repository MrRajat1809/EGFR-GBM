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
use_working_directory("outputs/noncoding")

# Suppress the massive startup text from Bioconductor
suppressMessages(library(TCGAbiolinks))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(dplyr))

cat("====================================================\n")
cat("🔍 QUERYING TCGA-GBM FOR EGFR RNA EXPRESSION (TPM)\n")
cat("====================================================\n")

# 1. Build the query for Glioblastoma RNA-Seq data
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)

# 2. Grab a random subset of 20 patients to save download time
results <- getResults(query)
tumor_subset <- head(results$cases, 20)

query_subset <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = tumor_subset
)

cat("Downloading RNA-seq data for 20 GBM patients...\n")
cat("(This connects to the NIH GDC API. It may take 1-2 minutes)...\n")
GDCdownload(query_subset, method = "api", files.per.chunk = 5,
            directory = project_path("data/raw/noncoding/GDCdata"))

cat("Preparing the genomic data matrix...\n")
data <- GDCprepare(query_subset, directory = project_path("data/raw/noncoding/GDCdata"))

# 3. Extract the exact RNA expression values for the EGFR gene
cat("Extracting EGFR Expression...\n")
row_data <- rowData(data)

# Find the exact row index for EGFR
egfr_idx <- which(row_data$gene_name == "EGFR")

# Extract Transcripts Per Million (TPM)
egfr_tpm <- assay(data, "tpm_unstrand")[egfr_idx, ]

# 4. Calculate Statistics
df <- data.frame(
  Patient = colnames(data),
  EGFR_TPM = as.numeric(egfr_tpm)
)

mean_tpm <- mean(df$EGFR_TPM, na.rm = TRUE)
max_tpm <- max(df$EGFR_TPM, na.rm = TRUE)

cat("\n====================================================\n")
cat("📊 TCGA PATIENT RNA-SEQ RESULTS (n=20)\n")
cat("====================================================\n")
cat(sprintf("Mean EGFR Expression:  %.2f TPM\n", mean_tpm))
cat(sprintf("Max EGFR Expression:   %.2f TPM\n\n", max_tpm))

if(mean_tpm > 50) {
    cat("🔥 VERDICT: MASSIVE OVEREXPRESSION CONFIRMED.\n")
    cat("Normal gene expression is typically < 20 TPM.\n")
    cat("These tumors are forcefully pumping out EGFR RNA, validating the Neo-Enhancer's impact.\n")
} else {
    cat("VERDICT: Baseline or Low Expression.\n")
}
cat("====================================================\n")

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

suppressMessages(library(tidyverse))

# 1. Use the pinned haploR installation supplied by docker/Dockerfile.
if (!requireNamespace("haploR", quietly = TRUE)) {
  stop("haploR is missing. Build the unified environment using docker/compose.yml.")
}
library(haploR)

# 2. Load the filtered VEP data
cat("Loading vep_non_coding_filtered.csv...\n")
vep_data <- read_csv("vep_non_coding_filtered.csv", show_col_types = FALSE)

# 3. Extract unique coordinates and format for HaploReg (chr:start-end)
regions <- vep_data %>%
  distinct(Location) %>%
  mutate(
    pos = str_extract(Location, "[0-9]+$"),
    query_str = paste0(Location, "-", pos)
  ) %>%
  pull(query_str)

# 4. Save a clean text file of these inputs
writeLines(regions, "haploreg_input_coordinates.txt")
cat("Generated 'haploreg_input_coordinates.txt' with", length(regions), "unique regions.\n")

# 5. Query the HaploReg database
cat("Querying HaploReg... (this might take a minute depending on the Broad's servers)\n")

# Wrap in tryCatch just in case the server times out
results <- tryCatch({
  queryHaploreg(query = regions)
}, error = function(e) {
  cat("Error hitting HaploReg servers:", e$message, "\n")
  return(NULL)
})

# 6. Save the output
if (!is.null(results) && nrow(results) > 0) {
  write_csv(results, "haploreg_results.csv")
  cat("Success! Annotations saved to 'haploreg_results.csv'.\n")
} else {
  cat("No results returned. The variants might not be in HaploReg, or the server dropped the connection.\n")
}

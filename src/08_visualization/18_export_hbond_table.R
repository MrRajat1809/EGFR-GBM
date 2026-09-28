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

library(dplyr)

# Custom function to safely read GROMACS .xvg files
read_xvg_hbond <- function(filepath, system_name) {
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^[#@]", lines)]
  df <- read.table(text = data_lines)
  df <- df[, 1:2]
  colnames(df) <- c("Time", "Hbonds")
  df$System <- system_name
  return(df)
}

# Reach into the GROMACS folders to grab the data
wt_hbond <- read_xvg_hbond("../gromacs/WT_Erlotinib/hbnum_200.xvg", "Wild Type")
v774m_hbond <- read_xvg_hbond("../gromacs/V774M_Erlotinib/hbnum_200.xvg", "V774M Mutant")
l861q_hbond <- read_xvg_hbond("../gromacs/L861Q_Erlotinib/hbnum_200.xvg", "L861Q Mutant")

# Combine all three datasets
master_hbond <- bind_rows(wt_hbond, v774m_hbond, l861q_hbond)

# Calculate summary statistics for the table
summary_table <- master_hbond %>%
  group_by(System) %>%
  summarise(
    `Average H-Bonds` = round(mean(Hbonds), 3),
    `Max H-Bonds observed` = max(Hbonds),
    # Calculate exactly what percentage of the simulation had at least 1 bond
    `% Time Bound (>= 1 H-Bond)` = paste0(round((sum(Hbonds >= 1) / n()) * 100, 1), "%")
  ) %>%
  # Arrange to match your previous graphs (WT first)
  arrange(factor(System, levels = c("Wild Type", "V774M Mutant", "L861Q Mutant")))

# Print the table right in the R console so you can see it
print(summary_table)

# Export to a clean CSV file that you can open in Excel or Word
write.csv(summary_table, "HBond_Summary_Table.csv", row.names = FALSE)
print("Success! Data exported to HBond_Summary_Table.csv")
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
use_working_directory(".")

library(ggplot2)
library(tidyr)
library(dplyr)

# --- 1. DATA PREP ---
data <- read.csv(project_path("outputs/standard_master_table/EGFR_Arm2_Master_Standardized.csv"))

tools_cols <- c("AlphaMissense_Standard", "VARITY_Standard", "SuSPect_Standard", "MetaLR_Standard", "MetaSVM_Standard", "MetaSNP_Standard", "PHDSNPg_Standard")
tools_clean <- c("AlphaMissense", "VARITY", "SuSPect", "MetaLR", "MetaSVM", "MetaSNP", "PHDSNPg")

data <- data %>% select(all_of(tools_cols))
colnames(data) <- tools_clean

summary_list <- lapply(tools_clean, function(t) {
  counts <- table(data[[t]])
  deleterious <- ifelse("Deleterious" %in% names(counts), counts["Deleterious"], 0)
  neutral <- ifelse("Benign" %in% names(counts), counts["Benign"], 0)
  data.frame(Tool = t, Deleterious = as.numeric(deleterious), Neutral = as.numeric(neutral))
})

summary_df <- do.call(rbind, summary_list)

plot_data <- summary_df %>%
  pivot_longer(cols = c("Deleterious", "Neutral"), names_to = "Classification", values_to = "Count") %>%
  group_by(Tool) %>%
  mutate(is_upper = Count == max(Count),
         vjust_pos = ifelse(is_upper, -1.2, 2.2)) %>%
  ungroup()

# --- 2. NEW AESTHETICS FOR ARM 2 ---
col_del <- "#E41A1C" 
col_neu <- "#4DAF4A" 
col_line <- "#2E2E2E" 

# --- 3. PLOTTING ---
p <- ggplot(plot_data, aes(x = Tool, y = Count)) +
  # Set global base size to 17
  theme_bw(base_size = 17) + 
  
  geom_segment(data = summary_df, 
               aes(x = Tool, xend = Tool, y = Neutral, yend = Deleterious), 
               color = col_line, linewidth = 1.2, inherit.aes = FALSE) +
  
  geom_point(aes(color = Classification), size = 6) +
  
  # Adjusted size to match ~17pt in mm
  geom_text(aes(label = Count, color = Classification, vjust = vjust_pos), 
            size = 6.0, show.legend = FALSE) +
  
  scale_color_manual(values = c("Deleterious" = col_del, "Neutral" = col_neu)) +
  
  # Reduced grid density: breaks every 10 instead of 5
  scale_y_continuous(breaks = seq(0, 60, by = 10), expand = expansion(mult = c(0.1, 0.2))) +
  
  theme(
    # Clean, minimal grid
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4), 
    panel.grid.minor = element_blank(), # Removed annoying minor grid lines
    panel.border = element_blank(),     
    
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_line(color = "black", linewidth = 0.6), 
    
    # Standardized size to 17
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 17, color = "black"),
    axis.text.y = element_text(size = 17, color = "black"),
    axis.title = element_text(size = 17),
    
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 17)
  ) +
  labs(x = NULL, y = "Number of Variants")

# --- 4. EXPORT ---
output_file <- project_path("outputs/figures/EGFR_Arm2_Disease_Tendency_Rangeplot.png")
ggsave(output_file, plot = p, width = 10, height = 7, dpi = 300, bg = "white")

cat("Success! Arm 2 Range plot saved to", output_file, "\n")
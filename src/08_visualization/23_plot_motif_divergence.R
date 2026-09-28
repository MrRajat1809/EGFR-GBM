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

# Load libraries pre-installed in r_gbm
suppressMessages(library(tidyverse))

cat("Loading JASPAR significant deltas...\n")
df <- read_csv("jaspar_significant_deltas.csv", show_col_types = FALSE)

# 1. Filter and clean the data to get the absolute top hits
# We want the Top 5 Gains (positive) and Top 5 Losses (negative)
top_gains <- df %>% 
  filter(delta > 0) %>% 
  arrange(desc(delta)) %>% 
  head(6) # Pulling 6 to ensure we get a clean mix of TFs

top_losses <- df %>% 
  filter(delta < 0) %>% 
  arrange(delta) %>% 
  head(6)

# Combine and create a clean factor for plotting
plot_data <- bind_rows(top_gains, top_losses) %>%
  mutate(
    # Clean up the TF names (remove MA.xxx.)
    TF_Name = str_replace(name, "MA[0-9]+\\.[0-9]+\\.", ""),
    # Create a Status column for coloring
    Status = ifelse(delta > 0, "Gain of Binding (Oncogenes)", "Loss of Binding (Tumor Suppressors)"),
    # Add the locus to the TF name for the label
    Label = paste0(TF_Name, " (", str_extract(variant, "[0-9]+"), ")")
  )

# Reorder the labels based on delta value so the plot flows nicely
plot_data$Label <- reorder(plot_data$Label, plot_data$delta)

# 2. Build the ggplot
cat("Generating publication-ready diverging plot...\n")
p <- ggplot(plot_data, aes(x = Label, y = delta, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  coord_flip() + # Flip to horizontal for readability
  scale_fill_manual(values = c("Gain of Binding (Oncogenes)" = "#e74c3c", 
                               "Loss of Binding (Tumor Suppressors)" = "#3498db")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30", margin = margin(b=15)),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t=10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Neo-Enhancer TF Hijacking in GBM",
    subtitle = "In-Silico Binding Affinity Deltas (JASPAR)",
    y = "Relative Binding Affinity Shift (%)"
  ) +
  # Add a vertical line at 0
  geom_hline(yintercept = 0, color = "black", size = 0.8)

# 3. Save the plot using the rsvg dependency available in r_gbm
ggsave("Motif_Divergence_Plot.pdf", plot = p, width = 10, height = 7, dpi = 300)
ggsave("Motif_Divergence_Plot.png", plot = p, width = 10, height = 7, dpi = 300)

cat("Success! Saved as Motif_Divergence_Plot.pdf and .png\n")
